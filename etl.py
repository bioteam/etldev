#!/usr/bin/env python
import argparse
import yaml
import sys
import re
import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import csv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="ETL config file")
    parser.add_argument(
        "-d", "--dictionary", help="file containing data dictionary"
    )
    parser.add_argument("-i", "--input", help="file data to be ETL'd")
    args = parser.parse_args()
    return args


class ETLdbGap:
    def __init__(self, config):
        self.config = config
        self._data_dictionary = []
        self._map_phenotype_to_concept = []
        self._data = []
        self._variables = {}

    def read_data_dictionary(self, dictfile):
        try:
            self.config["dictformat"]
        except KeyError:
            format = "areds"  # This is the default
        else:
            format = self.config["dictformat"]
        if format == "areds2":
            self.read_areds2_data_dictionary(dictfile)
        else:
            self.read_areds_data_dictionary(dictfile)

    #
    # The AREDS dictionaries may have leading comment lines
    #
    def read_areds_data_dictionary(self, dictfile):
        with open(dictfile, encoding="latin-1") as csvfile:
            ptr = csvfile.tell()
            line = csvfile.readline()
            # Skip leading comment lines
            while line.startswith("!#!"):
                ptr = csvfile.tell()
                line = csvfile.readline()
            csvfile.seek(ptr)
            reader = csv.DictReader(csvfile)
            for row in reader:
                trimmed_row = {}
                for key, value in row.items():
                    trimmed_row[
                        key
                    ] = value.strip()  # Trimming leading/trailing whitespace
                self._data_dictionary.append(trimmed_row)

            # for row in reader:
            # Skip empty lines
            # if list(row.values())[0]:
            #   self._data_dictionary.append(row)

    #
    # The AREDS2 dictionaries have enumerated values separated by commas
    #
    def read_areds2_data_dictionary(self, dictfile):
        with open(dictfile, "r", encoding="utf-8-sig") as csvfile:
            ptr = csvfile.tell()
            #
            # Read headers explicitly and let DictReader make a list of
            # the enumerated values
            #
            line = csvfile.readline()
            names = line.split(sep=",")[:16]

            reader = csv.DictReader(
                csvfile, fieldnames=names, restkey="VALUES"
            )
            for row in reader:
                # Skip empty lines
                if list(row.values())[0]:
                    self._data_dictionary.append(row)

    def write_concepts(self, conceptsfile):
        split_data = []
        for row in self._data_dictionary:
            varname = row[self.config["varname"]]
            if varname == self.config["patientid"]:
                continue
            vartype = row[self.config["typename"]]
            if type(row[self.config["enumname"]]) is list:
                if "encoded value" in vartype:
                    values = list(filter(len, row[self.config["enumname"]]))
                else:
                    values = row[self.config["enumname"]][0]
            else:
                values = re.split(
                    self.config["separator"], row[self.config["enumname"]]
                )

            if self.config["description"] != "":
                description = row[self.config["description"]]
                path = "/".join(["", self.config["pathroot"], description, ""])
            else:
                path = "/".join(["", self.config["pathroot"], varname, ""])
            # Check what i2b2 type of variable it is: integer, float,
            # string or large-string
            if len(values) <= 1:
                if vartype.lower() == "num" or vartype.lower() == "integer":
                    i2b2vartype = "integer"
                elif (
                    vartype.lower() == "decimal" or vartype.lower() == "float"
                ):
                    i2b2vartype = "float"
                else:
                    i2b2vartype = "string"

                i2b2code = varname
                dbgap_code_id = -1
                conceptpath = path
                split_data.append((conceptpath, i2b2code, i2b2vartype))
                self._map_phenotype_to_concept.append(
                    (
                        conceptpath,
                        i2b2code,
                        i2b2vartype,
                        dbgap_code_id,
                        varname,
                    )
                )
            else:
                for i, value in enumerate(values):
                    value.replace('"', "")
                    # Limit the split to 1 as some descriptions contain "="
                    clin_name = value.split("=", 1)
                    # Remove spaces
                    clin_name = [x.strip(" ") for x in clin_name]
                    dbgap_code_id = ""
                    i2b2conceptlabel = ""
                    i2b2concept = ""
                    if len(clin_name) > 1:
                        dbgap_code_id = clin_name[0].replace('"', "").lstrip()
                        i2b2concept = clin_name[1].replace('"', "")
                        i2b2concept = i2b2concept.replace(",", " or ")
                        i2b2conceptlabel = "".join(
                            filter(str.isalnum, i2b2concept)
                        )
                    if (
                        i2b2conceptlabel == dbgap_code_id
                    ):  # TODO: Ensure i2b2code is unique in the ontology
                        varcode = "".join(filter(str.isalnum, i2b2concept))
                    else:
                        varcode = (
                            "".join(filter(str.isalnum, i2b2concept))
                            + dbgap_code_id
                        )
                    varname4i2b2 = "".join(varname.split())
                    if (len(varname4i2b2) + len(varcode)) > 50:
                        truncate = 50 - len(varname4i2b2)
                        i2b2code = varname4i2b2 + varcode[-truncate:]
                    else:
                        i2b2code = varname4i2b2 + varcode
                    conceptpath = path + i2b2concept
                    split_data.append((conceptpath, i2b2code, "assertion"))
                    self._map_phenotype_to_concept.append(
                        (
                            conceptpath,
                            i2b2code,
                            "assertion",
                            dbgap_code_id,
                            varname,
                        )
                    )
        with open(conceptsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["path", "code", "type"])
            for row in split_data:
                writer.writerow(row)

    def read_facts(self, phenocsvfile):
        with open(phenocsvfile, "r", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            self._data = list(reader)
        i = 0
        while i < len(self._data[0]):
            self._variables[self._data[0][i]] = i
            i += 1

    def add_time(self, visitdateformat, beginDate, timediff):
        startdate = beginDate
        if visitdateformat == 2:  # days
            startdate = beginDate + datetime.timedelta(days=timediff)
        elif visitdateformat == 3:  # months
            startdate = beginDate + relativedelta(
                months=int(timediff)
            )  # relativedelta does not support Non-integer years and months which are ambiguous and not currently supported.
        elif visitdateformat == 4:  # years
            startdate = beginDate + datetime.timedelta(days=timediff * 365)
        return startdate

    #
    # Timestamps for individual facts.
    #
    # Mode 0: all times are self.config["basedate"]
    # Mode 1a: all times are self.config["timevar"]["default"] + "basedate"
    # Mode 1b: facts are tied to certain time variables via explicit config
    # Mode 2: facts are tied to certain time variables via regex
    # Mode 3: a list of possible additional time points, represented as deltas
    # Mode 4: a list of possible additional time points, relative to "basedate"
    # Mode 5: time is defined by a visit variable that is parsable using a regex, e.g. F04, for 4 monrths visit
    def fact_time(self, i, j):
        visitbaselinedate = self.config["basedate"]
        visitdateformat = int(self.config["dateformat"])
        beginDate = datetime.datetime.strptime(visitbaselinedate, "%d/%m/%Y")
        if int(self.config["datemode"]) == 0:
            return beginDate
        varname = list(self._variables.keys())[
            list(self._variables.values()).index(j)
        ]
        if (self.config["datemode"]) == 1:
            defaulttimevar = self.config["timevar"]["default"]
            # Is there a specific time variable for this variable?
            try:
                self.config["timevar"][varname]
            except KeyError:
                timevar = defaulttimevar
            else:
                timevar = self.config["timevar"][varname]
            if self._data[i][self._variables[timevar]].isalnum():
                timediff = float(self._data[i][self._variables[timevar]])
            else:
                timediff = float(
                    self._data[i][self._variables[defaulttimevar]]
                )

            startdate = self.add_time(visitdateformat, beginDate, timediff)
            return startdate.strftime("%Y-%m-%d")
        elif (self.config["datemode"]) == 2:
            startdates = []
            for tv in self.config["timevar"]:
                tvre = self.config["timevar"][tv]
                timediff = float(self._data[i][self._variables[tv]])
                startdate = self.add_time(visitdateformat, beginDate, timediff)
                if re.search(tvre, self._data[0][j]):
                    return startdate.strftime("%Y-%m-%d")
                startdates.append(startdate)

            print(f"No matching time for: {self._data[0][j]!r}")
            if len(startdates):
                laststartdate = max(startdates)
                return laststartdate.strftime("%Y-%m-%d")
            else:
                return -1
        elif (self.config["datemode"]) == 3:
            defaulttimevar = self.config["timevar"]["default"]
            addltimevars = self.config["additionaltimevar"]
            timediff = float(self._data[i][self._variables[defaulttimevar]])
            addltime = 0
            for tv in addltimevars:
                if self._data[i][self._variables[tv]].isnumeric():
                    timeval = int(self._data[i][self._variables[tv]])
                    addltime = max(addltime, timeval)
            additionaldatedifftimeunits = int(
                self.config["additionaldatedifftimeunits"]
            )
            startdate = self.add_time(visitdateformat, beginDate, timediff)
            startdate = self.add_time(
                additionaldatedifftimeunits, startdate, addltime
            )
            return startdate.strftime("%Y-%m-%d")
        elif (self.config["datemode"]) == 4:
            defaulttimevar = self.config["timevar"]["default"]
            addltimevars = self.config["additionaltimevar"]
            timediff = float(self._data[i][self._variables[defaulttimevar]])
            addltime = 0
            for tv in addltimevars:
                if self._data[i][self._variables[tv]].isnumeric():
                    timeval = int(self._data[i][self._variables[tv]])
                    addltime = max(addltime, timeval)
            timediff = max(timediff, addltime)
            startdate = self.add_time(visitdateformat, beginDate, timediff)
            return startdate.strftime("%Y-%m-%d")
        elif (self.config["datemode"]) == 5:
            datevar = self.config["timevar"]["default"]
            timediff = 0
            visitval = self._data[i][self._variables[datevar]]
            match = re.search(r"\d+", visitval)
            if match:
                timediff = int(match.group())
            startdate = self.add_time(visitdateformat, beginDate, timediff)
            return startdate.strftime("%Y-%m-%d")

    def collect_facts(self):
        facts = []
        code = ""
        value = ""
        dt_string = ""
        for i in range(len(self._data)):
            if i > 0:
                # Patient ID
                mrn = self._data[i][self._variables[self.config["patientid"]]]

                # Loop through cells in row
                # Should skip time variables
                if self.config["datemode"] == 2:
                    skiplist = list(self.config["timevar"].keys())
                else:
                    skiplist = list(self.config["timevar"].values())
                for j, value in enumerate(self._data[i]):
                    varname = list(self._variables.keys())[
                        list(self._variables.values()).index(j)
                    ]
                    if (
                        j == self._variables[self.config["patientid"]]
                        or self._data[i][j].strip() == ""
                        or varname in skiplist
                    ):
                        continue
                    else:
                        if self._data[i][j].strip() == "":
                            continue

                        code = "-"
                        for x in range(len(self._map_phenotype_to_concept)):
                            # Check if it is an enumerated value,
                            # then only add code
                            if (
                                self._map_phenotype_to_concept[x][3] == value
                                and self._map_phenotype_to_concept[x][4]
                                == self._data[0][j]
                            ):
                                code = self._map_phenotype_to_concept[x][1]
                                value = ""

                        if code == "-":
                            code = self._data[0][j]
                            value = self._data[i][j]

                        dt_string = self.fact_time(i, j)
                        facts.append((mrn, str(dt_string), code, value))

        return facts

    def write_facts(self, factsfile):
        facts = self.collect_facts()

        with open(factsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mrn", "start-date", "code", "value"])
            for row in facts:
                writer.writerow(row)


#
# Command-line arguments. Could add dictionary and input to config, but
# there are instances where the same config will work with different inputs
#
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="ETL config file")
    parser.add_argument(
        "-d", "--dictionary", help="file containing data dictionary"
    )
    parser.add_argument("-i", "--input", help="file data to be ETL'd")
    args = parser.parse_args()
    return args


#
# Configuration file with parsing specifications
#
def load_conf(config_file):
    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"The config file '{config_file}' does not exist.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error reading the config file: {e}")
        sys.exit(1)
    return config


def main():
    inputs = parse_args()
    etl_conf = load_conf(inputs.config)
    etl = ETLdbGap(etl_conf)
    etl.read_data_dictionary(inputs.dictionary)
    conceptsfile = "concepts_" + etl_conf["filebase"] + ".csv"
    etl.write_concepts(conceptsfile)
    etl.read_facts(inputs.input)
    factsfile = "facts_" + etl_conf["filebase"] + ".csv"
    etl.write_facts(factsfile)


if __name__ == "__main__":
    main()
