#!/usr/bin/env python
import argparse
import yaml
import sys
import re
import datetime
from datetime import timedelta
import csv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="ETL config file", default="etl.yml")
    parser.add_argument("-d", "--dictionary", help="file containing data dictionary")
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
            values = re.split(self.config["separator"], row[self.config["enumname"]])
            path = "/".join(["", self.config["pathroot"], varname, ""])
            # Check what i2b2 type of variable it is: integer, float,
            # string or large-string
            if len(values) == 1:
                if vartype.lower() == "num" or vartype.lower() == "integer":
                    i2b2vartype = "integer"
                elif vartype.lower() == "decimal" or vartype.lower() == "float":
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
                        i2b2conceptlabel = "".join(filter(str.isalnum, i2b2concept))
                    if (
                        i2b2conceptlabel == dbgap_code_id
                    ):  # TODO: Ensure i2b2code is unique in the ontology
                        varcode = "".join(filter(str.isalnum, i2b2concept))
                    else:
                        varcode = (
                            "".join(filter(str.isalnum, i2b2concept)) + dbgap_code_id
                        )
                    varname4i2b2 = "".join(varname.split())
                    if (len(varname4i2b2) + len(varcode)) > 50:
                        truncate = 50 - len(vavarname4i2b2rname)
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

    def write_adverse_facts(self, factsfile):
        mrn = ""
        startdate = datetime.datetime.now()
        code = ""
        value = ""

        # Collect facts for facts.csv
        facts = []
        datecolignore = -1
        code = ""
        value = ""
        dt_string = ""
        visitdatevarname = self.config["datevar"]
        visitdateformat = int(self.config["dateformat"])
        visitbaselinedate = self.config["basedate"]
        additionaldatediffvarname = self.config["additionaldatediffvarname"]
        additionaldatedifftimeunits = int(self.config["additionaldatedifftimeunits"])
        for i in range(len(self._data)):
            if i > 0:
                # Patient ID
                mrn = self._data[i][self._variables[self.config["patientid"]]]
                # print ("patient"+str(mrn))
                # Visit Date
                if visitdateformat == 1:  # Use date from column visitdatevarname
                    startdate = self._data[i][self._variables[visitdatevarname]]
                    datecolignore = self._variables[visitdatevarname]
                else:
                    timediff2 = 0
                    timediff = float(self._data[i][self._variables[visitdatevarname]])
                    datecolignore = self._variables[visitdatevarname]
                    # carry out conversion between string
                    # to datetime object
                    if additionaldatediffvarname != "NA":
                        if self._data[i][
                            self._variables[additionaldatediffvarname]
                        ].isnumeric():
                            timediff2 = float(
                                self._data[i][
                                    self._variables[additionaldatediffvarname]
                                ]
                            )

                    beginDate = datetime.datetime.strptime(
                        visitbaselinedate, "%d/%m/%Y"
                    )

                    if visitdateformat == 2:  # use visitbaselinedate and add Days
                        startdate = beginDate + datetime.timedelta(days=int(timediff))
                    elif visitdateformat == 3:  # use visitbaselinedate and add Months
                        startdate = beginDate + datetime.timedelta(
                            days=int(timediff * 12)
                        )
                    elif visitdateformat == 4:  # use visitbaselinedate and add Years
                        startdate = beginDate + datetime.timedelta(
                            days=int(timediff * 365)
                        )

                    if (
                        additionaldatedifftimeunits > 1
                        and additionaldatedifftimeunits < 5
                    ):
                        if (
                            additionaldatedifftimeunits == 2
                        ):  # use visitbaselinedate and add Days
                            startdate = startdate + datetime.timedelta(
                                days=int(timediff2)
                            )
                        elif (
                            additionaldatedifftimeunits == 3
                        ):  # use visitbaselinedate and add Months
                            startdate = startdate + datetime.timedelta(
                                days=int(timediff2 * 12)
                            )
                        elif (
                            additionaldatedifftimeunits == 4
                        ):  # use visitbaselinedate and add Years
                            startdate = startdate + datetime.timedelta(
                                days=int(timediff2 * 365)
                            )

                    dt_string = startdate.strftime("%Y-%m-%d")
                # Loop through cells in row
                for j, value in enumerate(self._data[i]):
                    if (
                        j == self._variables[self.config["patientid"]]
                        or j == datecolignore
                    ):
                        continue
                    else:
                        if self._data[i][j].strip() == "":
                            continue

                        code = "-"
                        for x in range(len(self._map_phenotype_to_concept)):
                            if x > 0:
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
                        facts.append((mrn, str(dt_string), code, value))
                        # print (mrn + ", " + str(dt_string) + ", " + str(code) + ", " + str (value) )

        # Write facts file

        # factsfile = 'facts_' + self.config['filebase'] + '.csv'
        with open(factsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mrn", "start-date", "code", "value"])
            for row in facts:
                writer.writerow(row)

    def write_fundus_facts(self, factsfile):
        mrn = ""
        startdate = datetime.datetime.now()
        code = ""
        value = ""

        # Collect facts for facts.csv
        facts = []
        datecolignore = -1
        code = ""
        value = ""
        dt_string = ""
        righteyedate = self.config["rightdate"]
        righteyere = r"\b(?:RE\w+|\w+RE)\b"
        lefteyedate = self.config["leftdate"]
        lefteyere = r"\b(?:LE\w+|\w+LE)\b"
        visitdateformat = int(self.config["dateformat"])
        visitbaselinedate = self.config["basedate"]
        for i in range(len(self._data)):
            if i > 0:
                # Patient ID
                mrn = self._data[i][self._variables[self.config["patientid"]]]

                beginDate = datetime.datetime.strptime(visitbaselinedate, "%d/%m/%Y")
                righttimediff = float(self._data[i][self._variables[righteyedate]])
                lefttimediff = float(self._data[i][self._variables[lefteyedate]])
                rightstartdate = beginDate + datetime.timedelta(
                    days=int(righttimediff * 365)
                )
                leftstartdate = beginDate + datetime.timedelta(
                    days=int(lefttimediff * 365)
                )
                laststartdate = max(rightstartdate, leftstartdate)

                # Loop through cells in row
                for j, value in enumerate(self._data[i]):
                    if (
                        j == self._variables[self.config["patientid"]]
                        or j == self._variables[righteyedate]
                        or j == self._variables[lefteyedate]
                    ):
                        continue
                    else:
                        if self._data[i][j].strip() == "":
                            continue

                        code = "-"
                        for x in range(len(self._map_phenotype_to_concept)):
                            if x > 0:
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
                        if re.search(righteyere, self._data[0][j]):
                            dt_string = rightstartdate.strftime("%Y-%m-%d")
                        elif re.search(lefteyere, self._data[0][j]):
                            dt_string = leftstartdate.strftime("%Y-%m-%d")
                        else:
                            dt_string = laststartdate.strftime("%Y-%m-%d")

                        facts.append((mrn, str(dt_string), code, value))
                        # print (mrn + ", " + str(dt_string) + ", " + str(code) + ", " + str (value) )

        # Write facts file

        # factsfile = 'facts_' + self.config['filebase'] + '.csv'
        with open(factsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mrn", "start-date", "code", "value"])
            for row in facts:
                writer.writerow(row)

    def write_mortality_facts(self, factsfile):
        mrn = ""
        startdate = datetime.datetime.now()
        code = ""
        value = ""

        # Collect facts for facts.csv
        facts = []
        datecolignore = -1
        code = ""
        value = ""
        dt_string = ""
        visitdatevarname = self.config["datevar"]
        visitdateformat = int(self.config["dateformat"])
        visitbaselinedate = self.config["basedate"]
        for i in range(len(self._data)):
            if i > 0:
                # Patient ID
                mrn = self._data[i][self._variables[self.config["patientid"]]]
                # Visit Date
                timediff = float(self._data[i][self._variables[visitdatevarname]])
                beginDate = datetime.datetime.strptime(visitbaselinedate, "%d/%m/%Y")
                startdate = beginDate + datetime.timedelta(days=int(timediff * 365))
                dt_string = startdate.strftime("%Y-%m-%d")
                # Loop through cells in row
                for j, value in enumerate(self._data[i]):
                    if (
                        j == self._variables[self.config["patientid"]]
                        or j == self._variables[self.config["datevar"]]
                    ):
                        continue
                    else:
                        if self._data[i][j].strip() == "":
                            continue

                        code = "-"
                        for x in range(len(self._map_phenotype_to_concept)):
                            if x > 0:
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
                        facts.append((mrn, str(dt_string), code, value))
                        # print (mrn + ", " + str(dt_string) + ", " + str(code) + ", " + str (value) )

        # Write facts file

        # factsfile = 'facts_' + self.config['filebase'] + '.csv'
        with open(factsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mrn", "start-date", "code", "value"])
            for row in facts:
                writer.writerow(row)

    def write_accord_key_facts(self, factsfile):
        mrn = ""
        startdate = datetime.datetime.now()
        code = ""
        value = ""

        # Collect facts for facts.csv
        facts = []
        datecolignore = -1
        code = ""
        value = ""
        dt_string = ""
        visitdatevarname = self.config["datevar"]
        visitdateformat = int(self.config["dateformat"])
        visitbaselinedate = self.config["basedate"]

        for i in range(len(self._data)):
            if i > 0:
                # Patient ID
                mrn = self._data[i][self._variables[self.config["patientid"]]]
                # print ("patient"+str(mrn))
                # Visit Date
                if visitdateformat == 1:  # Use date from column visitdatevarname
                    startdate = self._data[i][self._variables[visitdatevarname]]
                    datecolignore = self._variables[visitdatevarname]
                elif visitdateformat == 0:
                    startdate = datetime.datetime.strptime(
                        visitbaselinedate, "%d/%m/%Y"
                    )
                    dt_string = startdate.strftime("%Y-%m-%d")
                else:
                    timediff2 = 0
                    timediff = float(self._data[i][self._variables[visitdatevarname]])
                    datecolignore = self._variables[visitdatevarname]
                    # carry out conversion between string
                    # to datetime object

                    beginDate = datetime.datetime.strptime(
                        visitbaselinedate, "%d/%m/%Y"
                    )

                    if visitdateformat == 2:  # use visitbaselinedate and add Days
                        startdate = beginDate + datetime.timedelta(days=int(timediff))
                    elif visitdateformat == 3:  # use visitbaselinedate and add Months
                        startdate = beginDate + datetime.timedelta(
                            days=int(timediff * 12)
                        )
                    elif visitdateformat == 4:  # use visitbaselinedate and add Years
                        startdate = beginDate + datetime.timedelta(
                            days=int(timediff * 365)
                        )

                    dt_string = startdate.strftime("%Y-%m-%d")
                # Loop through cells in row
                for j, value in enumerate(self._data[i]):
                    if (
                        j == self._variables[self.config["patientid"]]
                        or j == datecolignore
                    ):
                        continue
                    else:
                        if self._data[i][j].strip() == "":
                            continue

                        code = "-"
                        for x in range(len(self._map_phenotype_to_concept)):
                            if x > 0:
                                # Check if it is an enumerated value, then only add code
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
                        facts.append((mrn, str(dt_string), code, value))
                        # print (mrn + ", " + str(dt_string) + ", " + str(code) + ", " + str (value) )

        # Write facts file

        # factsfile = 'facts_' + self.config['filebase'] + '.csv'
        with open(factsfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mrn", "start-date", "code", "value"])
            for row in facts:
                writer.writerow(row)

    def write_facts(self, factsfile):
        if self.config["filebase"] == "adverse":
            self.write_adverse_facts(factsfile)
        elif self.config["filebase"] == "fundus":
            self.write_fundus_facts(factsfile)
        elif self.config["filebase"] == "mortality":
            self.write_mortality_facts(factsfile)
        elif self.config["filebase"] == "accord_key":
            self.write_accord_key_facts(factsfile)


#
# Command-line arguments. Could add dictionary and input to config, but
# there are instances where the same config will work with different inputs
#
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="ETL config file", default="etl.yml")
    parser.add_argument("-d", "--dictionary", help="file containing data dictionary")
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
