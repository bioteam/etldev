#!/usr/bin/env python
import argparse
import yaml
import datetime
from datetime import timedelta
import csv

def parse_args(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="ETL config file", default = "etl.yml")
    parser.add_argument("-d", "--dictionary", help="file containing data dictionary")
    parser.add_argument("-i", "--input", help="file data to be ETL'd")
    args = parser.parse_args()
    return args

def load_conf(config_file):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
        etl_conf = config[0]["etl"]
    return etl_conf

def parse_dictionary(dictfile):
    with open(dictfile, encoding='latin-1') as csvfile:
        ptr = csvfile.tell()
        line = csvfile.readline()
        # Skip leading comment lines
        while line.startswith("!#!"):
            ptr = csvfile.tell()
            line = csvfile.readline()
        csvfile.seek(ptr)
        reader = csv.DictReader(csvfile)
        elem = []
        for row in reader:
            # Skip empty lines
            if list(row.values())[0]:
                elem.append(row)
    return elem
  
def etl_concepts(etl_conf,dd):
    split_data = []
    map_phenotype_to_concept = []
    for row in dd:
        varname = row[etl_conf['varname']]
        if varname == etl_conf['patientid']:
            continue
        vartype = row[etl_conf['typename']]
        values = row[etl_conf['enumname']].split(etl_conf['separator'])
        path = "/".join(['',etl_conf['pathroot'],varname,''])
        # Check what i2b2 type of variable it is: integer, float, string or large-string
        if len(values) == 1:
            if vartype.lower() == "num" or vartype.lower() == "integer" :
                i2b2vartype= "integer"
            elif vartype.lower() == "decimal" or vartype.lower() == "float" :
                i2b2vartype= "float"
            else:
                i2b2vartype= "string"

            i2b2code = varname
            dbgap_code_id = -1
            conceptpath = path
            split_data.append((conceptpath, i2b2code,i2b2vartype))
            map_phenotype_to_concept.append((conceptpath, i2b2code,i2b2vartype, dbgap_code_id, varname))
        for i, value in enumerate(values):
            value.replace('"','')      
            # Limit the split to 1 as some descriptions contain "="
            clin_name= value.strip().split('=',1)
            dbgap_code_id = ''
            i2b2conceptlabel=''
            i2b2concept=''
            if len(clin_name) > 1:
                dbgap_code_id = clin_name[0].replace('"','').lstrip()
                i2b2concept = clin_name[1].replace('"','') 
                i2b2concept = i2b2concept.replace(',',' or ') 
                i2b2conceptlabel= ''.join(filter(str.isalnum, i2b2concept))
            if (i2b2conceptlabel == dbgap_code_id): # TODO: Ensure i2b2code is unique in the ontology
                i2b2code = varname + ''.join(filter(str.isalnum, i2b2concept))   
            else:   
                i2b2code= varname + ''.join(filter(str.isalnum, i2b2concept)) + dbgap_code_id        
            if len(i2b2code)>50: # Ensure i2b2code is no longer than 50 characters, max : 50, can be any length
                leftposition=len(i2b2code)-50
                i2b2code=i2b2code[leftposition:]
            conceptpath = path + i2b2concept
            split_data.append((conceptpath, i2b2code,"assertion"))
            map_phenotype_to_concept.append((conceptpath, i2b2code,"assertion", dbgap_code_id, varname))
    conceptsfile = 'concepts_' + etl_conf['filebase'] + '.csv'
    with open(conceptsfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['path','code','type'])
        for row in split_data:
            writer.writerow(row)

    return map_phenotype_to_concept

def etl_facts(etl_conf, map_phenotype_to_concept, phenocsvfile):
    with open(phenocsvfile, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
    variables = {}
    i = 0
    while i < len(data[0]):
        variables[data[0][i]] = i
        i += 1
    
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
    visitdatevarname = etl_conf['datevar']
    visitdateformat = int(etl_conf['dateformat'])
    visitbaselinedate = etl_conf['basedate']
    additionaldatediffvarname = etl_conf['additionaldatediffvarname']
    additionaldatedifftimeunits = int(etl_conf['additionaldatedifftimeunits'])
    for i in range(len(data)):
        if i > 0:
            # Patient ID
            mrn =  data[i][variables[etl_conf['patientid']]]
            #print ("patient"+str(mrn))
            # Visit Date
            if visitdateformat == 1: # Use date from column visitdatevarname
                startdate = data[i][variables[visitdatevarname]]
                datecolignore = variables[visitdatevarname]
            else:
                timediff2=0
                timediff=float(data[i][variables[visitdatevarname]])
                datecolignore = variables[visitdatevarname]
                # carry out conversion between string 
                # to datetime object
                if additionaldatediffvarname != "NA":
                    if data[i][variables[additionaldatediffvarname]].isnumeric():
                        timediff2 = float(data[i][variables[additionaldatediffvarname]])

                beginDate = datetime.datetime.strptime(visitbaselinedate, "%d/%m/%Y")

                if visitdateformat == 2: # use visitbaselinedate and add Days
                    startdate = beginDate  + datetime.timedelta(days=int(timediff)) 
                elif visitdateformat == 3: # use visitbaselinedate and add Months  
                    startdate = beginDate  + datetime.timedelta(days=int(timediff * 12)) 
                elif visitdateformat == 4: # use visitbaselinedate and add Years
                    startdate = beginDate  + datetime.timedelta(days=int(timediff * 365) )

                if additionaldatedifftimeunits > 1 and additionaldatedifftimeunits < 5:
                    if additionaldatedifftimeunits == 2: # use visitbaselinedate and add Days
                        startdate = startdate  + datetime.timedelta(days=int(timediff2)) 
                    elif additionaldatedifftimeunits == 3: # use visitbaselinedate and add Months  
                        startdate = startdate  + datetime.timedelta(days=int(timediff2 * 12)) 
                    elif additionaldatedifftimeunits == 4: # use visitbaselinedate and add Years
                        startdate = startdate  + datetime.timedelta(days=int(timediff2 * 365) )
                
                dt_string = startdate.strftime("%Y-%m-%d")
            # Loop through cells in row
            for j, value in enumerate(data[i]):
                if j == variables[etl_conf['patientid']] or j == datecolignore:
                    continue
                else:
                    if data [i][j].strip() == "" :
                        continue
                
                    code = "-"
                    for x in range(len(map_phenotype_to_concept)):
                        if x > 0 :
                        # Check if it is an enumerated value, then only add code
                            if map_phenotype_to_concept[x][3] == value and map_phenotype_to_concept[x][4] == data[0][j] :
                                code = map_phenotype_to_concept[x][1]
                                value = ""

                    if code =="-" :
                        code = data[0][j]
                        value = data [i][j]
                    facts.append((mrn,str(dt_string), code, value))
                    #print (mrn + ", " + str(dt_string) + ", " + str(code) + ", " + str (value) )
     
# Write facts file

    factsfile = 'facts_' + etl_conf['filebase'] + '.csv'
    with open(factsfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['mrn', 'start-date', 'code', 'value'])
        for row in facts:
            writer.writerow(row)


def main(): 
    inputs = parse_args()
    etl_conf = load_conf(inputs.config)
    dd = parse_dictionary(inputs.dictionary)
    map2concepts = etl_concepts(etl_conf,dd)
    etl_facts(etl_conf,map2concepts,inputs.input)

if __name__ == '__main__':
    main()

