# Transform dbGaP files to i2b2 files

`etl.py` transforms dbGaP data dictionary and data files to the i2b2 file format of fact and concept files for i2b2-etl-docker.   
Run with 3 arguments: Configuration file ( "-c", "--config"), input file ("-i", "--input"), and data dictionary ("-d", "--dictionary") file    
 
Example:   
```python etl.py -c AREDS_followup.yml -i areds_followup.txt -d dbGaP_Data_Dictionary_AREDS_followup_x.csv```

# Basic YML Configuration file:

## dictformat
dbGap dictionary format. Data dictionaries within the same study have similar format. Currently 2 studies are supported as examples:   
 dictformat: 'areds2'   
 dictformat: 'areds'   
Data files have to be obtained from dbGaP via a DUA.  

 ## separator
Separator token. Provide it in single quotes, e.g. ‘;’  

## enumname
Column name in data dictionary containing enumerated values 

## typename 
Column name in data dictionary that defines the type of variable such as string, int, enumerated  

## varname 
Column name in data dictionary with variable names for the ontology   

## patientid 
Column name in data dictionary with the Patient id 

## basedate
The baseline date or enrollment date for a clinical trial participation  

## dateformat 
Used for calculating time differences between dates, e.g. between recrutiment/baseline date and visit date.    
If dateformat not 0 or 1, select variable name containing the date difference from baseline date in days, months, years.    
Options:   
0 = no date - when using baseline date or datemode-Mode 6 
1 = date  
2 = days from baseline  
3 = months from baseline  
4 = years from baseline   
5 = 0.5 years from baseline   
dateformat is ignored for datemode-Mode 6 using 
 ```visitdatefile``` coupled with ["timevar"]["default"]. Example: Areds_followup.yml   

## datemode
Mode 0: all times are relative to "basedate". Example: AREDS2_dem.yml   
Mode 1: DEPRECATED - replaced by Mode 5 for float and Mode 7 for integer
Mode 2: facts are tied to certain time variables via regex, e.g. left and right eye exam dates. Example: AREDS_fundus.yml  
Mode 3: DEPRECATED - a list of possible additional time points, represented as deltas   
Mode 4: DEPRECATED - a list of possible additional time points, relative to "basedate"   
Mode 5: visit is calculated by a time difference as a decimal (parsable within a string). Example: AREDS_adverse.yml, AREDS_mortality.yml  
Mode 6: Use Visno file coupled with ["timevar"]["default"]. Does not use dateformat. Example AREDS2_rcf.yml   
Mode 7: visit is calculated by a time difference as an integer (parsable within a string), e.g. F04, for 4 months visit. Example ACCORD_f34.yml   

## pathroot 
Name of the ontology root for concepts.csv  

## filebase 
Name for output files usually correlating with study name  

## visitdatefile [Optional]
CSV file mapping visit numbers to visit dates 
Example: 'AREDS2_VISNO_Dates.csv'

## timevar [Optional]
When different dates are attached to diferent variables within the same file
### default
This variable may be linked to   
```datemode: 6 ```  
with a  
```visitdatefile```  
For an example, see AREDS2_rcf.yml

## demographics_file [Optional]
CSV file mapping dbGaP demographics as defined in data dicitonary into i2b2 demographic codes.
The i2b2code should be of length 50 and contain only characters and numbers. It is recommended to run etl.py without demographics_file in the yml file to consult and use the codes in concepts.csv before running it with this configuration.
Example: 'Dem-AREDS2.csv'

## codeprefix [Optional] 
A prefix to use for the i2b2 codes to separate the facts from this data file from other data files using the same concepts.
Used most often when testing Bring-Your-Own-Data by uploading the same file multiple times.
 
  
\# Additional Dates  
* additionaldatediffvarname &emsp;# NA or Name of variable when a second date variable contains difference to the first. E.g. Adverse event date  
* additionaldatedifftimeunits &emsp;# Units for the difference in time of the additional time variable (example diff between first occurrence  and adverse event resolution)  
                            &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;# If NA, write 0, 2=days from baseline, 3= months from baseline, 4= years from baseline  

