# etldev

`etl.py` transforms dbGaP data to the i2b2 data format.   
Run with 3 arguments: Configuration file, input file, and data dictionary file

# Basic YML Configuration file:

* separator &nbsp;&emsp;# Separator token. Provide it in single quotes, e.g. ‘;’  
* enumname &emsp;# Column name in data dictionary containing enumerated values  
* typename &nbsp;&emsp;# Column name in data dictionary defining  the type of variable such as string, int, enumerated  
* varname &emsp;# Column name in data dictionary with variable names for the ontology   
* patientid &emsp;# Patient id name in data dictionary file  
* basedate &emsp;# The baseline or enrollment date for a clinical trial participation  
* dateformat &emsp;# Options: 1=date, 2=days from baseline, 3= months from baseline, 4= years from baseline  
* datevar &emsp;# If visitdateformat not 1 : variable name containing the difference in date for the visit  
* pathroot &emsp;# Name of the ontology root for concepts.csv  
* filebase &emsp;# Name for output files usually correlating with study name  
  
\# Additional Dates  
* additionaldatediffvarname &emsp;# NA or Name of variable when a second date variable contains difference to the first. E.g. Adverse event date  
* additionaldatedifftimeunits &emsp;# Units for the difference in time of the additional time variable (example diff between first occurrence  and adverse event resolution)  
                            &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;# If NA, write 0, 2=days from baseline, 3= months from baseline, 4= years from baseline  

