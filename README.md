# FRAPAnalysis
MATLAB code for analysis of FRAP data

## Additional Source Code
MATLAB code in bfmatlab and FRAP_analysis-master are both required for FRAP_Analysis_v1.m. Make sure to update the path locations in lines 10 and 11 of FRAP_Analysis_v1.m to match the local directories for your workstation.

## Data Collection Settings
Prior to using the code, make sure to update lines 50 and 51 cooresponding to the bleachFrame and frameLength to match your experimental conditions.

## Running FRAP_Analysis_v1.m
### Compute Fractional Recovery
To compute fractional recovery curves for every movie in a directory, supply the directory location on line 13 for the "outdir" variable. Multiple directories can be specified by listing directory in an array (ex: {folder1, folder2, ...}). Additionally confirm that "ToDo" on line 39 is set to "0".

### Perform Curve Fitting
Set the "ToDo" variable to "1" on line 39 and ensure that the file selection scheme for the "myExt: variable in line 45 matches the files you would like to pool for analysis. Additionally, update the "compiled_Csv" variable on line 48 where the outputs from two-exponential fits will be recorded.

## Additional Support
Direct any issues or comments to jferrie@berkeley.edu
