# Script Output 
NS_results.txt contains the python script output. This output is in a tabular format that was later used for statistical analysis using R and RStudio.
## Setting up Data for R Analysis
Editing of text files must be done before the NS output can be used in R
```
$ vim NS_results.txt
# first go to the dashed bar under the headers and type dd
# then go to N/S and change to N_S
:%s/\s\+/\t/g #to create tab format for R
:wq
```
