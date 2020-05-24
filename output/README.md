# Script Output 
"NS_results.txt" contains the python script output. 
## Setting up Data for R Analysis
Editing of text files must be done before the NS output can be used in R
```
$ vim NS_results.txt
# first go to the dashed bar under the headers and type dd
# then go to N/S and change to N_S
:%s/\s\+/\t/g #to create tab format for R
:wq
```
# Rstudio Output 
Output should include 4 Figures: "Figure1.png", "Figure2.png", "Figure3.png", "Figure4.png".

Statical test results were compiled into 2 tables and added to output: "Table.1.txt", "Table.2.txt"
