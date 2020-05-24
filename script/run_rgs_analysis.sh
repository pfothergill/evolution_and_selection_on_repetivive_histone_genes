#!/bin/bash

# remove old debugger file so a new one can be created (might want to create copies after each program run to keep this data)
if [ -e ../output/debugger.txt ]; then
    rm -f ../output/debugger.txt
fi

# setup the fields for the temp_results.txt file. Using this temp file is a way to obtain the results
# and then manipulate the data for a "pretty print", saving that "pretty print" to the final file
# and then removing this temp file!
echo "Lineage   Ploidy  Histone  N  S  N/S" > ../output/temp_NS_results.txt
echo "-------   ------  -------  ---------------------- ---------------------- ---------------------" >> ../output/temp_NS_results.txt

# iterate through each file in the data/deviate_histones_raw directory, setting each filename equal to $f and
# running the current $f through the analyze_raw_files.py python script
for f in ../data/deviate_histones_raw/*; do python3 histone_analyzer.py $(basename $f); done

# arrange the files in an easier to read format and remove the temp file (temp file has the exact same
# data as the file results.txt, it's just not as organized/aligned as the new results.txt)
column -t ../output/temp_NS_results.txt > ../output/NS_results.txt && rm -f ../output/temp_NS_results.txt

# DONE :)
