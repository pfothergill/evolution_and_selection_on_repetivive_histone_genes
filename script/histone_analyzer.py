#!/usr/bin/python3

import sys
from math import log


# This is useful as a command line arguement that allows us to better understand what call we are on
# For example, when Yellow2_Histone_H2A.raw is ran through this python script, it will set
# filename = "../data/deviate_histones_raw/Yellow2_Histone_H2A.raw" and from there:
# sample_id = Yellow2
# placeholder = Histone (IT WILL ALWAYS EQUAL "Histone" AS HISTONE IS THE MIDDLE ITEM IN "Yellow2_Histone_H2A.raw"
# when the above is split at each "_" symbol )
# histone = H2A (the histone.split(".") splits the original histone (equal to H2A.raw) into 'H2A' and 'raw' and then
# sets histone equal to the zeroth item, in this case, H2A)
filename = "../data/deviate_histones_raw/" + sys.argv[1]
sample_id, placeholder, histone = sys.argv[1].split('_')
histone = histone.split('.')[0]

# Dictionary for codons. Each dictory key is a sequence and the value for the key is
# in [possible synonymous changes, "one letter code"] format
CodonDict = {
    "ATG":[0, "M"], "TGG":[0, "W"], "ACT":[1, "T"], "ACC":[1, "T"], "ACA":[1, "T"], "ACG":[1, "T"], "TCT":[1, "S"], 
    "TCC":[1, "S"], "TCA":[1, "S"], "TCG":[1, "S"], "CTT":[1, "L"], "CTA":[1, "L"], "CTC":[1, "L"], "CTG":[1, "L"],
    "CCT":[1, "P"], "CCA":[1, "P"], "CCC":[1, "P"], "CCG":[1, "P"], "CGT":[1, "R"], "CGC":[1, "R"], "GTT":[1, "V"],
    "GTC":[1, "V"], "GTA":[1, "V"], "GTC":[1, "V"], "GTG":[1, "V"], "GCT":[1, "A"], "GCA":[1, "A"], "GCC":[1, "A"], 
    "GCG":[1, "A"], "GGT":[1, "G"], "GGC":[1, "G"], "GGA":[1, "G"], "GGG":[1, "G"], 
    "AAT":[1/3, "N"], "AAC":[1/3, "N"], "AAA":[1/3, "K"], "AAG":[1/3, "K"], "AGT":[1/3, "S"], "AGC":[1/3, "S"], "TTT":[1/3, "F"],
    "TTC":[1/3, "F"], "TAT":[1/3, "Y"], "TAC":[1/3, "Y"], "TGT":[1/3, "C"], "TGC":[1/3, "C"], "TGA":[1/3, "X"], "CAT":[1/3, "H"], 
    "CAC":[1/3, "H"], "CAA":[1/3, "Q"], "CAG":[1/3, "Q"], "GAT":[1/3, "D"], "GAA":[1/3, "E"], "GAC":[1/3, "D"], "GAG":[1/3, "E"], 
    "ATT":[2/3, "I"], "ATC":[2/3, "I"], "ATA":[2/3, "I"], "AGA":[2/3, "R"], "AGG":[2/3, "R"], "TTA":[2/3, "L"], "TTG":[2/3, "L"], 
    "TAA":[2/3, "X"], "TAG":[2/3, "X"], "CGA":[4/3, "R"], "CGG":[4/3, "R"]}

# Global values used for sum_for_avg_s and sum_for_avg_n. These will be reset to 0 during each call to this python script and only change
# during a syn (sum_for_avg_s) or non-syn (sum_for_avg_n) nucleotide change
#for non-synonymous
sum_for_avg_n = 0
sum_K = 0

# for synonymous
sum_for_avg_s = 0
sum_L = 0

# Function to computationally find sum_for_avg_s or sum_for_avg_n for each given codon and adds those respective values to their
# global sum_for_avg_s or sum_for_avg_n value. This code only runs when the below code finds a base that isn't refbase but is greater 
# than the 1x value divided by 2 (normal rounding where x.5 rounds up). This 1/2 1x parameter to find changes among sequences was 
# chosen to not over look possible changes that were 1 or 2 less then the 1x value. For debugging purposes, this part of
# the code also writes to a debugger_file (debugger_file.txt) located in the output directory of this git repo
def run_codon(codon_list, orig_codon_sequence, tmp_info, number_of_changes):
    
    global sum_for_avg_n
    global sum_L

    global sum_for_avg_s
    global sum_K
    
    with open("../output/debugger_file.txt", "a") as debugger_txt:
        for x in range(3):
            debugger_txt.write(str(codon_list[x])+'\n')
        debugger_txt.write("tmp_info: {}\n".format(tmp_info))
        debugger_txt.write("Number of Changes: {}\n".format(number_of_changes))
        debugger_txt.write("Original sequence: {}\n".format(orig_codon_sequence))
        debugger_txt.write("Original aa: {}\n\n".format(CodonDict[orig_codon_sequence][1]))
        changed_codons = []
        for z in range(len(tmp_info)):
            for y in range(len(tmp_info[z][0])):
                tmp_codon_sequence = ""
                change = tmp_info[z][0][y]
                for x in range(3):
                    if tmp_info[z][1] == x:
                        tmp_codon_sequence += change
                    else:
                        tmp_codon_sequence += orig_codon_sequence[x]
                changed_codons.append(tmp_codon_sequence)
        for x in range(len(changed_codons)):
            debugger_txt.write("Changed sequence: {}\n".format(changed_codons[x]))
            debugger_txt.write("Changed aa: {}\n".format(CodonDict[changed_codons[x]][1]))
            if CodonDict[orig_codon_sequence][1] == CodonDict[changed_codons[x]][1]:
                debugger_txt.write("Type of change: Synonymous\n")
                sum_for_avg_s += number_of_changes[x]
                sum_L += CodonDict[orig_codon_sequence][0]
            else:
                debugger_txt.write("Type of change: Non-synonymous\n")
                sum_for_avg_n += number_of_changes[x]
                sum_K += (3-CodonDict[orig_codon_sequence][0])
            debugger_txt.write('\n')

# open the Lineage1xAndCN.txt file as to get information used in later calculations
with open("../data/Lineage1xAndCN.txt", "r") as info_file:
    information = []
    for line in info_file:
        line = line.split()
        if line[0] == sample_id:
            ploidy = int(line[1][0])
            one_x = float(line[2])
            if histone == 'H2A':
                copy_number = int(line[3])
            elif histone == 'H2B':
                copy_number = int(line[4])
            elif histone == 'H3':
                copy_number = int(line[5])
            elif histone == 'H4':
                copy_number = int(line[6])
            break

# iterate through the inputted filename, 3 lines at a time and at the end, call run_codon() with 
# all the relevant information to find sum_for_avg_s and sum_for_avg_n for each codon that changes
with open(filename, "r") as raw_file:
    codon_list = []
    position = 0
    codon_sequence = ""
    need_to_run_test = False
    tmp_sequence = []
    tmp_base = ""
    number_of_changes = []
    nucleotide_counter=0
    for line in raw_file:
        line = line.split()
        if line[0][0] != "#":
            nucleotide_counter += 1
            TEfam = line[1]
            raw_file_sample_id = line[2]
            refbase = line[3]
            A = int(line[4])
            C = int(line[5])
            G = int(line[6])
            T = int(line[7])
            if refbase == 'A':
                if (C > round(one_x/2)):
                    tmp_base += 'C'
                    number_of_changes.append(C/one_x)
                if (G > round(one_x/2)):
                    tmp_base += 'G'
                    number_of_changes.append(G/one_x)
                if (T > round(one_x/2)):
                    tmp_base += 'T'
                    number_of_changes.append(T/one_x)
                if tmp_base != "":
                    need_to_run_test = True
                    info = tmp_base, position
                    tmp_sequence.append(info)
                    tmp_base = ""
            elif refbase == 'C':
                if (A > round(one_x/2)):
                    tmp_base += 'A'
                    number_of_changes.append(A/one_x)
                if (G > round(one_x/2)):
                    tmp_base += 'G'
                    number_of_changes.append(G/one_x)
                if (T > round(one_x/2)):
                    tmp_base += 'T'
                    number_of_changes.append(T/one_x)
                if tmp_base != "":
                    need_to_run_test = True
                    info = tmp_base, position
                    tmp_sequence.append(info)
                    tmp_base = ""
            elif refbase == 'G':
                if (A > round(one_x/2)):
                    tmp_base += 'A'
                    number_of_changes.append(A/one_x)
                if (C > round(one_x/2)):
                    tmp_base += 'C'
                    number_of_changes.append(C/one_x)
                if (T > round(one_x/2)):
                    tmp_base += 'T'
                    number_of_changes.append(T/one_x)
                if tmp_base != "":
                    need_to_run_test = True
                    info = tmp_base, position
                    tmp_sequence.append(info)
                    tmp_base = ""
            elif refbase == 'T':
                if (A > round(one_x/2)):
                    tmp_base += 'A'
                    number_of_changes.append(A/one_x)
                if (C > round(one_x/2)):
                    tmp_base += 'C'
                    number_of_changes.append(C/one_x)
                if (G > round(one_x/2)):
                    tmp_base += 'G'
                    number_of_changes.append(G/one_x)
                if tmp_base != "":
                    need_to_run_test = True
                    info = tmp_base, position
                    tmp_sequence.append(info)
                    tmp_base = ""
            codon_list.append(line)
            codon_sequence += line[3]
            if (position == 2):
                if (need_to_run_test == True):
                    run_codon(codon_list, codon_sequence, tmp_sequence, number_of_changes)
                codon_list = []
                number_of_changes = []
                position = 0
                codon_sequence = ""
                need_to_run_test = False
                tmp_sequence = []
                tmp_base = ""
            else:
                position += 1

N = (sum_for_avg_n/nucleotide_counter) / (sum_L * copy_number)
S = (sum_for_avg_s/nucleotide_counter) / (sum_K * copy_number)

# write the results to the temp_results.txt file in output directory. This is the file that is written to
# because later we will take that file and do some fancy things to align the columns properly and the
# output will be named something else and this file will be deleted
with open("../output/temp_NS_results.txt", "a") as NS_results_file:  
    NS_results_file.write("{}   {}   {}   {}   {}   {}\n".format(sample_id, ploidy, histone, N, S, N/S))

