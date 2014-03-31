""" cufflink_evaluator.py
Usage: python cufflink_evaluator.py -d cufflink_directory -a expected_count_A -b expected_count_B 
Input:  -d input/output directory (select_pseudogene_128_v2/10X_100L_1A_30X_100L_1A)
		-a expected number of fragment in sample A
		-b expected number of fragment in sample B
Output: line of cufflink output combining count and fpkm {XLOC_id \t Gene_id \t Overlap \t Count_A \t Count_B \t FPKM_A \t FPKM_B \t Expected_A \t Expected_B}
Function: 1. combine the information from cufflink result with expected data

Date: 2014-03-30
Author: Chelsea Ju
"""

GENES = {}
XLOC = {}


import sys, re, os, subprocess, random, argparse

"""
	Function : build gene gene id
"""
def build_name(file):
	fh = open(file, 'rb')
	for line in fh:
		line = line.rstrip()
		data = line.split("\t")
		xloc = data[0]
		genes = data[3].split(",")

		for g in genes:
			GENES[g] = (xloc)
	fh.close()

"""
	Function : build observation
"""
def build_observation(file):

	fh = open(file, 'rb')
	for line in fh:
		line = line.rstrip()
		data = line.split("\t")
		xname = data[0]
		data.append("0")
		data.append("0")
		XLOC{xname} = data
	fh.close()

"""
	Function : read in expected count for each gene
"""
def build_expectation(fileA, fileB):

	fh = open(fileA, 'rb')
	for line in fh:
		(name, count) = line




def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    file_A = options.fileA
    file_B = options.fileB
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    name_file = dir + "cufflink_genes.txt"
    result = dir + "cufflink_results.txt"
    outfile = dir + "cufflink_evaluation.txt"

    build_name(name_file)
    build_observation(result)
    build_expectation(file_A, file_B)


#    evaluation = combine_expectation_observation(a_exp, b_exp, result)

#    export_data(evaluation, outfile)

    print ""
    print "Writing Gene List to File : %s" %(outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='cufflink_evaluator.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="cufflink directory of input and output files", required = True)
    parser.add_argument("-a", "--sampleA", dest="fileA", type=str, help="expected number of count from sample A", required = True)
    parser.add_argument("-b", "--sampleB", dest="fileB", type=str, help="expected number of count from sample B", required = True)

    main(parser)
