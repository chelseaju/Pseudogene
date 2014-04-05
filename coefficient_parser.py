""" locus_parser.py
Usage: python coefficient_parser.py -d directory -t type 
Input:  -d input/output directory (select_pseudogene_128_v2)
		-t type = gene or isoform
Output: line of coefficient with genomic locus {Gene_id \t Chr \t Start \t End \t Count\n}
Function: Parse out the genomic location for each coefficient
Date: 2014-04-01
Author: Chelsea Ju
"""

import sys, re, os, subprocess, random, argparse


## local computer
ENSG = "/home/chelseaju/Database/Pseudogene/ENSG_Pseudogene_74.bed"
ENST = "/home/chelseaju/Database/Pseudogene/ENST_Pseudogene_74.bed"

LOCUS={}


"""
	Function : read in the genomic location for each gene
"""
def build_locus(bedfile):

	fh = open(bedfile, 'rb')

	for line in fh:
		data = line.split("\t")
		chromosome = data[0]
		start = data[1]
		end = data[2]
		name = data[3]

		LOCUS[name] = (chromosome, start, end)

	fh.close()

"""
	Function : combine the locus with coefficients
"""
def combine_locus_coefficient(coefficientfile):

	combine_data = []

	fh = open(coefficientfile, 'rb')
	for line in fh:
		line = line.rstrip()
		data = line.split("\t")
		name = data[0]
		coef = data[1]

		unknown = re.match("Unknown_", name)
		if(unknown):
			(prefix, chromosome, start, end) = name.split("_")
			combine_data.append((name, chromosome, start, end, coef))

		else:
			(chromosome, start, end) = LOCUS[name]
			combine_data.append((name, chromosome, start, end, coef))
	fh.close()

	return combine_data

"""
	Function : export data
"""
def export_data(data, outfile):
	
	fh = open(outfile, 'w')
	for line in data:
		fh.write("%s\t%s\t%s\t%s\t%s\n" %(line[1], line[2], line[3], line[0], line[4]))
	fh.close()

def main(parser):

    options = parser.parse_args()
    dir = options.dir
    type = options.type
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    if(type == "ENSG"):
    	bed_file = ENSG
    else:
    	bed_file = ENST
    coefficient_file = dir + type + "_glmnet_coefficient.xls"
    outfile = dir + type + "_locus_coefficient.xls"

    build_locus(bed_file)
    locus_coefficient = combine_locus_coefficient(coefficient_file)

    export_data(locus_coefficient, outfile)

    print ""
    print "Writing Coefficient to File : %s" %(outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='coefficient_parser.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-t", "--type", dest="type", type=str, help="gene or ENSG", required = True)
    main(parser)
