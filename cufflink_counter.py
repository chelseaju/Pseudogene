""" cufflink_counter.py
Usage: python cufflink_counter.py -d cufflink_directory 
Input:  -d input/output directory (select_pseudogene_128_v2/10X_100L_1A_30X_100L_1A)
Output: line of cufflink output combining count and fpkm {XLOC_id \t Gene_id \t Overlap \t Count_A \t Count_B \t FPKM_A \t FPKM_B}
Function: 1. combine the information from genes.fpkm_tracking and genes.count_tracking with gene ids

Date: 2014-03-30
Author: Chelsea Ju
"""

import sys, re, os, subprocess, random, argparse

GENES = {}

"""
	Function : read in the file created from cufflink_identifier.py, build the gene table
"""
def build_gene_table(gene_table):

	fh = open(gene_table, "rb")

	for line in fh:
		(xloc, gene, coverage, alternatives) = line.split("\t")
		GENES[xloc] = (gene, coverage)

	fh.close()

"""
	Function : combine the fpkm and gene count results
"""
def combine_fpkm_count(gene_fpkm, gene_count):

	result = {}

	fpkm = open(gene_fpkm, 'rb')
	count = open(gene_count, 'rb')

	line = fpkm.readline()
	for line in fpkm:
		data = line.split("\t")
		xloc = data[0]
		fpkm_a = data[9]
		fpkm_b = data[13]

		if(GENES.has_key(xloc)):
			(gene_name, gene_overlap) = GENES[xloc]
		else:
			gene_name = xloc
			gene_overlap = "1"

		result[xloc] = [gene_name, gene_overlap, fpkm_a, fpkm_b]

	fpkm.close()

	line = count.readline()
	for line in count:
		data = line.split("\t")
		xloc = data[0]
		count_a = data[1]
		count_b = data[6]

		result[xloc].append(count_a)
		result[xloc].append(count_b)

	count.close()

	return result

"""
	Function : export result to file
"""
def export_result(result, outfile):
	fh = open(outfile, 'w')

	fh.write("XLOC\tGene_Name\tOverlap\tFPKM_A\tFPKM_B\tCount_A\tCount_B\n")
	for loc in result.keys():
		r = result[loc]

		fh.write(loc)
		fh.write("\t")
		fh.write("\t".join(r))
		fh.write("\n")

	fh.close()


def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

	gene_fpkm = dir + "genes.fpkm_tracking"
	gene_count = dir + "genes.count_tracking"
	gene_table = dir + "cufflink_genes.txt"

	outfile = dir + "cufflink_results.txt"

	build_gene_table(gene_table)
	result = combine_fpkm_count(gene_fpkm, gene_count)

	export_result(result, outfile)


    print ""
    print "Writing Gene List to File : %s" %(outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='cufflink_counter.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="cufflink directory of input and output files", required = True)

    main(parser)
