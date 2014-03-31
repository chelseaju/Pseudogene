""" cufflink_evaluator_v2.py
Usage: python cufflink_evaluator_v2.py -d cufflink_directory -t type -a expected_count_A -b expected_count_B 
Input:  -d input/output directory (select_pseudogene_128_v2/10X_100L_1A_30X_100L_1A)
		-a expected number of fragment in sample A
		-b expected number of fragment in sample B
		-t type = gene or isoform
Output: line of cufflink output combining count and fpkm {XLOC_id \t Gene_id \t Count_A \t Count_B \t FPKM_A \t FPKM_B \t Expected_A \t Expected_B}
Function: 1. combine the information from cufflink result with expected data

Date: 2014-03-31
Author: Chelsea Ju
"""

XLOC = {}
EXPECTATION = {}


import sys, re, os, subprocess, random, argparse

"""
	Function : read in expected count for each gene
"""
def build_expectation(file, order):

	fh = open(file, 'rb')
	for line in fh:
		line = line.rstrip()
		(name, count) = line.split("\t")

		if(EXPECTATION.has_key(name)):
			EXPECTATION[name][order] = count
		else:
			count_data = ["0","0"]
			count_data[order] = count
			EXPECTATION[name] = count_data

	fh.close()

"""
	Function : combine the observation and expected value
		if the id is not in the list, expected values are zeros
		if the expected ids are not identified by cufflink, they are added to the bottom of the list
"""
def combine_expectation_observation(result):

	evaluation_data = []
	seen_ids = {}

	# mark all expected ids
	for k in EXPECTATION.keys():
		seen_ids[k] = False

	fh = open(result, 'rb')
	for line in fh:
		data = line.split("\t")
		t_id = data[0]
		x_id = data[1]
		gene_name = data[2]
		locus = data[3]
		a_fpkm = data[7]
		b_fpkm = data[8]
		log = data[9]
		a_count = data[15]
		b_count = data[20]

		expect = ["0","0"]

		if(EXPECTATION.has_key(gene_name)):
			expect = EXPECTATION[gene_name]
			seen_ids[gene_name] = True

		output_data = [t_id, x_id, gene_name, locus, a_fpkm, b_fpkm, log, a_count, b_count] + expect
		evaluation_data.append(output_data)
	fh.close()

	for k in sorted(seen_ids, key=seen_ids.get):

		if(seen_ids[k] == False):
			data = ["-",k,"-","-","-","-","-","-","-"] + EXPECTATION[k]
			evaluation_data.append(data)
		else:
			break


	return evaluation_data

"""
	Function : output data
"""
def export_data(result, file):
	fh = open(file, 'w')
	fh.write("T_LOC\tXLOC\tGENENAME\tLOCUS\tFPKM_A\tFPKM_B\tLOG\tCOUNT_A\tCOUNT_B\tEXPECTED_A\tEXPECTED_B\n")
	for r in result:
		fh.write("\t".join(r))
		fh.write("\n")

	fh.close()

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    file_A = options.fileA
    file_B = options.fileB
    block = options.block
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    fpkm = dir + block + "_exp.diff"
    count = dir + block + "s.count_tracking"
    fpkm_tmp = dir + "tmp1.txt"
    count_tmp = dir + "tmp2.txt"
    result = dir + block + "_results.txt"
    outfile = dir + block + "_evaluation.txt"

    ## merge the fpkm and count data first
    os.system("sort %s -t '\t' > %s" %(fpkm, fpkm_tmp))
    os.system("sort %s -t '\t' > %s" %(count, count_tmp))
    os.system("join -t '\t' %s %s > %s" %(fpkm_tmp, count_tmp, result))
    os.system('rm %s %s' %(fpkm_tmp, count_tmp))

    build_expectation(file_A, 0)
    build_expectation(file_B, 1)

    evaluation = combine_expectation_observation(result)
    export_data(evaluation, outfile)

    print ""
    print "Writing Gene List to File : %s" %(outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='cufflink_evaluator.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="cufflink directory of input and output files", required = True)
    parser.add_argument("-a", "--sampleA", dest="fileA", type=str, help="expected number of count from sample A", required = True)
    parser.add_argument("-b", "--sampleB", dest="fileB", type=str, help="expected number of count from sample B", required = True)
    parser.add_argument("-t", "--blocktype", dest="block", type=str, help="gene or isoform", required = True)

    main(parser)
