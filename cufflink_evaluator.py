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

XLOC = {}
EXPECTATION = {}


import sys, re, os, subprocess, random, argparse

"""
	Function : build gene gene id
"""
def build_name(file):
	fh = open(file, 'rb')
	for line in fh:
		line = line.rstrip()
		data = line.split("\t")
		x = data[0]
		alternative_ids = data[3]
		XLOC[x] = alternative_ids.split(",")
	fh.close()

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

def combine_expectation_observation(result):

	evaluation_data = []
	seen_ids = {}

	# mark all expected ids
	for k in EXPECTATION.keys():
		seen_ids[k] = False

	fh = open(result, 'rb')
	line = fh.readline()
	for line in fh:
		line = line.rstrip()
		data = line.split("\t")
		xid = data[0]
		xgene = data[1]

		if(EXPECTATION.has_key(xgene)):
			data = data + [xgene] + EXPECTATION[xgene]
			seen_ids[xgene] = True
		else:
			found = False
			if(XLOC.has_key(xid)):
				for a_gene in XLOC[xid]:
					if(EXPECTATION.has_key(a_gene)):
						data = data + [a_gene] + EXPECTATION[a_gene]
						seen_ids[a_gene] = True
						found = True
						break

			if(not found):
				data = data + ["-", "0","0"]

		evaluation_data.append(data)

	fh.close()

	for k in sorted(seen_ids, key=seen_ids.get):
		if(seen_ids[k] == False):
			data = ["-",k,"-","-","-","-","-",k] + EXPECTATION[k]
			evaluation_data.append(data)
		else:
			break


	return evaluation_data

def export_data(result, file):
	fh = open(file, 'w')
	fh.write("XLOC\tGENENAME\tOVERLAP\tFPKM_A\tFPKM_B\tCOUNT_A\tCOUNT_B\tEXPECTED\tEXPECTED_A\tEXPECTED_B\n")
	for r in result:
		fh.write("\t".join(r))
		fh.write("\n")

	fh.close()

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

    main(parser)
