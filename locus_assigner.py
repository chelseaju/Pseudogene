"""
Usage: python locus_assigner.py -y expected_abundance -x distribution_matrix -o outfile
Input: -y expected abundance predicted by lasso
	-x distribution matrix reconstructed by A matrix
Output: -o a list of gene ids with their abundance and genomic locus, separated by tab
Function: link the gene, expected gene count, and genomic location to one file
Author: Chelsea Ju
Date: 2014-05-12
"""

import sys, os, re, argparse, random, datetime

# LAB
DB = "/home/chelseaju/Database/Pseudogene/ParentENSG_Pseudogene_74.bed"

# HOFFMAN
DB = "/u/home/c/chelseaj/project/database/Pseudogene/ParentENSG_Pseudogene_74.bed"

# MAC
DB = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ParentENSG_Pseudogene_74.bed"

ABUNDANCE = {}


"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


"""
	Function : parse the first line of X matrix
"""
def parseX(xfile):
	fh = open(xfile, 'rb')
	wh = open("temp_file.txt", 'wb')

	line = fh.readline()

	for name in line.rstrip().split("\t"):		
		gene_name = name.split("_")

		if(len(gene_name) > 1):
			ABUNDANCE[name] = (0, gene_name[1], gene_name[2], gene_name[3])
			wh.write("%s\n" %(name))
		else:
			ABUNDANCE[gene_name[0]] = (0, "", "", "")

	fh.close()
	wh.close()

"""
	Function : parse Y abundance
"""
def parseY(yfile):
	fh = open(yfile, 'rb')
	for line in fh:
		data = line.rstrip().split("\t")
		if(len(data) == 2):
			ABUNDANCE[data[0]][0] = data[1]
	fh.close()

def main(parser):
	options = parser.parse_args()
	yfile = options.yvector
	xfile = options.xmatrix
	outfile = options.output

	parseX(xfile)
	parseY(yfile)

	for x in ABUNDANCE.keys():
		print x, ABUNDANCE[x]



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='locus_assigner.py')
    parser.add_argument("-y", "--yvector", dest="yvector", type=str, help="expected abundance", required = True)
    parser.add_argument("-x", "--xmatrix", dest="xmatrix", type=str, help="distribution matrix", required = True)
    parser.add_argument("-o", "--output", dest="output", type=str, help="gene ids and their abundance and genomic locus", required = True)

    main(parser)

