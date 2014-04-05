""" alignment_resolver.py
Usage: python alignment_resolver.py -d directory -e expected_count
Input:  -b bamfile
        -e expected_count
Output: a list of multiple alignment with its best destination
Function: 1. retrieve the fragments for each locus in lasso coefficient
		  2. retain all the uniquely mapped fragments
		  3. rank the multireads by their quality score
		  4. if the number of fragments has reached the expectation, stop
		  5. otherwise, take the next multireads that haven't been assigned yet 
		  6. repeat 4 and 5

Date: 2014-04-04
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

ASSIGNMENT={}

"""
	Function : 
		1. given a region, query the reads in this region
		2. keep all the uniquely mapped reads
		3. rank the multireads by their quality scores
		4. if the number of fragments has reach the expectation, stop
		5. otherwise, take the next multireads that haven't been assigned yet
		6. assigned this multiread to this region by putthing it into hash
		7. repeat 5-7
"""
def query_reads(chromosome, start, end, count, bamFH):

	for read in bamFH.fetch(chromosome, int(start), int(end)):
		print read.pos



"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
"""
def resolver(bamfile, expected):

	fh = open(expected, 'rb')
	bamFH = pysam.Samfile( bamfile, "rb" )
    
	for line in fh:
		(name, count, chromosome, start, end) = line.rstrip().split(" ")
		query_reads(chromosome, start, end, count, bamFH)


	bamFH.close()
	fh.close() 


def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    expected_file = options.expectation

    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    bamfile = dir + "accepted_hits_sorted.bam"
    output = dir + "multireads_assignment.txt"

    multireads = resolver(bamfile, expected_file)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='alignment_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-e", "--expectation", dest="expectation", type=str, default=80, help="expected fragment count by lasso", required = True)
    main(parser)
