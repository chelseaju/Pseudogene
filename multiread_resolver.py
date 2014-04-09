""" multiread_resolver.py
Usage: python multiread_resolver.py -d directory
Input:  -b bamfile
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

	actual_count = 0	
	rank_multiread = []

	for read in bamFH.fetch(chromosome, int(start), int(end)):
		name = read.qname

		# scoring scheme:
		# proper_pair = +50 for paired end
		# mapq = + mapq

		score = 0
		if(ASSIGNMENT.has_key(name)):
			# properly paired = mate is mapped, same chromosome
			if(read.is_paired and not read.mate_is_unmapped and read.rnext == read.tid):
				score += 50

			score += read.mapq

			# position = +1 base
			rank_multiread.append((name, read.tid, read.flag, read.tid, read.pos+1))

		else:
			actual_count += 1

	# actual_count = fragment count
	# needs to multiple 2 for count, margin = 1
	print name, actual_count, count		


"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
"""
def resolver(bamfile, expected):

	print expected
	fh = open(expected, 'rb')
	bamFH = pysam.Samfile( bamfile, "rb" )
	for line in fh:
		(name, count, chromosome, start, end) = line.rstrip().split("\t")

		if(name == "PGOHUM00000258652"):
			query_reads(chromosome, start, end, round(float(count)), bamFH)

	bamFH.close()
	fh.close() 

"""
	Function : read in a list of multireads name
"""
def import_multiread(file):
	fh = open(file, 'r')
	for line in fh:
		line = line.strip()
		ASSIGNMENT[line] = ""
	fh.close()

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir

    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    bamfile = dir + "accepted_hits_sorted.bam"
    multiread_file = dir + "multireads.txt"
    output = dir + "multireads_assignment.txt"
    expected_file = dir + "locus_expectation.txt"

    import_multiread(multiread_file)
    resolver(bamfile, expected_file)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    main(parser)
