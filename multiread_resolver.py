""" multiread_resolver.py
Usage: python multiread_resolver.py -d directory -p is_paired
Input:  -d the directory and prefix of output files ex: select_pseudogene_128_v1/10X_101L_1A/tophat_out
		-p are reads paired? [True|False] 
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
		if(ASSIGNMENT.has_key(name) and ASSIGNMENT[name] == ""):
			# properly paired = mate is mapped, same chromosome
			if(read.is_paired and not read.mate_is_unmapped and read.rnext == read.tid):
				score += 50

			score += read.mapq

			# position = +1 base
			rank_multiread.append((name, score))

		else:
			actual_count += 1

	rank_multiread = sorted(rank_multiread, key=lambda k: -k[1])
	multiread_index = 0
	while(actual_count < count and multiread_index < len(rank_multiread)):
		(select_name, select_score) = rank_multiread[multiread_index]
		ASSIGNMENT[select_name] = (chromosome, start, end)
		multiread_index += 1
		actual_count += 1

"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
"""
def resolver(bamfile, expected, pairend):

	fh = open(expected, 'rb')
	bamFH = pysam.Samfile( bamfile, "rb" )
	for line in fh:
		(name, count, chromosome, start, end) = line.rstrip().split("\t")
		count = float(count)

		if(pairend):
			count = count * 2
		query_reads(chromosome, start, end, int(round(count)), bamFH)

	bamFH.close()
	fh.close() 

"""
	Function : read in a list of multireads name
"""
def import_multiread(file):
	fh = open(file, 'rb')
	for line in fh:
		line = line.strip()
		ASSIGNMENT[line] = ""
	fh.close()

"""
	Function : export multireads assignment
"""
def export_multiread(file):
	fh = open(file, 'w')
	for key in ASSIGNMENT.keys():
		if(ASSIGNMENT[key] != ""):
			(chr, start, end) = ASSIGNMENT[key]
			fh.write("%s\t%s\t%s\t%s\n" %(key, chr, start, end))
		else:
			fh.write("%s\n" %(key))
	fh.close()

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    pairend = options.pair

    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    bamfile = dir + "accepted_hits_sorted.bam"
    multiread_file = dir + "multireads.txt"
    output = dir + "multireads_assignment.txt"
    expected_file = dir + "locus_expectation.txt"
    assignment_file = dir + "multireads_assignment.txt"

    import_multiread(multiread_file)
    resolver(bamfile, expected_file, pairend)
    export_multiread(assignment_file)

    print ""
    print "Writing Multiread Assignment to File : %s" %(assignment_file)    

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)
