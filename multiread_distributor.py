""" multiread_distributor.py
Usage: python multiread_distributor.py -d directory -p is_paired
Input:  -d the directory and prefix of output files ex: select_pseudogene_128_v1/10X_101L_1A/tophat_out
		-p are reads paired? [True|False] 
Output: 
Function: 1. retrieve the fragments for each locus in lasso coefficient
		  2. retain all the uniquely mapped fragments
		  

Date: 2014-04-04
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

DISTRIBUTION={}
EXPECTATION={}

"""
	Function : 
		1. given a region, query the reads in this region
		2. keep all the uniquely mapped reads
		3. for multireads, record the best alignment for the given positions
		3. return the value of max((expected_count - actual_count), 0) 
"""
def query_reads(chromosome, start, end, count, gene, bamFH):

	actual_count = 0
	for read in bamFH.fetch(chromosome, int(start), int(end)):
		name = read.qname

		if(read.positions[0]+1 >= int(start) and read.positions[0] +1 <= int(end)):

			# scoring scheme:
			# proper_pair = +50 for paired end
			# mapq = + mapq

			score = 0
			if(DISTRIBUTION.has_key(name)):	
				# properly paired = mate is mapped, same chromosome
				if(read.is_paired and not read.mate_is_unmapped and read.rnext == read.tid):
					score += 50
				score += read.mapq

				if(DISTRIBUTION[name].has_key(gene)):
					(assign_chr, assign_start, assign_end, assign_score) = DISTRIBUTION[name][gene]
					if(score > assign_score):
						DISTRIBUTION[name][gene] = (chromosome, int(read.positions[0]+1), read.pnext+1, score)
				else:
					DISTRIBUTION[name][gene] = (chromosome, int(read.positions[0]+1), read.pnext+1, score)
			else:
				actual_count += 1

#	rank_multiread = sorted(rank_multiread, key=lambda k: -k[1])
#	print chromosome, start, end, count, actual_count
#	multiread_index = 0
#	while(actual_count < count and multiread_index < len(rank_multiread)):
#		(select_name, select_score) = rank_multiread[multiread_index]
#		ASSIGNMENT[select_name] = (chromosome, start, end)
#		multiread_index += 1
#		actual_count += 1

	return max(count - actual_count, 0)


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
		missing_count = query_reads(chromosome, start, end, int(round(count)), name, bamFH)

		if(missing_count > 0):
			EXPECTATION[name] = missing_count

	bamFH.close()
	fh.close() 

"""
	Function : read in a list of multireads name
"""
def import_multiread(file):
	fh = open(file, 'rb')
	for line in fh:
		line = line.strip()
		DISTRIBUTION[line] = {}
	fh.close()

"""
	Function : export multireads assignment
"""
def export_multiread(file):
	fh = open(file, 'w')
	for key in DISTRIBUTION.keys():
		if(len(DISTRIBUTION[key]) == 1):
			gene_name = DISTRIBUTION[key].keys()[0]
			(chromosome, start, pair_start, score) = DISTRIBUTION[key][gene_name]
			if(EXPECTATION.has_key(gene_name)):
				if(pair_start == 0):
					EXPECTATION[gene_name] -=1
				else:
					EXPECTATION[gene_name] -=2

			fh.write("%s\t%s\t%s\t%s\n" %(key, chromosome, start, pair_start))

		else:
			print key, DISTRIBUTION[key]


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
   
    parser = argparse.ArgumentParser(prog='multiread_disributor.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)
