""" uniqueread_assigner.py
Usage: python uniqueread_assigner.py -d directory -p is_paired
Input:  -d the directory of input, output files ex: select_pseudogene_128_v1/10X_101L_1A/tophat_out
		-p are reads paired? [True|False] 
Output: uniquely mapped read with its best destination (in sam file)
	and remaining reads for each gene
Function: 1. retrieve the fragments for each locus
		  2. rank all the uniquely mapped fragments and keep them up to the given threshold

Date: 2014-04-04
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse, datetime


ASSIGNMENT={}


"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))

"""
	Function : keep potential reads up to "count"
		return the number of missing reads, and leftover reads
"""
def assign_reads(potential_reads, count, paired):

	previous_name = ""
	actual_count = 0
	leftover = []

	for pr in potential_reads:
		(name, score, read) = pr

		# for pair end reads
		if(paired and name == previous_name):
			continue
		elif(actual_count < count):
			previous_name = name
			actual_count += 1
		else:
			leftover.append(read)
	return (count - actual_count, leftover)

"""
	Function : 
		1. given a region, query the reads in this region
		2. rank all the uniquely mapped read by their quality scores
"""
def query_reads(chromosome, start, end, bamfh):

	unique_reads = []

	for read in bamfh.fetch(chromosome, start, end):
		name = read.qname
		read_start = read.positions[0] + 1
		read_end = read.positions[-1] + 1

		# find unique read
		if(not ASSIGNMENT.has_key(name)):
			# check if the read is completely within the range +/- 80bp
			if(read_start >= (start - 80) and read_start <= (end + 80) and read_end >= (start - 80) and read_end <= (end + 80)):

				# scoring scheme:
				# proper_pair = +50 for paired end
				# mapq = + mapq

				score = 0
				# properly paired = mate is mapped, same chromosome
				if(read.is_paired and not read.mate_is_unmapped and read.rnext == read.tid):
					score += 50

				score += read.mapq

				unique_reads.append((name, score, read))

	sorted_reads = sorted(unique_reads, key=lambda x:(x[1], x[0]), reverse = True)

	return sorted_reads

"""
	Function : break the given line by tab-delimited to retrieve the genomic location
"""
def get_locus(line):
	data = re.split('\t|_', line.rstrip())

	(name, chromosome, start, end, strand, count) = ('', '', 0, 0, '', 0.0)	
	if(len(data) > 3 and data[0] != ''):
		chromosome = data[1]
		start = int(data[2])
		end = int(data[3])

		if(data[0] == "Unknown"):
			name = data[0] + "_" + data[1] + "_" + data[2] + "_" + data[3]
		else:
			name = data[0]
			strand = data[4]

		if(len(data) > 5 and data[5] != ''):
			count = float(data[5])

	return(name, chromosome, start, end, strand, count)

"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
"""
def retrieve_uniqueread(bamfile, expected_file, post_expected_file, pairend, dir):

	wh = open(post_expected_file, 'wb')
	fh = open(expected_file, 'rb')
	bamfh = pysam.Samfile(bamfile, 'rb')

	line = fh.readline()
	for line in fh:
		(name, chromosome, start, end, strand, count) = get_locus(line)
		count = int(round(count))

		potential_reads = query_reads(chromosome, int(start), int(end), bamfh)
		(missing, leftover_reads) = assign_reads(potential_reads, count, pairend)

		# update count after uniqueread assigner
		if(missing > 0):
			wh.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(name, missing, chromosome, start, end, strand))


	 	# output leftover unique reads
	 	if(len(leftover_reads) > 0):
			leftover = dir + "correction/" + name + ".bam"
			leftoverfh = pysam.Samfile(leftover, 'wb', template=bamfh)
			for r in leftover_reads:
				leftoverfh.write(r)
			leftoverfh.close()

	bamfh.close()
	fh.close() 
	wh.close()
	echo("Writing Missing Counts to %s" %(post_expected_file))

"""
	Function : read in a list of multireads name
"""
def import_multiread(file):
	fh = open(file, 'rb')
	for line in fh:
		line = line.strip()
		ASSIGNMENT[line] = ""
	fh.close()

def main(parser):

    options = parser.parse_args()
    directory = options.dir
    pairend = options.pair

    ## check dir
    if(directory[-1] != "/"):
    	directory += "/"

    bamfile = directory + "accepted_hits_sorted.bam"
    multiread_file = directory + "correction/multireads.txt"
    region_file = directory + "correction/"

    expected_file = directory + "correction/expected_count.txt"
    post_expected_file =  directory + "correction/post_expected_count.txt" # count after uniquread assigner
    

    import_multiread(multiread_file)
    retrieve_uniqueread(bamfile, expected_file, post_expected_file, pairend, directory)
    echo("Writing Leftover Reads to %s" %(region_file))


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)


