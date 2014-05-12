""" unique_read_assigner.py
Usage: python unique_read_assigner.py -d directory -p is_paired
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
	Function : keep potential reads up to "count", and write it to outbamfh
		the rest is written to leftoverfh
		return the number of missing reads
"""
def assign_reads(potential_reads, count, paired, outbamfh, leftoverfh):

	previous_name = ""
	actual_count = 0

	for pr in potential_reads:
		(name, score, read) = pr

		# for pair end reads
		if(paired and name == previous_name):
			outbamfh.write(read)

		elif(actual_count < count):
			previous_name = name
			outbamfh.write(read)
			actual_count += 1

		else:
			leftoverfh.write(read)

	return (count - actual_count)

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

	sorted_reads = sorted(unique_reads, key=lambda x:(x[1], x[0]))

	return sorted_reads

"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
"""
def retrieve_uniqueread(bamfh, expected_file, outbamfh, pairend, dir):

	fh = open(expected_file, 'rb')
	for line in fh:
		(name, count, chromosome, start, end) = line.rstrip().split("\t")
		count = float(count)

		leftover = dir + "correction/" + name + ".bam"
		leftoverfh = pysam.Samfile(leftover, 'wb', template=bamfh)

		if(name == "ENSG00000197223"):
			potential_reads = query_reads(chromosome, int(start), int(end), bamfh)
			missing = assign_reads(potential_reads, int(round(count)), pairend, outbamfh, leftoverfh)
			print missing

		leftoverfh.close()
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

def main(parser):

    options = parser.parse_args()
    directory = options.dir
    pairend = options.pair

    ## check dir
    if(directory[-1] != "/"):
    	directory += "/"

    bamfile = directory + "accepted_hits_sorted.bam"
    multiread_file = directory + "correction/multireads.txt"
    outbam = directory + "corrected_hits.bam"
    expected_file = directory + "correction/expected_count.txt"
    
    bamFH = pysam.Samfile(bamfile, 'rb')
    outbamFH = pysam.Samfile(outbam, "wb", template=bamFH)

    import_multiread(multiread_file)
    retrieve_uniqueread(bamFH, expected_file, outbamFH, pairend, directory)

    bamFH.close()
    outbamFH.close()

    echo("Writing Leftover Reads to %s" %(directory+"/corection/"))
    echo("Writing New Assignment to %s" %(outbam))


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)
