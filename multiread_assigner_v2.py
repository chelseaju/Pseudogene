""" multiread_assigner.py
Usage: python multiread_assigner.py -d directory -p is_paired
Input:  -d the directory of input, output files ex: select_pseudogene_128_v1/10X_101L_1A/tophat_out
		-p are reads paired? [True|False] 
Output:  1. update multiread list with the new gene name it assigned to
		2. output reolved_hits.bam for a list of selected Multireads
		3. output remove_list to flag a list of reads to be removed
		4. output post_post_expected_count.txt for the remaining counts
Function: 1. Rank the list based on the number of missing reads (large to small)
		  2. Iterate through the list, and assign the multireads fall in this region to itself
		  3. Update the multiread list if the read has already been assigned
		  4. Write the newly assigned multiread to resolved_hits.bam
		  5. Repeate 2-4 until no more multireads to be assigned

Date: 2014-05-12
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse, datetime


ASSIGNMENT={}
COUNT={}

"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


"""
	Function : keep potential reads up to "count", update ASSIGNMENT
		output the assigned multireads to resolvedbam
		update the remove list (add the newly assigned multireads for removal)
"""
def assign_reads(geneName, potential_reads, count, paired, resolvedbam, removefh):

	previous_name = ""
	actual_count = 0	

	for pr in potential_reads:
		(name, score, read) = pr

		# for pair end reads
		if(paired and name == previous_name):
			resolvedbam.write(read)
			ASSIGNMENT[name] = geneName
		elif(actual_count < count):
			previous_name = name
			actual_count += 1
			resolvedbam.write(read)
			removefh.write("%s\n" %(name))
			ASSIGNMENT[name] = geneName

	return (count - actual_count)



"""
	Function : 
		1. given a region, query the reads in this region
		2. rank all the multireads by their quality scores
"""
def query_reads(chromosome, start, end, bamfh):

	multi_reads = []

	for read in bamfh.fetch(chromosome, start, end):
		name = read.qname
		read_start = read.positions[0] + 1
		read_end = read.positions[-1] + 1

		# find multi read
		if(ASSIGNMENT.has_key(name) and ASSIGNMENT[name] == ""):

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

				multi_reads.append((name, score, read))

	sorted_reads = sorted(multi_reads, key=lambda x:(x[1], x[0]))

	return sorted_reads



"""
	Function : iterat throught the gene list with expected count
		query the fragments within the region
		to avoid assigning all multireads to one region, "proportion" can be adjusted
		1 is optimal so far
"""
def retrieve_multiread(bamfh, resolvedfh, pairend, proportion, removefh):


	for data in sorted(COUNT.items(), key=lambda x: x[1][0]):
		name = data[0]
		(expected, chromosome, start, end, strand) = data[1]

		if(expected > 0):
			potential_reads = query_reads(chromosome, int(start), int(end), bamfh)

			missing = assign_reads(name, potential_reads, int(round(expected*proportion)), pairend, resolvedfh, removefh)

			remains = expected - int(round(expected*proportion)) + missing
			COUNT[name] = (remains, chromosome, start, end, strand) # update count
"""
	Function : read in the expected_count.txt
"""
def import_count(file):
	fh = open(file, 'rb')
	for line in fh:
		(name, count, chromosome, start, end, strand) = line.rstrip().split("\t")
		COUNT[name] = (int(count), chromosome, start, end, strand)
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
	Function : update expected_count
"""
def update_count(file):
	fh = open(file, 'w')
	for key in COUNT.keys():
		(count, chromosome, start, end, strand) = COUNT[key]

		if(count > 0):
			fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(key, count, chromosome, start, end, strand))

	fh.close()
	echo("Writing Updated Read Count to %s" %(file))

"""
	Function : update multireads with the assigned region
"""
def update_multiread(file):
	fh = open(file, 'w')
	for k in ASSIGNMENT.keys():
		fh.write("%s\t%s\n" %(k, ASSIGNMENT[k]))
	fh.close()

def main(parser):

    options = parser.parse_args()
    directory = options.dir
    pairend = options.pair

    ## check dir
    if(directory[-1] != "/"):
    	directory += "/"

    bamfile = directory + "accepted_hits_sorted.bam"
    resolvedbam = directory + "correction/resolved_hits.bam"
    multiread_file = directory + "correction/multireads.txt"
    remove_list = directory + "correction/removelist.txt"

    expected_file =  directory + "correction/post_expected_count.txt" # count after uniquread assigner
    post_expected_file =  directory + "correction/post_post_expected_count.txt" # count after uniquread assigner
    
    bamfh = pysam.Samfile(bamfile, 'rb')
    resolvedfh = pysam.Samfile(resolvedbam, 'wb', template = bamfh)
    removefh = open(remove_list, 'w')

    # load data into memory
    import_multiread(multiread_file)
    import_count(expected_file)
    
    # Assign multireads
    pre_count = 0
    total_count = sum([c[1][0] for c in COUNT.items()])

    counter = 1
    while(total_count > 0 and total_count != pre_count):    	
    	echo("Updating Multireads Count, Cycle %d" %(counter))
    	pre_count = total_count
    	retrieve_multiread(bamfh, resolvedfh, pairend, 1, removefh)
    	total_count = sum([c[1][0] for c in COUNT.items()])
    	counter += 1

    # Update Multiread Assignment, and Read Counts
    update_multiread(multiread_file)
    update_count(post_expected_file)

    bamfh.close()
    resolvedfh.close()
    removefh.close()
    echo("Writing Selected Multireads to %s" %(resolvedbam))
    echo("Writing Update Multireads at %s" %(multiread_file))
    echo("Writing Delete Multireads to %s" %(remove_list))



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver_v2.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)
