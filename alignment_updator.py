"""
Usage: python alignment_updator.py -d directory -p paired
Input: -d the directory and prefix of output files
Output: corrected_hits.bam
Function: Go through the bam file, accepted_hits.bam
	1. remove the hits if it is in removed_list.txt
	2. add the resolved_hits to file
	3. result is written to corrected_hits.bam
Author: Chelsea Ju
Date: 2014-05-20
"""

import pysam, argparse, datetime, os, pysam

"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


"""
	Function : go through each hit in bamfile
		if the read exists in delete_list, do not include in the new bam file
		add the resolved hits from resolved_files to the new bam file
"""
def correct_reads(old_bam, new_bam, delete_list, resolved_files):

	old_bam_fh = pysam.Samfile(old_bam, 'rb')
	new_bam_fh = pysam.Samfile(new_bam, 'wb', template=old_bam_fh)

	## remove reads
	echo("Removing Reads")
	for read in old_bam_fh:
		if(not delete_list.has_key(read.qname)):
			new_bam_fh.write(read)

	old_bam_fh.close()

	## add resolved reads
	echo("Adding Resolved Reads")
	for f in resolved_files:
		fh = pysam.Samfile(f, 'rb')
		for new_read in fh:
			new_bam_fh.write(new_read)
		fh.close()

	new_bam_fh.close()

	echo("Writing new bam file to %s" %(new_bam))

"""
	Function : retreive a list of bam file name 
"""
def get_resolved_file(path):

	resolved_list = []
	for f in os.listdir(path):
		if(f.endswith("resolved_hits.bam")):
			resolved_list.append(path + "/" + f)

	return resolved_list

"""
	Function : load the list to memory
"""
def import_remove_reads(remove_file):
	remove_list = {}

	remove_fh = open(remove_file, 'rb')
	for line in remove_fh:
		remove_list[line.rstrip()] = 0

	return remove_list

def main(parser):

	options = parser.parse_args()
	directory = options.dir

	if(directory[-1] != "/"):
		directory += "/"

	old_bam = directory + "accepted_hits.bam"
	new_bam = directory + "corrected_hits.bam"

	remove_list = directory + "correction/removelist.txt"
	resolved_list = get_resolved_file(directory + "correction/")

	reads_to_be_removed = import_remove_reads(remove_list)

	# read correction
	correct_reads(old_bam, new_bam, reads_to_be_removed, resolved_list)
	
if __name__ == "__main__":   

	parser = argparse.ArgumentParser(prog='alignment_updator.py')
	parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
	main(parser)
