""" uniqueread_reassigner.py
Usage: python uniqueread_reassigner.py -d directory -p is_paired
Input:  -d the directory of input, output files ex: select_pseudogene_128_v1/10X_101L_1A/tophat_out
		-p are reads paired? [True|False] 
Output: retrieve the unique read from close sequence similarity region, and reassign the positions
Function: 

Date: 2014-05-13
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse, datetime, subprocess
import blast_parser

from Bio.Blast.Applications import NcbiblastnCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


## GENOME FILE

# LAB
GENOME="/home/chelseaju/Database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"
# HOFFMAN
#GENOME = "/u/home/c/chelseaj/project/database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"

BLAST_THRESHOLD = 0.000001
COUNT={}

"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))

"""
	Function : convert the modified blast result to sam
"""
def blast2sam(oldread, blast, tid, start):

	blast = sorted(blast, key = lambda x: x[0])

	blast_index = 0
	pair_index = -1

	if(len(blast) > 1):
		if(oldread.is_paired and oldread.is_read2):
			blast_index = 1
			pair_index = 0
		else:
			pair_index = 1

	read = pysam.AlignedRead()
	read.qname = oldread.qname
	read.seq = oldread.seq
	read.pos = start + blast[blast_index][1] -1
	read.cigar = blast[blast_index][4]
	read.rname = tid
	read.mapq = oldread.mapq
	read.mrnm = tid if pair_index != -1 else -1
	read.mpos = start + blast[pair_index][1] - 1 if pair_index != -1 else -1 # position of mate
	read.tlen = 0
	read.qual = oldread.qual
	
	# flags
	read.flag = oldread.flag

	if(len(blast) > 1):
		read.is_proper_pair = True
		read.mate_is_unmapped = False
		read.mate_is_reverse = blast[pair_index][5]
	elif(read.is_paired):
		read.is_proper_pair = False
		read.mate_is_unmapped = True

	read.is_reverse = blast[blast_index][5]

	# build insert size
	if(pair_index != -1):
		length =  max(blast[pair_index][2], blast[blast_index][2]) - min(blast[pair_index][1], blast[blast_index][1])
		read.tlen = length if blast[blast_index][1] < blast[pair_index][1] else length*(-1)

#		print (blast[blast_index][1], blast[blast_index][2], blast[pair_index][1], blast[pair_index][2], 
#			length if blast[blast_index][1] < blast[pair_index][1] else length*(-1))

	return read


"""
	Function : select the first hits up to "count"
"""
def select_blast_hits(count, blast_hits, pairend):

	selection = {}
	actual_count = 0

	for h in blast_hits:
		name_splitter = re.match(r"(.*)\/(\d)$", h[0])
		if(name_splitter):
			name = name_splitter.group(1)
		else:
			name = h[0]

		# for pair end reads
		if(pairend and selection.has_key(name)):
			selection[name].append(h)

		elif(actual_count < count):
			selection[name] = [h]
			actual_count += 1

	return(count - actual_count, selection)

"""
	Function : 1. put blast_hits into resolved_bam up to query_count
			2. remove them from region_file_prefix
			3. flag them for removal in removed_bam 
"""
def resolve_blast_hits(selected_hits, region_file_prefix, chromosome, start, name):

	region_bam = region_file_prefix + ".bam"
	region_bam_reassign = region_file_prefix + "_reassign.bam"
	region_bam_new = region_file_prefix + "_new.bam"
	region_bam_remove = region_file_prefix + "_remove.bam"

	region_bam_fh = pysam.Samfile(region_bam, 'rb')
	region_bam_fh_reassign = pysam.Samfile(region_bam_reassign, 'wb', template = region_bam_fh)
	region_bam_fh_new = pysam.Samfile(region_bam_new, 'wb', template = region_bam_fh)
	region_bam_fh_remove = pysam.Samfile(region_bam_remove, 'wb', template = region_bam_fh)

	for read in region_bam_fh:
		if(selected_hits.has_key(read.qname) and len(selected_hits[read.qname]) > 0):
			region_bam_fh_remove.write(read)

			new_read = blast2sam(read, selected_hits[read.qname], region_bam_fh.gettid(chromosome), int(start))
			region_bam_fh_reassign.write(new_read)
# for debug
# 			print selected_hits[read.qname], name, chromosome, start, int(start) + selected_hits[read.qname][0][1], int(start) + selected_hits[read.qname][0][2]
		else:
			region_bam_fh_new.write(read)

	region_bam_fh.close()
	region_bam_fh_new.close()
	region_bam_fh_reassign.close()
	region_bam_fh_remove.close()

	# replace the region bam file by the new region bam file
#	os.system('mv %s %s' %(region_bam_new, region_bam))


"""
	Function : query left over reads from the given gene name, write the sequence to fasta
"""
def query_reads(region_file_prefix):

	query_fasta = region_file_prefix + ".fasta"
	query_bam= region_file_prefix + ".bam"
	query_fh = pysam.Samfile(query_bam, 'rb')
	fasta_fh = open(query_fasta, 'wb')
	reads = []

	for r in query_fh:
		reads.append(r)

		name = r.qname
		if(r.is_paired and r.is_read1):
			name += "/1"
		elif(r.is_paired and r.is_read2):
			name += "/2"
		else:
			name += "/0"

		fasta_fh.write(">%s\n%s\n" %(name, r.seq))

	fasta_fh.close()
	return (query_fasta, reads)


def reassign_reads(fasta_file, matrix_file, region_file, pairend):

	for name in COUNT.keys():
		(count, chr, start, end) = COUNT[name]
		fasta = subprocess.check_output(["grep", "-A1", name, fasta_file])
		(title, sequence, empty) = fasta.split("\n")

		# write reference fasta to file
		reference_fasta_file = region_file + "ref_" + name + ".fasta"
		ref_fasta_fh = open(reference_fasta_file, 'w')
		ref_fasta_fh.write(">%s\n%s\n" %(title, sequence))
		ref_fasta_fh.close()

		## grep distribution matrix
		matrix_info = subprocess.check_output(["grep", name, matrix_file])
		(header, distribution, empty) = matrix_info.split("\n")

		distribution_info = filter(lambda x:float(x[1]) > 0, sorted(zip(header.split("\t"), distribution.split("\t")[1:]), key=lambda x: x[1], reverse = True))

		for (query_name, query_count) in distribution_info:
			if(round(float(query_count)) > 0 and query_name != name and COUNT[name][0] > 0):
				region_file_prefix = region_file + query_name
				blast_output_file = region_file_prefix + "_" + name + ".xml"
				(query_fasta_files, read_obj) = query_reads(region_file_prefix)

				#call blastn
				NcbiblastnCommandline(query=query_fasta_files, subject=reference_fasta_file, outfmt=5, out=blast_output_file)()
				blast_hits = blast_parser.parse_blast_result(blast_output_file)
				(resolve_count, selected_hits) = select_blast_hits(count, blast_hits, pairend)
				resolve_blast_hits(selected_hits, region_file_prefix, chr, start, name)

				# upate counter
				COUNT[name] = (resolve_count, chr, start, end)


"""
	Function : convert the target genes to bed format
"""
def count_2_bed(bedfile):
	fh = open(bedfile, 'wb')
	for k in COUNT.keys():
		fh.write("%s\t%s\t%s\t%s\n" %(COUNT[k][1], COUNT[k][2], COUNT[k][3], k))
	fh.close()

	echo("Writing target genes to bed file %s" %(bedfile))


"""
	Function : read in the expected_count.txt
"""
def import_count(file):
	fh = open(file, 'rb')
	for line in fh:
		(name, count, chromosome, start, end) = line.rstrip().split("\t")
		COUNT[name] = (int(count), chromosome, start, end)
	fh.close()



def main(parser):

    options = parser.parse_args()
    directory = options.dir
    pairend = options.pair

    ## check dir
    if(directory[-1] != "/"):
    	directory += "/"

    region_file = directory + "correction/"
    post_expected_file =  directory + "correction/post_expected_count.txt" # count after uniquread assigner

    target_genes_bed = directory + "correction/target_genes.bed"
    target_genes_fasta = directory + "correction/target_genes.fasta"

    distribution_matrix_file = directory + "correction/distribution_matrix.txt"

    ## read in data
    import_count(post_expected_file)
    count_2_bed(target_genes_bed)
    os.system('bedtools getfasta -name -fi %s -bed %s -fo %s' %(GENOME, target_genes_bed, target_genes_fasta))
    echo("Writing target genes to fasta file %s" %(target_genes_fasta))


    current_count = 0
    total_count =  sum([COUNT[k][0] for k in COUNT.keys()])

    while(total_count != current_count):
	    current_count = sum([COUNT[k][0] for k in COUNT.keys()])
	    reassign_reads(target_genes_fasta, distribution_matrix_file, region_file, pairend)
	    total_count =  sum([COUNT[k][0] for k in COUNT.keys()])

	 

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)


