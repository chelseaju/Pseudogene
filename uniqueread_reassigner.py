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



def blast_sequences(query_file, reference_file, out_file, pairend):
	output = NcbiblastnCommandline(query=query_file, subject=reference_file, outfmt=5, out=out_file)()
	blast_result_record = NCBIXML.parse(open(out_file))

	## Query Name = record.query
	## Query Length = record.query_len
	## Query Match Start / End = hsp.query_start / hsp.query_end
	## Subject Match Start / End = hsp.sbjct_start / hsp.sbjct_end


	for record in blast_result_record:
		for align in record.alignments:
			for hsp in align.hsps:
				print record.query, hsp.query_start, hsp.query_end




"""
	Function : query left over reads from the given gene name, write the sequence to fasta
"""
def query_reads(query_name, region_file):

	query_fasta = region_file + query_name + ".fasta"
	query_bam= region_file + query_name + ".bam"
	query_fh = pysam.Samfile(query_bam, 'rb')
	fasta_fh = open(query_fasta, 'wb')
	reads = []

	for r in query_fh:
		reads.append(r)
		fasta_fh.write(">%s\n%s\n" %(r.qname, r.seq))

	fasta_fh.close()
	return (query_fasta, reads)


def reassign_reads(fasta_file, matrix_file, region_file, removed_bam, resolved_bam, pairend):

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
			if(round(float(query_count)) > 0 and query_name != name):
				blast_output_file = region_file + query_name + "_" + name + ".xml"
				(query_fasta_files, read_obj) = query_reads(query_name, region_file)
				blast_sequences(query_fasta_files, reference_fasta_file, blast_output_file, pairend)




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
    removed_bam = directory + "removed_hits.bam"
    resolved_bam = directory + "resolved_hits.bam"

    target_genes_bed = directory + "correction/target_genes.bed"
    target_genes_fasta = directory + "correction/target_genes.fasta"

    distribution_matrix_file = directory + "correction/distribution_matrix.txt"

    ## read in data
    import_count(post_expected_file)
    count_2_bed(target_genes_bed)
    os.system('bedtools getfasta -name -fi %s -bed %s -fo %s' %(GENOME, target_genes_bed, target_genes_fasta))
    echo("Writing target genes to fasta file %s" %(target_genes_fasta))

    reassign_reads(target_genes_fasta, distribution_matrix_file, region_file, removed_bam, resolved_bam, pairend)

#    retrieve_uniqueread(bamfile, expected_file, post_expected_file, pairend, directory)
 

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_resolver.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="are reads paired? [True|False]", required = True)
    main(parser)


