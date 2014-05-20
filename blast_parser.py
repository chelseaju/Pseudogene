"""
	Parse Blast Result
"""

import sys, re, pysam, os, random, argparse, datetime, subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

"""
	Function : return a list of blast result, with the bast cigar score
"""
def parse_blast_result(blast_outfile):
	blast_result_record = NCBIXML.parse(open(blast_outfile))
	result_list = [] # each recored contains (query_name, best_start, best_end, score, cigar)


	for record in blast_result_record:
		if(len(record.alignments) > 0):

			(start, end, score, cigar, rc) = parse_record(record)
			result_list.append((record.query, start, end, score, cigar, rc))

	return sorted(result_list, key= lambda x: (x[3], x[0]), reverse = True)


def parse_record(record):
	
	align = record.alignments[0]

	potential_hits = []

	for hsp in align.hsps:
#			print record.query, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end
		cigar = get_cigar(hsp.sbjct, hsp.match, hsp.query)
		potential_hits.append((hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, cigar))

	# sort hits by start
	potential_hits = sorted(potential_hits, key = lambda x:x[0])
	potential_hits = remove_duplicate_hits(potential_hits)

	(start, end, score, cigar, rc) = concatenate_hits(potential_hits, record.query_length, record.query)
	return (start, end, score, cigar, rc)

def concatenate_hits(hits, length, name):

	previous_q_start = 0
	previous_q_end = 0

	previous_s_start = 0
	previous_s_end = 0

	best_start = 0
	best_end = 0

	cigar_in_string = ""

	for h in hits:
		(qstart, qend, sstart, send, cigar) = h

		# check overlap
		if(previous_q_end != 0 and qstart <= previous_q_end):
			overlap = previous_q_end  - qstart + 1
			qstart = qstart + overlap

			if(sstart < send):
				sstart = sstart + overlap
			else:
				sstart = sstart - overlap

			cigar = cigar[overlap:]

		# add junction 
		if(previous_s_end != 0 and abs(sstart - previous_s_end) != 1):
#			print name, previous_s_end, sstart,  "add deletion"
			cigar_in_string += "%dN" %(abs(sstart - previous_s_end - 1))

		# add soft clip
		if(previous_q_start == 0 and qstart != 1):
			cigar_in_string = "%dS" %(qstart - 1) + cigar_in_string


		cigar_in_string += cigar_to_string(cigar)

		# update start and end
		if(best_start == 0):
			best_start = sstart

		best_end = send

		previous_q_start = qstart
		previous_q_end = qend
		previous_s_start = sstart
		previous_s_end = send

	## add soft clip at the end
	if(previous_q_end != length):
		cigar_in_string += "%dS" %(length - previous_q_end)

#	print name, best_start, best_end, cigar_in_string
	(score, cigar_in_bam) = cigar_to_bam(cigar_in_string)

	if(best_start > best_end):
		cigar_in_bam.reverse()

	reverse_complement = (best_start > best_end)

	return (min(best_start, best_end), max(best_start, best_end), score, cigar_in_bam, reverse_complement)

def remove_duplicate_hits(hits):

	previous_start = 0
	previous_end = 0
	previous_length = 0
	flag_for_remove = []

	for i in xrange(len(hits)):
		h = hits[i]

		# overlap
		if(previous_start > 0 and h[0] < previous_start + 10):
			if(h[1] - h[0]) > previous_length: #replace
				previous_start = h[0]
				previous_length = h[1] - h[0]
				flag_for_remove.append((previous_start, previous_end))
			else:
				flag_for_remove.append((h[0], h[1]))

		previous_start = h[0]
		previous_end = h[1]

	hits = filter(lambda x: (x[0],x[1]) not in flag_for_remove, hits)

	return hits

"""
	Function : scan through cigar string, remove duplication, convert it to bam format

	M	BAM_CMATCH	0
	I	BAM_CINS	1
	D	BAM_CDEL	2
	N	BAM_CREF_SKIP	3
	S	BAM_CSOFT_CLIP	4
	H	BAM_CHARD_CLIP	5
	P	BAM_CPAD	6
	=	BAM_CEQUAL	7
	X	BAM_CDIFF	8
"""
def cigar_to_bam(cigar_in_string):


	lookup_table = { "M": 0, "I":1, "D":2, "N":3, "S":4, "H":5, "P":6, "=":7, "X":8}

	score = 0
	cigar_array = []
	previous_chr = ""
	previous_count = ""

	cigar_splitter = re.compile(r'(\d+\w)').findall(cigar_in_string)

	for s in cigar_splitter:
		current_chr = s[-1]
		if(current_chr == previous_chr):
			previous_count +=  int(s[0:-1])
		else:
			if(previous_chr != ""):
				cigar_array.append((lookup_table[previous_chr], previous_count))

			previous_count = int(s[0:-1])
			previous_chr = current_chr

		if(current_chr == "M"):
			score += int(s[0:-1])

	if(previous_chr != ""):
		cigar_array.append((lookup_table[previous_chr], previous_count))

	return(score, cigar_array)



"""
	Function : convert cigar to string
"""
def cigar_to_string(c):

	pre_label = ""
	count = 0
	cigar_string = ""
	for i in c:
		if(i == pre_label):
			count += 1
		else:
			if(pre_label != ""):
				cigar_string += "%d%s" %(count, pre_label)
			pre_label = i
			count = 1

	if(pre_label != ""):
		cigar_string += "%d%s" %(count, pre_label)

	return cigar_string



"""
	Function : build cigar score based on blast result
"""
def get_cigar(sbjct, match, query):

	cigar = ["."]*len(match)
	for i in xrange(len(match)):
		# match or mismatch
		if(match[i] == "|" or (sbjct[i] != "-" and query[i] != "-")):
			cigar[i] = "M"
		elif(sbjct[i] == "-"):
			cigar[i] = "I"
		else:
			cigar[i] = "D"

	return cigar
