""" pseudogene_pseudo_parser.py 
Usage: python pseudogene_pseudo_parser.py -i pseudogene_file -o output_prefix
Input: pseudogene file download from pseudogene database
Output: pseudo gene sequences in fasta format, and location into BED file
        a list of parent protein ensemble ids
Function: parse out the sequence for pseudogenes and parse out parent protein ids. 
Author: Chelsea Ju
Date: 2013-07-19
Note: Modify from pseudogene_parent_parser
"""

import sys, getopt, re
from time import sleep

def usage():
    options = [ ("-i pseudogene_file, --infile = pseudogene file", ""),
                ("-o output_prefix, --outfile = prefix of output files", "")]
    
    print "\nUsage: python pseudogene_parent_parser.py -i infile -o outfile_prefix\nOptions:"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def parser(in_fh, parent_list_fh, pseudo_bed_fh, pseudo_seq_fh):
    parents = {}
   
    next(in_fh) #remove the title line   

    ## pseudogenes
    for rec in in_fh:
        info = rec.split('\t')
        pseudo_id = info[0]
        pseudo_chr = "chr" + info[1]
        pseudo_start = int(info[2])
        pseudo_end = int(info[3])
        pseudo_strand = info[4]
        parent_id = info[5]
        pseudo_exons = info[18]
        pseudo_seq = info[21]
        
        pseudo_block_count = 0
        pseudo_block_size = []
        pseudo_block_start = []
        
        exons = pseudo_exons.split("],")
        for e in exons:
            e = re.sub(r'\[|\]|\s', '', e)
            (e_start, e_end) = e.split(",")
            pseudo_block_count += 1
            pseudo_block_start.append(str(int(e_start) - int(pseudo_start)))
            pseudo_block_size.append(str(int(e_end) - int(e_start) + 1))
        
        pseudo_start = int(pseudo_start) - 1
        pseudo_end = int(pseudo_end) - 1
        pseudo_bed_fh.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" %(pseudo_chr, pseudo_start, pseudo_end, pseudo_id, 0, pseudo_strand, pseudo_start, pseudo_start, 0, pseudo_block_count, (",".join(pseudo_block_size)), (",".join(pseudo_block_start))))
        pseudo_seq_fh.write(">%s\n%s\n" %(pseudo_id, pseudo_seq))
        parents[parent_id] = ""

    
    ## write parents to file
    for p in parents:
        parent_list_fh.write("%s\n" %(p))


def main():
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["infile=", "outprefix="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = str(a)
        elif o in ("-o", "--outfile_prefix"):
            outfile = str(a)
        else:
            assert False, "unhandled option"

    try:
        infile != None
        outfile != None

    except:
        print "\nMissing input file or output prefix"
        usage()
        sys.exit(2)

    parent_list = outfile + "_parents.peptide"
    pseudo_seq = outfile + "_pseudo.fa"
    pseudo_bed = outfile + "_pseudo.bed"

    in_fh = open(infile, 'rb')
    parent_list_fh = open(parent_list, 'w')
    pseudo_seq_fh = open(pseudo_seq, 'w')
    pseudo_bed_fh = open(pseudo_bed, 'w')

    parser(in_fh, parent_list_fh, pseudo_bed_fh, pseudo_seq_fh)
    
    in_fh.close()
    parent_list_fh.close()
    pseudo_seq_fh.close()
    pseudo_bed_fh.close()
    print "Write to file %s" %(parent_list)
    print "Write to file %s" %(pseudo_seq)
    print "Write to file %s" %(pseudo_bed)
    print "DONE\n\n"


if __name__ == "__main__":
    main()


