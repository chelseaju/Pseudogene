""" check_overlap_pseudo_parents.py
Usage: python check_overlap_pseudo_parents.py -p parents.bed -s pseudo.bed -o out_prefix
Input: BED format files for parents and pseudogene information
Output: overlapped area, parent ids and pseudogene ids
Function: check for overlapping area for the list of parents an pseudogenes 
Author: Chelsea Ju
Date: 2013-07-31
"""

import sys, getopt, re
from gene import *

#global vairable
PARENTLIST = []


def usage():
    options = [ ("-p parents.bed, --parent = parent information in bed format", ""),
                ("-s pseudo.bed, --pseudo = pseudogene information in bed format", ""),
                ("-o out_file, --output = prefix of the output file")]
    
    print "Usage: python check_overlap_pseudo_parents.py -p parents.bed -s pseudo.bed -o out_file\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"



def build_parent_database(parent_fh):
    global PARENTLIST
    
    for x in parent_fh:
        PARENTLIST.append(from_bed_to_gene(x))

    
def scan_pseudogenes(pseudo_fh, complete_fh, short_fh):
    global PARENTLIST
    
    for x in pseudo_fh:
        pseudogene = from_bed_to_gene(x)
        result = filter(lambda gene: gene.chr == pseudogene.chr and
                                     gene.strand == pseudogene.strand and (
                                     (gene.start < pseudogene.start and gene.end > pseudogene.end) or
                                     (pseudogene.start <= gene.start and pseudogene.end > gene.start) or
                                     (pseudogene.end >= gene.end and pseudogene.start < gene.end and pseudogene.start > gene.start)) 
                        , PARENTLIST)
        ### note: 
        ##     let -|--|- be pseudogene, -*--*- be parent gene
        ##        gene.start < pseudogene.start and gene.end > pseudogene.end test for:
        ##            --*--|--|--*--
        ##        pseudogene.start <= gene.start and pseudogene.end > gene.start test for:
        ##            --|--*--|--*--; --|--*--*|; --|--*--*--|--; |*--|--*--; |*--*|; |*--*--|
        ##        pseudogene.end >= gene.end and pseudogene.start < gene.end and pseudogene.start > gene.start) test for:
        ##            *--|--*--|--;  *--|--*|-- 
        
        if len(result) > 0:
            short_fh.write("%s\n" %(pseudogene.name))
            complete_fh.write("%s\t%s\t%s\n" %(pseudogene.name, pseudogene, result))


def from_bed_to_gene(line):
    data = line.split("\t")
    return Gene(data[3], data[0], data[5], data[1], data[2])

def main():
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:s:o:", ["parent=", "pseudo=", "outfile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-p", "--parent"):
            parent = str(a)
        elif o in ("-s", "--pseudo"):
            pseudo = str(a)
        elif o in ("-o", "--outfile"):
            outfile = str(a)
        else:
            assert False, "unhandled option"

    try:
        parent != None
        pseudo != None
        outfile != None

    except:
        print "\nMissing input files"
        usage()
        sys.exit(2)

    out_complete = outfile + "_overlapped_pseudo_parents.txt"
    out_short = outfile + "_overlapped_pseudogene.txt"
    parent_fh = open(parent, 'rb')
    pseudo_fh = open(pseudo, 'rb')

    out_complete_fh = open(out_complete, 'w')
    out_short_fh = open(out_short, 'w')
        
    build_parent_database(parent_fh)
    scan_pseudogenes(pseudo_fh, out_complete_fh, out_short_fh)
   
    parent_fh.close()
    pseudo_fh.close()
    out_complete_fh.close()
    out_short_fh.close()
    
    print "Write to file %s" %(out_complete)
    print "Write to file %s" %(out_short)
    print "DONE\n\n"


if __name__ == "__main__":
    main()
