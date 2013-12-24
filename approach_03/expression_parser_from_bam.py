""" expression_parser_from_bam.py 
Usage: python expression_parser_from_bam.py -r relationship -p parent.bed -s pseudo.bed -e accepted_hit.bam -o outprefix
Input: gene information in bed format, and mappable read in bam format
Output: gene with its expression in fragment count, and FPKM
Author: Chelsea Ju
Function: Read in the gene information, and the expression data. Return a file with expression of parents and pseudo genes in number of fragments and FPKM.
Date: 2013-07-31
Last Modify: 2013-08-06
"""

import sys, getopt, re, pysam, os
from gene import *

PARENT = dict()
PSEUDO = dict()
EXPRESSION = []
def usage():
    options = [ ("-r relationship, --relationship = file describe pseudogene and parent genes relation", ""),
                ("-p parents.bed --parents = parents.bed", ""),
                ("-s pseudo.bed --pseudo = pseudo.bed", ""),
                ("-e accepted_hit.bam, --expression = bam file", ""),
                ("-o outprefix, --outprefix = prefix of output file name", "")
    ]
    
    print "Usage: python expression_parser_from_bam.py -r relationship -p parent.bed -s pseudo.bed -e accepted_hit.bam -o outprefix"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"
    

def retrieve_expression(out_fh, total_fragment):
    global GENELIST
    
    for g in GENELIST:
        out_fh.write("%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\n" %(g.name, g.chr, g.strand, g.start, g.end, g.size, g.get_fragment(), g.get_fpkm(total_fragment)))


    

def extract_relationship(relationship_fh, out_fh):
    global EXPRESSION
    global PSEUDO
    global PARENT
    
    fragment_count = 0
    output_data = []
    line = relationship_fh.readline()
    for line in relationship_fh:
        data = line.split('\t')
        parent_obj = PARENT[data[5]]
        pseudo_obj = PSEUDO[data[0]]
        
        parent_fragments = filter(lambda read: read.chr == parent_obj.chr and read.start >= parent_obj.start and read.end <= parent_obj.end, EXPRESSION)
        pseudo_fragments = filter(lambda read: read.chr == pseudo_obj.chr and read.start >= pseudo_obj.start and read.end <= pseudo_obj.end, EXPRESSION)

        overlapped = set(parent_fragments).intersection(set(pseudo_fragments))

        parent_obj.add_primary_fragment(len(parent_fragments) - len(overlapped))
        parent_obj.add_secondary_fragment(len(overlapped))
        pseudo_obj.add_primary_fragment(len(pseudo_fragments) - len(overlapped))
        pseudo_obj.add_secondary_fragment(len(overlapped))
        
        fragment_count += len(parent_fragments) + len(pseudo_fragments) - len(overlapped)

        output_data.append((data[0], data[5], data[20]))
    

    ## write to file
    out_fh.write('pseudo_name\tpseudo_chr\tpseudo_strand\tpseudo_start\tpseudo_end\tpseudo_size\tpseudo_fragment\tpseudo_fpkm\tpseudo_primary_fragment\tpseudo_secondary_fragment\t')
    out_fh.write('parent_name\tparent_chr\tparent_strand\tparent_start\tparent_end\tparent_size\tparent_fragment\tparent_fpkm\tparent_primary_fragment\tparent_secondary_fragment\tpseudo_type\n')

    for (a,b,c) in output_data:
        a_obj = PSEUDO[a]
        b_obj = PARENT[b]

        out_fh.write('%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t' %(a_obj.name, a_obj.chr, a_obj.strand, a_obj.start, a_obj.end, a_obj.size, a_obj.get_fragment(0), a_obj.get_fpkm(fragment_count, 0), a_obj.get_fragment(1), a_obj.get_fragment(2)))
        out_fh.write('%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t' %(b_obj.name, b_obj.chr, b_obj.strand, b_obj.start, b_obj.end, b_obj.size, b_obj.get_fragment(0), b_obj.get_fpkm(fragment_count, 0), b_obj.get_fragment(1), b_obj.get_fragment(2)))
        out_fh.write('%s\n' %(c))




def build_expression(expression_fh):
    global EXPRESSION
    for read in expression_fh:
        chr = expression_fh.getrname(read.tid)
        strand = "."
        # paired-end read counts as one fragment
        if(read.is_proper_pair and read.is_read1):
            if(read.is_reverse):
                start = read.mpos
            else:
                start = read.pos
            
            end = start + abs(read.isize)
        # if the mate is not properly paired, count it separately
        elif(not read.is_proper_pair):
            start = read.pos
            end = read.aend
        else:
            continue
        
        EXPRESSION.append(Gene(read.qname, chr, strand, start, end))
        

def build_genelist(geneinfo_fh, type):    
    global PARENT
    global PSEUDO

    if(type == 1):
        GENELIST = PARENT
    else:
        GENELIST = PSEUDO
        
    for x in geneinfo_fh:
        gene = from_bed_to_gene(x)
        match = re.match('^(.*)_(.*)', gene.name)
        if(match):
            gname = match.group(1)
        else:
            gname = gene.name
                
        GENELIST[gname] = gene

def from_bed_to_gene(line):
    data = line.split("\t")
    g = Gene(data[3], data[0], data[5], data[1], data[2])
    g.set_size(data[10])
    return g
    
def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:p:s:e:o:", ["relationship=", "parents=", "pseudo=", "expression=", "outprefix=" ])
    
    except getopt.GetoptError, err:
    
    # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-p", "--parents"):
            parent_bed = str(a)
        elif o in ("-s", "--pseudo"):
            pseudo_bed = str(a)
        elif o in ("-r", "--relationship"):
            relationship_file = str(a)
        elif o in ("-o", "--outprefix"):
            outfile = str(a)
        elif o in ("-e", "--expression"):
            expression = str(a)
        else:
        	assert False, "unhandled option"
    
    try:
        parent_bed != None
        pseudo_bed != None
        relationship_file != None
        expression != None
        outfile != None
    
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)

    mapped_outfile = outfile + ".expression"

    parent_fh = open(parent_bed, 'rb')
    pseudo_fh = open(pseudo_bed, 'rb')
    relationship_fh = open(relationship_file, 'rb')
    expression_fh = pysam.Samfile(expression, 'rb')
    mapped_fh = open(mapped_outfile, 'w')
 
    # build gene information from bed files
    build_genelist(parent_fh, 1)
    build_genelist(pseudo_fh, 2)
    build_expression(expression_fh)
    extract_relationship(relationship_fh, mapped_fh)
    
    
    parent_fh.close()
    pseudo_fh.close()
    relationship_fh.close()
    expression_fh.close()
    mapped_fh.close()
    
    print "Write to file %s" %(mapped_outfile)
    print "DONE\n\n"



if __name__ == "__main__":
    main()





