""" compute_expected_fpkm.py 
Usage: python compute_expected_fpkm.py -p parent.bed -s pseudo.bed -r read.fq -o outprefix
Input: compute the expected fpkm from reads
Output: gene with its expected expression in fragment count, and FPKM
Author: Chelsea Ju
Function: Read in the gene information, and the read data. Return a file with expression of parents and pseudo genes in number of fragments and FPKM.
Date: 2013-08-14
"""


import sys, getopt, re, pysam, os, random
from gene import *

GENELIST = {}

def usage():
    options = [ ("-p parents.bed --parents = parents.bed", ""),
                ("-s pseudo.bed --pseudo = pseudo.bed", ""),
                ("-r reads.fq, --read = read file", ""),
                ("-o outprefix, --outprefix = prefix of output file name", "")
                
    ]
    
    print "Usage: python compute_expected_fpkm.py -p parent.bed -s pseudo.bed -r read.fq -o outprefix"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"



def export_fragment(out_fh, count):
    global GENELIST
    
    out_fh.write('name\tchr\tstrand\tstart\tend\tsize\tfragment\tfpkm\n')

    for name in sorted(GENELIST):
        gene = GENELIST[name]
        match = re.match('^(.*)_(.*)', gene.name)
        if(match):
            gname = match.group(1)
        else:
            gname = gene.name
        
        out_fh.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\n' %(gname, gene.chr, gene.strand, gene.start, gene.end, gene.size, gene.get_fragment(1), gene.get_fpkm(count,1)))


def count_reads(read_fh):
    global GENELIST

    total_count = 0
    for r in read_fh:
        m = re.search('^@(ENSP\d+_ENST\d+)', r)
        if(m != None):
            gname = m.group(1)
            gene_obj = GENELIST[gname]
            gene_obj.add_primary_fragment(1)
            total_count += 1
    
    return total_count

def build_genelist(geneinfo_fh):    
    global GENELIST
        
    for x in geneinfo_fh:
        gene = from_bed_to_gene(x)
        GENELIST[gene.name] = gene


def from_bed_to_gene(line):
    data = line.split("\t")
    g = Gene(data[3], data[0], data[5], data[1], data[2])
    g.set_size(data[10])
    return g
    

def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:s:r:o:", ["parents=", "pseudo=", "reads=", "outprefix=" ])
    
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
        elif o in ("-o", "--outprefix"):
            outfile = str(a)
        elif o in ("-r", "--reads"):
            readfile = str(a)
        else:
            assert False, "unhandled option"
    
    try:
        parent_bed != None
        pseudo_bed != None
        readfile != None
        outfile != None        
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)
        
    outfile = outfile + "_expected.fpkm"

    parent_fh = open(parent_bed, 'rb')
    pseudo_fh = open(pseudo_bed, 'rb')
    read_fh = open(readfile, 'rb')
    expected_fh = open(outfile, 'w')

 
    # build gene information from bed files
    build_genelist(parent_fh)
    build_genelist(pseudo_fh)
    fragment_count = count_reads(read_fh)
    export_fragment(expected_fh, fragment_count)
    
    
    parent_fh.close()
    pseudo_fh.close()
    expected_fh.close()
    
    print "Write to file %s" %(outfile)
    print "DONE\n\n"



if __name__ == "__main__":
    main()
