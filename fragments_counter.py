""" fragment_counter.py 
Usage: python fragment_counter.py -p parent.bed -s pseudo.bed -e accepted_hit.bam -o outprefix -d [0 = by_gene| 1 = by_region] -c [0 = ignore multiple mapped reads | 1 = randomly select one | 2 = over-count multiple reads ]
Input: gene information in bed format, mapped read in bam format. -d indicates weather the reads are mapped to gene or by "transcription locus" (ie region with overlapped genes).
       -c indicates the counting method. 0 to ignore all the reads mapped to multiple places, 1 randomly select a region, 2 counts all the places
Output: a n(n+1) matrix that described the fragments and fpkm distribution between each gene and others, where n = |parent + pseudogene|
Author: Chelsea Ju
Function: Read in the gene information, and the expression data. Compute the distribution of mapped fragment of each transcript/locus. 
Date: 2013-07-31
Last Modify: 2013-08-14  (modify from expression_parser_from_bam.py)
"""

import sys, getopt, re, pysam, os, random
from gene import *

GENELIST = []  ## also record the expected number of fragments
DISTRIBUTION = []
GENELIST_INDEX = dict()
GENELIST_MEMBER = dict()
EXPRESSION = dict()

def usage():
    options = [ ("-p parents.bed --parents = parents.bed", ""),
                ("-s pseudo.bed --pseudo = pseudo.bed", ""),
                ("-e accepted_hit.bam, --expression = bam file", ""),
                ("-o outprefix, --outprefix = prefix of output file name", ""),
                ("-d gene_definition, --definition = [0 = by_gene (default)| 1 = by_region]", ""),
                ("-c count_method, --count = [0 = ignore multiple mapping reads | 1 = randomly select a place | 2 = count all the places]", "")

    ]
    
    print "Usage: python fragment_counter.py -p parent.bed -s pseudo.bed -e accepted_hit.bam -o outprefix -d definition -c count_method"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"

def export_matrix(fpkm_fh, fragment_fh, count):
    global GENELIST
    global GENELIST_INDEX
    global DISTRIBUTION

    title  = "\t".join([x.name for x in GENELIST])
    fpkm_fh.write("\t%s\t%s\n" % (title, "Others"))
    fragment_fh.write("\t%s\t%s\n" % (title, "Others"))

    size = len(GENELIST)
    for r_index in range(len(DISTRIBUTION)):
        current_gene = GENELIST[r_index]
        fragment_fh.write("%s\t%d\t%s\n" %(current_gene.name, current_gene.get_fragment(1), "\t".join(str(x) for x in DISTRIBUTION[r_index])))
        
        fpkm_array = [0]*(size+1)
        for c_index in range(size):
            f = DISTRIBUTION[r_index][c_index]
            if(f > 0):
                primary_gene = GENELIST[c_index]
                primary_gene.set_secondary_fragment(f)
                fpkm_array[c_index] = primary_gene.get_fpkm(count,2)
        fpkm_array[c_index] = DISTRIBUTION[r_index][c_index]
#        print "%s\t%s\t%s\n" %(current_gene.name, current_gene.get_fpkm(count,1), "\t".join(str(x) for x in fpkm_array))
	fpkm_fh.write("%s\t%f\t%s\n" %(current_gene.name, current_gene.get_fpkm(count,1), "\t".join(str(x) for x in fpkm_array)))
        


"""
    Mapping type = 
        0: do not count the reads that mapped to multiple places
        1: randomly select a place for the read
        2: over-count the reads that mapped to multiple places [default]

    Region definition = 
        0: by transcript
        1: by locus
"""
def map_fragments(mapping_type, region_definition):
    global EXPRESSION
    global GENELIST
    global GENELIST_INDEX
    global DISTRIBUTION
    global GENELIST_MEMBER
      
    print "\nMapping Fragments: mapping type = %d....." %(mapping_type)
    total_fragment = 0
    
    for(read_name, mapped_list) in EXPRESSION.items():
        ## compute expected number of fragments for each transcript
        match = re.match("(ENSP\d+_ENST\d+)", read_name)
        rname = match.group(1)

        if(region_definition == 0):        
            index = GENELIST_INDEX[rname]
            GENELIST[index].add_primary_fragment(1)
        else:
            index = GENELIST_MEMBER[rname]            
            GENELIST[index].add_primary_fragment(1)
                
        to_be_extract = []
        if(mapping_type == 0 and len(mapped_list) > 1):  ## skip multiple reads
            pass
        elif(mapping_type == 1): ## random select one
            to_be_extract = [random.choice(mapped_list)]
        else:
            to_be_extract = mapped_list

        all_mapped_gene = []
        for (chr, start, end, quality) in to_be_extract:
            total_fragment += 1
            mapped_gene = filter(lambda gene: gene.start <= start and gene.end >= end and gene.chr == chr, GENELIST)
            if(len(mapped_gene) == 0):
                all_mapped_gene += ["Others"]
            else:
                all_mapped_gene += mapped_gene
       
        # count the gene once only from mapped read           
        for g in set(all_mapped_gene):
            if(g == "Others"):
                DISTRIBUTION[index][-1] =  DISTRIBUTION[index][-1] + 1
            else:
                i = GENELIST_INDEX[g.name]
                DISTRIBUTION[index][i] = DISTRIBUTION[index][i] + 1

    return total_fragment
    

## iterate through the bam file, and discard the mapping with low quality
def build_expression(bam_file):
    global EXPRESSION
    
    # sort the bam file based on query name
    sort_bam = bam_file+".sorted"
    pysam.sort("-n",bam_file, sort_bam)

    print "\nWriting Sorted bam file : %s.bam" %(sort_bam)

    ## iterate through sorted bam file
    bam_fh = pysam.Samfile(sort_bam+".bam")
    count = 0
    for read in bam_fh:
        count = count + 1
        chr = bam_fh.getrname(read.tid)
        strand = "."
        name = read.qname
        
        # paired-end read counts as one fragment
        if(read.is_proper_pair):
            mate = bam_fh.next()
            read1 = read
            if(read.is_read2):
                read1 = mate
            
            if(read1.is_reverse):
                start = read1.mpos
            else:
                start = read1.pos
                
            end = start + abs(read1.isize)
            
            quality = min(read.mapq, mate.mapq)

        else:
            start = read.pos
            end = read.aend
            quality = read.mapq
 
        if(EXPRESSION.has_key(name)):
            (exist_chr, exist_start, exist_end, exist_quality) = EXPRESSION[name][0]
            if(quality > exist_quality):
                EXPRESSION[name] = [(chr, start, end, quality)]
            elif(quality == exist_quality):
                EXPRESSION[name].append((chr, start, end, quality))
        else:
            EXPRESSION[name] = [(chr, start, end, quality)]
    bam_fh.close()
    

def export_genelist_member(locus_fh):
    global GENELIST_MEMBER
    global GENELIST
    
    for g in GENELIST:
        index = GENELIST_INDEX[g.name]
        transcripts = filter(lambda (member_id, member_index): member_index== index, GENELIST_MEMBER.items())
        locus_fh.write("%s\t%s\n" %(g.name, ",".join([x for (x,y) in transcripts])))
    

    
def build_matrix(n,m):
    global DISTRIBUTION
    for i in range(n):
        DISTRIBUTION.append([0]*(m))

def build_genelist_by_locus(geneinfo_fh):
    global GENELIST
    global GENELIST_INDEX
    global GENELIST_MEMBER

    for x in geneinfo_fh:
        gene = from_bed_to_gene(x)
        
        # search for existing region
        exist_region = filter(lambda region: gene.chr == region.chr and ((gene.end <= region.end and gene.end >= region.start) or (gene.start <= region.end and gene.start >= region.start) or (gene.start < region.start and gene.end > region.end)), GENELIST)
        if(len(exist_region) > 0):
            region = exist_region[0]
            start = min(gene.start, region.start)
            end = max(gene.end, region.end)
        
            region.set_start(start)
            region.set_end(end)
            region.set_size(str(end - start + 1))
            
            GENELIST_MEMBER[gene.name] = GENELIST_INDEX[region.name]
        
        else:
            name = "R" + str(len(GENELIST) + 1)
            GENELIST_INDEX[name] = len(GENELIST)
            GENELIST_MEMBER[gene.name] = len(GENELIST)
            gene.set_name(name)
            gene.set_size(str(gene.end - gene.start + 1))
            GENELIST.append(gene)
     
    return len(GENELIST)
    
def build_genelist_by_transcript(geneinfo_fh):    
    global GENELIST
    global GENELIST_INDEX
    current_count = len(GENELIST)
    count = 0
    for x in geneinfo_fh:
        gene = from_bed_to_gene(x)
        GENELIST.append(gene)        
        GENELIST_INDEX[gene.name] = count + current_count
        count += 1
    return count
        
def from_bed_to_gene(line):
    data = line.split("\t")
    g = Gene(data[3], data[0], data[5], data[1], data[2])
    g.set_size(data[10])
    return g
    

def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:s:e:o:d:c:", ["parents=", "pseudo=", "expression=", "outprefix=", "definition="])
    
    except getopt.GetoptError, err:
    
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    definition = 0
    for o, a in opts:
        if o in ("-p", "--parents"):
            parent_bed = str(a)
        elif o in ("-s", "--pseudo"):
            pseudo_bed = str(a)
        elif o in ("-o", "--outprefix"):
            outfile = str(a)
        elif o in ("-e", "--expression"):
            expression = str(a)
        elif o in ("-d", "--definition"):
            definition = int(a)
        elif o in ("-c", "--count_method"):
            count_method = int(a)

        else:
            assert False, "unhandled option"
    
    try:
        parent_bed != None
        pseudo_bed != None
        expression != None
        outfile != None
        count_method != None
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)
    
    fpkm_outfile = outfile + ".fpkm"
    fragment_outfile = outfile + ".fragment"

    parent_fh = open(parent_bed, 'rb')
    pseudo_fh = open(pseudo_bed, 'rb')
    fpkm_fh = open(fpkm_outfile, 'w')
    fragment_fh = open(fragment_outfile, 'w')


    build_expression(expression)
 
    # build gene information from bed files
    if(definition == 0):
        parent_size = build_genelist_by_transcript(parent_fh)
        pseudo_size = build_genelist_by_transcript(pseudo_fh)

        # analysis + data output: reset the distribution matrix for each analysis
        build_matrix(parent_size, parent_size+pseudo_size+1) # build n*(n+1) matrix
        fragment_count = map_fragments(count_method, definition)
        export_matrix(fpkm_fh, fragment_fh, fragment_count)

    else:
        total_size = build_genelist_by_locus(parent_fh)
        total_size = build_genelist_by_locus(pseudo_fh)
        
        # analysis + data output: reset the distribution matrix for each analysis
        build_matrix(total_size, total_size+1) # build n*(n+1) matrix
        fragment_count = map_fragments(count_method, definition)
        export_matrix(fpkm_fh, fragment_fh, fragment_count)

        # output locus information
        locus_information = outfile + "_locus.info"
        locus_fh = open(locus_information, 'w')
        export_genelist_member(locus_fh)
        locus_fh.close()
        print "Write to file %s" %(locus_information)
    
    parent_fh.close()
    pseudo_fh.close()
    fpkm_fh.close()
    fragment_fh.close()
 
    
    print "Write to file %s" %(fpkm_outfile)
    print "Write to file %s" %(fragment_outfile)
    print "DONE\n\n"



if __name__ == "__main__":
    main()
