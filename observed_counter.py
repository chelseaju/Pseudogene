"""
Usage: python observed_counter.py -i accepted_hit.bam -c chromosome_name -d directory -t dataType
Input: -i bam file with mapped reads, -d the directory and prefix of output files, -c chromosome name -t dataType[genes or transcripts]
Output: A list of gene (or transcript) with the number of reads. {name \t number_of_read}
Function: Scan through the bam file and count the number of unique reads in a given region
       
Date: 2014-01-07
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse


"""
    using samtool to sort the bam file
    equivalent to samtool sort bam sort.bam
    
    after the sorting, create an index file on the sort.bam
"""      
def sort_bam_file(bam_file, dir):

    ## get the file prefix
    prefix = ""
    prefix_match = re.match(r"(.*).bam", bam_file)

    try:
        prefix = prefix_match.group(1)
    except:
        print "Existing: Invalid bam file -i %s" %(bam_file)
        sys.exit(2)
        

    # sort the bam file
    bam_input = dir + bam_file
    sort_bam = dir +  prefix + "_sorted"
    pysam.sort(bam_input, sort_bam)
    sort_bam = sort_bam + ".bam"
    
    # index the sort bam file
    pysam.index(sort_bam)

    print ""
    print "Writing Sorted Bam File : %s" %(sort_bam)
    print "Writing Index Sorted Bam File : %s.bai" %(sort_bam)
    
    return sort_bam


"""
    Function: iterate through each gene (from gene_file)
        and count the number of unique reads that covers this gene region
        unique gene is defined by the name of the read
"""
def observed_read_counter(sorted_bam, gene_file):

    count_array = []
    distribution_array = []
    gene_fh = open(gene_file, 'rb')
    sorted_bam_fh = pysam.Samfile(sorted_bam)
    
    for gene in gene_fh:
        (id, chr, start, end) = gene.split("\t")
        
        unique_reads = {}
        
        # fetch the read in that region
        for read in sorted_bam_fh.fetch(chr, int(start)-1, int(end)-1): # start and end are 1-based
            unique_reads[read.qname] = ""
        
        count_array.append((id, len(unique_reads)))
    sorted_bam_fh.close()
    gene_fh.close()

    return count_array

"""
    Write the observed count to file
"""
def export_array(array, outfile):
    out_fh = open(outfile, "w")

    for (id, count) in array:
        out_fh.write("%s\t%d\n" %(id, count))
    
    out_fh.close()
    print ""
    print "Writing Observed Count to File : %s" %(outfile)    


def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    input = options.bamFile
    dir = options.dir
    dataType = options.data_type
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
        
    ## check sorted bam file and sorted bam file index
    prefix = ""
    prefix_match = re.match(r"(.*).bam", input)

    try:
        prefix = prefix_match.group(1)
    except:
        print "Existing: Invalid bam file -i %s" %(bam_file)
        sys.exit(2)

    sorted_input = sort_bam_file(input, dir)

    ## input file
    input_gene_file = dir + "mapping/" + chromosome_name + "_" + dataType + ".txt"    
  
    ## start counting the read
    counts_array = observed_read_counter(sorted_input, input_gene_file)
      
    ## output data
    outfile = dir + "mapping/" + chromosome_name + "_" + dataType + "_observed.txt"
    export_array(counts_array, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='observed_counter.py')
    parser.add_argument("-i", "--input", dest="bamFile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)