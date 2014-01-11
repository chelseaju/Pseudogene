"""
Usage: python observed_separator.py -i accepted_hit.bam -c chromosome_name -d directory
Input: -i bam file with mapped reads, -d the directory and prefix of output files, -c chromosome name
Output: a file with three columns {original_gene \t mapped_region \t read_count_of_mapped_region} 
Function: Similar to observed_counter.py, it gathers all the reads for newly identified gene regions,
    and splits them according to origin of the reads. The origin is identified by the prefix of read name.
    Each line in the output file corresponds to one origin and one mapped region, and the same origin may have multiple mapped regions. 

Date: 2014-01-08
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
        and separate reads that cover this gene region based on the origin of the reads
"""
def observed_read_separator(sorted_bam, gene_file):

    distribution_array = []
    gene_fh = open(gene_file, 'rb')
    sorted_bam_fh = pysam.Samfile(sorted_bam)
    
    for gene in gene_fh:
        (id, chr, start, end) = gene.split("\t")
        
        unique_origin = {}
        
        # fetch the read in that region
        for read in sorted_bam_fh.fetch(chr, int(start)-1, int(end)-1): # start and end are 1-based
            name = read.qname

            # extract the origin
            prefix_match = re.match(r"(.*?):.*", name)
            if(prefix_match):
                prefix = prefix_match.group(1)
                if(unique_origin.has_key(prefix)):
                    tmp_hash = unique_origin[prefix]
                    tmp_hash[name] = ""
                else:
                    unique_origin[prefix] = {}

        # copy the information from unique_origin to distribution_array           
        for k, v in unique_origin.items():
            distribution_array.append((k, id, len(v)))
        
    sorted_bam_fh.close()
    gene_fh.close()

    return distribution_array

"""
    Write the data to file
"""
def export_array(array, outfile):
    out_fh = open(outfile, "w")

    for (origin, mapped, count) in array:
        out_fh.write("%s\t%s\t%d\n" %(origin, mapped, count))
    
    out_fh.close()
    print ""
    print "Writing Read Distribution to File : %s" %(outfile)    


def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    input = options.bamFile
    dir = options.dir
    
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
    input_gene_file = dir + "mapping/" + chromosome_name + "_genes.txt"
  
    ## start counting the read
    distribution = observed_read_separator(sorted_input, input_gene_file)
          
    ## output data
    outfile = dir + "mapping/" + chromosome_name + "_distribution.txt"
    export_array(distribution, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='observed_separator.py')
    parser.add_argument("-i", "--input", dest="bamFile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)

    main(parser)