"""
Usage: python observed_separator_v7.py -c chromosome_name -d directory -t dataType
Input: -d the directory and prefix of output files, -c chromosome name -t dataType[genes or transcripts]
Output: a file with three columns {original_gene \t mapped_region \t read_count_of_mapped_region} 
Function: Similar to observed_counter.py, it gathers all the reads for newly identified gene regions,
    and splits them according to origin of the reads. The origin is identified by the prefix of read name.
    Each line in the output file corresponds to one origin and one mapped region, and the same origin may have multiple mapped regions. 

Date: 2014-01-08
Update: 2014-05-01 - change on line 55
Update: 2014-05-01 - change on line 59, add function to remove overlapped region
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse, datetime


"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))



"""
    Function: iterate through each gene (from gene_list)
        and separate reads that cover this gene region based on the origin of the reads
"""
def observed_read_separator(sorted_bam, gene_file, chromosome_name):

    distribution_array = []

    gene_fh = open(gene_file, 'rb')
    sorted_bam_fh = pysam.Samfile(sorted_bam)       

    for gene in gene_fh:
        (gene, mapped_region, chr, start, end ) = gene.rstrip().split("\t")
        unique_reads = {}

        # fetch the read from a region 
        # (note: a valid region should have enough evidence, we assume a region should be longer than 102 bp)
        if(int(end) - int(start) > 102):
            for read in sorted_bam_fh.fetch(chromosome_name, int(start)-1, int(end)-1): # start and end are 1-based, but bam is 0-based
                name = read.qname
                # extract the origin
                prefix_match = re.match(r"(.*?):.*", name)
                if(prefix_match):
                    prefix = prefix_match.group(1)

                    # count number of reads from gene mapped to specific region
                    if(prefix == gene):
                        ## valid count is considred as
                        # 1. the same read name is only counted once for this region
                        # 2. the read starting point has to fall within this region
                        if(not unique_reads.has_key(name) and read.positions[0] >= int(start) and read.positions[0] <= int(end)):
                            unique_reads[name] = 0
            
            distribution_array.append((gene, mapped_region, start, end, len(unique_reads)))
         
    gene_fh.close()
    sorted_bam_fh.close()

    return distribution_array

"""
    Write the data to file:
        First, remove the overlapped genes
"""
def export_array(array, outfile):

    # sort the array by gene, and start positions
    sort_array = sorted(array, key=lambda x: (x[0],x[2]))

    # write to output
    out_fh = open(outfile, "w")

    previous_origin = ""
    previous_mapped = ""
    previous_start = 0
    previous_count = 0

    for (origin, mapped, start, end, count) in sort_array:
        if(count > 2):
            # overlap case
            if(origin == previous_origin and int(start) >= previous_start and int(start) <= previous_end):

                # replace previous
                if(origin == mapped):
                    previous_origin = origin
                    previous_mapped = mapped
                    previous_start = int(start)
                    previous_end = int(end)
                    previous_count = int(count)

            else:
                if(previous_origin != ""):
                    out_fh.write("%s\t%s\t%d\n" %(previous_origin, previous_mapped, previous_count))

                previous_origin = origin
                previous_mapped = mapped
                previous_start = int(start)
                previous_end = int(end)
                previous_count = int(count)

    # output the last case
    if(previous_origin != ""):
        out_fh.write("%s\t%s\t%d\n" %(previous_origin, previous_mapped, previous_count))

    out_fh.close()

    echo("Writing Read Distribution to File : %s" %(outfile))


def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    dir = options.dir
    dataType = options.data_type
   
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
        
    sorted_input = dir + "accepted_hits_sorted.bam"

    ## input file
    input_gene_file = dir + "mapping/" + chromosome_name + "_" + dataType + ".txt"    
  
    ## start counting the read
    distribution = observed_read_separator(sorted_input, input_gene_file, chromosome_name)
          
    ## output data
    outfile = dir + "mapping/" + chromosome_name + "_" + dataType + "_distribution.txt"
    export_array(distribution, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='observed_separator_v7.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)


    main(parser)