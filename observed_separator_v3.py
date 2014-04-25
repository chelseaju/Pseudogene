"""
Usage: python observed_separator_v3.py -c chromosome_name -d directory -t dataType
Input: -d the directory and prefix of output files, -c chromosome name -t dataType[genes or transcripts]
Output: a file with three columns {original_gene \t mapped_region \t read_count_of_mapped_region} 
Function: Similar to observed_counter.py, it gathers all the reads for newly identified gene regions,
    and splits them according to origin of the reads. The origin is identified by the prefix of read name.
    Each line in the output file corresponds to one origin and one mapped region, and the same origin may have multiple mapped regions. 

Date: 2014-01-08
Last Update: 2014-04-18
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

"""
    Function: iterate through the list of candidate genes and filer out the ambigious mapped genes
        Ambigious gene is defined as:
        - having too few reads covering the region (ie less than 300 bp) 
        - have low coverage if it is an overlapped gene
"""
def gene_filter(input):

    final_gene_list = []
    gene_fh = open(input, 'rb')

    gene_list = [(line.strip().split("\t")) for line in gene_fh]
#    gene_list = sorted(gene_list, key = lambda gene_list: (gene_list[1], gene_list[2]))

#    previous_start = 0
#    previous_end = 0
#    previous_chr = 0

    for gene in gene_list:
        (ids, chr, start, end, region) = gene
        id_list = ids.split("/")
        region_list = region.split("/")
        best_coverage = 0
        best_id = ""

        # select the best coverage if the region contains overlapped genes
        for i in xrange(0,len(id_list)):
            gene_id = id_list[i]
            (gene_start, gene_end) = region_list[i].split("-")

            coverage = (min(float(gene_end), float(end)) - max(float(gene_start), float(start))) / max(float(gene_end) - float(gene_start), float(end) - float(start))
            print gene_id, coverage

            if(coverage > best_coverage):
                best_coverage = coverage
                best_id = gene_id

        # remove gene covered with only a few reads
        if(int(end) - int(start) > 110):
            final_gene_list.append((best_id, chr, start, end, best_coverage))
    
    gene_fh.close()
    return final_gene_list

"""
    Function: iterate through each gene (from gene_list)
        and separate reads that cover this gene region based on the origin of the reads
"""
def observed_read_separator(sorted_bam, gene_list, chromosome_name):

    distribution_array = []
    sorted_bam_fh = pysam.Samfile(sorted_bam)
    
    for gene in gene_list:        
        (ids, chr, start, end, quality) = gene    
        unique_origin = {}

        # fetch the read in that region
        for read in sorted_bam_fh.fetch(chromosome_name, int(start), int(end)): # start and end are 0-based
            name = read.qname
            # extract the origin
            prefix_match = re.match(r"(.*?):.*", name)
            if(prefix_match):
                prefix = prefix_match.group(1)
                if(unique_origin.has_key(prefix)):
                    tmp_hash = unique_origin[prefix]
                    tmp_hash[name] = ""
                else:
                    tmp_hash = {}
                    tmp_hash[name] = ""
                    unique_origin[prefix] = tmp_hash
   

         # copy the information from unique_origin to distribution_array           
        for k, v in unique_origin.items():
            distribution_array.append((k, ids, len(v)))
#            isoform_names = ids.split("/")
#            isoform_quality = quality.split("/")

#            if(k in isoform_names):
#                distribution_array.append((k, k, len(v)))
#            else:
#                best_index = isoform_quality.index(max(isoform_quality))
#                distribution_array.append((k, isoform_names[best_index], len(v)))
         
    sorted_bam_fh.close()

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
    dir = options.dir
    dataType = options.data_type
   
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
        
    sorted_input = dir + "accepted_hits_sorted.bam"

    ## input file
    input_gene_file = dir + "mapping/" + chromosome_name + "_" + dataType + ".txt"    
  
    ## filter out genes with little read evidence 
    genes_to_count = gene_filter(input_gene_file)

    ## start counting the read
    distribution = observed_read_separator(sorted_input, genes_to_count, chromosome_name)
          
    ## output data
    outfile = dir + "mapping/" + chromosome_name + "_" + dataType + "_distribution.txt"
    export_array(distribution, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='observed_separator_v2.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)


    main(parser)