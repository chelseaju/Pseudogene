"""
Usage: python expected_counter.py -i accepted_hit.bam -d directory
Input: -i bam file with mapped reads, -d the directory and prefix of output files
Output: A list of gene with the number of reads. {gene_name \t number_of_read}
Function: Assuming read with the same prefix (gene name) comes from the same gene.
        The script iterates through the bam file and count expected number of read for each gene.
Date: 2013-12-29
Author: Chelsea Ju
Note: this script is adopted from fragments_counter
"""

import sys, re, pysam, os, random, argparse

"""
    using samtool to sort the bam file by read name
    equivalent to samtool sort -n bam sort.bam
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
    sort_bam = dir +  prefix + "_sortedByName"
    pysam.sort('-n', bam_input, sort_bam)
    sort_bam = sort_bam + ".bam"
    
    print ""
    print "Writing Sorted Bam File : %s" %(sort_bam)
     
    return sort_bam


"""
    iterate through the bam file, which is sorted by the read name
    the file may contain multiple alignments and pair end mates, so the read is only counted when seen the first time
"""
def expected_read_counter(sorted_input):
    
    expected_hash = {}
    previous_read = ""
    sorted_bam_fh = pysam.Samfile(sorted_input)

    for read in sorted_bam_fh:
        name = read.qname
        
        # the first time seeing this read  
        # (note: reads contain multiple alignment and pair-end mates, but each read is counted once)
        if(name != previous_read):
            prefix_match = re.match(r"(.*?):.*", name)
            if(prefix_match):
                prefix = prefix_match.group(1)
    
                # if the read name exists
                if(expected_hash.has_key(prefix)):
                    expected_hash[prefix] += 1
                else:
                    expected_hash[prefix] = 1
        
        # update previous_read
        previous_read = name
    
    sorted_bam_fh.close()
    return expected_hash
        
"""
    write the expected counts to file
"""
def output_array(hash, outfile):
    out_fh = open(outfile, "w")

    for k, v in hash.items():
        out_fh.write("%s\t%s\n" %(k,str(v)))
    
    out_fh.close()
    print "Writing Expected Count to File : %s" %(outfile)    

        

def main(parser):
    
    options = parser.parse_args()
    input = options.bamFile
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    ## start counting the readd
    sorted_input = sort_bam_file(input, dir)
    counts_hash = expected_read_counter(sorted_input)
      
    ## output data
    outfile = dir + "expected_read_count.txt"
    output_array(counts_hash, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='expected_counter.py')
    parser.add_argument("-i", "--input", dest="bamFile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)

    main(parser)
