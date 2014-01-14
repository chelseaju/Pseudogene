"""
Usage: python bamfile_sorter.py -i accepted_hit.bam -d directory
Input: -i bam file with mapped reads, -d the directory and prefix of output files
Output: accepted_hit_sorted.bam accepted_hit_sorted.bam.bai accepted_hit_sortedByName.bam accepted_hit_sortedByName.bam.bai
Function: calling samtool to sort the bam file, and index the sorted bam file
    two different sorts are performed: one by read name, and one by position
Date: 2014-01-13
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

"""
    using samtool to sort the bam file
    equivalent to samtool sort bam sort.bam
    
    after the sorting, create an index file on the sort.bam
"""      
def sort_by_position(bam_file, dir):

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
    using samtool to sort the bam file by read name
    equivalent to samtool sort -n bam sort.bam
    index the sorted bam file
"""
def sort_by_name(bam_file, dir):

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


     

def main(parser):
    
    options = parser.parse_args()
    input = options.bamFile
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    ## input file
    sort_by_position(input, dir)
    sort_by_name(input, dir)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='bamfile_sorter.py')
    parser.add_argument("-i", "--input", dest="bamFile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)

    main(parser)

