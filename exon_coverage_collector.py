""" exon_coverage_counter
Usage: python exon_coverage_counter.py -c chromosome -d directory -m min_intron_size
Input:  -c chromosome name, also use to retrieve position file (ie chr22, position file = chr22_positionCount.txt)
        -d input/output directory 
        -m minimum intron size, default is 80 as suggested in "Junk in Your Genome: Intron Size and Distribution" 
        (http://sandwalk.blogspot.com/2008/02/junk-in-your-genome-intron-size-and.html)
Output: the header line is the chromosome name, follow by lines with {start \t end \t median readcount \n}
Function: 1. go through the position count file
          2. an exon region is considered as having consecutive positions with non-zero read count. 
          3. due to sequencing technology, a few nucleotides may not be covered, and results a small gap in the exon.
            thus, consecutive positions are relaxed to allow gaps smaller than the defined threshold (mininum intron size).
Date: 2014-01-01
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

"""
    Function: return the median value from the array
"""

def median(array):
    size = len(array)
    array = sorted(array)
    if(size % 2 != 0):
        return float(array[int(size/2)])
    else:
        return float((array[int(size/2)] + array[int(size/2)-1])/2)
    
    


"""
    Function: traverse the position count file (infile), where each line contains {position, readCount}
    if the differences between two positions are less than the given threshold (max_gap), 
    these positions are considered residing on the same exon
    
    output = {start, end, medianCount}

"""

def collect_exons(infile, max_gap):
    
    exons = []
    line = infile.readline()
    
    previous_position = 0
    start_position = 0
    counts = []
    
    for line in infile:
        line.rstrip()
        (position, count) = line.split("\t")
        position = int(position)
        count = int(count)
    
        # consecutive positions may have gap as they are not covered by reads
        # if the gap is less than the given threshold, it's less likely to be an intron
        # thus, a new exon started with a position that is far apart from the previous position (>= given threshold)
        if( (position - previous_position) >= max_gap):
            
            # an indicator to put the valid exon to array
            if(previous_position != 0):
                end_position = previous_position
            
                # exon count is calculated as the median value of all positions
                exons.append((start_position, end_position, median(counts))) 
                
            # reset the start position and counts
            start_position = position
            counts = []

        # accumulate count for the exon region
        counts.append(count)
        previous_position = position
        
    return exons

def export_exons(exons, outfh):
    
    for e in exons:
       outfh.write("%s\t%s\t%f\n" % (e[0], e[1], e[2]))
    
def main(parser):
    options = parser.parse_args()
    chromosome_name = options.chromosome
    dir = options.dir
    min_intron_size = options.intronSize
    
    # checking directory
    if(dir[-1] != "/"):
        dir += "/"
    
    # input file name
    input_file = dir + chromosome_name + "_positionCount.txt"
    output_file = dir + chromosome_name + "_exonCount.txt"
    
    input_fh = open(input_file, 'rb')
    output_fh = open(output_file, 'w')
    
    # collapse positions to be exons
    exons = collect_exons(input_fh, min_intron_size)    
    
    # output data
    export_exons(exons, output_fh)

    input_fh.close()
    output_fh.close()
    print "\nWriting Exon Counts to file : %s" %(output_file)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='exon_coverage_counter.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-m", "--min_intron", dest="intronSize", type=int, default=80, help="minimum intron size [OPTIONAL]")

    main(parser)
