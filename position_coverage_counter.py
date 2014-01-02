""" position_coverage_counter
Usage: python position_coverage_counter.py -c chromosome -i accepted_hit.bam -t threshold -d directory
Input: -c chromosome name, -d input/output directory -i bam file with mapped reads
Output: the header line is the chromosome name and size, follow by lines with {position \t readcount \n}
    note: the ouput only report the position with non-zero count   
Function: 1. convert the bam file to sam file
          2. traverse the sam file and select the mapping for a specific chromosome
          3. ignore the mapping whose quality falls below a certain threshold
          4. include multiple mapping (ie count all the positions the read maps to)
          5. note: this program also handles pair-end reads
Date: 2013-12-26
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

# chromosome size: from chr1 - chr23, X, Y, mitochondira
CHR_SIZE = [247199719, 242751149, 199446827, 191263063, 180837866,
            170896993, 158821424, 146274826, 140442298, 135374737,
            134452384, 132289534, 114127980, 106360585, 100338915,
            88822254, 78654742, 76117153, 63806651, 62435965,
            46944323, 49528953, 154913754, 57741652, 16569
            ]

"""
    convert the chromosome name to an index
    0 = chr1,
    22 = chrX,
    23 = ChrY,
    24 = Mitochondria
"""
def get_chr_index(name):
    
    name = name.lower()
    match_chr_number = re.match(r"\D*(\d+)", name)
    match_chr_X = re.match(r"\D*(x)", name)
    match_chr_Y = re.match(r"\D*(y)", name)
    match_chr_M = re.match(r"\D*(m)", name)
    
    if(match_chr_number):
        return int(match_chr_number.group(1)) - 1
    
    elif(match_chr_X):
        return 22

    elif(match_chr_Y):
        return 23

    elif(match_chr_M):
        return 24

    else:
        print ("unknown input %s", name)
        sys.exit(2)
    
"""
    reverse function of get_chr_index
    given an index, return the chr name    
"""
def get_chr_name(index):
    
    if(index < 22):
        return "chr" + str(index+1)
    elif(index == 22):
        return "chrX"
    elif(index == 23):
        return "chrY"
    elif(index == 24):
        return "chrM"
    else:
        print ("unknown input %s", str(index))
        sys.exit(2)
        
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
    The main function of the script:
    1. it fetches the reads of a given chromosome name from the sorted bam file
    2. for each read, it marks the positions the read covered
        if the read mapped to multiple places in the given chromosome, all positions are marked
        however, each position is only marked once by one read.
        for example, readA has three alignments (1-100, 1-40:101-160, 201-300),
                then position 1-160 and 201-300 are marked once (+1) by this read    

    note: because of the counting method, each read is considered as single read. 
        pair-end properties (such as "read mapped in proper pair") are not taken into consideration.
        in another word, accepted reads (all in accepted_hit.bam) are not discarded.
"""
def read_counter(sorted_input, chromosome_name, threshold):

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    counts_array = [0]*chr_size
    
    unique_read_name = {}

    samfile = pysam.Samfile( sorted_input, "rb" )

    for read in samfile.fetch(chromosome_name, 0, chr_size):
        read_name = read.qname
        
        # applies to those reads that passed the specified mapq threshold
        if(read.mapq > threshold):
            if(not unique_read_name.has_key(read_name)):
                unique_read_name[read_name] = []
    
            # collect a list of positions the reads aligned to, including the positions resulted in multiple alignments
            unique_read_name[read_name].extend(read.positions)

    # each read has a list of positions it mapped to
    # each position is updated only once by one read:
    #    the duplicated positions are removed
    #    update counts_array by the sorted list of positions
    for k,v in unique_read_name.items():
        sort_v = sorted(set(v))
        for sv in sort_v:
            counts_array[sv] +=1
    
    samfile.close() 
    return counts_array

"""
    write the read counts to file
"""
def output_array(array, chromosome_name, outfile):
    header = chromosome_name + "_size_" + str(len(array)) 
    total_count = sum(array)
    
    out_fh = open(outfile, 'w')
    out_fh.write("%s\n" % (header))
    
    if(total_count > 0):
        position = 0
        for a in array:
            position += 1 
            if(a > 0):
                out_fh.write("%d\t%d\n" %(position, a))

    out_fh.close()
    print "Writing Read Counts to file : %s" %(outfile)
    
def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    input = options.bamFile
    dir = options.dir
    threshold = options.threshold
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## start counting the reads
    sorted_input = sort_bam_file(input, dir)
    counts_array = read_counter(sorted_input, chromosome_name, threshold)
      
    ## output data
    outfile = dir + chromosome_name + "_positionCount.txt"
    output_array(counts_array, chromosome_name, outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='position_coverage_counter.py')
    parser.add_argument("-t", "--threshold", dest="threshold", type=int, default=0, help="mapQ cutoff [OPTIONAL]")
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-i", "--input", dest="bamFile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
