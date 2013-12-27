""" coverage_counter
Usage: python coverage_counter.py -c chromosome -i accepted_hit.bam -t threshold -o outputPrefix
Input: -c chromosome name, -i bam file with mapped reads, -o the directory and prefix of output files
Output: the first line is the chromosome name and size, follow by lines with {position \t readcount \n}
Function: 1. convert the bam file to sam file
          2. traverse the sam file and select the mapping for a specific chromosome
          3. ignore the mapping whose quality falls below a certain threshold
          4. include multiple mapping (ie count all the positions the read maps to)
          5. note: this program can also handle pair-end reads
Date: 2013-12-26
Author: Chelsea Ju
"""

import sys, getopt, re, pysam, os, random

# chromosome size: from chr1 - chr23, X, Y, mitochondira
CHR_SIZE = [247199719, 242751149, 199446827, 191263063, 180837866,
            170896993, 158821424, 146274826, 140442298, 135374737,
            134452384, 132289534, 114127980, 106360585, 100338915,
            88822254, 78654742, 76117153, 63806651, 62435965,
            46944323, 49528953, 154913754, 57741652, 16569
            ]


def usage():
    options = [  ("-c chromosome --chromosome = chromosome_name", ""),
                ("-i accepted_hit.bam --input = accepted_hit.bam", ""),
                ("-t threshold --threshold = quality threshold", ""),
                ("-o outputPrefix --output = prefix of output file", "")                
    ]
    
    print "Usage: python coverage_counter.py -c chromosome -i accepted_hit.bam -t threshold -o outputPrefix"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"



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
        
def sort_bam_file(bam_file, outputPrefix):

    # sort the bam file based on query name
    sort_bam = bam_file + ".sorted"
    pysam.sort("-n",bam_file, sort_bam)
    sort_bam = bam_file + ".sorted" + ".bam"

    pysam.index(sort_bam)

    print "\nWriting Sorted bam file : %s" %(sort_bam)
    print "\nWriting Index Sorted bam file : %s.bai" %(sort_bam)
    
    return sort_bam

"""
    the main function of the script
    it fetches the reads of a given chromosome name from the sorted bam file
    for each read, it marks the positions the read covered
    
"""

def read_counter(sorted_input, chromosome_name, threshold):

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    counts_array = [0]*chr_size

    samfile = pysam.Samfile( sorted_input, "rb" )

    for read in samfile.fetch(chromosome_name, 0, 10):
        print read

    return counts_array

    
def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:i:t:o:", ["chromosome=", "input=", "threshold=", "outputPrefix="])
    
    except getopt.GetoptError, err:
    
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    threshould = 0
    for o, a in opts:
        if o in ("-c", "--chromosome"):
            chromosome_name = str(a)
        elif o in ("-i", "--input"):
            input = str(a)
        elif o in ("-o", "--outputPreix"):
            output = str(a)
        elif o in ("-t", "--threshold"):
            threshold = int(a)

        else:
            assert False, "unhandled option"
    
    try:
        chromosome_name != None
        input != None
        output != None
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)

    
    ## start counting the reads
    sorted_input = sort_bam_file(input, output)
    counts_array = read_counter(sorted_input, chromosome_name, threshold)
    
    ## output data
    outfile = output + chromosome_name + "_count.txt"
#    output_array(counts_array, chromosome_name, outfile)
    
    
    
    
   
    


if __name__ == "__main__":
    main()
