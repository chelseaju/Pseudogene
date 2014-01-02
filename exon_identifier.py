""" exon_identifier
Usage: python exon_identifier.py -c chromosome -d directory -m min_intron_size -i bamfile
Input:  -i bamfile (ex accepted_hits.bam)
        -c chromosome name
        -d input/output directory 
        -m minimum intron size, default is 80 as suggested in "Junk in Your Genome: Intron Size and Distribution" 
        (http://sandwalk.blogspot.com/2008/02/junk-in-your-genome-intron-size-and.html)
Output: lines of exons with {start \t end \n}
Function: 1. sort and index the bam file
          2. fetch the bam file for a specific chromosome
          3. traverse each position to identify the coverage
          4. collect the exon by collapsing the consecutive positions with non-zero coverage.
             however, due to sequencing technology, a few nucleotides may not be covered, and results a small gap in the exon.
             thus, consecutive positions are relaxed to allow gaps smaller than the defined threshold (mininum intron size).

Date: 2014-01-01
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
    Function:
    1. create an array with the size of the chromosome
    2. traverse the bamfile with a given chromosome name, and obtain the coverage for each position
    3. if the position contains a non-zero coverage, mark the position in the array true
    4. return teh array
"""
def pileup_position(bam, chromosome_name):

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    counts_array = [False]*chr_size
    
    samfile = pysam.Samfile( bam, "rb" )
    
    for pileupcolumn in samfile.pileup(chromosome_name, 0, chr_size):
        if(pileupcolumn.n > 0):
            counts_array[pileupcolumn.pos] = True
 
    return counts_array
    
"""
    Function: traverse the position array. 
            if the position is marked true, the position is covered by read(s).
            exons is defined as a consecutive positions. 
            if the difference between two marked positions are less than the given threshold(max_gap) 
            these positions are considered residing on the same exon    
    output = [(start, end)]
"""
def collect_exon(position_array, max_gap):

    exons = []
    
    pre_pos = 0
    start_pos = 0
    
    for i in xrange(len(position_array)):
        if(position_array[i] == True):
            if((i - pre_pos) >= max_gap ):
                if(pre_pos != 0):
                    exons.append((start_pos, pre_pos))
                
                start_pos = i
            pre_pos = i
        
    return exons

"""
    Function: export exon to file, in the format of
                {start_pos \t end_pos}
    note - the output position is 1-based
"""
def export_exon(exon_array, outfile):

    outfh = open(outfile, "w")
    for e in exon_array:
        outfh.write("%d\t%d\n" %(e[0]+1, e[1]+1)) # the position in e is 0-based

    outfh.close()
    print "Writing Exon List to File : %s" %(outfile)

def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    input = options.bamfile
    dir = options.dir
    min_intron_size = options.intronsize
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## start counting the reads
    sorted_input = sort_bam_file(input, dir)
    position_array = pileup_position(sorted_input, chromosome_name)     


    ## collapse the positions into exons
    exon_array = collect_exon(position_array, min_intron_size)

    ## output data
    outfile = dir + chromosome_name + "_exons.txt"
    export_exon(exon_array, outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='exon_identifier.py')
    parser.add_argument("-i", "--input", dest="bamfile", type=str, help="bam file name, ex accepted_hit.bam", required = True)
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-m", "--min_intron", dest="intronsize", type=int, default=80, help="minimum intron size [OPTIONAL]")

    main(parser)
