""" exon_identifier
Usage: python exon_identifier.py -c chromosome -d directory -m min_intron_size
Input:  -c chromosome name
        -d input/output directory 
        -m minimum intron size, default is 80 as suggested in "Junk in Your Genome: Intron Size and Distribution" 
        (http://sandwalk.blogspot.com/2008/02/junk-in-your-genome-intron-size-and.html)
Output: lines of exons with {start \t end \n}
Function: 1. fetch the bam file for a specific chromosome
          2. traverse each position to identify the coverage
          3. collect the exon by collapsing the consecutive positions with non-zero coverage.
             however, due to sequencing technology, a few nucleotides may not be covered, and results a small gap in the exon.
             thus, consecutive positions are relaxed to allow gaps smaller than the defined threshold (mininum intron size).

Date: 2014-01-01
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse

# chromosome size: from chr1 - chr23, X, Y, mitochondira
CHR_SIZE = [249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 81195210, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560, 59373566, 16569
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
    Function:
    1. create an array with the size of the chromosome
    2. traverse the bamfile with a given chromosome name, and obtain the coverage for each position
    3. if the position contains a non-zero coverage, mark the position in the array true
    4. return teh array
"""   
def mark_position(bam, chromosome_name):

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    counts_array = [False]*chr_size
    
    samfile = pysam.Samfile( bam, "rb" )
    
    for read in samfile.fetch(chromosome_name, 0, chr_size):
        for pos in read.positions:
            counts_array[pos] = True
 
    return counts_array

"""
## Pileup approach, which ignores pair end reads that mapped to two different chromosomes
def pileup_position(bam, chromosome_name):

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    counts_array = [False]*chr_size
    
    samfile = pysam.Samfile( bam, "rb" )
    
    for pileupcolumn in samfile.pileup(chromosome_name, 17152953, 17153061):
        if(pileupcolumn.n > 0):
            counts_array[pileupcolumn.pos] = True
 
    return counts_array
""" 

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

            # consecutive positions may have gap as they are not covered by reads
            # if the gap is less than the given threshold, it's less likely to be an intron
            # thus, a new exon started with a position that is far apart from the previous position (>= given threshold)            
            if((i - pre_pos) >= max_gap ):

                # pre_pos !=0 serves as an indicator to put the valid exon to array
                # valid exon is considered to have at least 10bp
                if(pre_pos != 0 and (pre_pos - start_pos + 1) > 10):
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
    dir = options.dir
    min_intron_size = options.intronsize
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## start counting the reads
    sorted_input = dir + "accepted_hits_sorted.bam"
    position_array = mark_position(sorted_input, chromosome_name)     
#    position_array = pileup_position(sorted_input, chromosome_name)     


    ## collapse the positions into exons
    exon_array = collect_exon(position_array, min_intron_size)

    ## output data
    output_dir = dir + "mapping/"

    if(not os.path.exists(output_dir)):
        os.mkdir(output_dir)
        
    outfile = output_dir + chromosome_name + "_exons.txt"
    export_exon(exon_array, outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='exon_identifier.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-m", "--min_intron", dest="intronsize", type=int, default=80, help="minimum intron size [OPTIONAL]")

    main(parser)
