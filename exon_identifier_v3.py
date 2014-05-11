""" exon_identifier_v3
Usage: python exon_identifier_v3.py -c chromosome -d directory -m min_intron_size
Input:  -c chromosome name
        -d input/output directory 
        -m minimum intron size, default is 80 as suggested in "Junk in Your Genome: Intron Size and Distribution" 
        (http://sandwalk.blogspot.com/2008/02/junk-in-your-genome-intron-size-and.html)
        -x maximum intron size, default is 4000 as suggested in the same article
Output: a list of region the gene mapped to, in bed format
    where "name" is the gene name (instead of region name)
Function: 1. fetch the bam file for a specific chromosome
          2. traverse each read to identify the covered positions
          3. collect the exon by collapsing the consecutive positions with non-zero coverage.
             however, due to sequencing technology, a few nucleotides may not be covered, and results a small gap in the exon.
             thus, consecutive positions are relaxed to allow gaps smaller than the defined threshold (mininum intron size).

Date: 2014-01-01
Update: 2014-04-24 change output format and traverse approach
Author: Chelsea Ju
"""

import sys, re, pysam, os, random, argparse, datetime

# chromosome size: from chr1 - chr23, X, Y, mitochondira
CHR_SIZE = [249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 81195210, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560, 59373566, 16569
            ]



"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


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
    match_chr_MT = re.match(r"\D*(mt)", name)
    
    if(match_chr_number):
        return int(match_chr_number.group(1)) - 1
    
    elif(match_chr_X):
        return 22

    elif(match_chr_Y):
        return 23

    elif(match_chr_MT):
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
        return str(index+1)
    elif(index == 22):
        return "X"
    elif(index == 23):
        return "Y"
    elif(index == 24):
        return "MT"
    else:
        print ("unknown input %s", str(index))
        sys.exit(2)

"""
    Function:
    1. create an array with the size of the chromosome
    2. traverse the bamfile with a given chromosome name, and obtain the coverage for each position
    3. if the position contains a non-zero coverage, mark the position in the array true
    output = [(name, chr, start, end)]
"""   
def mark_position(bam, chromosome_name, max_intron):

    position_array = []
    samfile = pysam.Samfile( bam, "rb" )

    chr_size = CHR_SIZE[get_chr_index(chromosome_name)]
    for read in samfile.fetch(chromosome_name, 0, chr_size):
        read_name = read.qname

        prefix_match = re.match(r"(.*?):.*", read_name)
        if(prefix_match):
            prefix = prefix_match.group(1)

            start = int(read.positions[0]) + 1
            end = int(read.positions[-1]) + 1
            index = -2
            while(end - start > max_intron): # start and end position shouldn't span longer than expected intron size
                end = int(read.positions[index]) + 1
                index -= 1

            # traverse from end, to see which one is longer
            if(index != -2):

                second_start = int(read.positions[0]) + 1
                second_end = int(read.positions[-1]) + 1
                index = 1
                while(second_end - second_start > max_intron):
                    second_start = int(read.positions[index]) + 1
                    index += 1
                if(second_end - second_start > end - start):
                    start = second_start
                    end = second_end

            position_array.append((prefix, start, end))             # bamfile is zero based

    return position_array


"""
    Function: traverse the sorted position array (by read_name and starting position)
        if two entries overlap, or have a gap less than the given threshold 
        collapse them into a single entry (defined as an exon)
    output = [(name, chr, start, end)]
"""
def collect_exon(position_array, max_gap):

    exon_array = []
    position_array = sorted(position_array, key=lambda x: (x[0], x[2]))    

    previous_gene = ""
    previous_end = 0
    previous_start = 0

    for p in position_array:
        (name, start, end) = p

        # continuous region
        if(name == previous_gene and (start - max_gap) <= previous_end):
            previous_end = max(previous_end, end)

        # output information
        else:
            if(previous_gene != ""):
                exon_array.append((previous_gene, previous_start, previous_end))

            previous_gene = name
            previous_start = start
            previous_end = end

    # output the last gene
    if(previous_gene != ""):
        exon_array.append((previous_gene, previous_start, previous_end))

    return exon_array

"""
    Function: export exon to file, in bed format
    note - the output position is 1-based
"""
def export_exon(exon_array, chromosome_name, outfile):

    outfh = open(outfile, "w")
    for e in exon_array:
        (name, start, end) = e
        outfh.write("%s\t%d\t%d\t%s\n" %(chromosome_name, start, end, name))

    outfh.close()
    echo("Writing Exon List to File : %s" %(outfile))


def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    dir = options.dir
    min_intron_size = options.min_intron
    max_intron_size = options.max_intron
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## start counting the reads
    sorted_input = dir + "accepted_hits_sorted.bam"
    position_array = mark_position(sorted_input, chromosome_name, max_intron_size)     

    ## collapse the positions into exons
    exon_array = collect_exon(position_array, min_intron_size)

    ## output data
    output_dir = dir + "mapping/"

    if(not os.path.exists(output_dir)):
        os.mkdir(output_dir)
        
    outfile = output_dir + chromosome_name + "_exons.bed"
    export_exon(exon_array, chromosome_name, outfile)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='exon_identifier_v3.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)
    parser.add_argument("-m", "--min_intron", dest="min_intron", type=int, default=80, help="minimum intron size [OPTIONAL]")
    parser.add_argument("-x", "--max_intron", dest="max_intron", type=int, default=4000, help="maximum intron size [OPTIONAL]")

    main(parser)
