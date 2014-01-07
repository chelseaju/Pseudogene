""" gene_identifier
Usage: python gene_identifier.py -c chromosome -d directory 
Input:  -c chromosome name
        -d input/output directory 
Output: lines of gene region or exon region with {name \t start \t end \n}
Function: 1. Map each exon position to the database to find a corresponding gene or pseudogene
          2. If the exon does not have a match, return this exon
          4. Return the mapped boundary of each unmarked exon

Date: 2014-01-06
Author: Chelsea Ju
"""
import sys, re, os, random, argparse, sqlite3
import sqlite3

"""
    make the chromosome name readable for database
"""
def convert_chromosome_name(name):
    
    name = name.lower()
    match_chr_number = re.match(r"\D*(\d+)", name)
    match_chr_X = re.match(r"\D*(x)", name)
    match_chr_Y = re.match(r"\D*(y)", name)
    match_chr_M = re.match(r"\D*(m)", name)

    if(match_chr_number):
        return match_chr_number.group(1)
    
    elif(match_chr_X):
        return match_chr_number.group(1)

    elif(match_chr_Y):
        return match_chr_number.group(1)

    elif(match_chr_M):
        return "MT"

    else:
        print ("unknown input %s", name)
        sys.exit(2)
  


def map_exon_to_gene(input_file, chr):
    
    mapped_regions = {}
    
    ## connect to database
    conn = sqlite3.connect('pseudogene.db')
    
    ## read file
    input_fh = open(input_file, 'rb')

    for line in input_fh:
        (start, end) = line.split('\t')
        
        gene_search = "SELECT t.transcript_id, t.transcript_start, t.transcript_end  "+
                      "FROM ensembl_exon as e INNER JOIN ensembl_mapping as m ON e.exon_id = m.exon_id INNER JOIN ensembl_transcript as t ON t.transcript_id = m.transcript_id " + 
                      "WHERE t.chromosome_name = '%s' AND e.exon_chr_end > %d AND e.exon_chr_start < %d" % (chr, start, end)
        
        pseudogene_search = "SELECT p.id, p.start_coordinate, p.stop_coordinate " + 
                            "FROM pseudogene as p " +
                            "WHERE p.chromosome_name = '%s' AND p.stop_coordinate > %d AND p.start_coordinate < %d " %(chr, start, end)

        c = conn.cursor()


    
    input_fh.close()
    conn.close()
    


def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## input file
    input_file = dir + chromosome_name + "_exons.txt"
    
    chr = convert_chromosome_name(chromosome_name)
    mapped_regions = map_exon_to_gene(input_file,chr)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='gene_identifier.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
