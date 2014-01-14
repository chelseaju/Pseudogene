""" ENSP2ENSG
Usage: python ENSP2ENSG.py -i input -o output -d directory 
Input:  -i input file name
        -o output file name
        -d input/output directory 
Output: any id with ENSP is converted to its corresponding gene name ENSG
Function: 1. Scan through the file to look for key word ENSPXXXXX
          2. Query the database to find the corresponding ENSGXXXX

Date: 2014-01-13
Author: Chelsea Ju
"""
import sys, re, os, random, argparse, sqlite3
import sqlite3

#DB = '/u/home/c/chelseaj/database/PseudogeneDB/pseudogene.db'
DB = '/home/chelseaju/Database/PseudogeneDB/pseudogene.db'


"""
    Function: scan through infile, and change the ENSP keyword to ENSG through retrieving the database
"""
def convert_ids(infile, outfile):
    
    ## connect to database
    conn = sqlite3.connect(DB)
    
    infh = open(infile, 'rb')
    outfh = open(outfile, 'w')
    
    for line in infh:
        # break up the line into ENSP token
        ensp = re.split(r"(ENSP\d+)", line)

        for e in ensp:            
            # check if each token has an ENSP ID
            ensp_id = re.match(r"(ENSP\d+)", e)
            if(ensp_id):
                query_id = ensp_id.group(1)
                
                query = "SELECT m.gene_id FROM ensembl_mapping m INNER JOIN ensembl_transcript t on m.transcript_id = t.transcript_id WHERE t.protein_id = '%s'" % (query_id)
                
                c = conn.cursor()      
                c.execute(query)
                
                ensg_id = c.fetchone()[0]
                line = line.replace(query_id, ensg_id)
    
        outfh.write(line)
    
    conn.close()
    infh.close()
    outfh.close()
    
    print ""
    print "Writing Conversion to File : %s" %(outfile)
            
            
def main(parser):
    
    options = parser.parse_args()
    input = options.input
    output = options.output
    dir = options.dir

    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    ## input file
    input_file = dir + input
    
    ## output file
    output_file = dir + output

    convert_ids(input_file, output_file)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='ENSP2ENSG.py')
    parser.add_argument("-i", "--input", dest="input", type=str, help="input file", required = True)
    parser.add_argument("-o", "--output", dest="output", type=str, help="output file", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
