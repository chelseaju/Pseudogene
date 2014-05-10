"""
Usage: python unknown_unifier_v2.py -d list_of_directories -t dataType 
Input: -d list of working directory -t dataType[genes or transcripts]
Output: a list of unique unknown identifiers (unknown_ids.txt) 
Function: Read in the column names from [gene|transcript]_distribution.matrix. 
    The column ids are generated from [gene|transcript]_identifier.py. 
    The script iterates the ids from different analysis, and combine the IDs if they overlap.
Author: Chelsea Ju
Date: 2014-01-19
Update: 2014-03-05 change to infile editing
Update: 2014-05-09 change to output the unique identifiers
"""

import sys, re, os, random, argparse, datetime
from numpy  import *

IDs = {}

"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


"""
    Function: read in IDs from files
"""
def parse_id(file):
    fh = open(file, "rb")
    line = fh.readline()
    line = line.rstrip()
    for name in line.split("\t"):
        unknown_match = re.match("Unknown_", name)
        if(unknown_match):
            IDs[name] = ""
    fh.close()

"""
    Function: combine IDs within 80 bp apart
"""
def merge_ids():

    merger = []
    
    current_chr = "0"
    current_start = 0
    current_end = 0


    for e in sorted(IDs.keys()):
        (unknown, chr, start, end) = e.split("_")
        start = int(start)
        end = int(end)
        
        ## search overlap
        if(current_chr == chr and start >= current_start and start <= current_end):
            current_start = min(current_start, start)
            current_end = max(current_end, end)

        ## not overlap, output current chr
        else:
            if(current_chr != "0"):
                new_name = "Unknown_" + current_chr + "_" + str(current_start) + "_" + str(current_end)
                merger.append(new_name)

            current_chr = chr
            current_start = start
            current_end = end

    # export last one
    if(current_chr != "0"):
        new_name = "Unknown_" + current_chr + "_" + str(current_start) + "_" + str(current_end)
        merger.append(new_name)

    return merger

def export(ids, outfile):
    fh = open(outfile, 'w')
    for i in ids:
        fh.write("%s\n" %(i))
    fh.close()

    echo("File written to %s" %(outfile))

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    outdir = options.outdir
    dataType = options.data_type

    if(outdir[-1] != "/"):
        outdir += "/"
    outfile = outdir + dataType + "_unknown_ids.txt"

    # read in IDs
    for i in xrange(0,len(dir)):
        d = dir[i]
        if(d[-1] != "/"):
            d+= "/"
        
        dir[i] = d + dataType + "_distribution.matrix"
        parse_id(dir[i])

    # merge ids
    merger = merge_ids()

    export(merger, outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='unknown_unifier_v2.py')
    parser.add_argument("-o", "--output_dir", dest="outdir", type=str, help="output directory", required = True)
    parser.add_argument("-d", "--directory", nargs='+', dest="dir", type=str, help="list_of_directories", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)