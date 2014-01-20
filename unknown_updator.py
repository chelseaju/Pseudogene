"""
Usage: python unknown_updator -f file_to_update -d lookup_table_directory
Input: -f file to update -d directory for the lookup table
Function: in-line find and replace. find the first column keyword of the lookup table in a given file (file to be update), 
    and replace it with the second column keyword

Author: Chelsea Ju
Date: 2014-01-19

"""

import sys, re, os, random, argparse

"""
    Function: construct the lookup table by fetch in data
"""
def build_lookup_table(file):

    fh = open(file, 'rb')
    table = fh.readlines()
    fh.close()    
    return table

"""
    Function: iterate through the lookup table, and replace the first keyword by the second keyword in file
        output the replacement into the same file
"""
def find_and_replace(lookup_table, infile):
    
    fh = open(infile, 'rb')
    data = "".join(fh.readlines())
    fh.close()
    
    for line in lookup_table:
        line = line.rstrip()
        (first, second) = line.split("\t")
        
        data = data.replace(first, second)
    
    fh = open(infile, 'w')
    fh.write(data)
    fh.close()
    
    print ""
    print "Unknown IDs Upated: %s" % (infile)
        
def main(parser):
    
    options = parser.parse_args()
    lookup_dir = options.lookup
    infile = options.input
    
    # build lookup table
    if(lookup_dir[-1] != "/"):
        lookup_dir += "/"
    lookup_file = lookup_dir + "lookup_table.txt"
    lookup_table = build_lookup_table(lookup_file)
    
    find_and_replace(lookup_table, infile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='unknown_updator.py')
    parser.add_argument("-d", "--lookupDir", dest="lookup", type=str, help="lookup table directory", required = True)
    parser.add_argument("-f", "--file", dest="input", type=str, help="file to be updated", required = True)

    main(parser)