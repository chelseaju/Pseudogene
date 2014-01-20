"""
Usage: python unknown_unifier -d list_of_directories -t dataType -o outdir
Input: -d list of working directory -t dataType[genes or transcripts] -o output directory
Output: look up table - lookup_table.txt
Function: Read in the column names from [gene|transcript]_distribution.matrix. 
    The column ids are generated from [gene|transcript]_identifier.py. 
    The script iterates the ids from different analysis, and combine the IDs if they are within 80 bp apart.
Author: Chelsea Ju
Date: 2014-01-19
"""

import sys, re, os, random, argparse

IDs = []

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
            IDs.append(name)
    fh.close()


"""
    Function: combine IDs within 80 bp apart
"""
def merge_ids():

    merger = []
    
    current_chr = "0"
    current_start = 0
    current_end = 0

    new_name = ""
    members = []

    for e in sorted(set(IDs)):
        (unknown, chr, start, end) = e.split("_")
        start = int(start)
        end = int(end)
        
        if(current_chr!= chr):
            ## output current information: only interested in the IDs that have merged
            if(len(members) > 1):
                merger.append((new_name, members))

            # reset name and members
            members = [e]
            new_name = "Unknown_" + chr + "_" + str(start) + "_" + str(end)
        else:           
            ## overlap region or within 80bp apart
            if((current_end > start and current_start < end) or (abs(start - current_end) <= 80)):
               (unknown, name_chr, name_start, name_end) = new_name.split("_")
               name_start = int(name_start)
               name_end = int(name_end)
               new_name = "Unknown_" + chr + "_" +str(min(name_start, start)) + "_" + str(max(name_end, end))
               members.append(e)
            ## a new region, output previous information if more than 2 IDs are merged
            else:
                if(len(members) > 1):
                    merger.append((new_name, members))

                # reset name and members
                members = [e]
                new_name = "Unknown_" + chr + "_" + str(start) + "_" + str(end)
                    
                
        ## reset current information            
        current_chr = chr
        current_start = start
        current_end = end

    # export last one
    if(len(members) > 1):
        merger.append((new_name, members))

    return merger

"""
    Function: output data
"""
def export_data(data, file):
    fh = open(file, 'w')
    for (name, array) in data:
        for id in array:
            fh.write("%s\t%s\n" %(id, name))
    
    fh.close()
    print ""
    print "Writing Lookup Table: %s" %(file)
    
def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    outdir = options.outdir
    dataType = options.data_type

    # read in IDs
    for d in dir:
        if(d[-1] != "/"):
            d+= "/"
        
        matrix_file = d + dataType + "_distribution.matrix"
        parse_id(matrix_file)

    # merge ids
    merger = merge_ids()

    # output file    
    if(outdir[-1] != "/"):
        outdir += "/"
    outfile = outdir + dataType + "_lookup.txt"
    export_data(merger, outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='unknown_unifier.py')
    parser.add_argument("-o", "--output_dir", dest="outdir", type=str, help="output directory", required = True)
    parser.add_argument("-d", "--directory", nargs='+', dest="dir", type=str, help="list_of_directories", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)