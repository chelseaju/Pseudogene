"""
Usage: python unknown_unifier -d list_of_directories -t dataType 
Input: -d list of working directory -t dataType[genes or transcripts]
Output: infile editing 
Function: Read in the column names from [gene|transcript]_distribution.matrix. 
    The column ids are generated from [gene|transcript]_identifier.py. 
    The script iterates the ids from different analysis, and combine the IDs if they are within 80 bp apart.
Author: Chelsea Ju
Date: 2014-01-19
Update: 2014-03-05 change to infile editing
"""

import sys, re, os, random, argparse
from numpy  import *

IDs = {}

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

    merger = {}
    
    current_chr = "0"
    current_start = 0
    current_end = 0

    new_name = ""
    members = []

    for e in sorted(IDs.keys()):
        (unknown, chr, start, end) = e.split("_")
        start = int(start)
        end = int(end)
        
        if(current_chr!= chr):
            ## output current information: only interested in the IDs that have merged
            if(len(members) > 1):
                for m in members:
                    merger[m] = new_name

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
                    for m in members:
                        merger[m] = new_name

                # reset name and members
                members = [e]
                new_name = "Unknown_" + chr + "_" + str(start) + "_" + str(end)
        ## reset current information            
        current_chr = chr
        current_start = start
        current_end = end

    # export last one
    if(len(members) > 1):
        for m in members:
            merger[m] = new_name

    return merger

"""
    Function: output data
"""
def infile_editing(data, file):

    fh = open(file, 'rb')
    matrix = [ line.rstrip().split('\t') for line in fh ]
    fh.close()

    column_to_be_update = {}
    column_to_be_remove = []


    # go through the matrix
    for i in xrange(0, len(matrix[0])):
        id = matrix[0][i]
        if(data.has_key(id)):
            id = data[id]
            matrix[0][i] = id

            if(column_to_be_update.has_key(id)):
                column_to_be_update[id].append(i)
            else:
                column_to_be_update[id] = [i]

    for k in column_to_be_update.keys():
        column_index = column_to_be_update[k]
        if len(column_index) > 1 :
            first_index = column_index[0]

            # iterate through all the indecies to be updated
            for j in column_index[1:]:
                for row in xrange(1, len(matrix)):
#                    print matrix[row][0], matrix[row][first_index], matrix[row][j], first_index, j  ## for senity check 
                    matrix[row][first_index] = str(max(int(matrix[row][first_index]), int(matrix[row][j])))
                    column_to_be_remove.append(j)

    # output value
    fh = open(file, 'w')
    for row_data in matrix:
        row_data = [x for y,x in enumerate(row_data) if y not in column_to_be_remove]
        fh.write("\t".join(row_data))
        fh.write("\n")
    fh.close()

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    outdir = options.outdir
    dataType = options.data_type

    # read in IDs
    for i in xrange(0,len(dir)):
        d = dir[i]
        if(d[-1] != "/"):
            d+= "/"
        
        dir[i] = d + dataType + "_distribution.matrix"
        parse_id(dir[i])

    # merge ids
    merger = merge_ids()

    for d in dir:
        infile_editing(merger, d)
        print "File Upate: %s" %(d)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='unknown_unifier.py')
    parser.add_argument("-o", "--output_dir", dest="outdir", type=str, help="output directory", required = True)
    parser.add_argument("-d", "--directory", nargs='+', dest="dir", type=str, help="list_of_directories", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)