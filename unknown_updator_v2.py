"""
Usage: python unknown_updator_v2.py -d directory -i file_to_be_update -t dataType 
Input: -d toppest directory -i the distribution matrix file -t dataType[genes or transcripts]
Output: the upated distribution matrix
Function: Read in the column names from [gene|transcript]_distribution.matrix. 
    If the name is "Unknown", replace the column name with the unified name
    If there are more than one columns with the same unified name, merge the columns

Author: Chelsea Ju
Date: 2014-05-09
"""

import sys, re, os, random, argparse, datetime

IDs = []
MERGE_COLUMNS = {}

"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))

"""
    Function : parse the unified ids
"""
def parse_ids(infile):

    fh = open(infile, 'rb')
    for line in fh:
        line = line.rstrip()
        (unknown, chromosome, start, end) = line.split("_")
        IDs.append((line, chromosome, int(start), int(end)))
        MERGE_COLUMNS[line] = []
    fh.close()

"""
    Function : flag the columns to be merged
"""
def flag_columns(infile):

    fh = open(infile, 'rb')
    colnames = fh.readline().rstrip().split("\t")
    fh.close

    sort_IDs =  sorted(IDs, key=lambda i:(i[1], i[2]))
    for i in xrange(0,len(colnames)):
        col_info = colnames[i].split("_")
        if(len(col_info) > 2):
            (unknown, chromosome, start, end) = col_info
            mapped_id = filter(lambda x: (chromosome == x[1] and int(start) >= x[2] and int(start) <= x[3]), sort_IDs)

            MERGE_COLUMNS[mapped_id[0][0]].append(i+1) # make it 1-base to be parsed by R

def export_merge_column(file):
    fh = open(file, 'w')
    for k in MERGE_COLUMNS.keys():
        if(len(MERGE_COLUMNS[k]) > 0):
            fh.write("%s\t%s\n" %(k, "\t".join([str(s) for s in MERGE_COLUMNS[k]])))
    fh.close()

    echo("File written to %s" %(file))

def main(parser):
    
    options = parser.parse_args()
    outdir = options.outdir
    dataType = options.data_type
    distribution_file = options.distribution

    if(outdir[-1] != "/"):
        outdir += "/"

    unify_file = outdir + dataType + "_unknown_ids.txt"

    merge_file = distribution_file + ".mergelist"

    # read in unified IDs
    parse_ids(unify_file)

    # read in the matrix
    flag_columns(distribution_file)
    export_merge_column(merge_file)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='unknown_updator_v2.py')
    parser.add_argument("-d", "--output_dir", dest="outdir", type=str, help="output directory", required = True)
    parser.add_argument("-i", "--distribution_matrix", dest="distribution", type=str, help="distribution matrix to be updated", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)