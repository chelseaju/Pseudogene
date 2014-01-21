"""
Usage: python distribution_matrix.py -d directory -t dataType
Input: -d the working directory  -t dataType[genes or transcripts]
Output: genes_distribution.matrix or transcripts_distribution.matrix
Function: Reformat the data from distribution.eqn to matix

Author: Chelsea Ju
Date: 2014-01-09
"""
import sys, re, os, random, argparse, glob


ORIGIN_ID = {} # ID => Index
NEW_ID = {} 

"""
    Function: determine the size an ID of the matrix
        matrix size should be N x M, where N is the number of original regions, and M is the number of mapped
        usually M > N
        It builds ORIGIN_ID and NEW_ID, where keys are the unique ID for original regions and the mapped region respectively,
            values are the index of the key.
        For mapped region, the index refers to the dth position after the original regions (ie N + d)            
"""
def retrieve_ID(input):

    n = 0
    m = 0
    input_fh = open(input, 'rb')    
    for line in input_fh:
        line = line.rstrip()
        (lhs, rhs) = line.split(" =")
        
        # LHS
        (o_count, o_id) = lhs.split(" * ")
        (o_pid, o_tid) = o_id.split("_")
        ORIGIN_ID[o_pid] = n
        n += 1
        
        # RHS
        for element in rhs.split(" +"):
            # ignore the empty RHS
            if(len(element) > 0):
                (m_count, m_id) = element.split(" * ")
                if(not ORIGIN_ID.has_key(m_id) and not NEW_ID.has_key(m_id)):
                    NEW_ID[m_id] = m 
                    m += 1
    input_fh.close()

"""
    Function: go through the distribution.eqn again, and compute, output the matrix to file
"""
def construct_export_matrix(input, output):
    
    in_fh = open(input, 'rb')
    out_fh = open(output, 'w')
    
    n = len(ORIGIN_ID)
    m = len(NEW_ID)

    # print column name to file first
    for k in sorted(ORIGIN_ID, key=ORIGIN_ID.get):
        out_fh.write("%s\t" %(k))

    i = 0
    for k in sorted(NEW_ID, key = NEW_ID.get):
        i = i+1
        if(i < m):
            out_fh.write("%s\t" %(k))
        else:
            out_fh.write("%s\n" %(k))
    
    for line in in_fh:
        line = line.strip()
        (lhs, rhs) = line.split(" =")
        
        data = [0]*(n+m)

        # LHS
        (o_count, o_id) = lhs.split(" * ")
        (o_pid, o_tid) = o_id.split("_")

        # RHS
        for element in rhs.split(" +"):
            # ignore the empty RHS
            if(len(element) > 0):
                (m_count, m_id) = element.split(" * ")
                m_count = m_count.strip(" ")
                if(ORIGIN_ID.has_key(m_id)):
                    m_id_index = ORIGIN_ID[m_id]
                else:
                    m_id_index = n + NEW_ID[m_id]
                
                data[m_id_index] = m_count
        
        out_fh.write("%s\t" %(o_pid))
        out_fh.write("\t".join([str(d) for d in data]))
        out_fh.write("\n")
    
    in_fh.close()
    out_fh.close()

    print ""
    print "Writing Distribution Matrix: %s" %(output)
 

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    dataType = options.data_type
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    # distribution equation file
    input_file = dir + dataType + "_distribution.eqn"
    retrieve_ID(input_file)
    
    # output the matrix to file
    output_file = dir + dataType + "_distribution.matrix"
    construct_export_matrix(input_file, output_file)
    
    

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='distribution_matrix.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)
