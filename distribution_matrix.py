"""
Usage: python distribution_matrix.py -d directory -t dataType
Input: -d the working directory  -t dataType[genes or transcripts]
Output: genes_distribution.matrix or transcripts_distribution.matrix
Function: Reformat the data from distribution.eqn to matix

Author: Chelsea Ju
Date: 2014-01-09
"""
import sys, re, os, random, argparse, glob


IDs = {} # ID => Index

"""
    Function: determine the size an ID of the matrix
        matrix size should be N x M, where N is the number of original regions, and M is the number of mapped
        usually M > N
        It builds ID, where key is the unique ID for original regions and the mapped region,
            value is the index of the key.
"""
def retrieve_ID(input):

    index = 0
    
    # traverse the file twice, the first time to get the original IDs, second time to get the mapped region
    input_fh = open(input, 'rb')
    for line in input_fh:
        line = line.rstrip()
        (lhs, rhs) = line.split(" =")
        
        # LHS
        (o_count, o_id) = lhs.split(" * ")
        (o_pid, o_tid) = o_id.split("_")
        if(not IDs.has_key(o_pid)):
            IDs[o_pid] = index
            index += 1

    # go back to the beginning of the file
    input_fh.seek(0)   
    for line in input_fh:
        line = line.rstrip()
        (lhs, rhs) = line.split(" =")
                
        # RHS
        for element in rhs.split(" +"):
            # ignore the empty RHS
            if(len(element) > 0):
                (m_count, m_id) = element.split(" * ")
                if(not IDs.has_key(m_id)):
                    IDs[m_id] = index
                    index += 1
    input_fh.close()

"""
    Function: go through the distribution.eqn again, and compute, output the matrix to file
"""
def construct_export_matrix(input, output):
    
    in_fh = open(input, 'rb')
    out_fh = open(output, 'w')
        
    # print column name to file first
    for k in sorted(IDs, key=IDs.get):
        out_fh.write("%s\t" %(k))
    
    out_fh.write("\n")
    for line in in_fh:
        line = line.strip()
        (lhs, rhs) = line.split(" =")
        
        data = [0]*(len(IDs))

        # LHS
        (o_count, o_id) = lhs.split(" * ")
        (o_pid, o_tid) = o_id.split("_")

        # RHS
        for element in rhs.split(" +"):
            # ignore the empty RHS
            if(len(element) > 0):
                (m_count, m_id) = element.split(" * ")
                m_count = m_count.strip(" ")
                if(IDs.has_key(m_id)):
                    m_id_index = IDs[m_id]
                else:
                    print "ID %s not found " %(m_id)
                    exit()
                    
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
