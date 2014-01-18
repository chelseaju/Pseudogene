"""
Usage: python distribution_coefficient_merger.py -d distribution_file -c coefficient_file
Input: -d path for distribution file, -c path for coefficient file
Output: a file with three columns {expected gene \t distributed regions \t lasso coefficients} 
Function: Combine the lasso coefficient with the distribution equation. 
    The first column of the output is the LHS of distribution equation; the second column is the RHS of distribution equation.
    The third column is a copy of RHS in distribution equation, but replace the read count by lasso coefficient.

Date: 2014-01-17
Author: Chelsea Ju
"""

import sys, re, os, random, argparse

COEFFICIENT  = {}


"""
    Function: read in the lasso coefficient, and construct the table for it
"""
def build_coefficient_table(file):
    fh = open(file, 'rb')

    # remove headers
    fh.readline()
    fh.readline()
    fh.readline()

    for line in fh:
        line = line.rstrip()
        (id, coef) = line.split()
        COEFFICIENT[id] = float(coef)
    fh.close()


"""
    Function: merge the distribution equation with lasso coefficient
     output three columns:
         1. region with expected read count
         2. distributed region with read count
         3. distributed region with lasso coefficient
"""
def merge_coefficient(file):
    
    data = []
    fh = open(file, 'rb')

    for line in fh:
        line = line.rstrip()
        (LHS, RHS) = line.split(" = ")
        
        distribution = RHS.split(" + ")
        
        new_RHS = ""
        for d in distribution:
            d = d.replace(" +", "")
            (count, id) = d.split(" * ")
            coef = COEFFICIENT[id]
            new_RHS += str(coef) + " * " + id + " + "
        
        data.append((LHS, RHS, new_RHS))

    fh.close()    
    return data

"""
    Function : export data
"""
def export_data(data, file):
    
    fh = open(file, 'w')
    for (LHS, RHS, new_RHS) in data:
        fh.write("%s\t%s\t%s\n" %(LHS, RHS, new_RHS))

    fh.close()
    
    print ""
    print "Writing Distribution Coefficient: %s" %(file)
        
        
def main(parser):
    
    options = parser.parse_args()
    distribution_file = options.dis
    coefficient_file = options.coef
    
    ## parse output dir
    prefix = re.match("(.*)/(.*)_.*", distribution_file)
    if(prefix):
        out_dir = prefix.group(1)
        data_type = prefix.group(2)
    else:
        out_dir = "./"
        data_type = "genes"

    outfile = out_dir + "/" + data_type + "_distribution_coefficient.eqn"
    
    ## function starts here
    build_coefficient_table(coefficient_file)
    merge_data = merge_coefficient(distribution_file)
    export_data(merge_data, outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='distribution_coefficient_merger.py')
    parser.add_argument("-d", "--distribution", dest="dis", type=str, help="path of distribution files", required = True)
    parser.add_argument("-c", "--coefficient", dest="coef", type=str, help="path of coefficient files", required = True)

    main(parser)