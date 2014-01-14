"""
Usage: python distribution_equation.py -d directory
Input: -d the working directory 
Output: distribution.eqn
Function: collects all the data generated by observed_separator, and formulate the distribution using equation
    for example: 108 * ORIGIN_1 = 100 * REGION_1 + 8 * REGION_2

Author: Chelsea Ju
Date: 2014-01-08
"""

import sys, re, os, random, argparse, glob

# Each item contains an array of information. 
# The first element is the expected count, and the following elements are tuples of (count, mapped_region)
EQUATION_HASH = {}

"""
    Function: initiate the equation hash by building the LHS of the equation
        ie     X * ORIGIN = RHS
"""
def equation_LHS(input):
    
    input_fh = open(input, 'rb')
    
    for line in input_fh:
        line = line.strip()
        (id, count) = line.split("\t")
        EQUATION_HASH[id] = [int(count)]
    
    input_fh.close()


"""
    Function: build the RHS of the equation
        ie     LHS = X * REGION1 + Y * REGION2 + X * REGION3
"""
def equation_RHS(input):
    
    input_fh = open(input, 'rb')
    for line in input_fh:
        line = line.strip()
        (origin, mapped, count) = line.split("\t")
        EQUATION_HASH[origin].append((int(count), mapped))
    input_fh.close()

"""
    Function: write data to file
"""
def export_data(output):

    output_fh = open(output, 'w')
    
    for k in sorted(EQUATION_HASH.keys()):
        v = EQUATION_HASH[k]
        # output region with mapped reads
        if(len(v) > 0):
            output_fh.write("%d * %s = " %(v[0], k)) # LHS            
            for (count, name) in v[1:]:
                output_fh.write("%d * %s + " %(count, name))            
            output_fh.write("\n")
    
    output_fh.close()
    print ""
    print "Writing Distribution Equation: %s" %(output)
    
def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    # file for expected count        
    expected_count_file = dir + "expected_read_count.txt"
    equation_LHS(expected_count_file)
    
    # data directory
    mapping_dir = dir + "mapping/"
    filenames = glob.glob(mapping_dir +"*_distribution.txt") ## collect all the files generated from observed_separator.py

    for f in filenames:
        equation_RHS(f)

    # outfile
    outfile = dir + "distribution.eqn"
    export_data(outfile)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='distribution_equation.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)

    main(parser)