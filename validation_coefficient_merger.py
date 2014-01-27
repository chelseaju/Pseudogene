"""
Usage: python validation_coefficient_merger.py -d output_directory -v validation_file -c coefficient_file -t dataType
Input: -d output and validation file path, -v validation file, -c path for coefficient file -t dataType[genes or transcripts]
Output: concatenation of validation data with distributed regions and their lasso coefficients 
    {Gene Name \t Expected \t Predicted \t Absolute Error \t (Prediction - Expectation)/Expectation \t Region Name \t Region Count \t Region Coefficient} 
Function: Combine the lasso coefficient with the distribution equation to the validation file.
 
Date: 2014-01-26
Author: Chelsea Ju
"""

import sys, re, os, random, argparse

COEFFICIENT  = {}
EQUATIONS = {}


"""
    Function: read in the lasso coefficient, and construct the table for it
"""
def build_coefficient_table(file):
    fh = open(file, 'rb')

    for line in fh:
        line = line.rstrip()
        if(re.match("^ENS|PGOHUM", line)):
            (id, coef) = line.split()
            COEFFICIENT[id] = float(coef)
    fh.close()


"""
    Function: read in the distribution equation, and construct the table for it
"""
def build_equation_table(file):
    fh = open(file, 'rb')

    for line in fh:
        line = line.rstrip()
        (LHS, RHS) = line.split(" = ")
        (count, id) = LHS.split(" * ")
        (g_id, t_id) = id.split("_")
        EQUATIONS[g_id] = RHS
    fh.close()    


"""
    Function: extend the validation data with distribution equation and lasso coefficient
     adding additional column-sets:
         1. disgributed region name
         2. distributed region with read count
         3. distributed region with lasso coefficient
"""
def extend_validation(file):
    
    data = []
    fh = open(file, 'rb')

    line = fh.readline()
    for line in fh:
        line = line.rstrip()
        (id, expected, predicted, absolute_error, relative_error) = line.split("\t")
        id = id.replace("\"", "")
        equation = EQUATIONS[id]
        
        regions = equation.split(" + ")
        output = [id, expected, predicted, absolute_error, relative_error]
        for r in regions:
            (count, id) = r.split(" * ")
            id = id.replace(" +", "")
            coef = COEFFICIENT[id]
            output.extend([str(id), str(count), str(coef)])
            
        data.append(output)
    fh.close() 
    return data

"""
    Function : export data
"""
def export_data(data, file):
    
    fh = open(file, 'w')
    fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Expected", "Predicted", "Absolute Error", "Relative Error", "Region ID", "Region Count", "Region Coeef"))

    for d in data:
        fh.write("\t".join(d))
        fh.write("\n")    
    fh.close()
    
    print ""
    print "Extending Validation Information: %s" %(file)
        
        
def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    validation_file = options.valid
    coefficient_file = options.coef
    dataType = options.data_type

    
    if(dir[-1] != "/"):
        dir += "/"

    ## file
    validation_file = dir + validation_file
    equation_file = dir + dataType + "_distribution.eqn"
    output = validation_file[0:-4] + "_extension.txt"
    
        
    ## function starts here
    build_coefficient_table(coefficient_file)
    build_equation_table(equation_file)
    merge_data = extend_validation(validation_file)
    export_data(merge_data, output)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='validation_coefficient_merger.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="path of output and validation file", required = True)
    parser.add_argument("-v", "--validation", dest="valid", type=str, help="validation file", required = True)
    parser.add_argument("-c", "--coefficient", dest="coef", type=str, help="path of coefficient files", required = True)
    parser.add_argument("-t", "--dataType", dest="data_type", type=str, help="genes or transcripts", required = True)

    main(parser)