""" assign_random_expression.py 
Usage: python assign_random_expression.py -i input_data -o output_prefix -e expression_range
Input: gene in bed format, expression = [v1, v2, v3] for specific values or [min-max] for a range
Output: replicate the gene with randomly assigned expression value, output in bed format 
Author: Chelsea Ju
Function: Read in the gene with bed format, and replicate the information with given expression value.
Date: 2013-08-26
Last Modify: 2013-09-24 [change the input]
"""

import sys, getopt, re, random

def usage():
    options = [ ("-i input_gene --input = input gene in bed format", ""),
                ("-o output_prefix, --output = output file prefix", ""),
                ("-e expression_range, --expression = expression range","")
     ]
    
    print "Usage: assign_random_expression.py -i input_data -o output_prefix -e expression_range"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


"""
    expression = [V1, V2, V3, ...] use "," to separate specific values
    expression = [MIN-MAX] use "-" to indicate the range
"""
def expand_expression(expression):

    values = expression.split("-")  # for range
    
    if(len(values) == 1):
        values = expression.split(",")  # for specific values
    else:
        values = range(int(values[0]), int(values[1])+1)
    
    return [int(i) for i in values]


def assign_expression(in_fh, bed_fh, expr_fh, expression_values):

    for line in in_fh:
        name = line.split("\t")[3]
        level = random.choice(expression_values)

        expr_fh.write("%s %d\n" %(name, level))
        bed_fh.write(line*level)
        

def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:e:", ["input=", "output=", "min=", "max="])
    
    except getopt.GetoptError, err:
    
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-i", "--input"):
            input = str(a)
        elif o in ("-o", "--output"):
            output = str(a)
        elif o in ("-e", "--expession"):
            expression = str(a)

    try:
        input != None
        output != None
        expression != None
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)


    out_bed = output + ".bed"
    out_expr = output + ".expression"


    in_fh = open(input, 'rb')
    bed_fh = open(out_bed, 'w')
    expression_fh = open(out_expr, 'w')
    
    expression_range = expand_expression(expression)
    assign_expression(in_fh, bed_fh, expression_fh, expression_range)
 
    print expression_range
    
    in_fh.close()
    bed_fh.close()
    expression_fh.close()
    
    print "Write to file %s" %(out_bed)
    print "Write to file %s" %(out_expr)

    print "DONE\n\n"



if __name__ == "__main__":
    main()

