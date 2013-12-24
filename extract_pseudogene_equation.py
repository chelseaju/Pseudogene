""" extract_pseudogene_equation.py 
Usage: python extract_pseudogene_equation.py -i input_eqn -o output
Input: distribution matrix equation
Output: a tab-delimited file contains [name, expected_value, itself_count, pseudogene_count, other_count]
Author: Chelsea Ju
Function: Read in the distribution equation, extract the equations containing pseudogenes
Date: 2013-08-19
"""

import sys, getopt, re

def usage():
    options = [ ("-i input_eqn --input = input matrix", ""),
                ("-o output, --output = name of output file", "")
    ]
    
    print "Usage: extract_pseudogene_distribution.py -i input_eqn -o output"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def extract_pseudogene_distribution(in_fh, out_fh):
    
    pseudo = 0

    out_fh.write("Name\tExpected\tItself\tPseudogenes\tOthers\n")
    for line in in_fh:
        line = line.rstrip()
        line = line.replace(' ','')
        if(line == "Itself+Pseudogenes+Others"):
            pseudo = 1
        if(line == "Itself"):
            pseudo = 0
        
        # decipher equation
        if(pseudo == 1):
            equal_sign = line.split("=")
            if(len(equal_sign) > 1):
                lhs = equal_sign[0] # expected value
                rhs = equal_sign[1] # composition

                (gene, expect) = value_name_pair(lhs)
                itself_value = 0
                pseudo_value = 0
                other_value = 0
                
                rhs_split = rhs.split("+")
                
                for r in rhs_split:
                    (name, value) = value_name_pair(r)
                    if(name == gene):
                        itself_value += value
                    elif(name == "Other" or check_gene_type(name) == 1):
                        other_value += value
                    elif(check_gene_type(name)==2):
                        pseudo_value += value
            
                out_fh.write("%s\t%f\t%f\t%f\t%f\n" %(gene, expect, itself_value, pseudo_value, other_value))

                    

def value_name_pair(s):
    data = s.split("*")
    return((data[1], float(data[0])))

"""
    Return Value:
        1 = parent
        2 = pseudo
"""

def check_gene_type(name):
    parent = re.match('^ENSP', name)
    if(parent):
        return 1

    pseudo = re.match('^PGOHUM', name)    
    if(pseudo):
        return 2

    return 0

def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["input=", "output="])
    
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

    try:
        input != None
        output != None
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)


    in_fh = open(input, 'rb')
    out_fh = open(output, 'w')
    
    extract_pseudogene_distribution(in_fh, out_fh)

    in_fh.close()
    out_fh.close()
    
    print "Write to file %s" %(output)
    print "DONE\n\n"



if __name__ == "__main__":
    main()
