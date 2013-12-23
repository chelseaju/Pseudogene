""" extract_distribution_matrix.py 
Usage: python extract_distribution_matrix.py -i input_matrix -o output_filename
Input: expression distribution matrix (generated from fragments_counter) 
Output: a file with formulas described the expression distribution
Author: Chelsea Ju
Function: Read in the distribution matrix, and output them in formulas
Date: 2013-08-19
"""

import sys, getopt, re

def usage():
    options = [ ("-i input_matrix --input = input matrix", ""),
                ("-o output, --output = output file name", "")
    ]
    
    print "Usage: extract_distribution_matrix.py -i input_matrix -o output_filename"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def extract_data(in_fh, out_fh):

    genes = in_fh.readline().rstrip().split("\t")

    itself_only = []
    itself_others = []
    itself_pseudogene = []
    itself_other_parents = []
    others_only = []

    for line in in_fh:
        data = line.split("\t")
        expression = data[2:]
        expected = float(data[1])
        current_gene = data[0]
                
        if(sum(float(x) for x in expression) > 0 and expected > 0):
            formula = [expression_to_string(expected, current_gene)]

            pseudogene = 0
            others = 0
            other_parents = 0
            itself = 0
            
            for index, value in enumerate(expression):
                if(float(value) > 0):
                    mapped_gene = genes[index+1]
                    formula.append(expression_to_string(value, mapped_gene))

                    match_pseudo = re.match('PGOHUM', mapped_gene)
                    match_parent = re.match('ENSP', mapped_gene)
                    match_current = re.match(mapped_gene, current_gene)
                    match_others = re.match('Others', mapped_gene)

                    
                    if(match_pseudo != None):
                        pseudogene = 1
                    
                    if(match_current == None and match_parent != None): # mapped to other parents, but not itself
                        other_parents = 1
                    
                    if(match_others != None):
                        others = 1
                    
                    if(match_current != None):
                        itself = 1
            
            if pseudogene == 1:
                itself_pseudogene.append("%s = %s" % (formula[0], " + ".join(formula[1:])))
            elif other_parents == 1:                        
                itself_other_parents.append("%s = %s" % (formula[0], " + ".join(formula[1:])))
            elif others == 1 and itself == 1:
                itself_others.append("%s = %s" % (formula[0], " + ".join(formula[1:])))
            elif itself == 1 and others == 0:
                itself_only.append("%s = %s" % (formula[0], " + ".join(formula[1:])))
            elif others == 1 and itself == 0:
                others_only.append("%s = %s" % (formula[0], " + ".join(formula[1:])))

    out_fh.write("Itself + Other Parents + Others\n")
    out_fh.write("\n".join(itself_other_parents))
    out_fh.write("\n\nItself + Others\n")
    out_fh.write("\n".join(itself_others))
    out_fh.write("\n\nItself + Pseudogenes + Others\n")
    out_fh.write("\n".join(itself_pseudogene))
    out_fh.write("\n\nItself \n")
    out_fh.write("\n".join(itself_only))
    out_fh.write("\n\nOthers Only \n")
    out_fh.write("\n".join(others_only))
            
            
def expression_to_string(coef, name):    
    return "%.2f * %s" %(float(coef), strip_gene_name(name))


def strip_gene_name(name):
    
    m = re.match('(ENSP.*)_ENST', name)
    if(m):
        return m.group(1)
    else:
        return name

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
    
    extract_data(in_fh, out_fh)

    in_fh.close()
    out_fh.close()
    
    print "Write to file %s" %(output)
    print "DONE\n\n"



if __name__ == "__main__":
    main()


