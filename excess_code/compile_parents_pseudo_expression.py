""" compile_parents_pseudo_expression.py
Usage: python compare_parents_pseudo_expression.py -r parent_pseudo.txt -p parent.expression -s pseudo.expression -o outfile_prefix 
Input: parent expression file, pseudo expression file, and parent vs pseudo relationship
Output: a file that combine the parent and pseudo gene expressions, and a file with linear models statistics
Author: Chelsea Ju
Function: Compile the fragment mapped data of pseudogene and parent gene
Date: 2013-08-03
"""


import sys, getopt, re, os

EXPRESSION = dict()

def usage():
    options = [ ("-r parent_pseudo.txt, --relation = file describe parent pseudo relationship", ""),
                ("-p parent.expression, --parent = parent expression file", ""),
                ("-s pseudo.expression, --pseudo = pseudo expression file", ""),
                ("-o outprefix, --outprefix = prefix of output file name", "")
    ]
    
    print "Usage: python compare_parents_pseudo_expression.py -r parent_pseudo.txt -p parent.expression -s pseudo.expression -o outfile_prefix"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def retrieve_expression(relation_fh, out_expression_fh):
    global EXPRESSION
    
    for x in relation_fh:
        data = x.split("\t")
        pseudo_name = data[0]
        parent_name = data[5]
        pseudo_type = data[20]
        
        if(EXPRESSION.has_key(pseudo_name) and EXPRESSION.has_key(parent_name)):
            out_expression_fh.write("%s\t%s\t%s\t%s\t%s\n" %(pseudo_name, EXPRESSION[pseudo_name], parent_name, EXPRESSION[parent_name], pseudo_type))
    
def build_expression_db(expression_fh):    
    global EXPRESSION
    
    for x in expression_fh:
        x = x.rstrip()
        data = x.split("\t")
        name = data[0]
        fragment = float(data[6])
        fpkm = float(data[7])

        match = re.search("(.*)_(.*)", name)
        if(match != None):
            name = match.group(1)
        EXPRESSION[name] = "%s\t%s" % (str(fragment), str(fpkm))
    

def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:p:s:o:", ["relation=", "parent=", "pseudo=", "outprefix=" ])
    
    except getopt.GetoptError, err:
    
    # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-r", "--relation"):
            relation_file = str(a)
        elif o in ("-o", "--outprefix"):
            outfile = str(a)
        elif o in ("-p", "--parent"):
            parent_file = str(a)
        elif o in ("-s", "--pseudo"):
            pseudo_file = str(a)

        else:
        	assert False, "unhandled option"
    
    try:
        relation_file != None
        parent_file != None
        pseudo_file != None
        outfile != None
    
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)

    out_expression = outfile + ".expression"

    relation_fh = open(relation_file, 'rb')
    parent_fh = open(parent_file, 'rb')
    pseudo_fh = open(pseudo_file, 'rb')

    out_expression_fh = open(out_expression, 'w')

    ## function starts here
    build_expression_db(parent_fh)
    build_expression_db(pseudo_fh)
    retrieve_expression(relation_fh, out_expression_fh)

    relation_fh.close()
    parent_fh.close()
    pseudo_fh.close()
    out_expression_fh.close()

    print "Write to file %s" %(out_expression)
    print "DONE\n\n"



if __name__ == "__main__":
    main()



