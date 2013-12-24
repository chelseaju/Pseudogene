""" pseudogene_expression_specifier.py 
Usage: python pseudogene_expression_specifier.py -i infile -o outfile
Argument: input = pseudogene expression output by pseudogene_expression_parser.py
          output = absolute path of the output file
Author: Chelsea Ju
Date: 2013-05-07
Function: Read in the pseudogene expression. Return a file with the best expression of the pseudogenes and the parent genes.
"""

import httplib2, sys, getopt, json
import sys
import getopt
import re
from Bio.Seq import Seq
from time import sleep

def usage():
    options = [ ("-i pseudogene_expression, --infile = pseudogene expression file", "expression file from pseudogene_exprssion_parser.py"),
                ("-o outfile, --outfile = output file name", "output file")
            ]
    
    print "Usage: python pseudogene_expression_specifier.py -i pseudogene_expression -o outfile"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"

"""
Function: extract the input file
"""
def retrieve_data(in_fh, out_fh):

    out_fh.write("Pseudo_Gene\tPseudo_Class\tPseudo_Chr\tPseudo_Start\tPseudo_End\tPseudo_Read_Start\tPseudo_Read_End\tPseudo_Strand\tPseudo_FPKM\tPseudo_Start_Diff\tPseudo_Notation\tPseudo_End_Diff\tPseudo_Seq\tParent_Gene\tParent_Chr\tParent_Start\tParent_End\tParent_Read_Start\tParent_Read_End\tParent_Strand\tParent_FPKM\tParent_Start_Diff\tParent_Notation\tParent_End_Diff\tParent_Seq\n")
    count = 0
    for line in in_fh.readlines():
        (pseudo_id, pseudo_chr, pseudo_start, pseudo_end, pseudo_strand, pseudo_type, pseudo_fpkm,
         parent_id, parent_chr, parent_start, parent_end, parent_strand, parent_fpkm) = line.split("\t")
         
        pseudo_expression = extract_fpkm(pseudo_chr, pseudo_strand, pseudo_start, pseudo_end, pseudo_fpkm)
        parent_expression = extract_fpkm(parent_chr, parent_strand, parent_start, parent_end, parent_fpkm)
 
        for pseudo in pseudo_expression:
            alignment_score = 0
            for parent in parent_expression:
                (pu_start, pu_end, pu_seq) = pseudo_expression[pseudo].split("||")
                (p_start, p_end, p_seq) = parent_expression[parent].split("||")
                
                pseudo_notation = get_boundary_notation(pseudo_start, pseudo_end, pu_start, pu_end)
                parent_notation = get_boundary_notation(parent_start, parent_end, p_start, p_end)
            
                out_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pseudo_id, pseudo_type, pseudo_chr, pseudo_start, pseudo_end, pu_start, pu_end, pseudo_strand, pseudo, pseudo_notation, pu_seq, parent_id, parent_chr, parent_start, parent_end, p_start, p_end, parent_strand, parent, parent_notation, p_seq))
        
        count +=1
        if(count % 100 == 0):
            sleep(1)
                   
 
"""
 Function: return a list of best representing fpkm from the given information       
"""        
def extract_fpkm(chr, strand, start, end, fpkm):

    expression = {}
    
    if(fpkm == ""):
        expression["0"] = "%s||%s" %(start, end)
    
    else:
        fpkms = fpkm.split(";")
        fpkms.pop()
        
        for e in fpkms:
            (r_start, r_end, r_expr) = e.split("||")
            position = "%s||%s" %(r_start, r_end)
            
            if(not expression.has_key(str(r_expr))):
                expression[str(r_expr)] = position
            else:
                record = expression[str(r_expr)]
                (db_start, db_end) = record.split("||")
                
                if((abs(int(r_start)-int(start)) + abs(int(r_end) - int(end))) < (abs(int(db_start) - int(start)) + abs(int(db_end) - int(end)))):
                    expression[str(r_expr)] = position

    # retrieve sequence from ensemble API  
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"

    for e in expression:
        if(e != 0):
            (r_start, r_end) = expression[e].split("||")  
            ext = "/sequence/region/human/%s:%s..%s:1?" %(chr,r_start,r_end)
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        
        if not resp.status == 200:
            print "Invalid response: ", resp.status
        else:
            decoded = json.loads(content)
            seq = decoded['seq']
            if(strand == "-1"):
                seq = Seq(seq).reverse_complement().tostring()
            expression[e] = "%s||%s||%s" %(r_start, r_end, seq)
        
    return expression
    
def get_boundary_notation(gene_start, gene_end, read_start, read_end):
    
    notation = ['','','','']
    
    if(gene_start == read_start):
        notation[0] = "|*"
    if(gene_start < read_start):
        notation[0] = "|"
        notation[1] = "*"
    if(gene_start > read_start):
        notation[0] = "*"
        notation[1] = "|"
    if(gene_end == read_end):
        notation[3] = "*|"
    if(gene_end > read_end):
        notation[2] = "*"
        notation[3] = "|"
    if(gene_end < read_end):
        notation[2] = "|"
        notation[3] = "*"
    
    start_diff = int(gene_start) - int(read_start)
    end_diff = int(gene_end) - int(read_end)
    
    symbol =  "--".join(notation)
    
    return "%s\t%s\t%s" %(str(start_diff), symbol, str(end_diff))
    
def main():

    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["infile=", "outfile=" ])
    
    except getopt.GetoptError, err:
    
    # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = str(a)
        elif o in ("-o", "--outfile"):
            outfile = str(a)
        else:
            assert False, "unhandled option"
    
    try:
        infile != None
        outfile != None
    
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)
    
    in_fh = open(infile, 'rb')
    out_fh = open(outfile, 'w')
    
    retrieve_data(in_fh, out_fh)
    
    
    out_fh.close()
    in_fh.close()
    
    print "Write to file %s" %(outfile)
    print "DONE\n\n"



if __name__ == "__main__":
    main()




