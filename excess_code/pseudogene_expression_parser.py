""" pseudogene_expression_parser.py 
Usage: python pseudogene_expression_parser.py -p pseudogene -e expression -o outfile
Author: Chelsea Ju
Date: 2013-05-07
Function: Read in the pseudogene information, and the expression data. Return a file with expression of the pseudogenes and the parent genes.
"""

import sys
import getopt
import re
import sqlite3

def usage():
    options = [ ("-p pseudogene, --pseudogene = pseudogene file", "pseudogene and parent gene information"),
                ("-e expression, --expression = expression file", "output from cufflinks"),
                ("-o outfile, --outfile = output file name", "output file")
    ]
    
    print "Usage: python pseudogene_expression_parser.py -p pseudogene -e expression -o outfile"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"
    


"""
Function: read in the expression data, and build a temporary database
"""
def build_expression_database(expression_fh, db_conn):
    
    print "\nBuilding the temporary expression database ........"
    
    # clear database
    drop_stat = "DROP TABLE IF EXISTS expression"
    db_conn.execute(drop_stat)
    
    # create database
    create_stat = '''CREATE TABLE expression 
    		 (chr char(255), start int, end int, strand int, fpkm float);'''
    index_chr_stat = '''CREATE INDEX chr_idx ON expression(chr);'''
    index_start_stat = '''CREATE INDEX start_idx ON expression(start);'''
    index_end_stat = '''CREATE INDEX end_idx ON expression(end);'''
    
    db_conn.execute(create_stat)
    db_conn.execute(index_chr_stat)
    db_conn.execute(index_start_stat)
    db_conn.execute(index_end_stat)
    
    
    # avoid repeate record, a hash table is used to mark the positions
    position = {}
    
    ## read in data
    for rec in expression_fh:
        (chromosome, program, seq_type, start, end, length, strand, frame, others) = rec.split('\t')
        
        if(strand == "+"):
            strand = 1
        else:
            strand = -1
        
        pos = "%s-%s-%s-%s" % (chromosome, strand, start, end)
        
        fpkm_search = re.match(r".+ FPKM \"(.+)\"; frac", others)
        
        if(fpkm_search != None):
            fpkm = float(fpkm_search.group(1))
        else:
            fpkm = 0
        
        position[pos] = (chromosome, int(start), int(end), strand, fpkm) 
    
    
    # insert statements into sqlite db
    db_conn.executemany('INSERT INTO expression VALUES (?,?,?,?,?)', position.values())
    db_conn.commit()

"""
Function: retrieve the fpkm of the closest transcript for given location
"""
def query_fpkm(query_chr, query_start, query_end, query_strand, db_conn):

    statement = """ SELECT * FROM expression WHERE chr = '%s' AND strand = '%s' AND (
    		(start = %d AND end = %d) OR
    		(start = %d AND end > %d) OR
    		(start = %d AND end < %d) OR
    		(start < %d AND end = %d) OR
    		(start > %d AND end = %d) OR
    		(start < %d AND end > %d) OR
    		(start < %d AND end > %d AND end < %d) OR
    		(start > %d AND start < %d AND end < %d) OR
    		(start > %d AND start < %d AND end > %d) ); 
    	    """ % (query_chr, query_strand, query_start, query_end,
    		   query_start, query_end, query_start, query_end,
    		   query_start, query_end, query_start, query_end,
    		   query_start, query_end,   
    		   query_start, query_start, query_end,
    		   query_start, query_end, query_end,
    		   query_start, query_end, query_end)
    cur = db_conn.cursor()
    cur.execute(statement)
    
    result = ""
    rows = cur.fetchall()
    
    for(chromosome, start, end, strand, fpkm)  in rows:
        output = "%s||%s||%f;" %(start, end, fpkm)
        result += output
    	
    
    return result


def retrieve_expression(pseudo_fh, out_fh, db_conn):

    print "\nRetrieving expressions ............................................"
    
    title = pseudo_fh.readline()
    for rec in pseudo_fh.readlines():
        rec = rec.rstrip()	
        (pseudo_id, pseudo_chr, pseudo_start, pseudo_end, pseudo_strand,pseudo_class, parent_gene, parent_chr, parent_strand, parent_start, parent_end) = rec.split('\t')
        pseudo_fpkm = query_fpkm(pseudo_chr, int(pseudo_start), int(pseudo_end), pseudo_strand, db_conn)
        if(pseudo_fpkm != "" and parent_chr !="NA"):
        
            parent_fpkm = query_fpkm(parent_chr, int(parent_start), int(parent_end), parent_strand, db_conn)
            
            output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pseudo_id, pseudo_chr, pseudo_start, pseudo_end, pseudo_strand, pseudo_class, pseudo_fpkm, parent_gene, parent_chr, parent_start, parent_end, parent_strand, parent_fpkm)
            
            out_fh.write(output)	
	
def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:e:o:", ["pseudogene=", "expression=", "outfile=" ])
    
    except getopt.GetoptError, err:
    
    # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-p", "--pseudogene"):
            pseudogene = str(a)
        elif o in ("-o", "--outfile"):
            outfile = str(a)
        elif o in ("-e", "--expression"):
            expression = str(a)
        else:
        	assert False, "unhandled option"
    
    try:
        pseudogene != None
        expression != None
        outfile != None
    
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)
    
    
    pseudo_fh = open(pseudogene, 'rb')
    expression_fh = open(expression, 'rb')
    out_fh = open(outfile, 'w')

    ## build database
    db_name = "cufflink_expression.db"
    conn = sqlite3.connect(db_name)
    
    build_expression_database(expression_fh, conn)
    #	print query_fpkm("13", 31903638, 31904805, "+", conn)
    retrieve_expression(pseudo_fh, out_fh, conn)
    conn.close()
    
    
    out_fh.close()
    pseudo_fh.close()
    expression_fh.close()
    
    print "Write to file %s" %(outfile)
    print "DONE\n\n"



if __name__ == "__main__":
    main()




