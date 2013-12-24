""" sequence_parser.py 
Usage: python sequence_parse.py -i pseudogene_file -o output_prefix
Author: Chelsea Ju
Date: 2013-06-18
Function: read in the pseudogene file, and retreive the sequence for parents and pseudogenes
"""

import httplib2, sys, getopt, json
from time import sleep

def usage():
    options = [ ("-i pseudogene_file, --infile = pseudogene file", "pseudogene annotation file"),
                ("-o output_prefix, --outfile = output file prefix", "output filename prefix")]
    
    print "Usage: python parse_pseudogene.py -i infile -o outfile_prefix\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"

"""
Function: parse out the psuedogene sequence, and retrieve all the varient sequences of the parent
"""
def sequence_parser(in_fh, parent_fh, pseudo_fh):
    next(in_fh) #remove the title line
    
    for rec in in_fh:
        info = rec.split('\t')
        pseudo_id = info[0]
        parent_id = info[8]
        pseudo_seq = info[21]
        
        #accessing the suerver
        print "Accessing %s" %(parent_id)
        http = httplib2.Http(".cache")
        server = "http://beta.rest.ensembl.org"
        ext = "/feature/id/"+parent_id+"?feature=cds"
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    
        if not resp.status == 200:
            print "Invalid response: ", resp.status
 
        decoded = json.loads(content)

        parent_seq = []
        for d in decoded:


def main():
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["infile=", "outprefix="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = str(a)
        elif o in ("-o", "--outfile_prefix"):
            outfile = str(a)
        else:
            assert False, "unhandled option"

    try:
        infile != None
        outfile != None

    except:
        print "Missing input file or output prefix"
        usage()
        sys.exit(2)

    parent_file = outfile+"_parents.fa"
    pseudo_file = outfile+"_pseudo.fa"    


    in_fh = open(infile, 'rb')
    parent_fh = open(parent_file, 'w')
    pseudo_fh = open(pseudo_file, 'w')
    
    sequence_parser(in_fh, parent_fh, pseudo_fh)
    
    in_fh.close()
    parent_fh.close()
    pseudo_fh.close()
    
    print "Write to file %s" %(parent_file)
    print "Write to file %s" %(pseudo_file)
    print "DONE\n\n"


if __name__ == "__main__":
    main()


