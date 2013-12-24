""" pseudogene_parser.py 
Usage: python parse_pseudogene.py -i pseudogene_file -o output_name
Author: Chelsea Ju
Date: 2013-05-06
Function: read in the pseudogene file, and parse out the positions and parents information for pseudogenes 
"""

import httplib2, sys, getopt, json
from time import sleep


def usage():
        options = [ ("-i pseudogene_file, --infile = pseudogene file", "pseudogene annotation file"),
		    ("-o output_name, --outfile = output file", "output file name")]

        print "Usage: python parse_pseudogene.py -i infile -o outfile\nOptions"
        for (o,d) in options:
                print "\t" + o + "\t" + d
        print "\n"


"""
Function: retrieve the ensembl database for given gene
"""
def retrieve_ensembl_gene(gene):

	print "Accessing %s" % (gene)
	http = httplib2.Http(".cache")
 
	server = "http://beta.rest.ensembl.org"
	ext = "/feature/id/"+gene+"?feature=gene"
	resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
 
	if not resp.status == 200:
  		print "Invalid response: ", resp.status
  		return ("NA","NA","NA","NA","NA")
 
	decoded = json.loads(content)
	
	name = str(decoded[0]['external_name'])
	chromosome = str(decoded[0]['seq_region_name'])
	start = decoded[0]['start']
	end = decoded[0]['end']
	strand = str(decoded[0]['strand'])

	if(strand == "1"):
		strand = "+"
	else:
		strand = "-"

	return (name, chromosome, strand, start, end)


"""
Function: iterates through the file handler and write to the output file handler
"""
def parser(in_h, out_h):
	
	 # Pseudogene format: <ID> <CHR> <START> <END> <STRAND> <PARENT_PROTEIN> <PARTEIN_START> <PROTEIN_END> <PARENT_GENE> <FRACTION> <INSERTION> <DELETIONS> <SHIFTS> <STOPS> <EVALUE> <POLYA> <DISABLEMENTS> <EXONS> <INTRONS> <CLASS> <SEQUENCE> <LINK>
       	next(in_h)
	count = 0
	for rec in in_h:
		info = rec.split('\t')
		pseudo_id = info[0]
		pseudo_chr = info[1]
		pseudo_start = info[2]
		pseudo_end = info[3]
		pseudo_strand = info[4]
		parent_gene = info[8]
		pseudo_type = info[20]	

		(parent_name, parent_chr, parent_strand, parent_start, parent_end) = retrieve_ensembl_gene(parent_gene)

		out_h.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(pseudo_id, pseudo_chr, pseudo_start, pseudo_end, pseudo_strand, pseudo_type, parent_gene, parent_chr, parent_strand, parent_start, parent_end))

		count +=1

		## sleep one second for every 100 retrieval to avoid site access blockage
		if(count % 100 == 0):
			sleep(1)	



def main():
        # read in user's input
        try:
                opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["infile=", "outfile="])
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
                input != None

        except:
                print "Missing input file"
                usage()
                sys.exit(2)

	if(outfile == None):
		outfile = "pseudo_genes_data.txt"


        in_fh = open(infile, 'rb')
        out_fh = open(outfile, 'w')

        parser(in_fh, out_fh);

        out_fh.close()
        in_fh.close()
	
	print "Write to file %s" %(outfile)
	print "DONE\n\n"


if __name__ == "__main__":
    main()

