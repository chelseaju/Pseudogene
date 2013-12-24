""" pseudogene_parent_retriever.py 
Usage: python pseudogene_parent_retriever.py -i ensemble_protien_list -o output_prefix
Input: a list of ensemble protein ids generated from pseudogene_pseudo_parser
Output: parent gene information in BED format
Function: retrieve parent protein information from ensemble API Rest. 
Author: Chelsea Ju
Date: 2013-07-19
Note: Modify from pseudogene_parent_parser
"""

import sys, getopt, httplib2, json, re
from time import sleep

def usage():
    options = [ ("-i ensemble_protein_ids, --infile = ensemble proteins", ""),
                ("-o output_prefix, --outfile = prefix of output file", "")]
    
    print "Usage: python pseudogene_parent_retriever.py -i ensemble_protien_list -o output_prefix\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def ensemble_retriever(in_fh, parent_bed_fh):
    # retrieve sequence from ensemble API  
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    access_counter = 0
    
    for ensembl in in_fh:
        ensembl = ensembl.rstrip()
        print ensembl
        ## retrieve gene information (+ or - strand)
        ext = "/sequence/id/"+ensembl+"?"
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        access_counter = reset_access(access_counter)

        protein_access = 0
        while (not resp.status == 200) and (protein_access < 10):
            sleep(2)
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
            protein_access +=1

        if( not resp.status == 200):
            print "Invalid ID: %s\n" %(ensembl)
        else:
            decoded = json.loads(content)
            sequence = decoded["seq"]
            seq_size = len(sequence) + 1

            # if the sequence length is valid, retrieve genomic sequence
            gene_access = 0           
            if(seq_size > 1):
                ext2 = "/map/translation/"+ensembl+"/1.."+str(seq_size)+"?"
                resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
                access_counter = reset_access(access_counter)

                while not (resp2.status == 200) and (gene_access < 10):
                    sleep(2)
                    resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
                    gene_access += 1

                if(not resp.status == 200):
                    print "Unable to retrieve server+ext2"
                    
                else:
                    decoded2 = json.loads(content2)
                    
                    for e in decoded2['mappings']:
                        if(e['strand'] == 0):
                            decoded2['mappings'].remove(e)
                    
                    exon_count = len(decoded2['mappings'])
                    block_size = []
                    block_start = []
                    
                    if(exon_count > 0):
                        strand = decoded2['mappings'][0]['strand']
                        chr = "chr" + decoded2['mappings'][0]['seq_region_name']
                        
                        # start and end position depends on the chromosome orientation
                        if(strand == -1):
                            start = decoded2['mappings'][-1]['start']
                            end = decoded2['mappings'][0]['end']
                        elif(strand == 1):
                            start = decoded2['mappings'][0]['start']
                            end = decoded2['mappings'][-1]['end']
                         
                        for e in decoded2['mappings']:
                            e_size = e['end'] - e['start'] + 1
                            e_relative_start = e['start'] - start 
            
                            block_size.append(str(e_size))
                            block_start.append(str(e_relative_start))

                        # UCSC ensemble is 0 based
                        start = start - 1
                        end = end - 1
                                   
                        parent_bed_fh.write("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" %(chr, start, end, ensembl, 0, strand, start, start, 0, exon_count, ",".join(block_size), ",".join(block_start)))



def reset_access(count):
    if(count % 5 ==0):
        sleep(3)
    
    count += 1
    return count



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
        print "\nMissing input file or output prefix"
        usage()
        sys.exit(2)

    parent_bed = outfile+"_parents.bed"

    in_fh = open(infile, 'rb')
    parent_bed_fh = open(parent_bed, 'w')

    ensemble_retriever(in_fh, parent_bed_fh)
    
    in_fh.close()
    parent_bed_fh.close()
    print "Write to file %s" %(parent_bed)
    print "DONE\n\n"


if __name__ == "__main__":
    main()
