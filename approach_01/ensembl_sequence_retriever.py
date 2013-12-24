""" ensemble_sequence_retriever.py 
Usage: python ensemble_sequence_retriever.py -i input_file -o output_file
Author: Chelsea Ju
Date: 2013-06-26
Function: read in a list of ensemble ids, query the ensemble site and output the sequence of transcripts
"""

import httplib2, sys, getopt, json
from time import sleep

def usage():
    options = [ ("-i input_file, --infile = input file", "input file"),
                ("-o output_file, --outfile = output file", "output file")]
    
    print "Usage: python ensemble_sequence_retriever.py -i input -o output\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"

def sequence_parser(in_fh, out_fh):

    # retrieve sequence from ensemble API  
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    access_counter = 0
    
    for ensemble in in_fh:
        ensemble = ensemble.rstrip()
        
        ## retrieve gene information (+ or - strand)
        ext = "/feature/id/"+ensemble+"?feature=gene"
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        access_counter = reset_access(access_counter)

        gene_access = 0
        while (not resp.status == 200) and (gene_access < 10):
            sleep(2)
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
            gene_access +=1

        if( not resp.status == 200):
	        print "Invalid ID: %s\n" %(ensemble)
        else:
            decoded = json.loads(content)
            strand = 0
            for g in decoded:
                id = g['ID']
                if(id == ensemble):
                    strand = g['strand']
                    break

            if(strand == 0):
                print "Warning, no valid information for %s" % (ensemble)

            ## query cds (transcrips)
            ext2 = "/feature/id/"+ensemble+"?feature=cds"
            resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
            access_counter = reset_access(access_counter)

            while not resp2.status == 200:
                sleep(2)
                resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
                print server+ext2

            decoded2 = json.loads(content2)
            seq = ""
            varient = 0
            for d in decoded2:
                start = d['start']
            	end = d['end']
            	strand2 = d['strand']
            	region = d['seq_region_name']
            	rank = d['rank']

                if(strand2 == strand):
                    ## retrieve sequence
                    ext3 = "/sequence/region/human/%s:%s..%s:%s?" %(region,start,end, strand)
                    resp3, content3 = http.request(server+ext3, method="GET", headers={"Content-Type":"application/json"})
                    access_counter = reset_access(access_counter)
                
                    while not resp3.status == 200:
                        sleep(5)
                        resp3, content3 = http.request(server+ext3, method="GET", headers={"Content-Type":"application/json"})
                        print server+ext3

                    seq_decoded = json.loads(content3)                    
                    if(rank == 1 and seq != ""):
                        out_fh.write(">%s_%s\n%s\n" %(ensemble,varient,seq))
                
                    if(rank == 1):
                        seq = seq_decoded['seq']
                        varient +=1
                    else:
                        seq += seq_decoded['seq']                         
            
            if(seq != ""):
                out_fh.write(">%s_%s\n%s\n" %(ensemble, varient, seq))


def reset_access(count):
    if(count % 5 ==0):
        sleep(3)
    
    count += 1
    return count


def main():
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["infile=", "outfile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        str(err) # will print something like "option -a not recognized"
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
        print "Missing input file or output file"
        usage()
        sys.exit(2)


    in_fh = open(infile, 'rb')
    out_fh = open(outfile, 'w')   
 
    sequence_parser(in_fh, out_fh)
    
    in_fh.close()
    out_fh.close()
 
    print "Write to file %s" %(outfile)
    print "DONE\n\n"


if __name__ == "__main__":
    main()



