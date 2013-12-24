""" pseudogene_parent_parser.py 
Usage: python pseudogene_parent_parser.py -i pseudogene_file -o output_prefix
Author: Chelsea Ju
Date: 2013-06-18
Function: read in the pseudogene file, parse out the sequence for pseudogenes and retrieve parents information. Pseudogenes and parents information are also outputted in BED format.
"""

import sys, getopt, httplib2, json, re
from time import sleep

def usage():
    options = [ ("-i pseudogene_file, --infile = pseudogene file", "pseudogene annotation file"),
                ("-o output_prefix, --outfile = output file prefix", "output filename prefix")]
    
    print "Usage: python pseudogene_parent_parser.py -i infile -o outfile_prefix\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"

"""
Function: parse out the psuedogene sequence, and the parent id
"""
def sequence_parser(in_fh, parent_bed_fh, pseudo_bed_fh, pseudo_seq_fh):
    parents = {}

    # retrieve sequence from ensemble API  
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    access_counter = 0
    
    next(in_fh) #remove the title line   

    ## pseudogenes
    for rec in in_fh:
        info = rec.split('\t')
        pseudo_id = info[0]
        pseudo_chr = "chr" + info[1]
        pseudo_start = int(info[2])
        pseudo_end = int(info[3])
        pseudo_strand = info[4]
        parent_id = info[8]
        pseudo_exons = info[18]
        pseudo_seq = info[21]
        
        pseudo_block_count = 0
        pseudo_block_size = []
        pseudo_block_start = []
        
        exons = pseudo_exons.split("],")
        for e in exons:
            e = re.sub(r'\[|\]|\s', '', e)
            (e_start, e_end) = e.split(",")
            pseudo_block_count += 1
            pseudo_block_start.append(str(int(e_start) - int(pseudo_start)))
            pseudo_block_size.append(str(int(e_end) - int(e_start) + 1))
        
        pseudo_bed_fh.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" %(pseudo_chr, pseudo_start, pseudo_end, pseudo_id, 0, pseudo_strand, pseudo_start, pseudo_start, 0, pseudo_block_count, (",".join(pseudo_block_size)), (",".join(pseudo_block_start))))
        pseudo_seq_fh.write(">%s\n%s\n" %(pseudo_id, pseudo_seq))
        parents[parent_id] = ""

    ## parents
    for p in parents:

        ext = "/feature/id/"+p+"?feature=gene"
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
                if(id == p):
                    strand = g['strand']
                    start = g['start']
                    end = g['end']
                    chr = "chr" + g['seq_region_name']
                    break

            if(strand == 0):
                print "Warning, no valid information for %s" % (ensemble)
            
            else:                
                ## query cds (transcrips)
                ext2 = "/feature/id/"+p+"?feature=cds"
                resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
                access_counter = reset_access(access_counter)
    
                while not resp2.status == 200:
                    sleep(2)
                    resp2, content2 = http.request(server+ext2, method="GET", headers={"Content-Type":"application/json"})
                    print server+ext2
    
                decoded2 = json.loads(content2)
                variant = 0
                exons_start = []
                exons_size = []
                exon_count = 0

                if(strand == 1):
                    p_strand = "+"
                else:
                    p_strand = "-"
                p_start = str(int(start) - 1)
                p_end = str(int(end) - 1)

                for d in decoded2:
                    e_start = d['start']
                    e_end = d['end']
                    e_strand = d['strand']
                    rank = d['rank']
    
                    if(e_strand == strand):
                        
                        if(rank == 1 and len(exons_start) !=0):
                            parent_bed_fh.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" %(chr, p_start, p_end, p+"_"+str(variant), 0, p_strand, p_start, p_start, 0, exon_count, (",".join(exons_size)), (",".join(exons_start))))
                        
                        if(rank == 1):
                            exons_start = []
                            exons_size = [] 
                            exon_count = 0                           
                            variant += 1

                        exon_count += 1
                        exons_start.append(str(int(e_start)-int(start)))
                        exons_size.append(str(int(e_end)-int(e_start)+1))
            
                if(len(exons_start)!=0):
                    parent_bed_fh.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" %(chr, p_start, p_end, p+"_"+str(variant), 0, p_strand, p_start, p_start, 0, exon_count, (",".join(exons_size)), (",".join(exons_start))))

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
        print "Missing input file or output prefix"
        usage()
        sys.exit(2)

    parent_bed = outfile+"_parents.bed"
    pseudo_bed = outfile+"_pseudo.bed"
    pseudo_seq = outfile+"_pseudo.fa"

    in_fh = open(infile, 'rb')
    parent_bed_fh = open(parent_bed, 'w')
    pseudo_bed_fh = open(pseudo_bed, 'w')
    pseudo_seq_fh = open(pseudo_seq, 'w')
    
    sequence_parser(in_fh, parent_bed_fh, pseudo_bed_fh, pseudo_seq_fh)
    
    in_fh.close()
    parent_bed_fh.close()
    pseudo_bed_fh.close()
    pseudo_seq_fh.close()
    print "Write to file %s" %(parent_bed)
    print "Write to file %s" %(pseudo_bed)
    print "Write to file %s" %(pseudo_seq)
    print "DONE\n\n"


if __name__ == "__main__":
    main()


