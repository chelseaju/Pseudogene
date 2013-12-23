""" pseudogene_parent_bed_retriever.py 
Usage: python pseudogene_parent_bed_retriever.py -i ensemble_transcript_sequence -o output_prefix
Input: a fasta file of ensemble transcript sequence generated from pseudogene_pseudo_parser.R
Output: parent mRNA information in BED format
Function: retrieve parent transcript information from ensemble API Rest. 
Author: Chelsea Ju
Date: 2013-07-19
"""

import sys, getopt, httplib2, json, re, itertools
from time import sleep

def usage():
    options = [ ("-i ensemble_transcript_ids, --infile = ensemble transcripts", ""),
                ("-o output_prefix, --outfile = prefix of output file", "")]
    
    print "Usage: python pseudogene_parent_bed_retriever.py -i ensemble_transcript_list -o output_prefix\nOptions"
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def ensemble_retriever(in_fh, parent_bed_fh):
    # retrieve sequence from ensemble API  
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    access_counter = 0

    for id, seq in itertools.izip_longest(in_fh, in_fh, fillvalue=''):
        id = id.rstrip()
        size = len(seq) -1
        (peptide, transcript) = id.split("_")

        ext = "/map/cdna/"+transcript+"/1.."+str(size)
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        access_counter = reset_access(access_counter)

        est_access = 0
        while (not resp.status == 200) and (est_access < 10):
            sleep(2)
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
            est_access +=1

        if( not resp.status == 200):
            print "Invalid ID: %s\n" %(transcript)
        else:
            decoded = json.loads(content)

            if(decoded['mappings'][0]['strand'] == 0):
                decoded['mappings'].pop(0)
            if(decoded['mappings'][-1]['strand'] == 0):
                decoded['mappings'].pop()
            
            ## exons order returned by ensembl: first coding exon to the last coding exon
            ## BED format: sequence positions always stay with forward strand. 
            ## For reverse sequence, it retrieves all exon sequences in forward strand, concatenates the sequences, and reverse transcription.
            #### It means..  the start position for forward strand is located on the first exon, but on the last exon for reverse strand
            ####             the end position for forward strand is located on the last exon, but on the first exon for forward strand
            ####             when write in BED format, exons order stay the same for forward strand, but needs to be reversed for the reverse strand
            
            chr = "chr" + decoded['mappings'][0]['seq_region_name']
            strand = decoded['mappings'][0]['strand']
            if(strand == 1):
                strand = "+"
                start = decoded['mappings'][0]['start']
                end = decoded['mappings'][-1]['end']
            else:
                strand = "-"
                start = decoded['mappings'][-1]['start']
                end = decoded['mappings'][0]['end']
 
 
            exon_starts = []
            exon_sizes = []
            exon_count = 0
            exon_total_size = 0
            
            for e in decoded['mappings']:
                exon_count += 1
                exon_starts.append(e['start'])
                e_size = e['end'] - e['start'] + 1
                exon_total_size += e_size
                exon_sizes.append(str(e_size))

            if(exon_total_size < size):
                extra_size = size - exon_total_size
                exon_sizes.insert(0,str(extra_size))
                exon_count += 1
                
                if(strand == "-"):
                    exon_starts.insert(0,(end + 1))
                    end = end + extra_size
                else:
                    exon_starts.insert(0, (start - extra_size))
                    start = start - extra_size
                    print extra_size

            exon_starts = [str(x - start) for x in exon_starts]
            
            # UCSC 0-based
            start = start - 1
            end = end - 1


            if(strand == "-"):
                exon_starts.reverse()
                exon_sizes.reverse()
            
#            print "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s" %(chr, start, end,id.lstrip(">"), 0, strand, start, start, 0, exon_count, ",".join(exon_sizes), ",".join(exon_starts))
            parent_bed_fh.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" %(chr, start, end,id.lstrip(">"), 0, strand, start, start, 0, len(decoded['mappings']), ",".join(exon_sizes), ",".join(exon_starts)))
            
            
            
            
"""
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

"""

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
