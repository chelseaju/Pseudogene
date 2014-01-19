""" transcript_identifier
Usage: python transcript_identifier.py -c chromosome -d directory 
Input:  -c chromosome name
        -d input/output directory 
Output: lines of transcript region or exon region with {name \t chr \t start \t end \n}
Function: 1. Map each exon position to the database to find a corresponding transcript or pseudogene
          2. If the exon does not have a match, return this exon
          4. Return the mapped boundary of each unmarked exon

Date: 2014-01-06
Author: Chelsea Ju
"""
import sys, re, os, random, argparse, sqlite3
import sqlite3

#DB = '/u/home/c/chelseaj/database/PseudogeneDB/pseudogene.db'
DB = '/home/chelseaju/Database/PseudogeneDB/pseudogene.db'


"""
    make the chromosome name readable for database
"""
def convert_chromosome_name(name):
    
    name = name.lower()
    match_chr_number = re.match(r"\D*(\d+)", name)
    match_chr_X = re.match(r"\D*(x)", name)
    match_chr_Y = re.match(r"\D*(y)", name)
    match_chr_M = re.match(r"\D*(m)", name)
    
    if(match_chr_number):
        return match_chr_number.group(1)
    
    elif(match_chr_X):
        return "X"

    elif(match_chr_Y):
        return "Y"

    elif(match_chr_M):
        return "MT"

    else:
        print ("unknown input %s", name)
        sys.exit(2)
  

"""
    Function: map the identified exon to known gene, and compute the read coverage range of the gene
        a valid mapping refers to covering at least 30% of the gene. if the coverage is less than 30%, it is considered a hypothetical gene
    
"""
def map_exon_to_gene(input_file, chr):
    
    mapped_regions = {}
    gene_list = []
    
    
    ## connect to database
    conn = sqlite3.connect(DB)
    
    ## read file
    input_fh = open(input_file, 'rb')

    for line in input_fh:
        line = line.rstrip()

        (start, end) = line.split('\t')

        gene_search = "SELECT t.transcript_id, t.protein_id, t.transcript_start, t.transcript_end FROM ensembl_exon as e INNER JOIN ensembl_mapping as m ON e.exon_id = m.exon_id INNER JOIN ensembl_transcript as t ON t.transcript_id = m.transcript_id WHERE t.chromosome_name = '%s' AND e.exon_chr_end > %d AND e.exon_chr_start < %d" % (str(chr), int(start), int(end))
        pseudogene_search = "SELECT p.id, p.start_coordinate, p.stop_coordinate FROM pseudogene as p WHERE p.chromosome = '%s' AND p.stop_coordinate > %d AND p.start_coordinate < %d " %(str(chr), int(start), int(end))

        c = conn.cursor()      
        c.execute(gene_search)
        found = False
        for r in c.fetchall():
            found = True
            if(r[1] != ""):
                id = r[1] + "::" + str(r[2]) + "::" + str(r[3])
                if(mapped_regions.has_key(id)):
                    (current_min, current_max) = mapped_regions[id]
                    mapped_regions[id] = (min(current_min, int(start)), max(current_max, int(end)))
                    
                else:
                    mapped_regions[id] = (int(start), int(end))
        
        # exon doesn't have any mapped parent gene, try the pseudogene        
        if(not found):
            c = conn.cursor()
            c.execute(pseudogene_search)
            
            for r in c.fetchall():
                found = True
                id = r[0] + "::" + str(r[1]) + "::" + str(r[2])
                
                if(mapped_regions.has_key(id)):
                    (current_min, current_max) = mapped_regions[id]
                    mapped_regions[id] = (min(current_min, int(start)), max(current_max, int(end)))
                    
                else:
                    mapped_regions[id] = (int(start), int(end))


        # exon doesn't have any mapped parent gene nor mapped psuedogene
        if(not found):
            id = "Unknown_" + chr +"_"+ str(start) + "_" + str(end)
            gene_list.append((id, int(start), int(end)))
     
    input_fh.close()
    conn.close()

    for k,v in mapped_regions.items():
        (gene_id, gene_start, gene_end) = k.split("::")
        (mapped_start, mapped_end) = v       
       
       ## a valid mapping requires covering at least 30% of the gene, 
       ## if it is less than 30%, return the exon range as a hypothetical gene
#        cover_length = min(int(gene_end), int(mapped_end)) - max(int(gene_start), int(mapped_start))
#        cover_percentage = float(cover_length) / float((int(gene_end) - int(gene_start)))
       
#        if(cover_percentage > 0.3):            
        gene_list.append((gene_id, mapped_start, mapped_end))
#        else:
#            new_id = str(chr) + "_" + str(mapped_start) + "_" + str(mapped_end)
#            gene_list.append((new_id, mapped_start, mapped_end))
    
    return gene_list

"""
    Function: print the list of gene to file
"""
def export_genes(gene_list, out_file, chr):
    
    out_fh = open(out_file, 'w')
    for (id, start, end) in gene_list:        
        out_fh.write("%s\t%s\t%d\t%d\n" %(id, chr, int(start), int(end)))
    
    out_fh.close()
    print ""
    print "Writing Gene List to File : %s" %(out_file)
       
def main(parser):
    
    options = parser.parse_args()
    chromosome_name = options.chromosome
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"
    
    ## input file
    input_file = dir + "mapping/" + chromosome_name + "_exons.txt"
    
    ## output file
    output_dir = dir + "mapping/"
    output_file = output_dir + chromosome_name + "_transcripts.txt"
    
    chr = convert_chromosome_name(chromosome_name)
    mapped_regions = map_exon_to_gene(input_file, chr)
    export_genes(mapped_regions, output_file, chromosome_name)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='gene_identifier.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
