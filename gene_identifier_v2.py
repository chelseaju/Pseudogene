""" gene_identifier_v2
Usage: python gene_identifier_v2.py -c chromosome -d directory 
Input:  -c chromosome name
        -d input/output directory 
Output: lines of gene region or unknown region with {name \t chr \t start \t end \n}
Function: 1. intersect a list of exons to ensembl gene database
          2. collapse the unmatch exons to create unknown region 

Date: 2014-03-02
Author: Chelsea Ju
Update from gene_identifier.py
"""
import sys, re, os, subprocess, random, argparse

ENSEMBL_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ensemblGene_pseudogene.bed"
#ENSEMBL_GENE = "/u/scratch/c/chelseaj/database/EnsemblTranscriptome/Genes/ensemblGene_pseudogene.bed"

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
    Function: map the identified exon to known gene
    
"""
def map_exon_to_gene(input_file):
        
    gene_list = {}
    unknown = []  
    unknown_collapse = []
    
    mapping = subprocess.check_output(["bedtools", "intersect", "-wb", "-loj", "-a", input_file, "-b", ENSEMBL_GENE])
    mapping_data = mapping.split("\n")
    
    for m in mapping_data:
        data = m.split("\t")
        if(len(data) > 2):
            mapped_gene_name = data[6]
            
            # found mapped gene
            if(mapped_gene_name != "-1" and mapped_gene_name != "."):
                gene_list[mapped_gene_name] = (data[3], data[4], data[5])
            
            else:
                unknown.append((data[0], data[1], data[2]))

    # collapse unknwon region
    unknwon = sorted(unknown)    
    if(len(unknown) > 0):
        unknown_collapse.append(unknwon[0])
    for u in unknown:
        previous = unknown_collapse[-1]
        if(int(u[1]) > int(previous[2]) and int(u[1]) - int(previous[2]) > 80):
            unknown_collapse.append(u)
        else:
            unknown_collapse[-1] = (previous[0], min(previous[1], u[1]), max(previous[2], u[2]))
 
    for u in unknown_collapse:
        name = "Unknown_" + str(u[0]) + "_" + str(u[1]) + "_" + str(u[2])
        gene_list[name] = u
            
    return gene_list

"""
    Function: print the list of gene to file
"""
def export_genes(gene_list, out_file):
    
    out_fh = open(out_file, 'w')
    for k in gene_list.keys():
        (chr, start, end) = gene_list[k]
        out_fh.write("%s\t%s\t%d\t%d\n" %(k, chr, int(start), int(end)))
    
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
    input_file = dir + "mapping/" + chromosome_name + "_exons.bed"
    
    ## output file
    output_dir = dir + "mapping/"
    output_file = output_dir + chromosome_name + "_genes.bed"
    
    chr = convert_chromosome_name(chromosome_name)
    mapped_regions = map_exon_to_gene(input_file)
    export_genes(mapped_regions, output_file)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='gene_identifier.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
