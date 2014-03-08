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

#ENSEMBL_GENE = "/home/chelseaju/Database/EnsemblTranscriptome/Genes/ensemblGene_pseudogene.bed"
#ENSEMBL_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ensemblGene_pseudogene.bed"
ENSEMBL_GENE = "/u/scratch/c/chelseaj/database/EnsemblTranscriptome/Genes/ensemblGene_pseudogene.bed"

PRIORITY_GENE = {}

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
    Function: List of gene to be prioritize, serve as selecting overlapped genes
"""
def read_priority_gene(filename):
    fh = open(filename, 'rb')

    for line in fh:
        (gid, count) = line.split("\t")
        (ensp, ensg) = gid.split("_")
        PRIORITY_GENE[ensg] = ""

    fh.close()
"""
    Function: map the identified exon to known gene    
"""
def map_exon_to_gene(input_file):

    final_list = []        
    gene_list = {}
    unknown = []  
    unknown_collapse = []

    if(os.stat(input_file)[6]!=0):
        mapping = subprocess.check_output(["bedtools", "intersect", "-wb", "-loj", "-a", input_file, "-b", ENSEMBL_GENE])
        mapping_data = mapping.split("\n")
        
        for m in mapping_data:
            data = m.split("\t")
            if(len(data) > 2):
                mapped_gene_name = data[6]
                
                # found mapped gene
                if(mapped_gene_name != "-1" and mapped_gene_name != "."):
                    if(gene_list.has_key(mapped_gene_name)):
                        exist_mapped_gene = gene_list[mapped_gene_name]
                        gene_list[mapped_gene_name] = (exist_mapped_gene[0], exist_mapped_gene[1], min(exist_mapped_gene[2], int(data[1])), max(exist_mapped_gene[2], int(data[2])))

                    else:
                        gene_list[mapped_gene_name] = (int(data[4]), int(data[5]), int(data[1]), int(data[2]))  # (exon_start, exon_end, ensembl_start, ensembl_end)
                else:
                    unknown.append((data[0], int(data[1]), int(data[2]))) #(chr, exon_start, exon_end)


        # remove overlap genes : some genes have overlapped positions. this step is to select the best matched gene
        # how to deal with prioritize list
        previous_gene = (0,0,0,0,0)  #(name, start, end, coverage, prioritize flag)
        for k in sorted(gene_list.items(), key = lambda x: x[1][0]):
            coverage = (min(float(k[1][1]), float(k[1][3])) - max(float(k[1][0]), float(k[1][2]))) / (float(k[1][1]) - float(k[1][0]))
            priority = PRIORITY_GENE.has_key(k[0])

            # check for overlap
            if(k[1][0] >= previous_gene[1] and k[1][0] <= previous_gene[2]):    # since the list is sorted, only check the starting position against previous stored record
                
                # select the prioritize gene
                if(priority and previous_gene[4] == 0):
                    previous_gene = (k[0], k[1][0], k[1][1], coverage, 1)

                # replace with best coverage gene
                elif ((not priority and previous_gene[4] == 0) or (priority and previous_gene[4] == 1)):
                    if(coverage > previous_gene[3]):
                        previous_gene = (k[0], k[1][0], k[1][1], coverage, int(priority))
            else:
                if(previous_gene[3] > 0.01):   # remove extremely low coverage gene
                    final_list.append((previous_gene[0], previous_gene[1], previous_gene[2]))
                previous_gene = (k[0], k[1][0], k[1][1], coverage, int(priority))
            
        # store the last record
        if(previous_gene[3] > 0.01):
            final_list.append((previous_gene[0], previous_gene[1], previous_gene[2]))
      
     
        # collapse unknown region
        unknwon = sorted(unknown)    
        if(len(unknown) > 0):
            unknown_collapse.append(unknwon[0])

        for u in unknown:
            previous = unknown_collapse[-1]

            # new region if the gap is bigger than 80
            if(u[1] > previous[2] and (u[1] - previous[2]) > 80):
                unknown_collapse.append(u)
            else:
                unknown_collapse[-1] = (u[0], min(previous[1], u[1]), max(previous[2], u[2]))
     
        for u in unknown_collapse:
            name = "Unknown_" + str(u[0]) + "_" + str(u[1]) + "_" + str(u[2])
            final_list.append((name, u[1], u[2]))
            
    return final_list

"""
    Function: print the list of gene to file
"""
def export_genes(gene_list, out_file, chromosome_name):
    
    out_fh = open(out_file, 'w')
    for (name, start, end) in gene_list:
        out_fh.write("%s\t%s\t%d\t%d\n" %(name, chromosome_name, int(start), int(end)))    
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
    output_file = output_dir + chromosome_name + "_genes.txt"
    
    chr = convert_chromosome_name(chromosome_name)
    read_priority_gene(dir + "ENSG_expected_read_count.txt")
    mapped_regions = map_exon_to_gene(input_file)
    export_genes(mapped_regions, output_file, chromosome_name)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='gene_identifier.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
