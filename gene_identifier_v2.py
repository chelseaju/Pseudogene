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


#ENSEMBL_GENE = "/home/chelseaju/Database/Ensembl/ENST_74.bed"
#PARENT_GENE = "/home/chelseaju/Database/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE = "/home/chelseaju/Database/Pseudogene/Pseudogene_74.bed"

ENSEMBL_GENE = "/u/home/c/chelseaj/project/database/Ensembl/ENSG_74.bed"
PARENT_GENE =  "/u/home/c/chelseaj/project/database/Pseudogene/Parent_ENST_74.bed"
PSEUDO_GENE =  "/u/home/c/chelseaj/project/database/Pseudogene/Pseudogene_74.bed"

#ENSEMBL_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ENST_74.bed"
#PARENT_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Pseudogene_74.bed"


"""
    make the chromosome name readable for database
"""
def convert_chromosome_name(name):
    
    name = name.lower()
    match_chr_number = re.match(r"\D*(\d+)", name)
    match_chr_X = re.match(r"\D*(x)", name)
    match_chr_Y = re.match(r"\D*(y)", name)
    match_chr_MT = re.match(r"\D*(mt)", name)
    
    if(match_chr_number):
        return match_chr_number.group(1)
    
    elif(match_chr_X):
        return "X"

    elif(match_chr_Y):
        return "Y"

    elif(match_chr_MT):
        return "MT"

    else:
        print ("unknown input %s", name)
        sys.exit(2)
  
"""
    Function: map the identified exon to known gene    
"""
def map_exon_to_gene(input_file, database):

    final_list = []        
    gene_list = {}
    unknown = []  

    if(os.stat(input_file)[6]!=0):
        mapping = subprocess.check_output(["bedtools", "intersect", "-wb", "-loj", "-a", input_file, "-b", database])
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
                    unknown.append("%s\t%s\t%s\n" %(data[0], data[1], data[2])) # keep the unmapped read from bam


        # remove overlap genes : some genes have overlapped positions. this step is to select the best matched gene
        previous_gene = (0,0,0,0)  #(name, start, end, coverage)
        for k in sorted(gene_list.items(), key = lambda x: x[1][0]):
            coverage = (min(float(k[1][1]), float(k[1][3])) - max(float(k[1][0]), float(k[1][2]))) / (float(k[1][1]) - float(k[1][0]))

            # check for overlap
            if(k[1][0] >= previous_gene[1] and k[1][0] <= previous_gene[2]):    # since the list is sorted, only check the starting position against previous stored record
                
                # replace with best coverage gene
                if(coverage > previous_gene[3]):
                    previous_gene = (k[0], k[1][0], k[1][1], coverage)
            else:
                if(previous_gene[3] > 0):   # remove extremely low coverage gene or zero entries
                    final_list.append((previous_gene[0], previous_gene[1], previous_gene[2]))

                previous_gene = (k[0], k[1][0], k[1][1], coverage)
            
        # store the last record
        if(previous_gene[3] > 0):
           final_list.append((previous_gene[0], previous_gene[1], previous_gene[2]))
           
            
    return (final_list, unknown)


"""
    Function : Collapse unknown region
"""
def collapse_unknown_region(unknown):

    final_list = []
    unknown_collapse = (0,0,0,0)  #(chr, start, end, read count)

    # collapse unknown region
    unknwon = sorted(unknown)

    for u_data in unknown:
        u_data = u_data.rstrip()
        u = u_data.split("\t")

        # new region if the gap is bigger than 80
        if(int(u[1]) > int(unknown_collapse[2]) and (int(u[1]) - int(unknown_collapse[2]) > 80)):

            # export unknown region
            if(unknown_collapse[3] > 0):
                name = "Unknown_" + str(unknown_collapse[0]) + "_" + str(unknown_collapse[1]) + "_" + str(unknown_collapse[2])
                final_list.append((name, unknown_collapse[1], unknown_collapse[2]))

            # new region
            unknown_collapse = (u[0], u[1], u[2], 1)

        else:
            unknown_collapse = (u[0], min(u[1], unknown_collapse[1]), max(u[2], unknown_collapse[2]), unknown_collapse[3]+1)


    # export the last one
    if(unknown_collapse[3] > 0):
        name = "Unknown_" + str(unknown_collapse[0]) + "_" + str(unknown_collapse[1]) + "_" + str(unknown_collapse[2])
        final_list.append((name, unknown_collapse[1], unknown_collapse[2]))

    return final_list
"""
    Function : Export unmatched reads to temporoary file
"""
def export_temp_unknown(unknown_list, unknown_file):

    fh = open(unknown_file, 'w')
    for u in unknown_list:
        fh.write(u)
    fh.close()

    print ""
    print "Write to Temporary File : %s" %(unknown_file)

"""
    Function: print the list of gene to file
"""
def export_genes(gene_list, out_file, chromosome_name, append):
    
    if(append):
        out_fh = open(out_file, 'a')
    else:
        out_fh = open(out_file, 'w')

    for (name, start, end) in gene_list:
        out_fh.write("%s\t%s\t%d\t%d\n" %(name, chromosome_name, int(start), int(end)))    
    out_fh.close()

       
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
    output_temp = output_dir + chromosome_name + "_unknown.txt"
    
    chr = convert_chromosome_name(chromosome_name)

    # map to parent genes
    (mapped_parents, unknown) = map_exon_to_gene(input_file, PARENT_GENE)
    export_temp_unknown(unknown, output_temp)
    export_genes(mapped_parents, output_file, chromosome_name, False)

    # map to pseudo genes
    (mapped_pseudo, unknown) = map_exon_to_gene(output_temp, PSEUDO_GENE)
    export_temp_unknown(unknown, output_temp)
    export_genes(mapped_pseudo, output_file, chromosome_name, True)

    # map to all genes
    (mapped_ensembl, unknown) = map_exon_to_gene(output_temp, ENSEMBL_GENE)
    export_genes(mapped_ensembl, output_file, chromosome_name, True)

    # collapse unknown
    grouped_unknown = collapse_unknown_region(unknown)
    export_genes(grouped_unknown, output_file, chromosome_name, True)

    # remove temporary file
    os.system("rm %s" %(output_temp))

    print ""
    print "Remove Temporary File : %s" %(output_temp)

    print ""
    print "Writing Gene List to File : %s" %(output_file)



if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='gene_identifier_v2.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
