""" transcript_identifier
Usage: python transcript_identifier_v2.py -c chromosome -d directory 
Input:  -c chromosome name
        -d input/output directory 
Output: 
Function: 1. Map each exon position to the database to find a corresponding transcript or pseudogene
          2. If the exon does not have a match, return this exon
          4. Return the mapped boundary of each unmarked exon

Date: 2014-04-24
Author: Chelsea Ju
"""
import sys, re, os, subprocess, random, argparse, datetime

# LAB
#DB = "/home/chelseaju/Database/Pseudogene/ENST_Pseudogene_74.bed"
#DB = "/home/chelseaju/Database/Pseudogene/ParentENST_Pseudogene_74.bed"
#ENSEMBL_GENE = "/home/chelseaju/Database/Ensembl/ENST_74.bed"
#PARENT_GENE = "/home/chelseaju/Database/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE = "/home/chelseaju/Database/Pseudogene/Pseudogene_74.bed"

# HOFFMAN
DB = "/u/home/c/chelseaj/project/database/Pseudogene/ParentENST_Pseudogene_74.bed"
ENSEMBL_GENE = "/u/home/c/chelseaj/project/database/Ensembl/ENST_74.bed"
#PARENT_GENE = "/u/home/c/chelseaj/project/database/Pseudogene/Parent_ENSG_74.bed"
#PSEUDO_GENE =  "/u/home/c/chelseaj/project/database/Pseudogene/Pseudogene_74.bed"

# MAC
#DB = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ParentENST_Pseudogene_74.bed"
#ENSEMBL_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ENST_74.bed"
#PARENT_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Pseudogene_74.bed"


"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))



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
    Function: map the identified exon to known transcript    
"""
def map_exon_to_gene(input_file, database):
    
    final_list = []        
    gene_list = {}
    unknown = []  

    if(os.stat(input_file)[6]!=0):
        mapping = subprocess.check_output(["bedtools", "intersect", "-wb", "-loj", "-a", input_file, "-b", database])
        mapping_data = mapping.split("\n")
    
        ## retrieve gene annotation based on overlapped region    
        for m in mapping_data:
            data = m.split("\t")
            # data = [a_chr, a_start, a_end, a_name, b_chr, b_start, b_end, b_name ....]

            if(len(data) > 2):
                read_name = data[3]
                mapped_name = data[7]

                # found mapped gene
                if(mapped_name != "-1" and mapped_name != "."):
                    mapped_gene_name = read_name + "_" + mapped_name
                    if(gene_list.has_key(mapped_gene_name)):
                        exist_mapped_gene = gene_list[mapped_gene_name]
                        gene_list[mapped_gene_name] = (exist_mapped_gene[0], exist_mapped_gene[1], exist_mapped_gene[2], min(exist_mapped_gene[3], int(data[1])), max(exist_mapped_gene[4], int(data[2])))
 
                    else:
                        gene_list[mapped_gene_name] = (read_name, int(data[5]), int(data[6]), int(data[1]), int(data[2]))  # (gene_name, ensembl_start, ensembl_end, exon_start, exon_end)
                else:
                    unknown.append("%s\t%s\t%s\t%s\n" %(data[0], data[1], data[2], data[3])) # keep the unmapped read from exon.bed

        ## remove overlap genes : some genes have overlapped positions. this step is to select the best matched gene
        # best matched gene is defined as 
        # 1. mapped_name matches with gene_name
        # 2. best coverage among all mapped genes

        previous_gene = ("","",0,0,0,0,0)  #(gene_name, mapped_name ensembl_start, ensembl_end, exon_start, exon_end, coverage)
        for k in sorted(gene_list.items(), key = lambda x: (x[1][0], x[1][1])):

            # k = (GENENAME_MAPPEDNAME => (gene_name, ensembl_start, ensembl_end, exon_start, exon_end))
            (gene_name, mapped_name) = k[0].split("_")
            coverage = (min(float(k[1][2]), float(k[1][4])) - max(float(k[1][1]), float(k[1][3]))) / max((float(k[1][2]) - float(k[1][1])), float(k[1][4]) - float(k[1][3]))

            # check for overlap - same gene name, overlapped positions
            if(gene_name == previous_gene[0] and k[1][1] >= previous_gene[2] and k[1][1] <= previous_gene[3]): # since the list is sorted, only check the starting position against previous stored record
                
                # replace previous gene only if the names do not match and coverage is better than previous gene
                if(previous_gene[0] != previous_gene[1] and (gene_name == mapped_name or coverage > previous_gene[6])):
                    previous_gene = (gene_name, mapped_name, min(previous_gene[2], k[1][1]), k[1][2], k[1][3], k[1][4], coverage)

            else:
                #output previous gene if not ""
                if(previous_gene[0] != ""):
#                    final_list.append((previous_gene[0], previous_gene[1], previous_gene[4], previous_gene[5]))
                    final_list.append((previous_gene[0], previous_gene[1], previous_gene[2], previous_gene[3]))

                previous_gene = (gene_name, mapped_name, k[1][1], k[1][2], k[1][3], k[1][4], coverage)


        # put last item into list
        if(previous_gene[0] != ""):
#            final_list.append((previous_gene[0], previous_gene[1], previous_gene[4], previous_gene[5]))
            final_list.append((previous_gene[0], previous_gene[1], previous_gene[2], previous_gene[3]))

    return (final_list, unknown)


"""
    Function : Collapse unknown region
"""
def collapse_unknown_region(unknown):

    final_list = []
    unknown_collapse = (0,0,0,0)  #(chr, start, end, gene_name)

    # collapse unknown region
    unknwon = sorted(unknown, key=lambda x:(x[3], x[1]))

    for u_data in unknown:
        u_data = u_data.rstrip()
        u = u_data.split("\t") # [chr, start, end, gene_name]

        # collapse region if two regions overlapped and are from the same gene
        if(u[3] == unknown_collapse[3] and int(u[1]) > int(unknown_collapse[1]) and int(u[1]) <= int(unknown_collapse[2])):
            unknown_collapse = (u[0], min(u[1], unknown_collapse[1]), max(u[2], unknown_collapse[2]), u[3])

        else:
            ## output collapsed region
            if(unknown_collapse[0] != 0):
                name = "Unknown_" + str(unknown_collapse[0]) + "_" + str(unknown_collapse[1]) + "_" + str(unknown_collapse[2])
                final_list.append((unknown_collapse[3], name, int(unknown_collapse[1]), int(unknown_collapse[2])))

            unknown_collapse = (u[0],u[1],u[2],u[3])

    # export the last one
    if(unknown_collapse[0] != 0):
        name = "Unknown_" + str(unknown_collapse[0]) + "_" + str(unknown_collapse[1]) + "_" + str(unknown_collapse[2])
        final_list.append((unknown_collapse[3], name, int(unknown_collapse[1]), int(unknown_collapse[2])))

    return final_list


"""
    Function : Export unmatched reads to temporoary file
"""
def export_temp_unknown(unknown_list, unknown_file):

    fh = open(unknown_file, 'w')
    for u in unknown_list:
        fh.write(u)
    fh.close()

    echo("Write to Temporary File")

"""
    Function: print the list of gene to file
"""
def export_genes(gene_list, out_file, chromosome_name, append):
    
    if(append):
        out_fh = open(out_file, 'a')
    else:
        out_fh = open(out_file, 'w')

    for (gene_name, mapped_name, start, end) in gene_list:
        out_fh.write("%s\t%s\t%s\t%d\t%d\n" %(gene_name, mapped_name, chromosome_name, start, end))    
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
    output_file = output_dir + chromosome_name + "_transcripts.txt"
    output_temp = output_dir + chromosome_name + "_unknown.txt"
    
    chr = convert_chromosome_name(chromosome_name)

    # map to parents and pseudogene
    (mapped_parents, unknown) = map_exon_to_gene(input_file, DB)
    export_temp_unknown(unknown, output_temp)
    export_genes(mapped_parents, output_file, chromosome_name, False)

    # map to all genes
    (mapped_ensembl, unknown) = map_exon_to_gene(output_temp, ENSEMBL_GENE)
    export_genes(mapped_ensembl, output_file, chromosome_name, True)


    # collapse unknown
    grouped_unknown = collapse_unknown_region(unknown)
    export_genes(grouped_unknown, output_file, chromosome_name, True)

    # remove temporary file
    os.system("rm %s" %(output_temp))
    echo("Remove Temporary File")
    echo("Write Gene to File : %s" %(output_file))


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='transcript_identifier_v2.py')
    parser.add_argument("-c", "--chromosome", dest="chromosome", type=str, help="chromosome name, ex chr1", required = True)
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input and output files", required = True)

    main(parser)
