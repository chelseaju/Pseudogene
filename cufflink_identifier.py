""" cufflink_identifier.py
Usage: python cufflink_identifier.py -d cufflink_directory 
Input:  -d input/output directory (select_pseudogene_128_v2/10X_100L_1A_30X_100L_1A)
Output: line of cufflink output matched to gene of interests {XLOC_XXX \t ENSG_XXX}
Function: 1. intersect a list of region identified by cufflinks with the annotation file
          2. output the corresponding gene name

Date: 2014-03-30
Author: Chelsea Ju
"""

import sys, re, os, subprocess, random, argparse

# LAB
DB = "/home/chelseaju/Database/Pseudogene/ParentENSG_Pseudogene_74.bed"
ENSEMBL_GENE="/home/chelseaju/Database/Ensembl/ENSG_74.bed"

# HOFFMAN
#ENSEMBL_GENE = "/u/home/c/chelseaj/project/database/Ensembl/ENSG_74.bed"
#PARENT_GENE =  "/u/home/c/chelseaj/project/database/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE =  "/u/home/c/chelseaj/project/database/Pseudogene/Pseudogene_74.bed"

# MAC
#ENSEMBL_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ENST_74.bed"
#PARENT_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Parent_ENST_74.bed"
#PSEUDO_GENE = "/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Pseudogene_74.bed"

"""
	Function: convert the cufflink ids to bed format
"""
def cufflink2bed(input, output):

	cufflinks = open(input, 'rb')
	bed = open(output, 'w')

	line = cufflinks.readline()
	for line in cufflinks:
		data = line.split("\t")
		xloc = data[0]
		(chromosome, position) = data[6].split(":")
		(start, end)  = position.split("-")

		bed.write("%s\t%s\t%s\t%s\n" %(chromosome, start, end, xloc))

	cufflinks.close()
	bed.close()
	

	print ""
	print "Write to %s\n" %(output)

"""
    Function: map the cufflink ids to known gene    
"""
def map_cufflink_to_gene(fpkm, database):

    gene_list = {}
    unknown = []  

    if(os.stat(fpkm)[6]!=0):
        mapping = subprocess.check_output(["bedtools", "intersect", "-wa", "-loj", "-a", fpkm, "-b", database])
        mapping_data = mapping.split("\n")

        for m in mapping_data:
            data = m.split("\t")
            if(len(data) > 2):
            	cufflink_name = data[3]
                mapped_gene_name = data[7]

                # found mapped gene
                if(mapped_gene_name != "-1" and mapped_gene_name != "."):

                	chrom = data[0]
                	c_start = int(data[1])
                	c_end = int(data[2])
                	g_start = int(data[5])
                	g_end = int(data[6])
                	coverage = float(min(c_end, g_end) - max(c_start, g_start)) / max(float(g_end - g_start), float(c_end - c_start))

                	if(gene_list.has_key(cufflink_name)):
                		exist_mapped_gene = gene_list[cufflink_name]
                		exist_gene_name = exist_mapped_gene[0]
                		exist_coverage = exist_mapped_gene[1]
                		exist_list = exist_mapped_gene[2]
                		exist_list.append(mapped_gene_name)

                		# select the gene with best overlap
                		if(coverage > exist_coverage):
                			exist_coverage = coverage
                			exist_gene_name = mapped_gene_name

                	else:
						gene_list[cufflink_name] = (mapped_gene_name, coverage, [mapped_gene_name]) # (ENSG___, coverage, [list of possible matches])

                else:
                    unknown.append("%s\t%s\t%s\t%s\n" %(chrom, c_start, c_end, cufflink_name)) # keep the unmapped cufflink id

    return (gene_list, unknown)

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
def export_genes(gene_list, out_file, append):
    
    if(append):
        out_fh = open(out_file, 'a')
    else:
        out_fh = open(out_file, 'w')

   	for k in gene_list.keys():
   		(gene_name, coverage, alternatives) = gene_list[k]
   		out_fh.write("%s\t%s\t%s\t%s\n" %(k, gene_name, coverage, ",".join(alternatives) ))

   	out_fh.close()



def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

	gene_fpkm = dir + "genes.fpkm_tracking"

	gene_bed = dir + "cufflinks.bed"
	output_temp = dir + "unknown.bed"

	outfile = dir + "cufflink_genes.txt"

	# convert the cufflinks IDs to bed format
	cufflink2bed(gene_fpkm, gene_bed)

	# map the cufflinks IDs (XLOC_) to Parent and Pseudogene
    (mapped_genes, unknown) = map_cufflink_to_gene(gene_bed, DB)
    export_temp_unknown(unknown, output_temp)
    export_genes(mapped_genes, outfile, False)

	# map the cufflinks IDs (XLOC_) to Ensembl
    (mapped_ensembl, unknown) = map_cufflink_to_gene(output_temp, ENSEMBL_GENE)
    export_genes(mapped_ensembl, outfile, True)

    # remove temporary file
    os.system("rm %s" %(output_temp))

    print ""
    print "Remove Temporary File : %s" %(output_temp)

    print ""
    print "Writing Gene List to File : %s" %(outfile)

if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='cufflink_identifier.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="cufflink directory of input and output files", required = True)

    main(parser)
