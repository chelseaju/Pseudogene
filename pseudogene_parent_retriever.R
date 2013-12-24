#####################  pseudogene_parent_retriever.py ############################# 
# Usage: R --no-save --slave < pseudogene_parent_retriever.R                      #
#          --args parent_list output_prefix                                       #
# Input: a list of parent protein (ensembl protein ids)                           #
# Output: mRNA sequence of the parent proteins, and a list of transcript ids      #
# Function: use biomaRt library to retrieve the ensembl transcript ids, and the   #
#         mRNA sequence                                                           # 
# Author: Chelsea Ju                                                              #
# Date: 2013-07-22                                                                #
###################################################################################

# load library
library(biomaRt)

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 2){
	stop(paste("Invalid Arguments\n",
			   "Usage: R--no-save --slave < pseudogene_parent_retriever.R --args parent_list output_prefix\n",
			   "\t parent_list = a list of ensembl protein ids\n",
			   "\t output_prefix = prefix for the output file\n"),
			   sep="");
}

in_file <- options[1];
out_prefix <- options[2];
fasta_file <- paste(out_prefix,"_parents.fa", sep="");
transcript_file <- paste(out_prefix, "_parents.transcript", sep="");

# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# read in ids
cat(paste("Reading IDs from ", in_file , "..............\n", sep=""));
ensembl_ids <- read.table(in_file, sep="\n");

# query biomart
results <- getBM(attributes = c("ensembl_transcript_id", "ensembl_peptide_id", "cdna"),
                 filters = "ensembl_peptide_id", values = ensembl_ids$V1,
                 mart = mart)

# export to file
apply(results, 1, function(x){ cat(paste(">",x[3],"_",x[2],"\n",x[1],"\n",sep=""), file = fasta_file, append = TRUE)
							}
	)

cat(paste("Write to file ", fasta_file, "\n", sep=""));

cat(paste("\n\nDone!!", sep=""));
				

                 
                 
                 
                 
                 
                 
                 
