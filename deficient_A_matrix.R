
# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 1){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < deficient_genes.R --args dir output threshold\n",
    "\t dir = directory for data input and output\n",
    sep=""));
}


dir <- options[1];
gene_file <- paste(dir, "/LassoValidation_rep54/problematic_genes.txt", sep="");
A_matrix_file <- paste(dir, "/LassoTraining_rep54/genes_matrix_A_Mean.txt", sep="");
validation_prefix <- c("30X_101L_4A", "30X_101L_6A", "30X_101L_8A", "30X_101L_R1A", "30X_101L_R2A", "30X_101L_R3A");


A_matrix <- read.table(A_matrix_file);
gene_list <- read.table(gene_file);

correlation_matrix <- matrix(0, nrow=nrow(gene_list), ncol=length(validation_prefix));
colnames(correlation_matrix) <- validation_prefix;
rownames(correlation_matrix) <- gene_list[,1];

for (i in 1:length(validation_prefix)){
	d <- validation_prefix[i];
	validation_file <- paste(dir, "/LassoValidation_rep54/", d, "_validation_x_matrix.txt", sep="");
	validation <- read.table(validation_file);

	for (j in 1:nrow(gene_list)){
		g <- gene_list[j,1];
		gene <- paste(d, g, sep="_");
		index <- which(rownames(validation) == gene);

		x_row <- validation[index,];
		a_row <- A_matrix[index,];

		correlation_matrix[j,i] <- cor(t(x_row), t(a_row));		
	}

}

outfile <- paste(dir, "/LassoValidation_rep54/problematic_genes_correlation.txt", sep="");
write.table(correlation_matrix, outfile, sep="\t");