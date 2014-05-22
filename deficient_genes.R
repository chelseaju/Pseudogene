

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < deficient_genes.R --args dir output threshold\n",
    "\t dir = directory for data input and output\n",
    "\t output = output file name\n",
    "\t threshold = similarity threshold\n",
    sep=""));
}


dir <- options[1];
output <- options[2];
threshold <- options[3];

output <- paste(dir , "/" , output, sep="");
write("", output);

file <- paste(dir, "/LassoTraining_rep54/genes_matrix_A_Mean.txt", sep="");
A_matrix <- read.table(file);

A_correlation <- cor(t(A_matrix));
group <- rep(0, ncol(A_correlation));

for (i in 1:nrow(A_correlation)){
	if(group[i] ==0){
		group[which(A_correlation[i, ] > threshold)] <- i;
	}
}

multi_groups <- unique(group[duplicated(group)]);

for (j in multi_groups){
	write(paste("Group ", j), output, append = T);
	write(t(colnames(A_correlation)[group==j]), output, append = T, sep="\t");
}

print(paste("Write to file", output))