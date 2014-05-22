########################## lasso_preparation_v2.R #########################
# Function: read in data (expected Y) and training matrix X from different#
#         directories. Organize data in two files for lasso training      #
# Usage: R --no-save < lasso_analysis.R --args dir datatype abund coverage#
# Arguments: dir = input directory (select_one_pseudogene_110)            #
#            datatype = [gene|transcript]                                 #
#            abundance = specific abundance (1A or 5A or 10A or all)      #
#            coverage = specific coverage (3X, 5X, 7X, 10X, 13X, 17X, 20X,#
#                                          23X, 27X, 30X, or all)         #
#            output = matrix_A_mean.txt                                   #
#                     matrix_A_median.txt                                 #
#                     matrix_A_variance.txt                               #
#                     variance_plot.png                                   #
# Author: Chelsea Ju                                                      #
# Date: 2014-04-10                                                        #
###########################################################################


library(plyr)

## self-defined function
read_distribution_matrix <- function(subdir, filetype){
	matrix_file <- paste(subdir, "tophat_out/", filetype, "_distribution.matrix", sep="");
	matrix_data <- read.table(matrix_file);
	matrix_data;
}

read_expected_data <- function(subdir, filetype, row_order, prefix){
	expected_file <- paste(subdir, "tophat_out/", filetype, "_expected_read_count.txt", sep="");
	expected_data <- read.table(expected_file, header = FALSE);
	c1 <- expected_data$V1;
	row_names <- unlist(strsplit(as.character(c1), "_")); ## modify the rownames: split ENSG00000100478_ENST00000542754 into "ENSG00000100478" "ENST00000542754"
#	index_vector <- c(1:length(row_names));
#	row_names <- row_names[which(index_vector %% 2 == 0)];
	
	rownames(expected_data) <- row_names;
    expected_data <- expected_data[row_order,];
    rownames(expected_data) <- paste(prefix, rownames(expected_data), sep="_")
	expected_data;
}


# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 4){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < data_examination_v2.R --args dir type\n",
    "\t dir = directory of input and output\n",
    "\t type = genes or transcripts\n",
    "\t abundance = 1A or 5A or 10A or all\n",
    "\t coverage = 3X or 5X or 7X or 10X or 13X or 17X or 20X or 23X or 27X or 30X or all\n",
    sep=""));
}

dir <- options[1];
type <- options[2];
abundance <- options[3];
coverage <- options[4];

## for percentage matrix and distribution matrix

depth <- c("5X", "7X", "10X", "13X", "15X", "17X", "20X", "23X", "27X", "30X");
expression <- c("4A", "6A", "8A", "R1A", "R2A", "R3A");


if(abundance == "all" && coverage == "all"){
	sample_size <- length(depth) * length(expression);
	directory <- as.vector(outer(depth, expression, paste, sep="_101L_"))

} else if(abundance == "all"){
	sample_size <- length(abundance);
	directory <- as.vector(outer(coverage, expression, paste, sep="_101L_"))

} else if(coverage == "all"){
	sample_size <- length(coverage);
	directory <- as.vector(outer(depth, abundance, paste, sep="_101L_"))
} else{
	stop(paste("Needs more matrices\n"))
}

rm(y);
rm(x);

## combind data
for (suffix in directory){

	d <- paste(dir, "/", suffix, "/", sep="")

	observed <- read_distribution_matrix(d, type);
	expected <- read_expected_data(d, type, rownames(observed), suffix )

	if(exists("y")){
		y <- rbind(y, expected);
	}else{
		y <- expected;
	}

	if(exists("x")){
		x <- rbind.fill(x, observed);
	}else{
		x <- observed;
	}

}

x_colname <- colnames(x);
x[is.na(x)] <- 0;

colnames(y) <- c("Gene_Name", "Read_Count");
rownames(x) <- rownames(y);


## normalize observation matrix against expected value
normed_x <- x / y$Read_Count;

## number of genes
gene_count <- nrow(x) / length(directory)

if(gene_count %% 1 != 0){
	stop("Gene Count Error");
}

## convert normed_x as a [sample x obs] matrix (in this case, 16 * 3840)
##  3840 has 128 genes x 304 regions
expand_x <- matrix(0, nrow = length(directory) , ncol = gene_count * ncol(x))

for(i in 1:length(directory)){
	start <- ((i - 1) * gene_count) + 1
	end <- start + gene_count - 1
	replicates <- unlist(normed_x[start:end,])
	expand_x[i,] <- replicates
}


norm_x_mean <- matrix(apply(expand_x, 2, mean), byrow = F, nrow = gene_count)
norm_x_median <- matrix(apply(expand_x, 2, median), byrow = F, nrow = gene_count)
norm_x_var <- matrix(apply(expand_x, 2, var), byrow = F, nrow = gene_count)

colnames(norm_x_mean) <- colnames(x)
colnames(norm_x_median) <- colnames(x)
colnames(norm_x_var) <- colnames(x)
rownames(norm_x_mean) <- rownames(x)[1:gene_count]
rownames(norm_x_median) <- rownames(x)[1:gene_count]
rownames(norm_x_var) <- rownames(x)[1:gene_count]

## output mean, var, median
output_subdir <- paste("observation_", coverage, "_", abundance, sep="")
dir.create(file.path(dir, output_subdir), showWarnings = FALSE)

mean_file <- paste(dir, "/", output_subdir, "/", type, "_matrix_A_Mean.txt", sep="")
median_file <- paste(dir, "/", output_subdir, "/", type, "_matrix_A_Median.txt", sep="")
var_file <- paste(dir, "/", output_subdir, "/", type, "_matrix_A_Variance.txt", sep="")

write.table(norm_x_mean, mean_file, sep="\t")
write.table(norm_x_median, median_file, sep="\t")
write.table(norm_x_var, var_file, sep="\t")

print(paste("Maximum variance = ", max(norm_x_var), "\n", sep=""));
print(paste("Minimum variance = ", min(norm_x_var), "\n", sep=""));
print(paste("Number of variance > 0.01 = ", length(which(norm_x_var > 0.01)), "\n", sep=""));


## plot histogram for those metabolites have high variance
large_variance <- which(norm_x_var > 0.01)

#png(paste(dir, "/", output_subdir, "/", type, "_variance_plot.png", sep=""), width=880, height=800)
#op <- par(mfrow = c(ceiling(sqrt(length(large_variance))), ceiling(sqrt(length(large_variance)))))

for (i in 1:length(large_variance)){
	index <- large_variance[i]

	r_index <- index %% gene_count
	c_index <- ceiling(index / gene_count)
	if(r_index == 0){
		r_index = gene_count
	}
	r_name <- colnames(x)[r_index]
	c_name <- colnames(x)[c_index]

	hist_title <- paste(r_name, "-", c_name, sep="")
	
	print(paste(hist_title, norm_x_var[index], sep="::"))

	hist(expand_x[, index], xlab = "Correlation", 
			main = hist_title, col = "lightgreen", cex.main = 1.5)
	curve(dnorm(x, mean=norm_x_mean[index], sd=sd(expand_x[, index])),
		 add=TRUE, col="darkblue", lwd=2) 

}
#dev.off()


