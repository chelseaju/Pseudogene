########################## lasso_preparation_v3.R #########################
# Function: read in data (expected Y) and training matrix X from different#
#         directories. Organize data in two files for lasso training      #
# Usage: R --no-save < lasso_analysis.R --args dir datatype rep           #
# Arguments: dir = input directory (select_one_pseudogene_110)            #
#            datatype = [gene|transcript]                                 #
#            rep = number of replicate used                               #
#            output = lasso_datatype_expected_Y.txt                       #
#                     lasso_datatype_matrix_X.txt                         #
# Author: Chelsea Ju                                                      #
# Date: 2014-05-04                                                        #
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
if(length(options) != 3){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_trainer.R --args dir type\n",
    "\t dir = directory of input and output\n",
    "\t type = genes or transcripts\n",
    "\t rep = number of replicates\n"),
    sep="");
}

dir <- options[1];
type <- options[2];
rep <- options[3];


if(as.integer(rep) > 54){
	stop(paste("Maximum replicates = 54"));
}

output_training <- paste("LassoTraining_rep", rep, sep="");
output_validation <- paste("LassoValidation_rep" , rep, sep="");

dir.create(file.path(dir, output_training), showWarnings = FALSE);
dir.create(file.path(dir, output_validation), showWarnings = FALSE);

## data preparation : training
depth <- c("5X", "7X", "10X", "13X", "15X", "17X", "20X", "23X", "27X");
expression <- c("4A", "6A", "8A", "R1A", "R2A", "R3A");

sample_size <- length(depth) * length(expression);
directory <- as.vector(outer(depth, expression, paste, sep="_101L_"));

select_dir <- directory[sample(1:sample_size, rep)];
rm(y);
rm(x);

## combind data
for (suffix in select_dir){

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
gene_count <- nrow(x) / length(select_dir);

if(gene_count %% 1 != 0){
	stop("Gene Count Error");
}

## convert normed_x as a [sample x obs] matrix (in this case, 16 * 3840)
##  3840 has 128 genes x 304 regions
expand_x <- matrix(0, nrow = length(directory) , ncol = gene_count * ncol(x))

for(i in 1:length(select_dir)){
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
mean_file <- paste(dir, "/", output_training, "/", type, "_matrix_A_Mean.txt", sep="");
median_file <- paste(dir, "/", output_training, "/", type, "_matrix_A_Median.txt", sep="");
var_file <- paste(dir, "/", output_training, "/", type, "_matrix_A_Variance.txt", sep="");

write.table(norm_x_mean, mean_file, sep="\t")
write.table(norm_x_median, median_file, sep="\t")
write.table(norm_x_var, var_file, sep="\t")

## output big matrix
lasso_x <- paste(dir, "/", output_training, "/", type, "_lasso_x_matrix.txt", sep="");
lasso_y <- paste(dir, "/", output_training, "/", type, "_lasso_y_vector.txt", sep="");

write.table(x, lasso_x, sep="\t");
write.table(y, lasso_y, sep="\t");

## data preparation : validation
validate_dir <- c("30X_101L_4A", "30X_101L_6A", "30X_101L_8A", "30X_101L_R1A", "30X_101L_R2A", "30X_101L_R3A");

## combind data
for (suffix in validate_dir){

	d <- paste(dir, "/", suffix, "/", sep="")

	observed <- read_distribution_matrix(d, type);
	expected <- read_expected_data(d, type, rownames(observed), suffix );
	observed <- rbind.fill(observed, x[1,]);
	observed[is.na(observed)] <- 0;
	observed <- observed[,colnames(x)];
	observed <- observed[-nrow(observed),];

	rownames(observed) <- rownames(expected);
	colnames(expected) <- c("Gene_Name", "Expected");

	## output validation data
	validation_x <- paste(dir, "/", output_validation, "/", suffix, "_validation_x_matrix.txt", sep="");
	validation_x_sum <- paste(dir, "/", output_validation, "/", suffix, "_validation_x_sum.txt", sep="");
	validation_y <- paste(dir, "/", output_validation, "/", suffix, "_validation_y_vector.txt", sep="");

	write.table(observed, validation_x, sep="\t");
	write.table(colSums(observed), validation_x_sum, sep="\t");
	write.table(expected, validation_y, sep="\t");

	print(paste("Writting Files to ", validation_x, sep=""));
	print(paste("Writting Files to ", validation_x_sum, sep=""));
	print(paste("Writting Files to ", validation_y, sep=""));

}


#output file \
flag_file <- paste(dir, "/finish_time.txt", sep="")
write.table("DONE", flag_file)







