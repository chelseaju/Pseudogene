########################## lasso_preparation_v2.R #########################
# Function: read in data (expected Y) and training matrix X from different#
#         directories. Organize data in two files for lasso training      #
# Usage: R --no-save < lasso_analysis.R --args dir datatype rep           #
# Arguments: dir = input directory (select_one_pseudogene_110)            #
#            datatype = [gene|transcript]                                 #
#            rep = number of replicate used                               #
#            output = lasso_datatype_expected_Y.txt                       #
#                     lasso_datatype_matrix_X.txt                         #
# Author: Chelsea Ju                                                      #
# Date: 2014-04-01                                                        #
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
if(length(options) != 2){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_trainer.R --args dir type\n",
    "\t dir = directory of input and output\n",
    "\t type = genes or transcripts\n",
    sep=""));
}

dir <- options[1];
type <- options[2];

## for percentage matrix and distribution matrix
dir_3XR1A <- paste(dir, "/3X_101L_R1A/", sep="")
dir_5XR1A <- paste(dir, "/5X_101L_R1A/", sep="")
dir_7XR1A <- paste(dir, "/7X_101L_R1A/", sep="")
dir_10XR1A <- paste(dir, "/10X_101L_R1A/", sep="")
dir_13XR1A <- paste(dir, "/13X_101L_R1A/", sep="")
dir_17XR1A <- paste(dir, "/17X_101L_R1A/", sep="")
dir_20XR1A <- paste(dir, "/20X_101L_R1A/", sep="")
dir_23XR1A <- paste(dir, "/23X_101L_R1A/", sep="")
dir_27XR1A <- paste(dir, "/27X_101L_R1A/", sep="")
dir_30XR1A <- paste(dir, "/30X_101L_R1A/", sep="")
dir_10X1A <- paste(dir, "/10X_101L_1A/", sep="") 
dir_10XR2A <- paste(dir, "/10X_101L_R2A/", sep="") 
dir_20X1A <- paste(dir, "/20X_101L_1A/", sep="") 
dir_20XR2A <- paste(dir, "/20X_101L_R2A/", sep="") 
dir_30X1A <- paste(dir, "/30X_101L_1A/", sep="") 
dir_30XR2A <- paste(dir, "/30X_101L_R2A/", sep="") 

# validation purpose
dir_10XR3A <- paste(dir, "/10X_101L_R3A/", sep="") 
dir_20XR3A <- paste(dir, "/20X_101L_R3A/", sep="") 
dir_30XR3A <- paste(dir, "/30X_101L_R3A/", sep="") 


observed_3XR1A <- read_distribution_matrix(dir_3XR1A, type);
observed_5XR1A <- read_distribution_matrix(dir_5XR1A, type);
observed_7XR1A <- read_distribution_matrix(dir_7XR1A, type);
observed_10XR1A <- read_distribution_matrix(dir_10XR1A, type);
observed_13XR1A <- read_distribution_matrix(dir_13XR1A, type);
observed_17XR1A <- read_distribution_matrix(dir_17XR1A, type);
observed_20XR1A <- read_distribution_matrix(dir_20XR1A, type);
observed_23XR1A <- read_distribution_matrix(dir_23XR1A, type);
observed_27XR1A <- read_distribution_matrix(dir_27XR1A, type);
observed_30XR1A <- read_distribution_matrix(dir_30XR1A, type);
observed_10X1A <- read_distribution_matrix(dir_10X1A, type);
observed_10XR2A <- read_distribution_matrix(dir_10XR2A, type);
observed_20X1A <- read_distribution_matrix(dir_20X1A, type);
observed_20XR2A <- read_distribution_matrix(dir_20XR2A, type);
observed_30X1A <- read_distribution_matrix(dir_30X1A, type);
observed_30XR2A <- read_distribution_matrix(dir_30XR2A, type);

# for validation
observed_10XR3A <- read_distribution_matrix(dir_10XR3A, type);
observed_20XR3A <- read_distribution_matrix(dir_20XR3A, type);
observed_30XR3A <- read_distribution_matrix(dir_30XR3A, type);


## expected vector
expected_3XR1A <- read_expected_data(dir_3XR1A, type, rownames(observed_3XR1A), "3XR1A");
expected_5XR1A <- read_expected_data(dir_5XR1A, type, rownames(observed_5XR1A), "5XR1A");
expected_7XR1A <- read_expected_data(dir_7XR1A, type, rownames(observed_7XR1A), "7XR1A");
expected_10XR1A <- read_expected_data(dir_10XR1A, type, rownames(observed_10XR1A), "10XR1A");
expected_13XR1A <- read_expected_data(dir_13XR1A, type, rownames(observed_13XR1A), "13XR1A");
expected_17XR1A <- read_expected_data(dir_17XR1A, type, rownames(observed_17XR1A), "17XR1A");
expected_20XR1A <- read_expected_data(dir_20XR1A, type, rownames(observed_20XR1A), "20XR1A");
expected_23XR1A <- read_expected_data(dir_23XR1A, type, rownames(observed_23XR1A), "23XR1A");
expected_27XR1A <- read_expected_data(dir_27XR1A, type, rownames(observed_27XR1A), "27XR1A");
expected_30XR1A <- read_expected_data(dir_30XR1A, type, rownames(observed_30XR1A), "30XR1A");
expected_10X1A<- read_expected_data(dir_10X1A, type, rownames(observed_10X1A), "10X1A");
expected_10XR2A <- read_expected_data(dir_10XR2A, type, rownames(observed_10XR2A), "10XR2A");
expected_20X1A <- read_expected_data(dir_20X1A, type, rownames(observed_20X1A), "20X1A");
expected_20XR2A <- read_expected_data(dir_20XR2A, type, rownames(observed_20XR2A), "20XR2A");
expected_30X1A <- read_expected_data(dir_30X1A, type, rownames(observed_30X1A), "30X1A");
expected_30XR2A <- read_expected_data(dir_30XR2A, type, rownames(observed_30XR2A), "30XR2A");


## for validation 
expected_10XR3A <- read_expected_data(dir_10XR3A, type, rownames(observed_10XR3A), "10XR3A");
expected_20XR3A <- read_expected_data(dir_20XR3A, type, rownames(observed_20XR3A), "20XR3A");
expected_30XR3A <- read_expected_data(dir_30XR3A, type, rownames(observed_20XR3A), "20XR3A");


y <- rbind(expected_3XR1A,expected_5XR1A,expected_7XR1A,expected_10XR1A,expected_13XR1A,expected_17XR1A,expected_20XR1A,expected_23XR1A,expected_27XR1A,expected_30XR1A,
	expected_10X1A,expected_10XR2A, expected_20X1A, expected_20XR2A, expected_30X1A, expected_30XR2A);

x <- rbind.fill(observed_3XR1A,observed_5XR1A,observed_7XR1A,observed_10XR1A,observed_13XR1A,observed_17XR1A,observed_20XR1A,observed_23XR1A,observed_27XR1A,observed_30XR1A,
        observed_10X1A,observed_10XR2A, observed_20X1A, observed_20XR2A, observed_30X1A, observed_30XR2A);

x_colname <- colnames(x);
x[is.na(x)] <- 0;

colnames(y) <- c("Gene_Name", "Read_Count");
rownames(x) <- rownames(y);


## normalize observation matrix against expected value
normed_x <- x / y$Read_Count;

## number of genes
gene_count <- nrow(x) / 16

## convert normed_x as a [sample x obs] matrix (in this case, 16 * 3840)
##  3840 has 128 genes x 304 regions
expand_x <- matrix(0, nrow = 16, ncol = gene_count * ncol(x))

for(i in 1:16){
	start <- ((i - 1) * gene_count) + 1
	end <- start + gene_count - 1
	replicates <- unlist(normed_x[start:end,])
	expand_x[i,] <- replicates
}


norm_x_mean <- matrix(apply(expand_x, 2, mean), byrow = F, nrow = gene_count)
norm_x_median <- matrix(apply(expand_x, 2, median), byrow = F, nrow = gene_count)
norm_x_var <- matrix(apply(expand_x, 2, var), byrow = F, nrow = gene_count)


## output mean, var, median
mean_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_Mean.txt", sep="")
median_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_Median.txt", sep="")
var_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_Variance.txt", sep="")

write.table(norm_x_mean, mean_file, sep="\t")
write.table(norm_x_median, median_file, sep="\t")
write.table(norm_x_var, var_file, sep="\t")

print(paste("Maximum variance = ", max(norm_x_var), "\n", sep=""));
print(paste("Minimum variance = ", min(norm_x_var), "\n", sep=""));
print(paste("Number of variance > 0.05 = ", length(which(norm_x_var > 0.05)), "\n", sep=""));


## plot histogram for those metabolites have high variance
large_variance <- which(norm_x_var > 0.05)

png(paste(dir,"/", "variance_plot.png", sep=""), width=880, height=800)
op <- par(mfrow = c(ceiling(sqrt(length(large_variance))), ceiling(sqrt(length(large_variance)))))

for (i in 1:length(large_variance)){
	index <- large_variance[i]

	c_name <- colnames(x)[ceiling(index / ncol(x))]
	r_name <- colnames(x)[index %% ncol(x)]

	hist_title <- paste(r_name, " Fragments Mapped to ", c_name, sep="")

	hist(expand_x[, index], xlab = "Correlation", 
			main = hist_title, col = "lightgreen", cex.main = 1)
	curve(dnorm(x, mean=norm_x_mean[index], sd=sd(expand_x[, index])),
		 add=TRUE, col="darkblue", lwd=2) 

}
dev.off()

## output mean, var, median
#filter_mean_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_MeanFilter.txt", sep="")
#filter_median_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_MedianFilter.txt", sep="")
#filter_var_file <- paste(dir,"/",output_training,"/", type, "_matrix_A_VarianceFilter.txt", sep="")

#write.table(filter_mean, filter_mean_file, sep="\t")
#write.table(filter_median, filter_median_file, sep="\t")
#write.table(filter_var, filter_var_file, sep="\t")


## output validation data
validate_10XR3A <- rbind.fill(observed_10XR3A, as.data.frame(diag_normalized_median));
validate_10XR3A[is.na(validate_10XR3A)] <- 0;
validate_10XR3A <- validate_10XR3A[1:gene_count,colnames(diag_normalized_median)];
rownames(validate_10XR3A) <- rownames(observed_10XR3A);

out_X_10XR3A_file <- paste(dir, "/",output_validation,"/", type, "_10XR3A_validation_X.txt", sep="");
out_Y_10XR3A_file <- paste(dir, "/",output_validation,"/", type, "_10XR3A_validation_Y.txt", sep="");
write.table(validate_10XR3A, out_X_10XR3A_file, sep="\t");
write.table(expected_10XR3A, out_Y_10XR3A_file, sep="\t");

print(paste("Written Data for Validation to ", out_X_10XR3A_file, sep=""));
print(paste("Written Data for Validation to ", out_Y_10XR3A_file, sep=""));

#filter_X_10XR3A <- validate_10XR3A[-filter_genes,-filter_genes]
#filter_X_10XR3A <- filter_X_10XR3A[complete.cases(filter_X_10XR3A)]
#filter_Y_10XR3A <- expected_10XR3A[-filter_genes,]

#out_filter_X_10XR3A_file <- paste(dir, "/",output_validation,"/", type, "_10XR3A_validation_X_filter.txt", sep="");
#out_filter_Y_10XR3A_file <- paste(dir, "/",output_validation,"/", type, "_10XR3A_validation_Y_filter.txt", sep="");
#write.table(filter_X_10XR3A, out_filter_X_10XR3A_file, sep="\t");
#write.table(filter_Y_10XR3A, out_filter_Y_10XR3A_file, sep="\t");

#print(paste("Written Data for Validation to ", out_filter_X_10XR3A_file, sep=""));
#print(paste("Written Data for Validation to ", out_filter_Y_10XR3A_file, sep=""));


## output validation data
validate_20XR3A <- rbind.fill(observed_20XR3A, as.data.frame(diag_normalized_median));
validate_20XR3A[is.na(validate_20XR3A)] <- 0;
validate_20XR3A <- validate_20XR3A[1:gene_count,colnames(diag_normalized_median)];
rownames(validate_20XR3A) <- rownames(observed_20XR3A);

out_X_20XR3A_file <- paste(dir, "/",output_validation,"/", type, "_20XR3A_validation_X.txt", sep="");
out_Y_20XR3A_file <- paste(dir, "/",output_validation,"/", type, "_20XR3A_validation_Y.txt", sep="");
write.table(validate_20XR3A, out_X_20XR3A_file, sep="\t");
write.table(expected_20XR3A, out_Y_20XR3A_file, sep="\t");

print(paste("Written Data for Validation to ", out_X_20XR3A_file, sep=""));
print(paste("Written Data for Validation to ", out_Y_20XR3A_file, sep=""));

#filter_X_20XR3A <- validate_20XR3A[-filter_genes,-filter_genes]
#filter_X_20XR3A <- filter_X_20XR3A[complete.cases(filter_X_20XR3A)]
#filter_Y_20XR3A <- expected_20XR3A[-filter_genes,]

#out_filter_X_20XR3A_file <- paste(dir, "/",output_validation,"/", type, "_20XR3A_validation_X_filter.txt", sep="");
#out_filter_Y_20XR3A_file <- paste(dir, "/",output_validation,"/", type, "_20XR3A_validation_Y_filter.txt", sep="");
#write.table(filter_X_20XR3A, out_filter_X_20XR3A_file, sep="\t");
#write.table(filter_Y_20XR3A, out_filter_Y_20XR3A_file, sep="\t");

#print(paste("Written Data for Validation to ", out_filter_X_20XR3A_file, sep=""));
#print(paste("Written Data for Validation to ", out_filter_Y_20XR3A_file, sep=""));



## output validation data
validate_30XR3A <- rbind.fill(observed_30XR3A, as.data.frame(diag_normalized_median));
validate_30XR3A[is.na(validate_30XR3A)] <- 0;
validate_30XR3A <- validate_30XR3A[1:gene_count,colnames(diag_normalized_median)];
rownames(validate_30XR3A) <- rownames(observed_30XR3A);

out_X_30XR3A_file <- paste(dir, "/",output_validation,"/", type, "_30XR3A_validation_X.txt", sep="");
out_Y_30XR3A_file <- paste(dir, "/",output_validation,"/", type, "_30XR3A_validation_Y.txt", sep="");
write.table(validate_30XR3A, out_X_30XR3A_file, sep="\t");
write.table(expected_30XR3A, out_Y_30XR3A_file, sep="\t");

print(paste("Written Data for Validation to ", out_X_30XR3A_file, sep=""));
print(paste("Written Data for Validation to ", out_Y_30XR3A_file, sep=""));

#filter_X_30XR3A <- validate_30XR3A[-filter_genes,-filter_genes]
#filter_X_30XR3A <- filter_X_30XR3A[complete.cases(filter_X_30XR3A)]
#filter_Y_30XR3A <- expected_30XR3A[-filter_genes,]

#out_filter_X_30XR3A_file <- paste(dir, "/",output_validation,"/", type, "_30XR3A_validation_X_filter.txt", sep="");
#out_filter_Y_30XR3A_file <- paste(dir, "/",output_validation,"/", type, "_30XR3A_validation_Y_filter.txt", sep="");
#write.table(filter_X_30XR3A, out_filter_X_30XR3A_file, sep="\t");
#write.table(filter_Y_30XR3A, out_filter_Y_30XR3A_file, sep="\t");

#print(paste("Written Data for Validation to ", out_filter_X_30XR3A_file, sep=""));
#print(paste("Written Data for Validation to ", out_filter_Y_30XR3A_file, sep=""));


## output X, Y for beta training
out_expected_file <- paste(dir, "/",output_training,"/", type, "_lasso_expected_Y.txt", sep="");
out_training_file <- paste(dir, "/",output_training,"/", type, "_lasso_matrix_X.txt", sep="");

write.table(final_x, out_training_file, sep="\t");
write.table(final_y, out_expected_file, sep="\t");

print(paste("Written the Expected Value (Y) to ", out_expected_file, sep=""));
print(paste("Written the Training Matrix (X) to ", out_training_file, sep=""));
