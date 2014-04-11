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
dir_3X1A <- paste(dir, "/3X_101L_1A/", sep="")
dir_3XR1A <- paste(dir, "/3X_101L_R1A/", sep="")
dir_3XR2A <- paste(dir, "/3X_101L_R2A/", sep="")
dir_3XR3A <- paste(dir, "/3X_101L_R3A/", sep="")

dir_5X1A <- paste(dir, "/5X_101L_1A/", sep="")
dir_5XR1A <- paste(dir, "/5X_101L_R1A/", sep="")
dir_5XR2A <- paste(dir, "/5X_101L_R2A/", sep="")
dir_5XR3A <- paste(dir, "/5X_101L_R3A/", sep="")

dir_7X1A <- paste(dir, "/7X_101L_1A/", sep="")
dir_7XR1A <- paste(dir, "/7X_101L_R1A/", sep="")
dir_7XR2A <- paste(dir, "/7X_101L_R2A/", sep="")
dir_7XR3A <- paste(dir, "/7X_101L_R3A/", sep="")

dir_10X1A <- paste(dir, "/10X_101L_1A/", sep="") 
dir_10XR1A <- paste(dir, "/10X_101L_R1A/", sep="")
dir_10XR2A <- paste(dir, "/10X_101L_R2A/", sep="") 
dir_10XR3A <- paste(dir, "/10X_101L_R3A/", sep="") 

dir_13X1A <- paste(dir, "/13X_101L_1A/", sep="")
dir_13XR1A <- paste(dir, "/13X_101L_R1A/", sep="")
dir_13XR2A <- paste(dir, "/13X_101L_R2A/", sep="")
dir_13XR3A <- paste(dir, "/13X_101L_R3A/", sep="")

dir_17X1A <- paste(dir, "/17X_101L_1A/", sep="")
dir_17XR1A <- paste(dir, "/17X_101L_R1A/", sep="")
dir_17XR2A <- paste(dir, "/17X_101L_R2A/", sep="")
dir_17XR3A <- paste(dir, "/17X_101L_R3A/", sep="")

dir_20X1A <- paste(dir, "/20X_101L_1A/", sep="") 
dir_20XR1A <- paste(dir, "/20X_101L_R1A/", sep="")
dir_20XR2A <- paste(dir, "/20X_101L_R2A/", sep="") 
dir_20XR3A <- paste(dir, "/20X_101L_R3A/", sep="")

dir_23X1A <- paste(dir, "/23X_101L_1A/", sep="")
dir_23XR1A <- paste(dir, "/23X_101L_R1A/", sep="")
dir_23XR2A <- paste(dir, "/23X_101L_R2A/", sep="")
dir_23XR3A <- paste(dir, "/23X_101L_R3A/", sep="")

dir_27X1A <- paste(dir, "/27X_101L_1A/", sep="")
dir_27XR1A <- paste(dir, "/27X_101L_R1A/", sep="")
dir_27XR2A <- paste(dir, "/27X_101L_R2A/", sep="")
dir_27XR3A <- paste(dir, "/27X_101L_R3A/", sep="")

dir_30X1A <- paste(dir, "/30X_101L_1A/", sep="") 
dir_30XR1A <- paste(dir, "/30X_101L_R1A/", sep="")
dir_30XR2A <- paste(dir, "/30X_101L_R2A/", sep="") 
dir_30XR3A <- paste(dir, "/30X_101L_R3A/", sep="") 


observed_3X1A <- read_distribution_matrix(dir_3X1A, type);
observed_3XR1A <- read_distribution_matrix(dir_3XR1A, type);
observed_3XR2A <- read_distribution_matrix(dir_3XR2A, type);
observed_3XR3A <- read_distribution_matrix(dir_3XR3A, type);

observed_5X1A <- read_distribution_matrix(dir_5X1A, type);
observed_5XR1A <- read_distribution_matrix(dir_5XR1A, type);
observed_5XR2A <- read_distribution_matrix(dir_5XR2A, type);
observed_5XR3A <- read_distribution_matrix(dir_5XR3A, type);

observed_7X1A <- read_distribution_matrix(dir_7X1A, type);
observed_7XR1A <- read_distribution_matrix(dir_7XR1A, type);
observed_7XR2A <- read_distribution_matrix(dir_7XR2A, type);
observed_7XR3A <- read_distribution_matrix(dir_7XR3A, type);

observed_10X1A <- read_distribution_matrix(dir_10X1A, type);
observed_10XR1A <- read_distribution_matrix(dir_10XR1A, type);
observed_10XR2A <- read_distribution_matrix(dir_10XR2A, type);
observed_10XR3A <- read_distribution_matrix(dir_10XR3A, type);

observed_13X1A <- read_distribution_matrix(dir_13X1A, type);
observed_13XR1A <- read_distribution_matrix(dir_13XR1A, type);
observed_13XR2A <- read_distribution_matrix(dir_13XR2A, type);
observed_13XR3A <- read_distribution_matrix(dir_13XR3A, type);

observed_17X1A <- read_distribution_matrix(dir_17X1A, type);
observed_17XR1A <- read_distribution_matrix(dir_17XR1A, type);
observed_17XR2A <- read_distribution_matrix(dir_17XR2A, type);
observed_17XR3A <- read_distribution_matrix(dir_17XR3A, type);

observed_20X1A <- read_distribution_matrix(dir_20X1A, type);
observed_20XR1A <- read_distribution_matrix(dir_20XR1A, type);
observed_20XR2A <- read_distribution_matrix(dir_20XR2A, type);
observed_20XR3A <- read_distribution_matrix(dir_20XR3A, type);

observed_23X1A <- read_distribution_matrix(dir_23X1A, type);
observed_23XR1A <- read_distribution_matrix(dir_23XR1A, type);
observed_23XR2A <- read_distribution_matrix(dir_23XR2A, type);
observed_23XR3A <- read_distribution_matrix(dir_23XR3A, type);

observed_27X1A <- read_distribution_matrix(dir_27X1A, type);
observed_27XR1A <- read_distribution_matrix(dir_27XR1A, type);
observed_27XR2A <- read_distribution_matrix(dir_27XR2A, type);
observed_27XR3A <- read_distribution_matrix(dir_27XR3A, type);

observed_30X1A <- read_distribution_matrix(dir_30X1A, type);
observed_30XR1A <- read_distribution_matrix(dir_30XR1A, type);
observed_30XR2A <- read_distribution_matrix(dir_30XR2A, type);
observed_30XR3A <- read_distribution_matrix(dir_30XR3A, type);


## expected vector
expected_3X1A <- read_expected_data(dir_3X1A, type, rownames(observed_3X1A), "3X1A");
expected_3XR1A <- read_expected_data(dir_3XR1A, type, rownames(observed_3XR1A), "3XR1A");
expected_3XR2A <- read_expected_data(dir_3XR2A, type, rownames(observed_3XR2A), "3XR2A");
expected_3XR3A <- read_expected_data(dir_3XR3A, type, rownames(observed_3XR3A), "3XR3A");

expected_5X1A <- read_expected_data(dir_5X1A, type, rownames(observed_5X1A), "5X1A");
expected_5XR1A <- read_expected_data(dir_5XR1A, type, rownames(observed_5XR1A), "5XR1A");
expected_5XR2A <- read_expected_data(dir_5XR2A, type, rownames(observed_5XR2A), "5XR2A");
expected_5XR3A <- read_expected_data(dir_5XR3A, type, rownames(observed_5XR3A), "5XR3A");

expected_7X1A <- read_expected_data(dir_7X1A, type, rownames(observed_7X1A), "7X1A");
expected_7XR1A <- read_expected_data(dir_7XR1A, type, rownames(observed_7XR1A), "7XR1A");
expected_7XR2A <- read_expected_data(dir_7XR2A, type, rownames(observed_7XR2A), "7XR2A");
expected_7XR3A <- read_expected_data(dir_7XR3A, type, rownames(observed_7XR3A), "7XR3A");

expected_10X1A <- read_expected_data(dir_10X1A, type, rownames(observed_10X1A), "10X1A");
expected_10XR1A <- read_expected_data(dir_10XR1A, type, rownames(observed_10XR1A), "10XR1A");
expected_10XR2A <- read_expected_data(dir_10XR2A, type, rownames(observed_10XR2A), "10XR2A");
expected_10XR3A <- read_expected_data(dir_10XR3A, type, rownames(observed_10XR3A), "10XR3A");

expected_13X1A <- read_expected_data(dir_13X1A, type, rownames(observed_13X1A), "13X1A");
expected_13XR1A <- read_expected_data(dir_13XR1A, type, rownames(observed_13XR1A), "13XR1A");
expected_13XR2A <- read_expected_data(dir_13XR2A, type, rownames(observed_13XR2A), "13XR2A");
expected_13XR3A <- read_expected_data(dir_13XR3A, type, rownames(observed_13XR3A), "13XR3A");

expected_17X1A <- read_expected_data(dir_17X1A, type, rownames(observed_17X1A), "17X1A");
expected_17XR1A <- read_expected_data(dir_17XR1A, type, rownames(observed_17XR1A), "17XR1A");
expected_17XR2A <- read_expected_data(dir_17XR2A, type, rownames(observed_17XR2A), "17XR2A");
expected_17XR3A <- read_expected_data(dir_17XR3A, type, rownames(observed_17XR3A), "17XR3A");

expected_20X1A <- read_expected_data(dir_20X1A, type, rownames(observed_20X1A), "20X1A");
expected_20XR1A <- read_expected_data(dir_20XR1A, type, rownames(observed_20XR1A), "20XR1A");
expected_20XR2A <- read_expected_data(dir_20XR2A, type, rownames(observed_20XR2A), "20XR2A");
expected_20XR3A <- read_expected_data(dir_20XR3A, type, rownames(observed_20XR3A), "20XR3A");

expected_23X1A <- read_expected_data(dir_23X1A, type, rownames(observed_23X1A), "23X1A");
expected_23XR1A <- read_expected_data(dir_23XR1A, type, rownames(observed_23XR1A), "23XR1A");
expected_23XR2A <- read_expected_data(dir_23XR2A, type, rownames(observed_23XR2A), "23XR2A");
expected_23XR3A <- read_expected_data(dir_23XR3A, type, rownames(observed_23XR3A), "23XR3A");

expected_27X1A <- read_expected_data(dir_27X1A, type, rownames(observed_27X1A), "27X1A");
expected_27XR1A <- read_expected_data(dir_27XR1A, type, rownames(observed_27XR1A), "27XR1A");
expected_27XR2A <- read_expected_data(dir_27XR2A, type, rownames(observed_27XR2A), "27XR2A");
expected_27XR3A <- read_expected_data(dir_27XR3A, type, rownames(observed_27XR3A), "27XR3A");

expected_30X1A <- read_expected_data(dir_30X1A, type, rownames(observed_30X1A), "30X1A");
expected_30XR1A <- read_expected_data(dir_30XR1A, type, rownames(observed_30XR1A), "30XR1A");
expected_30XR2A <- read_expected_data(dir_30XR2A, type, rownames(observed_30XR2A), "30XR2A");
expected_30XR3A <- read_expected_data(dir_30XR3A, type, rownames(observed_30XR3A), "30XR3A");


y <- rbind(expected_3X1A, expected_3XR1A, expected_3XR2A, expected_3XR3A,
	expected_5XR1A, expected_5XR2A, expected_5XR3A,
	expected_7X1A, expected_7XR1A, expected_7XR2A, expected_7XR3A,
	expected_10X1A, expected_10XR1A, expected_10XR2A, expected_10XR3A,
	expected_13X1A, expected_13XR1A, expected_13XR2A, expected_13XR3A,
	expected_17X1A, expected_17XR1A, expected_17XR2A, expected_17XR3A,
	expected_20X1A, expected_20XR1A, expected_20XR2A, expected_20XR3A,
	expected_23X1A, expected_23XR1A, expected_23XR2A, expected_23XR3A,
	expected_27X1A, expected_27XR1A, expected_27XR2A, expected_27XR3A,
	expected_30X1A, expected_30XR1A, expected_30XR2A, expected_30XR3A)

x <- rbind.fill(observed_3X1A, observed_3XR1A, observed_3XR2A, observed_3XR3A,
				observed_5XR1A, observed_5XR2A, observed_5XR3A,
				observed_7X1A, observed_7XR1A, observed_7XR2A, observed_7XR3A,
				observed_10X1A, observed_10XR1A, observed_10XR2A, observed_10XR3A,
				observed_13X1A, observed_13XR1A, observed_13XR2A, observed_13XR3A,
				observed_17X1A, observed_17XR1A, observed_17XR2A, observed_17XR3A,
				observed_20X1A, observed_20XR1A, observed_20XR2A, observed_20XR3A,
				observed_23X1A, observed_23XR1A, observed_23XR2A, observed_23XR3A,
				observed_27X1A, observed_27XR1A, observed_27XR2A, observed_27XR3A,
				observed_30X1A, observed_30XR1A, observed_30XR2A, observed_30XR3A)

x_colname <- colnames(x);
x[is.na(x)] <- 0;

colnames(y) <- c("Gene_Name", "Read_Count");
rownames(x) <- rownames(y);


## normalize observation matrix against expected value
normed_x <- x / y$Read_Count;

## number of genes
gene_count <- nrow(x) / 39

## convert normed_x as a [sample x obs] matrix (in this case, 16 * 3840)
##  3840 has 128 genes x 304 regions
expand_x <- matrix(0, nrow = 39, ncol = gene_count * ncol(x))

for(i in 1:39){
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
mean_file <- paste(dir,"/observation/", type, "_matrix_A_Mean.txt", sep="")
median_file <- paste(dir,"/observation/", type, "_matrix_A_Median.txt", sep="")
var_file <- paste(dir,"/observation/", type, "_matrix_A_Variance.txt", sep="")

write.table(norm_x_mean, mean_file, sep="\t")
write.table(norm_x_median, median_file, sep="\t")
write.table(norm_x_var, var_file, sep="\t")

print(paste("Maximum variance = ", max(norm_x_var), "\n", sep=""));
print(paste("Minimum variance = ", min(norm_x_var), "\n", sep=""));
print(paste("Number of variance > 0.05 = ", length(which(norm_x_var > 0.05)), "\n", sep=""));


## plot histogram for those metabolites have high variance
large_variance <- which(norm_x_var > 0.05)

png(paste(dir,"/observation/", "variance_plot.png", sep=""), width=880, height=800)
op <- par(mfrow = c(ceiling(sqrt(length(large_variance))), ceiling(sqrt(length(large_variance)))))

for (i in 1:length(large_variance)){
	index <- large_variance[i]

	c_name <- colnames(x)[ceiling(index / ncol(x))]
	r_name <- colnames(x)[index %% ncol(x)]

	hist_title <- paste(r_name, " Mapped to ", c_name, sep="")

	hist(expand_x[, index], xlab = "Correlation", 
			main = hist_title, col = "lightgreen", cex.main = 0.8)
	curve(dnorm(x, mean=norm_x_mean[index], sd=sd(expand_x[, index])),
		 add=TRUE, col="darkblue", lwd=2) 

}
dev.off()


