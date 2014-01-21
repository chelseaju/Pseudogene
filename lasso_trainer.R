############################# lasso_trainer.R #############################
# Function: read in the fragment experiment matrix,                       #
#         matrix. Use Lasso to train on the cov matrix.                   #
# Usage: R --no-save < lasso_analysis.R --args dir type                   #
# Arguments: input =  directory for input and output                      #
#            type = genes or transcripts                                  #
#            output = type_lasso_coefficient.txt                          #
# Author: Chelsea Ju                                                      #
# Date: 2014-01-14                                                        #
# Modify from lasso_analysis.R                                            #
###########################################################################

library(genlasso)
library(plyr)
library(glmnet)


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
	index_vector <- c(1:length(row_names));
	row_names <- row_names[which(index_vector %% 2 != 0)];
	
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
    "\t type = genes or transcripts\n"),
    sep="");
}

dir <- options[1];
type <- options[2];

dir_10X10A <- paste(dir, "/10X_100L_10A/", sep="") 
dir_10X30A <- paste(dir, "/10X_100L_30A/", sep="") 
dir_10X50A <- paste(dir, "/10X_100L_50A/", sep="") 

dir_20X10A <- paste(dir, "/20X_100L_10A/", sep="") 
dir_20X30A <- paste(dir, "/20X_100L_30A/", sep="") 
dir_20X50A <- paste(dir, "/20X_100L_50A/", sep="") 

observed_10X10A <- read_distribution_matrix(dir_10X10A, type);
observed_10X30A <- read_distribution_matrix(dir_10X30A, type);
observed_10X50A <- read_distribution_matrix(dir_10X50A, type);

observed_20X10A <- read_distribution_matrix(dir_20X10A, type);
observed_20X30A <- read_distribution_matrix(dir_20X30A, type);
observed_20X50A <- read_distribution_matrix(dir_20X50A, type);

expected_10X10A <- read_expected_data(dir_10X10A, type, rownames(observed_10X10A), "10X10A");
expected_10X30A <- read_expected_data(dir_10X30A, type, rownames(observed_10X30A), "10X30A");
expected_10X50A <- read_expected_data(dir_10X50A, type, rownames(observed_10X50A), "10X50A");

expected_20X10A <- read_expected_data(dir_20X10A, type, rownames(observed_20X10A), "20X10A");
expected_20X30A <- read_expected_data(dir_20X30A, type, rownames(observed_20X30A), "20X30A");
expected_20X50A <- read_expected_data(dir_20X50A, type, rownames(observed_20X50A), "20X50A");

y <- rbind(expected_10X10A, expected_10X30A, expected_10X50A, expected_20X10A, expected_20X30A, expected_20X50A); 
x <- rbind.fill(observed_10X10A, observed_10X30A, observed_10X50A, observed_20X10A, observed_20X30A, observed_20X50A);

D <- diag(1, ncol(x), ncol(x));

x_colname <- colnames(x);
x[is.na(x)] <- 0;

#out <- genlasso(y[,2], X=as.matrix(x), D=D);
#summary(out);

#beta <- coef(out, lambda=sqrt(nrow(x) * log(ncol(x))));
#beta_value <- beta$beta;
#rownames(beta_value) <- x_colname;

fit <- glmnet(as.matrix(x), y[,2]);
beta <- as.matrix(coef(fit, s = 0.01));

filename = paste(dir, "/", type, "_lasso_coefficient2.xls", sep="");
#write.table(paste("lambda", sqrt(nrow(x) * log(ncol(x))), sep="\t"), file=filename);
#write.table(beta_value[order(beta_value),], file=filename, append = T, quote = FALSE);

write.table(beta, file=filename, quote = FALSE);

print(paste("Written the Coefficient to ", filename, sep=""));










