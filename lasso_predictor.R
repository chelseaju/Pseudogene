########################### lasso_predictor.R #############################
# Function: use the trained coefficient generated from lasso_trainer.R    #
#   to predict the fragment count. two types of prediction is performed:  #
#   1. take in the matrix X to compute Y                                  #
#   2. combine the read count for a specific region (mimicing the real    #
#      scenario), and compute Y                                           #
# Usage: R --no-save < lasso_predictor.R --args dir subdir type           #
# Arguments: dir =  directory for coefficient                             #
#           subdir = directory for the predicting file and output         #
#            type = genes or transcripts                                  #
# Output = type_lasso_predictor_v1.txt type_lasso_predictor_v2.txt        #
# Author: Chelsea Ju                                                      #
# Date: 2014-01-15                                                        #
# Modify from lasso_analysis_validator.R                                  #
###########################################################################

library(plyr)

read_coefficient <- function(file){
    data <- read.table(file, skip = 3);
    colnames(data) <- c("Gene", "Coefficient");
    rownames(data) <- data$Gene;
    data;
}

read_coefficient2 <- function(file){
    data <- read.table(file, skip = 2);
    colnames(data) <- c("Gene", "Coefficient");
    rownames(data) <- data$Gene;
    data;
}

## read in the distribution, and expand the matrix based on the trained coefficient
## if the column is missing from the coefficient, add a new zero column to the matrix
read_distribution<- function(file, row_to_expend){
	matrix_data <- read.table(file);	
    extra_row <- row_to_expend[!row_to_expend %in% colnames(matrix_data)];
    extra_columns <- as.data.frame(matrix(rep(0, nrow(matrix_data) * length(extra_row)), ncol = length(extra_row)))
    colnames(extra_columns) <- extra_row;
    matrix_data <- cbind(matrix_data, extra_columns);
    matrix_data;
}



# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_predictor.R --args dir subdir type\n",
    "\t dir = directory of coefficient file\n",
    "\t subdir = directory of predicting file and output\n",
    "\t type = genes or transcripts\n"),
    sep="");
}

dir <- options[1];
subdir <- options[2];
type <- options[3];

coefficient_file <- paste(dir,"/", type, "_lasso_coefficient2.xls", sep="");
coefficient <- read_coefficient2(coefficient_file);

distribution_file <- paste(dir, "/", subdir, "/", type, "_distribution.matrix", sep="");
distribution <- read_distribution(distribution_file, rownames(coefficient));

# reorder the coefficient for the matrix
coefficient <- coefficient[colnames(distribution),]

# predict y from matrix
predict_y_v1 <- as.matrix(distribution) %*% as.matrix(coefficient$Coefficient);

# mimicing the real scenario, where we do not know how the reads are distributed
observation <- colSums(distribution);

predict_y_v2 <- t(as.matrix(observation)) %*% diag(c(coefficient$Coefficient), nrow(coefficient));
predict_y_v2 <- t(predict_y_v2);

rownames(predict_y_v2) <- colnames(distribution);


## write to file
output_v1 <- paste(dir, "/", subdir, "/", type, "_lasso_prediction_v1.txt", sep="");
output_v2 <- paste(dir, "/", subdir, "/", type, "_lasso_prediction_v2.txt", sep="");

write.table(predict_y_v1, file = output_v1, sep="\t", col.names = FALSE);
write.table(predict_y_v2, file = output_v2, sep="\t", col.names = FALSE);

print(paste("Writing the prediction to ", output_v1, sep=""));
print(paste("Writing the prediction to ", output_v2, sep=""));









