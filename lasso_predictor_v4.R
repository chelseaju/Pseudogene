########################### lasso_predictor_v2.R ##########################
# Function: use the trained coefficient generated from glmnet to predict  #
#   the fragment count. two types of prediction is performed:             #
#   1. take in the matrix X to compute Y                                  #
#   2. combine the read count for a specific region (mimicing the real    #
#      scenario), use alpha matrix to reconstruct the x matrix            #
# Usage: R --no-save < lasso_predictor_v2.R --args dir subdir type        #
# Arguments: dir =  directory                                             #
#            type = genes or transcripts or ENSG                          #
# Output = type_lasso_predictor_v1.txt type_lasso_predictor_v2.txt        #
# Author: Chelsea Ju                                                      #
# Date: 2014-03-28                                                        #
# Modify from lasso_analysis_validator.R                                  #
###########################################################################

library(plyr)

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_predictor.R --args dir subdir type lasso\n",
    "\t dir = directory \n",
    "\t subdir = subdirectory \n",
    "\t type = genes or transcripts\n"),
    sep="");
}

dir <- options[1];
subdir<- options[2];
type <- options[3];



## read in training X matrix to obtain the column and row names
training_X_file <- paste(dir,"/LassoTraining/", type, "_lasso_matrix_X.txt", sep="");
training_X <- read.table(training_X_file);
c_name <- colnames(training_X);
r_name <- rownames(training_X);

## read in betas
beta_file <-  paste(dir,"/LassoTraining/", type, "_B_Values.txt", sep="");
betas <- read.table(beta_file, header = F);
rownames(betas) <- c_name;

## read in alpha matrix
alpha_file <- paste(dir, "/LassoTraining/", type, "_A_Matrix.txt", sep="");
alphas <- read.table(alpha_file, header = F);
colnames(alphas) <- c_name;

## read in data
X_file <- paste(dir, "/LassoValidation/", type, "_", subdir, "_matrix_X.txt", sep="");

X_matrix <- read.table(X_file);

X_sum <- colSums(X_matrix); ## mimicing the real data

## prediction
mutual_betas <- intersect(colnames(X_matrix), rownames(betas))

filter_betas <- betas[mutual_betas,];
filter_alphas <- alphas[,mutual_betas];

filter_X <- X_matrix[,mutual_betas];
filter_X_sum <- colSums(filter_X);

predict_y_v1 <- as.matrix(filter_X) %*% as.matrix(filter_betas);

predict_x <- as.matrix(filter_alphas) %*% diag(filter_X_sum);
colnames(predict_x) <- colnames(filter_X);
rownames(predict_x) <- rownames(filter_X);
predict_y_v2 <- as.matrix(predict_x) %*% as.matrix(filter_betas);

## write to file
output_v1 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_prediction_v1.txt", sep="");
output_v2 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_prediction_v2.txt", sep="");

output_predictX <- paste(dir, "/LassoValidation/", type, "_", subdir, "_prediction_x.txt", sep="");
output_expectedX <- paste(dir, "/LassoValidation/", type, "_", subdir, "_expected_x.txt", sep="");

write.table(predict_y_v1, file = output_v1, sep="\t", col.names = FALSE);
write.table(predict_y_v2, file = output_v2, sep="\t", col.names = FALSE);

write.table(predict_x, file = output_predictX, sep="\t");
write.table(filter_X, file = output_expectedX, sep="\t");


print(paste("Writing the prediction to ", output_v1, sep=""));
print(paste("Writing the prediction to ", output_v2, sep=""));
print(paste("Writing the predicted X to ", output_predictX, sep=""));
print(paste("Writing the expected X to ", output_expectedX, sep=""));









