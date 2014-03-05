############################# lasso_trainer.R #############################
# Function: read in the expected value and training matrix from files     #
# 	use lasso to train the beta value (genlasso and glmnet)               #
# Usage: R --no-save < lasso_trainer.R --args dir type                    #
# Arguments: input =  directory for input and output                      #
#            type = genes or transcripts                                  #
#            output = lasso_type_coefficient.txt                          #
# Author: Chelsea Ju                                                      #
# Date: 2014-01-14                                                        #
# Modify from lasso_analysis.R                                            #
# Last Modify : 2014-03-05 separate data preparation                      #
###########################################################################

library(genlasso)
library(glmnet)

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

expected_file <- paste(dir, "/", type, "_", "lasso_expected_Y.txt", sep="");
training_file <- paste(dir, "/", type, "_", "lasso_matrix_X.txt", sep="");

x <- read.table(training_file, header = T);
y <- read.table(expected_file, header = T);


## genlasso
D <- diag(1, ncol(x), ncol(x));
out <- genlasso(y[,2], X=as.matrix(x), D=D);
summary(out);

beta <- coef(out, lambda=sqrt(nrow(x) * log(ncol(x))));
beta_value <- beta$beta;
rownames(beta_value) <- x_colname;

## glmnet
fit <- glmnet(as.matrix(x), y[,2]);
beta2 <- as.matrix(coef(fit, s = 0.01));

## output
filename_genlasso <- paste(dir, "/", type, "_genlasso_coefficient.xls", sep="");
filename_glmnet <- paste(dir, "/", type, "_glmnet_coefficient.xls", sep="");

write.table(paste("lambda", sqrt(nrow(x) * log(ncol(x))), sep="\t"), file=filename_genlasso);
write.table(beta_value[order(beta_value),], file=filename_genlasso, append = T, quote = FALSE);
print(paste("Written the Coefficient to ", filename_genlasso, sep=""));

write.table(beta2, file=filename_glmnet, quote = FALSE);

print(paste("Written the Coefficient to ", filename_glmnet, sep=""));










