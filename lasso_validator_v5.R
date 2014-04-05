################################## lasso_validator_v5.R ################################
# Function: Validate the prediction generated from lasso with real data                #
# Usage: R --no-save < lasso_validator.R --args dir X Y beta expected_X output         #
# Arguments: dir =  top directory for data                                             #
#            X = reconstructed X matrix                                                #
#            Y = expected Y value                                                      #
#       output = output name for file                                                  # 
# Author: Chelsea Ju                                                                   #
# Date: 2014-04-01                                                                     #
########################################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 6){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_validator_v5.R --args dir X Y beta expected_X output\n",
    "\t dir = top directory for data input and output\n",
    "\t X = reconstructed X matrix\n",
    "\t Y = expected Y values\n",
    "\t beta = beta file\n",
    "\t expected_X = expected X file\n",
    "\t output = output file\n",
    sep=""));
}

dir <- options[1];
X <- options[2];
Y <- options[3];
B <- options[4];
expected_X <- options[5];
output <- options[6];

v1_outfile <- paste(dir, "/LassoValidation_v2/v1_", output, sep="");
v2_outfile <- paste(dir, "/LassoValidation_v2/v2_", output, sep="");


## read in data
X_data <- read.table(paste(dir, "/LassoValidation_v2/", X, sep=""));
Y_data <- read.table(paste(dir, "/LassoValidation_v2/", Y, sep=""), header = TRUE);
B_data <- read.table(paste(dir, "/", B, sep=""));
eX_data <- read.table(paste(dir,"/LassoValidation_v2/", expected_X, sep=""), header = TRUE);
rownames(B_data) <- B_data$V1;

rownames(X_data) <- rownames(eX_data);
colnames(X_data) <- colnames(eX_data);


## X * Beta
# v1 raw matrix

## check for data names
name_verification <- colnames(eX_data) == rownames(B_data);
if(length(name_verification[name_verification == FALSE]) > 0){
#	stop(paste("Rownames of Beta do not match with Colnames of X"));
	
	B_data <- B_data[colnames(eX_data),]

}

estimated_Y <- as.matrix(eX_data) %*% as.matrix(B_data$V2);
estimated_Y <- estimated_Y[Y_data[,1],]; # reorganized Y

expected_Y <- Y_data[,2];

validation <- cbind(expected_Y, estimated_Y, (estimated_Y - expected_Y), abs(estimated_Y - expected_Y)/expected_Y);

colnames(validation) <- c("Expected", "Prediction", "Error", "Relative Error");

write.table(validation, v1_outfile, sep="\t");

print(paste("Written File to ", v1_outfile, sep=""));




# v2 reconstructed matrix

## check for data names
name_verification <- colnames(X_data) == rownames(B_data);
if(length(name_verification[name_verification == FALSE]) > 0){
	stop(paste("Rownames of Beta do not match with Colnames of X"));
}

estimated_Y <- as.matrix(X_data) %*% as.matrix(B_data$V2);
estimated_Y <- estimated_Y[Y_data[,1],]; # reorganized Y

expected_Y <- Y_data[,2];

validation <- cbind(expected_Y, estimated_Y, (estimated_Y - expected_Y), abs(estimated_Y - expected_Y)/expected_Y);

colnames(validation) <- c("Expected", "Prediction", "Error", "Relative Error");

write.table(validation, v2_outfile, sep="\t");

print(paste("Written File to ", v2_outfile, sep=""));








