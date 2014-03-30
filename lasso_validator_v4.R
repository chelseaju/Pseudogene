################################## lasso_validator_v4.R ################################
# Function: Validate the prediction generated from lasso_predictor.R with real data    #
# Usage: R --no-save < lasso_validator.R --args dir subdir type                        #
# Arguments: dir =  directory for expected data                                        #
#           subdir = directory for the prediction and output validation laasso         #
#     output = type_lasso_validator_v1.txt type_lasso_validator_v2.txt                 # 
# Author: Chelsea Ju                                                                   #
# Date: 2013-12-04                                                                     #
# Last Modify: 2013-12-22                                                              #
########################################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_predictor.R --args dir subdir type\n",
    "\t dir = directory for expected data\n",
    "\t subdir = directory for the prediction and output validation \n",
    "\t type = genes or transcripts\n"),
    sep="");
}

dir <- options[1];
subdir <- options[2];
type <- options[3];


Y_file <- paste(dir, "/LassoValidation/", type, "_", subdir, "_expected_Y.txt", sep="");
Y_vector <- read.table(Y_file, col.names=c("Gene", "Expected"));
rownames(Y_vector) <- Y_vector$Gene;


predicted_file_v1 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_prediction_v1.txt", sep="");
predicted_file_v2 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_prediction_v2.txt", sep="");

predicted_y_v1 <- read.table(predicted_file_v1, col.names = c("Gene","Prediction"));
predicted_y_v2 <- read.table(predicted_file_v2, col.names = c("Gene","Prediction"));
rownames(predicted_y_v1) <- predicted_y_v1$Gene;
rownames(predicted_y_v2) <- predicted_y_v2$Gene;

## validate analysis
Y <- Y_vector[rownames(predicted_y_v1),]
validation_v1 <- cbind(Y$Expected, predicted_y_v1$Prediction, abs(predicted_y_v1$Prediction - Y$Expected), abs(predicted_y_v1$Prediction - Y$Expected) / Y$Expected);

Y <- Y_vector[rownames(predicted_y_v2),]
validation_v2 <- cbind(Y$Expected, predicted_y_v2$Prediction, abs(predicted_y_v2$Prediction - Y$Expected), abs(predicted_y_v2$Prediction - Y$Expected) / Y$Expected);

colnames(validation_v1) <- c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");
colnames(validation_v2) <- c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");
rownames(validation_v1) <- rownames(predicted_y_v1);
rownames(validation_v2) <- rownames(predicted_y_v2);


## write to file

validation_file_v1 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_validation_v1.txt", sep="");
validation_file_v2 <- paste(dir, "/LassoValidation/", type, "_", subdir, "_validation_v2.txt", sep="");

write.table(validation_v1, file = validation_file_v1, sep="\t");
write.table(validation_v2, file = validation_file_v2, sep="\t");

print(paste("Writing the validation to ", validation_file_v1, sep=""));
print(paste("Writing the validation to ", validation_file_v2, sep=""));




