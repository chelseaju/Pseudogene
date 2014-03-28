################################## lasso_validator.R ###################################
# Function: Validate the prediction generated from lasso_predictor.R with real data    #
# Usage: R --no-save < lasso_validator.R --args dir subdir type                        #
# Arguments: dir =  directory for expected data                                        #
#           subdir = directory for the prediction and output validation laasso         #
#     output = type_lasso_validator_v1.txt type_lasso_validator_v2.txt                 # 
# Author: Chelsea Ju                                                                   #
# Date: 2013-12-04                                                                     #
# Last Modify: 2013-12-22                                                              #
########################################################################################


# self-defined function
read_expectation <- function(file){
    data <- read.table(file);
    c1 <- data$V1;
    row_names <- unlist(strsplit(as.character(c1), "_"));
    index_vector <- c(1:length(row_names));
	row_names <- row_names[which(index_vector %% 2 != 0)];
    rownames(data) <- row_names;
    colnames(data) <- c("Name", "Expected");
    data;
}

read_prediction <- function(file, row_order){
    data <- read.table(file);
    rownames(data) <- data[,1];
    data <- data[row_order,];
    colnames(data) <- c("Name", "Predicted");
    data;    
}

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 4){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < lasso_predictor.R --args dir subdir type\n",
    "\t dir = directory for expected data\n",
    "\t subdir = directory for the prediction and output validation \n",
    "\t type = genes or transcripts\n",
    "\t lasso = genlasso or glmnet\n"),
    sep="");
}

dir <- options[1];
subdir <- options[2];
type <- options[3];
lasso <- options[4];

expected_file <- paste(dir,"/", subdir, "/", type, "_expected_read_count.txt", sep="");
expectation <- read_expectation(expected_file);

predicted_file_v1 <- paste(dir,"/", subdir, "/", type, "_", lasso, "_training_2_prediction_v1.txt", sep="");
predicted_file_v2 <- paste(dir,"/", subdir, "/", type, "_", lasso, "_training_2_prediction_v2.txt", sep="");

prediction_v1 <- read_prediction(predicted_file_v1, rownames(expectation));
prediction_v2 <- read_prediction(predicted_file_v2, rownames(expectation));


## validate analysis
validation_v1 <- cbind(expectation$Expected, prediction_v1$Predicted, abs(prediction_v1$Predicted - expectation$Expected), abs(prediction_v1$Predicted - expectation$Expected) / expectation$Expected);
validation_v2 <- cbind(expectation$Expected, prediction_v2$Predicted, abs(prediction_v2$Predicted - expectation$Expected), abs(prediction_v2$Predicted - expectation$Expected) / expectation$Expected);

colnames(validation_v1) <- c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");
colnames(validation_v2) <- c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");
rownames(validation_v1) <- rownames(prediction_v1);
rownames(validation_v2) <- rownames(prediction_v2);


## write to file
output_v1 <- paste(dir, "/", subdir, "/", type, "_", lasso, "_training_2_validation_v1.txt", sep="");
output_v2 <- paste(dir, "/", subdir, "/", type, "_", lasso, "_training_2_validation_v2.txt", sep="");

write.table(validation_v1, file = output_v1, sep="\t");
write.table(validation_v2, file = output_v2, sep="\t");

print(paste("Writing the validation to ", output_v1, sep=""));
print(paste("Writing the validation to ", output_v2, sep=""));




