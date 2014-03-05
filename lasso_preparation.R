########################## lasso_preparation.R ############################
# Function: read in data (expected Y) and training matrix X from different#
#         directories. Organize data in two files for lasso training      #
# Usage: R --no-save < lasso_analysis.R --args dir datatype               #
# Arguments: dir = input directory (select_one_pseudogene_110)            #
#            datatype = [gene|transcript]                                 #
#            output = lasso_datatype_expected_Y.txt                       #
#                     lasso_datatype_matrix_X.txt                         #
# Author: Chelsea Ju                                                      #
# Date: 2014-03-05                                                        #
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
	index_vector <- c(1:length(row_names));
	row_names <- row_names[which(index_vector %% 2 == 0)];
	
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

x_colname <- colnames(x);
x[is.na(x)] <- 0;

colnames(y) <- c("Gene_Name", "Read_Count");
rownames(x) <- rownames(y);


# output data
out_expected_file <- paste(dir, "/", type, "_lasso_expected_Y.txt", sep="");
out_training_file <- paste(dir, "/", type, "_lasso_matrix_X.txt", sep="");

write.table(x, out_training_file, sep="\t");
write.table(y, out_expected_file, sep="\t");


print(paste("Written the Expected Value (Y) to ", out_expected_file, sep=""));
print(paste("Written the Training Matrix (X) to ", out_training_file, sep=""));
