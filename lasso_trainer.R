############################# lasso_trainer.R #############################
# Function: read in the fragment experiment matrix,                       #
#         matrix. Use Lasso to train on the cov matrix.                   #
# Usage: R --no-save < lasso_analysis.R --args dir type                   #
# Arguments: input =  directory for input and output                      #
#            type = genes or transcripts                                  #
#            output = 
# Author: Chelsea Ju                                                      #
# Date: 2014-01-14                                                        #
# Modify from lasso_analysis.R                                            #
###########################################################################

library(genlasso)
library(plyr)

## self-defined function
read_distribution_matrix <- function(subdir, filetype){
	matrix_file <- paste(subdir, "tophat_out/", filetype, "_distribution.matrix", sep="");
	matrix_data <- read.table(matrix_file);
	
	matrix_data;
}

read_expected_data <- function(subdir, filetype){
	expected_file <- paste(subdir, "tophat_out/", filetype, "_expected_read_count.txt", sep="");
	expected_data <- read.table(expected_file, header = FALSE);
	c1 <- expected_data$V1;
	row_names <- unlist(strsplit(as.character(c1), "_")); ## modify the rownames: split ENSG00000100478_ENST00000542754 into "ENSG00000100478" "ENST00000542754"
	index_vector <- c(1:length(row_names));
	row_names <- row_names[which(index_vector %% 2 != 0)];
	
	rownames(expected_data) <- row_names;
	expected_data <- expected_data[-1];
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

expected_10X10A <- read_expected_data(dir_10X10A, type);
expected_10X30A <- read_expected_data(dir_10X30A, type);
expected_10X50A <- read_expected_data(dir_10X50A, type);

expected_20X10A <- read_expected_data(dir_20X10A, type);
expected_20X30A <- read_expected_data(dir_20X30A, type);
expected_20X50A <- read_expected_data(dir_20X50A, type);


frag_10X10A_file <- "10X_100L_10A/human72_filter_transcript_all.fragment";
data_10X10A <- read_fragment_matrix(frag_10X10A_file);

frag_10X30A_file <- "10X_100L_30A/human72_filter_transcript_all.fragment";
data_10X30A <- read_fragment_matrix(frag_10X30A_file);

frag_10X50A_file <- "10X_100L_50A/human72_filter_transcript_all.fragment";
data_10X50A <- read_fragment_matrix(frag_10X50A_file);

frag_20X10A_file <- "20X_100L_10A/human72_filter_transcript_all.fragment";
data_20X10A <- read_fragment_matrix(frag_20X10A_file);

frag_20X30A_file <- "20X_100L_30A/human72_filter_transcript_all.fragment";
data_20X30A <- read_fragment_matrix(frag_20X30A_file);

frag_20X50A_file <- "20X_100L_50A/human72_filter_transcript_all.fragment";
data_20X50A <- read_fragment_matrix(frag_20X50A_file);


y <- c(data_10X10A[,1], data_10X30A[,1], data_10X50A[,1], data_20X10A[,1], data_20X30A[,1], data_20X50A[,1], data_30X10A[,1], data_30X30A[,1], data_30X50A[,1]);
x <- rbind.fill(data_10X10A[,2:ncol(data_10X10A)],data_10X30A[,2:ncol(data_10X30A)], data_10X50A[,2:ncol(data_10X50A)],
		data_20X10A[,2:ncol(data_20X10A)],data_20X30A[,2:ncol(data_20X30A)], data_20X50A[,2:ncol(data_20X50A)],
		data_30X10A[,2:ncol(data_30X10A)],data_30X30A[,2:ncol(data_30X30A)], data_30X50A[,2:ncol(data_30X50A)]) ;

D <- diag(1, ncol(x));

x_colname <- colnames(x);
x_rowname <- c(paste("10X10A",rownames(data_10X10A), sep="_"), paste("10X30A",rownames(data_10X30A), sep="_"), paste("10X50A",rownames(data_10X50A), sep="_"),
				paste("20X10A",rownames(data_20X10A), sep="_"), paste("20X30A",rownames(data_20X30A), sep="_"), paste("20X50A",rownames(data_20X50A), sep="_"),
				paste("30X10A",rownames(data_30X10A), sep="_"), paste("30X30A",rownames(data_30X30A), sep="_"), paste("30X50A",rownames(data_30X50A), sep="_"));
x[is.na(x)] <- 0;

out <- genlasso(y, X=as.matrix(x), D=D);
summary(out);

beta <- coef(out, lambda=sqrt(nrow(x) * log(ncol(x))))
beta_value <- beta$beta;
rowname(beta_value) <- x_colname;

write.table(paste("lambda", sqrt(nrow(x) * log(ncol(x))), sep="\t"), file="lasso_coefficient.xls")
write.table(beta_value[order(beta_value),], file="lasso_coefficient.xls", append = T);


