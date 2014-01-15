############################# lasso_analysis.R ############################
# Function: read in the fragment experiment matrix,                       #
#         matrix. Use Lasso to train on the cov matrix.                   #
# Usage: R --no-save < lasso_analysis.R --args input expected             #
# Arguments: input =  gene_experiment_fragment_matrix.txt                 #
#            expected = gene_experiment_expected_matrix.txt               #
#            output = 
# Author: Chelsea Ju                                                      #
# Date: 2013-12-04                                                        #
###########################################################################

library(genlasso)
library(plyr)


read_fragment_matrix <- function(filename){

	data <- read.table(filename, fill = T);
	rowname <- data$V1;
	rowname <- c("", as.character(rowname[2:length(rowname)]));
	colname <- t(data[1,1:(ncol(data)-2)]);
	colname <- c("Expected", colname);
	
	# assign rowname
	rownames(data) <- rowname;
	data <- data[,-1];
	
	# assign colname
	colnames(data) <- colname;
	data <- data[-1,];

	# filter data: remove all zero rows
	data <- data[data$Expected != 0,];
	
	# filter data: remove all zero columns
	tmp_data <- as.matrix(data);
	matrix <- matrix(as.numeric(tmp_data), ncol = ncol(tmp_data), nrow = nrow(tmp_data)); ## convert to numeric value
	tmp_data <- matrix[, !(colSums(matrix) == 0)];
	
	tmp_data <- as.data.frame(tmp_data);
	colnames(tmp_data) <- colname[!(colSums(matrix) == 0)];
	rownames(tmp_data) <- rownames(data);
	
	data <- tmp_data;	
	data;
	
}

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

frag_30X10A_file <- "30X_100L_10A/human72_filter_transcript_all.fragment";
data_30X10A <- read_fragment_matrix(frag_30X10A_file);

frag_30X30A_file <- "30X_100L_30A/human72_filter_transcript_all.fragment";
data_30X30A <- read_fragment_matrix(frag_30X30A_file);

frag_30X50A_file <- "30X_100L_50A/human72_filter_transcript_all.fragment";
data_30X50A <- read_fragment_matrix(frag_30X50A_file);


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


