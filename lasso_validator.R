################################## lasso_validator.R ###################################
# Function: Y = Xb                                                                     #
#         use the trained coefficient "b" to predict the fragment counts Y'.           #
#          is compared with the expected value.                                        #
# Usage: R --no-save < lasso_analysis.R --args input expected                        #
# Arguments: input =  gene_experiment_fragment_matrix.txt                            #
#            expected = gene_experiment_expected_matrix.txt                          #
#            output = lasso_coefficient.txt				                             # 
# Author: Chelsea Ju                                                                 #
# Date: 2013-12-04                                                                   #
# Last Modify: 2013-12-22                                                            #
######################################################################################

library(genlasso)
library(plyr)

setwd("/Users/Chelsea/Bioinformatics/WangLab/Pseudogene/select_one_pseudogene_110");

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


y <- as.data.frame(c(data_10X10A[,1], data_10X30A[,1], data_10X50A[,1], data_20X10A[,1], data_20X30A[,1], data_20X50A[,1], data_30X10A[,1], data_30X30A[,1], data_30X50A[,1]));
rownames(y) <- c(paste("10X10A_",rownames(data_10X10A), sep=""), paste("10X30A_",rownames(data_10X30A), sep=""), paste("10X50A_",rownames(data_10X50A), sep=""),
                 paste("20X10A_",rownames(data_20X10A), sep=""), paste("20X30A_",rownames(data_20X30A), sep=""), paste("20X50A_",rownames(data_20X50A), sep=""),
                 paste("30X10A_",rownames(data_30X10A), sep=""), paste("30X30A_",rownames(data_30X30A), sep=""), paste("30X50A_",rownames(data_30X50A), sep=""))

x <- rbind.fill(data_10X10A[,2:ncol(data_10X10A)],data_10X30A[,2:ncol(data_10X30A)], data_10X50A[,2:ncol(data_10X50A)],
                data_20X10A[,2:ncol(data_20X10A)],data_20X30A[,2:ncol(data_20X30A)], data_20X50A[,2:ncol(data_20X50A)],
                data_30X10A[,2:ncol(data_30X10A)],data_30X30A[,2:ncol(data_30X30A)], data_30X50A[,2:ncol(data_30X50A)]) ;


## take out the last three experiment for validataion (30X..)
trained_y <- y[1:660,];
validate_y <- y[661:990,];

trained_x <- x[1:660,];
validate_x <- x[661:990,];

D <- diag(1, ncol(trained_x));

trained_x_colname <- colnames(trained_x);
trained_x_rowname <- c(paste("10X10A",rownames(data_10X10A), sep="_"), paste("10X30A",rownames(data_10X30A), sep="_"), paste("10X50A",rownames(data_10X50A), sep="_"),
				paste("20X10A",rownames(data_20X10A), sep="_"), paste("20X30A",rownames(data_20X30A), sep="_"), paste("20X50A",rownames(data_20X50A), sep="_"));
trained_x[is.na(trained_x)] <- 0;

validate_x_colname <- colnames(validate_x);
validate_x_rowname <- c(paste("30X10A",rownames(data_30X10A), sep="_"), paste("30X30A",rownames(data_30X30A), sep="_"), paste("30X50A",rownames(data_30X50A), sep="_"));
validate_x[is.na(validate_x)] <- 0;


## run lasso
out <- genlasso(trained_y, X=as.matrix(trained_x), D=D);
summary(out);

beta <- coef(out, lambda=sqrt(nrow(trained_x) * log(ncol(trained_x))))
beta_value <- beta$beta;
rownames(beta_value) <- trained_x_colname;

write.table(paste("lambda", sqrt(nrow(trained_x) * log(ncol(trained_x))), sep="\t"), file="lasso_coefficient.txt")
write.table(beta_value[order(beta_value),], file="lasso_coefficient.txt", append = T);


## validate analysis
predict_y <- as.matrix(validate_x) %*% as.matrix(beta_value);

validate_y <- cbind(validate_y, predict_y, abs(predict_y - validate_y), abs(predict_y - validate_y)/validate_y);
rownames(validate_y) <- validate_x_rowname;

colnames(validate_y) <- c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");

write.table(validate_y, file="lasso_prediction.xls", sep="\t");

## validate using the real scenario: not knowing where the reads comes from (sum up column)
validate_x_10A <-  colSums(validate_x[1:110,]);
validate_x_30A <-  colSums(validate_x[111:220,]);
validate_x_50A <-  colSums(validate_x[221:330,]);

predict_y_10A <- t(as.matrix(validate_x_10A)) %*% diag(c(beta_value), length(beta_value));
predict_y_30A <- t(as.matrix(validate_x_30A)) %*% diag(c(beta_value), length(beta_value));
predict_y_50A <- t(as.matrix(validate_x_50A)) %*% diag(c(beta_value), length(beta_value));

predict_y_10A <- t(predict_y_10A);
predict_y_30A <- t(predict_y_30A);
predict_y_50A <- t(predict_y_50A);

rownames(predict_y_10A) <- paste("30X10A_",colnames(validate_x[1:110,]), sep="");
rownames(predict_y_30A) <- paste("30X30A_",colnames(validate_x[111:220,]), sep="");
rownames(predict_y_50A) <- paste("30X50A_",colnames(validate_x[221:330,]), sep="");

predict_y_10A<- cbind(c(y[rownames(predict_y_10A),]), predict_y_10A);
predict_y_10A[is.na(predict_y_10A)] <- 0;
predict_y_10A <- cbind(predict_y_10A, abs(predict_y_10A[,2] - predict_y_10A[,1]),  abs(predict_y_10A[,2] - predict_y_10A[,1])/predict_y_10A[,1]);

predict_y_30A<- cbind(c(y[rownames(predict_y_30A),]), predict_y_30A);
predict_y_30A[is.na(predict_y_30A)] <- 0;
predict_y_30A <- cbind(predict_y_30A, abs(predict_y_30A[,2] - predict_y_30A[,1]),  abs(predict_y_30A[,2] - predict_y_30A[,1])/predict_y_30A[,1]);

predict_y_50A<- cbind(c(y[rownames(predict_y_50A),]), predict_y_50A);
predict_y_50A[is.na(predict_y_50A)] <- 0;
predict_y_50A <- cbind(predict_y_50A, abs(predict_y_50A[,2] - predict_y_50A[,1]),  abs(predict_y_50A[,2] - predict_y_50A[,1])/predict_y_50A[,1]);

predict_y_all <- rbind(predict_y_10A, predict_y_30A, predict_y_50A);
colnames(predict_y_all) <-  c("Expected", "Predicted", "Absolute Error", "(Prediction - Expectation)/Expectation");

write.table(predict_y_all, file="lasso_prediction_v2.xls", sep="\t");

