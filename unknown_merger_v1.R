########################## unknown_merger_v1.R ############################
# Function: read in distribution matrix X and a list of column to merge   #
# Usage: R --no-save < unknown_merger_v1.R --args dir datatype            #
# Arguments: dir = input directory                                        #
#	(select_one_pseudogene_110/30X_101L_10A/tophat_out)                   #
#            datatype = [gene|transcript]                                 #
# Author: Chelsea Ju                                                      #
# Date: 2014-05-10                                                        #
###########################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 2){
    stop(paste("Invalid Arguments\n",
    "Usage: R--no-save --slave < unknown_merger.R --args dir type\n",
    "\t dir = directory of input and output\n",
    "\t type = genes or transcripts\n",
    sep=""));
}

dir <- options[1];
type <- options[2];

column_file <- paste(dir,"/", type, "_distribution.matrix.mergelist", sep="");
old_matrix_file <- paste(dir,"/", type, "_distribution.backup", sep="");
new_matrix_file <- paste(dir,"/", type, "_distribution.matrix", sep="");

column_index <- read.table(column_file, sep="\n");
old_matrix <- read.table(old_matrix_file, header = T);


# update column
for (i in 1:nrow(column_index)){
	positions <- unlist(strsplit(as.character(column_index[i,]), "\t"));
	name <- positions[1];
	sum_column <- rep(0, nrow(old_matrix));

	for (j in 2:length(positions)){
		c_index <- as.numeric(positions[j]);
		sum_column <- sum_column + old_matrix[,c_index];
		old_matrix[,c_index] <- rep(-1, nrow(old_matrix)); # flag for removal
	}

	# update old_matrix
	old_matrix[,as.numeric(positions[2])] <- sum_column;
	colnames(old_matrix)[as.numeric(positions[2])] <- name;
}

# remove flagged column
remove_column_index <- which(colSums(old_matrix) < 0);
if(length(remove_column_index)){
	old_matrix <- old_matrix[,-remove_column_index];
}

# write to file
write.table(old_matrix, new_matrix_file, sep="\t");

