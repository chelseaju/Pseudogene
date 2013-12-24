######################## model_expression_value.R ################################## 
# Usage: R --no-save --slave < model_expression_value.R --args INPUT OUT_PREFIX    #
# Input: expression (fragments and FPKM) of pseudogene and its parent              #
# Output: Linear regression of parents and pseudo fragments, and plots             #
# Function: use linear regression model to estimate the fragments mapped to parents#
#           and pseudogene                                                         # 
# Author: Chelsea Ju                                                               #
# Date: 2013-08-04                                                                 #
####################################################################################

library(ggplot2)

# define function
# source from http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph

fpkm_eqn = function(data){
    m = lm(parent_fpkm ~ pseudo_fpkm, data);
 	eq <- substitute(italic(y) == italic(a) + (italic(b)) %.% italic(x)*","~~italic(r)^2~"="~italic(r2), 
         			 list(a = format(coef(m)[1], digits = 2), 
              			  b = format(coef(m)[2], digits = 2), 
             			  r2 = format(summary(m)$r.squared, digits = 3)
             			  )
             		)
	as.character(as.expression(eq));                 
}

fragment_eqn = function(data){
    m = lm(parent_fragment ~ pseudo_fragment, data);
 	eq <- substitute(italic(y) == italic(a) + (italic(b)) %.% italic(x)*","~~italic(r)^2~"="~italic(r2), 
         			 list(a = format(coef(m)[1], digits = 2), 
              			  b = format(coef(m)[2], digits = 2), 
             			  r2 = format(summary(m)$r.squared, digits = 3)
             			  )
             		)
	as.character(as.expression(eq));                 
}


# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 2){
	stop(paste("Invalid Arguments\n",
			   "Usage: R--no-save --slave < model_expression_value.R --args expression output_prefix\n",
			   "\t expression = file that combine the exprssion value of pseudogenes and parent genes\n",
			   "\t output_prefix = prefix for the output file\n"),
			   sep="");
}


in_file <- options[1];
out_prefix <- options[2];
information <- paste(out_prefix, "_expression_fact.txt", sep="")
plot_processed_fpkm <- paste(out_prefix, "_fpkm_processed.png", sep="");
plot_duplicated_fpkm <- paste(out_prefix, "_fpkm_duplicated.png", sep="");
plot_ambiguous_fpkm <- paste(out_prefix, "_fpkm_ambiguous.png", sep="");
plot_all_fpkm <- paste(out_prefix, "_fpkm_all.png", sep="");


plot_processed_fragment <- paste(out_prefix, "_fragment_processed.png", sep="");
plot_duplicated_fragment <- paste(out_prefix, "_fragment_duplicated.png", sep="");
plot_ambiguous_fragment <- paste(out_prefix, "_fragment_ambiguous.png", sep="");
plot_all_fragment <- paste(out_prefix, "_fragment_all.png", sep="");


## read in data
data <- read.table(in_file, header = TRUE);
process_data <- data[data$pseudo_type == "Processed",]
duplicate_data <- data[data$pseudo_type == "Duplicated",]
ambiguous_data <- data[data$pseudo_type == "Ambiguous",]

data_nozero <- data[data$pseudo_fpkm > 0,]
process_data_nozero <- process_data[process_data$pseudo_fpkm > 0,]
duplicate_data_nozero <- duplicate_data[duplicate_data$pseudo_fpkm > 0,]
ambiguous_data_nozero <- ambiguous_data[ambiguous_data$pseudo_fpkm > 0,]

## plot range
#fpkm_xrange <- c(0,1.5**ceiling(log(max(data$pseudo_fpkm),1.5)))
#fpkm_yrange <- c(0,1.5**ceiling(log(max(data$parent_fpkm),1.5)))
#fragment_xrange <- c(0,1.5**ceiling(log(max(data$pseudo_fragment),1.5)))
#fragment_yrange <- c(0,1.5**ceiling(log(max(data$parent_fragment),1.5)))


## text position
#fpkm_xtext <- fpkm_xrange[2] / 2
#fpkm_ytext <- floor(fpkm_yrange[2])
#fragment_xtext <- fragment_xrange[2] / 2
#fragment_ytext <- floor(fragment_yrange[2])


## plot graph

all_fpkm <- ggplot(data = data, aes(x = pseudo_fpkm, y = parent_fpkm)) +
		            geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        		    geom_point() +
        		    xlab("Pseudogene FPKM") +
        		    ylab("Parent Genes FPKM") +
					ggtitle("FPKM of All Psuedogene") +
        		    geom_text(aes(x = (max(data$pseudo_fpkm) / 2), y = (max(data$parent_fpkm)-1000), label = fpkm_eqn(data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")

all_fragment <- ggplot(data =data, aes(x = pseudo_fragment, y = parent_fragment)) +
				  	  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
				  	  geom_point() +
				  	  xlab("Pseudogene Fragments") +
				  	  ylab("Parent Genes Fragments") +
					  ggtitle("Fragments of All Psuedogene") +
				  	  geom_text(aes(x = (max(data$pseudo_fragment) / 2), y = (max(data$parent_fragment)-1000), label = fragment_eqn(data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")


processed_fpkm <- ggplot(data = process_data, aes(x = pseudo_fpkm, y = parent_fpkm)) +
		            geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        		    geom_point() +
        		    xlab("Pseudogene FPKM") +
        		    ylab("Parent Genes FPKM") +
					ggtitle("FPKM of Processed Psuedogene") +
        		    geom_text(aes(x = (max(process_data$pseudo_fpkm) / 2), y = (max(process_data$parent_fpkm)-1000), label = fpkm_eqn(process_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")
 
processed_fragment <- ggplot(data = process_data, aes(x = pseudo_fragment, y = parent_fragment)) +
				  	  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
				  	  geom_point() +
				  	  xlab("Pseudogene Fragments") +
				  	  ylab("Parent Genes Fragments") +
					  ggtitle("Fragments of Processed Psuedogene") +
				  	  geom_text(aes(x = (max(process_data$pseudo_fragment) / 2), y = (max(process_data$parent_fragment)-1000), label = fragment_eqn(process_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")


duplicated_fpkm <- ggplot(data = duplicate_data, aes(x = pseudo_fpkm, y = parent_fpkm)) +
  		            geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
          		    geom_point() +
          		    xlab("Pseudogene FPKM") +
          		    ylab("Parent Genes FPKM") +
 					ggtitle("FPKM of Duplicated Psuedogene") +
          		    geom_text(aes(x = (max(duplicate_data$pseudo_fpkm) / 2), y = (max(duplicate_data$parent_fpkm)-1000), label = fpkm_eqn(duplicate_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")
					
duplicated_fragment <- ggplot(data = duplicate_data, aes(x = pseudo_fragment, y = parent_fragment)) +
 				  	  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
 				  	  geom_point() +
				  	  xlab("Pseudogene Fragments") +
 				  	  ylab("Parent Genes Fragments") +
 	  				  ggtitle("Fragments of Duplicated Psuedogene") +
 	        		  geom_text(aes(x = (max(duplicate_data$pseudo_fragment) / 2), y = (max(duplicate_data$parent_fragment)-1000), label = fpkm_eqn(duplicate_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")

        	
ambiguous_fpkm <- ggplot(data = ambiguous_data, aes(x = pseudo_fpkm, y = parent_fpkm)) +
		            geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        		    geom_point() +
        		    xlab("Pseudogene FPKM") +
        		    ylab("Parent Genes FPKM") + 
        		    ggtitle("FPKM of Ambiguous Psuedogene") +
        		    geom_text(aes(x = (max(ambiguous_data$pseudo_fpkm) / 2), y = (max(ambiguous_data$parent_fpkm)-1000), label = fpkm_eqn(ambiguous_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red")



ambiguous_fragment <- ggplot(data = ambiguous_data, aes(x = pseudo_fragment, y = parent_fragment)) +
					geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        		    geom_point() +
        		    xlab("Pseudogene Fragments") +
        		    ylab("Parent Genes Fragments") + 
        		    ggtitle("Fragments of Ambiguous Psuedogene") +
        		    geom_text(aes(x = (max(ambiguous_data$pseudo_fragment) / 2), y = (max(ambiguous_data$parent_fragment)-1000), label = fragment_eqn(ambiguous_data)), parse = TRUE, family="Times", size=4, lineheight=.8, color="red") 
 
# print to file

png(plot_all_fpkm, width=800, height=800)
all_fpkm
dev.off()

png(plot_all_fragment, width=800, height=800)
all_fragment
dev.off()


png(plot_processed_fpkm, width=800, height=800)
processed_fpkm
dev.off()

png(plot_processed_fragment, width=800, height=800)
processed_fragment
dev.off()


png(plot_duplicated_fpkm, width=800, height=800)
duplicated_fpkm
dev.off()

png(plot_duplicated_fragment, width=800, height=800)
duplicated_fragment
dev.off()


png(plot_ambiguous_fpkm, width=800, height=800)
ambiguous_fpkm
dev.off()

png(plot_ambiguous_fragment, width=800, height=800)
ambiguous_fragment
dev.off()

write.table("Genes with Zero Fragment:\n", file = information, quote = F, col.names = F, row.names = F)
write.table(data[data$parent_fragment ==0,], file = information, append = T, quote = F)
write.table("\n\n", file = information, quote = F, col.names = F, row.names = F, append = T)
write.table("Fragments Mapped to Both Pseudogene and Parent Genes:", file = information, quote = F, col.names = F, row.names = F, append = T)
write.table(data[data$parent_secondary_fragment >0,], file = information, append = T, quote = F)
write.table("\n\n", file = information, quote = F, col.names = F, row.names = F, append = T)


print(paste("Write to ", plot_all_fpkm, sep=""))
print(paste("Write to ", plot_processed_fpkm, sep=""))
print(paste("Write to ", plot_duplicated_fpkm, sep=""))
print(paste("Write to ", plot_ambiguous_fpkm, sep=""))
print(paste("Write to ", plot_all_fragment, sep=""))
print(paste("Write to ", plot_processed_fragment, sep=""))
print(paste("Write to ", plot_duplicated_fragment, sep=""))
print(paste("Write to ", plot_ambiguous_fragment, sep=""))

print(paste("Write to ", information, sep=""))
