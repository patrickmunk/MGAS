#!/usr/bin/env Rscript
#This script is used to rarify a count matrix and perform multi-dimensional scaling analysis

#Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL,
              help="read count matrix file name", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="file name prefix [default= %default]", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Determine the prefix to be used for output files
if (is.null(opt$prefix)){
        prefix = cat(gsub(" ","_",date()))
        } else {
                prefix <- opt$prefix
        }

##Check if matrix is supplied
if (is.null(opt$matrix)){
  print_help(opt_parser)
  stop("A matrix with read counts needs to be supplied", call.=FALSE)
}
df = read.table(opt$matrix, header=TRUE, sep="\t")


################# PROGRAM START #################

#Make rarefaction curve
#Calculate step size for rarefaction. Smallest of 1/4 of minimum sample, 1/100 of max or 1000.
step_size <- round(min(c((min(colSums(df)/4)),(max(colSums(df)/100)),1000)), digits=0)

library(vegan)

#Produce rarefaction plot
filename <- paste(prefix,"rarecurve.pdf",sep="_")
pdf(file=filename, height=10, width=10)
rarecurve(t(df),step=step_size,col="blue",ylab="Sample richness",xlab="Subsampled counts")
dev.off()

#Produce common ecology diversity and richness indeces for samples
df <- t(df)
div_df <- data.frame(shannon=diversity(df, index="shannon"), simpson=diversity(df, index="simpson"), fisher.alpha=fisher.alpha(df), richness=specnumber(df), rarefy_min_count=rarefy(df, sample=min(rowSums(df))))
div_df <- div_df[order(div_df$shannon, decreasing = T),]

#Export table with diversity measures
filename <- paste(prefix,"diversity.txt",sep="_")
write.table(div_df, file=filename, row.names=TRUE, quote=F, sep="\t")
