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

#Set parameters
min_depth <- 10000

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
df = read.delim(opt$matrix, header=TRUE, sep="\t")

#Choose to only include samples with more than x read counts
df <- df[,colSums(df)>min_depth]

#Rarefies matrix to the level of the sample with least data
library(vegan)
df <- t(rrarefy(t(df),sample=min(colSums(df))))

#Export rarefied data matrix (obs: low samples excluded)
filename <- paste(prefix,"rarefy_matrix.txt",sep="_")
write.table(df, file=filename, row.names=TRUE, quote=F, sep="\t")

#Remove rows with too little data (median = 0)
df <- df[apply(df,1,median)>0,]

#Perform non-metric multidimensional scaling
data_NMDS=metaMDS(df, k=3, trymax=50)

#Make a stressplot
filename <- paste(prefix,"stressplot.pdf",sep="_")
pdf(file=filename)
stressplot(data_NMDS)
dev.off()

#Plot the NMDS
filename <- paste(prefix,"NMDS.pdf",sep="_")
pdf(file=filename)
plot(data_NMDS)
dev.off()

#Find most abundant features to label
df_sort <- df[order(rowSums(df),decreasing=T),]
plot_names <- rownames(df_sort)[1:min(c(25,nrow(df_sort)))]
plot_names <- rownames(df) %in% plot_names

#Plot more NMDS
filename <- paste(prefix,"NMDS_anno.pdf",sep="_")
pdf(file=filename)
ordiplot(data_NMDS)
orditorp(data_NMDS,display="species",col="red",air=0.01)
orditorp(data_NMDS,display="sites",col="blue",air=0.01, cex=0.5, select=plot_names)
dev.off()
