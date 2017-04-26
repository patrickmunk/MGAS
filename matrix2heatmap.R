#!/usr/bin/env Rscript
#This script is used to produce a pretty heatmap from a matrix of relative abundances

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
df = read.table(opt$matrix, header=TRUE, sep="\t", check.names=F, row.names=1)


################# PROGRAM START #################

library(pheatmap)

#Select subset of data for plotting
df_sort <- df[order(rowSums(df),decreasing=T),]
df <- df[apply(df,1,median)>0,]
plot_names <- rownames(df_sort)[1:min(c(50,nrow(df_sort)))]
df <- df[rownames(df) %in% plot_names,]

#Calculate appropriate dimensions for exported figure
hm_height <- round((nrow(df)/8)+4, digits=0)
hm_width <- round((ncol(df)/5)+4, digits=0)

#Produce heatmaps
filename <- paste(prefix,"heatmap.pdf",sep="_")
pdf(file=filename, height=hm_height, width=hm_width, onefile=FALSE)
pheatmap(df, margins=c(8,8), scale="none", treeheight_row = 100, treeheight_col = 100)
dev.off()

filename <- paste(prefix,"log10_heatmap.pdf",sep="_")
pdf(file=filename, height=hm_height, width=hm_width, onefile=FALSE)
pheatmap(log10(df+1), margins=c(8,8), scale="none", treeheight_row = 100, treeheight_col = 100)
dev.off()

filename <- paste(prefix,"log10_corclust_heatmap.pdf",sep="_")
pdf(file=filename, height=hm_height, width=hm_width, onefile=FALSE)
pheatmap(log10(df+1), margins=c(8,8), scale="none", treeheight_row = 100, treeheight_col = 100, clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()

filename <- paste(prefix,"rowscale_heatmap.pdf",sep="_")
pdf(file=filename, height=hm_height, width=hm_width, onefile=FALSE)
pheatmap(log10(df+1), margins=c(8,8), scale="row", treeheight_row = 100, treeheight_col = 100)
dev.off()

