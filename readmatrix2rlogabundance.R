#!/usr/bin/env Rscript
#This script is used to produce produce relative abundances from a matrix, using the DESeq2 rlog transformation
#It uses the blind setting (transformation is blind to experimental design). This works poorly if there are strong experimental effects

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
df = read.table(opt$matrix, header=TRUE)


################# PROGRAM START #################

library(DESeq2)

#Use DESeq2 for rlog transformation
abundance_matrix <- df
metadata <- data.frame(sample=cbind(colnames(abundance_matrix)),fakedata=c(1:ncol(abundance_matrix)),row.names=colnames(abundance_matrix))

dds <- DESeqDataSetFromMatrix(countData = abundance_matrix, colData = metadata, design =~1)
rlog_dds <- rlog(dds,blind=T)
relabundance <- assay(rlog_dds)

#Export rlog-transformed relative abundances
filename <- paste(prefix,"relabun_rlog.txt",sep="_")
write.table(relabundance, file=filename, sep="\t", quote=F)