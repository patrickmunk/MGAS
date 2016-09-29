#!/usr/bin/env Rscript
#This script is used to produce relative abundance matrices from raw read counts

#Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL,
              help="read count matrix file name", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="file name prefix [default= %default]", metavar="character"),
        make_option(c("-s", "--sizefactor"), type="character", default=NULL,
              help="File with read counts per sample (one column) in the same order as samples are rows in the matrix. If none is supplied, the matrix column sum will be used instead.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##Check if matrix is supplied
if (is.null(opt$matrix)){
  print_help(opt_parser)
  stop("A matrix with read counts needs to be supplied.n", call.=FALSE)
}
df = read.table(opt$matrix, header=TRUE)

#Get the size factors for the samples (either supplied or calculated from matrix)
if (is.null(opt$size_factor)){
        size_factor <- colSums(df)
        } else {
                size_factor = read.table(opt$sizefactor, header=TRUE)
                size_factor <- size_factor[,1]
        }

#Calculate relative abundance from the matrix and size factor
df_rpm <- sweep(df,2,size_factor,"/")*10^6
df_perc <- sweep(df,2,size_factor,"/")*10^2

#Determine the prefix to be used for output files
if (is.null(opt$prefix)){
        prefix = cat(gsub(" ","_",date()))
        } else {
                prefix <- opt$prefix
        }

#Export matrix with normalized counts (RPM)
filename <- paste(prefix,"abundances_RPM.txt",sep="_")
write.table(df_rpm, file=filename, row.names=TRUE, quote=F, sep="\t")

#Export matrix with normalized counts (percent)
filename <- paste(prefix,"abundances_percent.txt",sep="_")
write.table(df_perc, file=filename, row.names=TRUE, quote=F, sep="\t")

#Produce some summary statistics
sumdf <- data.frame(reference=rownames(data.frame(df_rpm)),sd=apply(df_rpm,1,sd), mean=apply(df_rpm,1,mean), var=apply(df_rpm,1,var), median=apply(df_rpm,1,median))
write.table(sumdf, file=paste(prefix,"reference_summary.txt",sep="_"), row.names=TRUE, quote=F, sep="\t")
