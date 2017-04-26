#!/usr/bin/env Rscript
#Takes a matrix with either raw counts or relative abundances and aggregates the rows into clusters
#Please make sure it makes sense to add abundances together (e.g. not rlog relative abundance)

#Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL,
              help="read count matrix file name", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="file name prefix [default= %default]", metavar="character"),
        make_option(c("-c", "--clustfile"), type="character", default=NULL,
              help="File with matrix reference names (1st col) and cluster names (2nd col) used to aggregate the matrix", metavar="character")
);

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
  stop("A matrix with abundances needs to be supplied", call.=FALSE)
}
df = read.delim(opt$matrix, header=TRUE, check.names = F, sep="\t", row.names=1)

##Check if cluster file is supplied
if (is.null(opt$matrix)){
  print_help(opt_parser)
  stop("A cluster file with reference and cluster names needs to be supplied", call.=FALSE)
}
clustfile = read.delim(opt$clustfile, header=TRUE, check.names = F, sep="\t")

################# PROGRAM START #################

library(data.table)

#Read in an abundance matrix (either raw or relative)
abundance_matrix <- df
cluster_ref <- clustfile
colnames(cluster_ref)[1:2] <- c("ref_name","clust_name")

#Remove any features without data
abundance_matrix <- abundance_matrix[rowSums(abundance_matrix)>0,]

#Determine if there are features in the data not included in the cluster file
not_in_cluster <- setdiff(rownames(abundance_matrix),cluster_ref$ref_name)
if (length(setdiff(rownames(abundance_matrix),cluster_ref$ref_name))>0) {
  print(c("Warning! The following references are in your data but not in the cluster file! These will not be aggregated and will be appended as-is to your output table.",not_in_cluster))
  not_in_cluster <- abundance_matrix[not_in_cluster,]
  not_in_cluster <- cbind(clust_name=rownames(not_in_cluster),not_in_cluster)
}

#Aggregate columns
abundance_matrix$ref_name <- rownames(abundance_matrix)
abundance_matrix <- merge(abundance_matrix, cluster_ref, by="ref_name")
abundance_matrix <- data.table(abundance_matrix[,2:ncol(abundance_matrix)])
abundance_matrix <- melt(abundance_matrix, id.vars = "clust_name")
abundance_matrix <- data.frame(dcast.data.table(abundance_matrix, clust_name ~ variable, value.var = "value", fun.aggregate = sum, na.rm = TRUE), check.names=FALSE)

#Write out table of aggregated abundances
filename <- paste(prefix,"FPKM_table.txt",sep="_")
write.table(abundance_matrix, file=filename, quote=F, sep="\t", row.names = F)

#Write out table of aggregated abundances with non-aggregated features appended
abundance_matrix_unclust <- rbind(abundance_matrix, not_in_cluster)
filename <- paste(prefix,"FPKM_table_all.txt",sep="_")
write.table(abundance_matrix_unclust, file=filename, quote=F, sep="\t", row.names = F)
