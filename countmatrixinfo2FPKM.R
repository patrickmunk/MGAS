#!/usr/bin/env Rscript
#Takes a raw abundance matrix (am) and sample data (sd) with sequencing depths and calculates FPKM values
#The 1st-3rd col of the am is feature name, description and feature size
#1st column of sd is sample name, 2nd is read count

#Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL,
              help="read count matrix file name", metavar="character"),
        make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="file name prefix [default= %default]", metavar="character"),
        make_option(c("-s", "--sample_data"), type="character", default=NULL,
              help="File with sample names (1st col) and sequencing depth (read counts, 2nd col)", metavar="character")
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


################# PROGRAM START #################

df = read.delim(opt$matrix, header=TRUE, check.names = F, sep="\t")

#Get the sample data with size factors (if not supplied, will be normalized to matrix columns)
if (is.null(opt$sample_data)){
        metadata <- data.frame(sample_name=colnames(df), read_count=colSums(df))
        } else {
                metadata = read.delim(opt$sample_data, header=TRUE, check.names = F)
        }

#Set some parameters - will implement user-specified later
insert_size <- 300
min_alignment <- 50

#Read in matrix file a la Oksana and save descriptions to a 2nd df
abundance_matrix <- df
rownames(abundance_matrix) <- abundance_matrix[,1]
feature_description <- abundance_matrix[,1:3]
abundance_matrix <- abundance_matrix[,4:ncol(abundance_matrix)]

#If read counts are considered factors, change them to numeric
indx <- sapply(abundance_matrix, is.factor)
abundance_matrix[indx] <- lapply(abundance_matrix[indx], function(x) as.numeric(as.character(x)))

#If when divided by two, no values become fractions, assume paired-end mapping
if ((all(abundance_matrix/2)%%1==0) == TRUE) {
  print("It looks like PE-mapping was used. Converting to #fragments")
  abundance_matrix <- abundance_matrix/2
  mapping_ends <- 2
}

#Make sure that metadata and matrix samples are in same order
if (all(metadata[,1] == colnames(abundance_matrix))) {
  print("Samples are in the same order in matrix and metadata. Good.")
} else {
  print("Samples in metadata and matrix are out of sync. Will try to fix it.")
  #Make sure the objects contain the same samples - subset if required
  metadata <- metadata[metadata[,1] %in% colnames(abundance_matrix),]
  abundance_matrix <- abundance_matrix[,colnames(abundance_matrix) %in% metadata[,1]]
  #Ensure the order is the same in the two objects
  metadata <- metadata[sort(metadata[,1]),]
  abundance_matrix <- abundance_matrix[,sort(colnames(abundance_matrix)),]
}

#Normalize each column with the total read count of a given sample (Reads per million (RPM))
if (mapping_ends == 2) {
  metadata[,2] <- metadata[,2]/2
}
abundance_matrix_rpm <- sweep(abundance_matrix, 2, metadata[,2], `/`)*10^6

#Normalize each gene row with the length of the reference gene minus the average fragment size
feature_description$hit_probability <- as.numeric(as.character(feature_description[,3])) - (insert_size-(mapping_ends*min_alignment))
abundance_matrix_rpm <- sweep(abundance_matrix_rpm, 1, feature_description$hit_probability, `/`)*1000

#Write out FPKM values
write.table(round(abundance_matrix_rpm, digits=4), file=paste(prefix, "FPKM_matrix.txt", sep="_"), quote=F, sep="\t")
