#### mergeCountsMetaphlan.R
# Nov 5 2021
# Description: This piece of code is meant to replicate the metaphlan python
# script merge_metaphlan_tables.py, but instead of merging and labelling
# the relative abundances, it will merge the estimated_number_of_reads_from_the_clade
# column that is outputted when metaphlan is ran with the -t rel_ab_w_read_stats
# flag

library(optparse)

option_list <- list(
  make_option(c("-d", "--dir"), default="",
              help="Directory which houses the individual metaphlan outputs"),
  make_option(c("-o", "--out"), default="output.tsv",
              help = "output file name [default \"%default\"]"),
  make_option(c("-p", "--pattern"), default="*s1.tsv$",
              help = "Pattern to detect all outputs of metaphlan tsv files [default \"%default\"]")
)
opt <- parse_args(OptionParser(option_list=option_list))
basedir <- getwd()
setwd(opt$dir)

# Read in all the metaphlan output files
files <- list.files(pattern=opt$pattern) #pattern <- "*s1.tsv$"
metas <- lapply(files, function(f){
  metaphlan <- read.table(f, sep="\t", comment.char='', skip=4, header=T,
                          check.names = F, stringsAsFactors = F)
  colnames(metaphlan)[1] <- gsub("^#", "", colnames(metaphlan)[1])
  metaphlan[,c('clade_name', 'clade_taxid', 'estimated_number_of_reads_from_the_clade')]
})

# Merge the individual files by the clade name and relabel the columns by sample names
meta_df <- Reduce(function(x,y) merge(x,y,by=c('clade_name', 'clade_taxid'), all=T), metas)
colnames(meta_df)[-c(1:2)] <- gsub("(.s1)?.tsv$", "", files)

#output
write.table(meta_df, file=file.path(basedir, opt$out), sep="\t",
            col.names = T, row.names = F, quote = F)
