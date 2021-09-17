library(ggplot2)
data_dir <- "~/Projects/mcgaha/teresa_metagenomics/16s/diversity/data/"
setwd(data_dir)

plotBox <- function(eve_df, value){
  if(any(grepl("Undetermined", eve_df$id))){
    eve_df <- eve_df[-grep("Undetermined", eve_df$id),]
  }
  eve_df$group <- gsub("-.*", "", eve_df$group.day)
  eve_df$day <- gsub("^.*-", "", eve_df$group.day)
  eve_df$day <- factor(eve_df$day, levels=as.character(sort(as.integer(unique(eve_df$day)))))
  ggp <- ggplot(data=eve_df, aes_string(x='day', y=value, fill='group')) +
    geom_boxplot(position=position_dodge(0.8)) +
    theme_minimal()
  tukey <- TukeyHSD(aov(data=eve_df, as.formula(paste0(value, "~group.day"))))$group.day
  tukey <- tukey[order(tukey[,'p adj']),]
  return(list("gg"=ggp, "stat"=tukey))
}

## Pielous Evennes
evef <- read.table(file.path(data_dir, "evenness_f.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)
ever <- read.table(file.path(data_dir, "evenness_r.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)

gg_evef <- plotBox(evef, value='pielou_evenness')
gg_ever <- plotBox(ever, value='pielou_evenness')
print(gg_evef$gg)
View(round(gg_evef$stat,3))
print(gg_ever$gg)
View(round(gg_ever$stat,3))

## Faiths Phylogenetic Diversity
fpdf <- read.table(file.path(data_dir, "faithpd_f.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)
fpdr <- read.table(file.path(data_dir, "faithpd_r.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)

gg_fpdf <- plotBox(fpdf, value='faith_pd')
gg_fpdr <- plotBox(fpdr, value='faith_pd')
print(gg_fpdf$gg)
View(round(gg_fpdf$stat,3))
print(gg_fpdr$gg)
View(round(gg_fpdr$stat,3))
