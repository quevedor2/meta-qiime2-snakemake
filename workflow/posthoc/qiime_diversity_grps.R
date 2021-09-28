library(ggplot2)
library(cowplot)
PDIR <- "~/Projects/mcgaha/teresa_metagenomics/16s/"
qc_dir <- file.path(PDIR, "qc/se/data")
data_dir <- file.path(PDIR, "diversity/data/")
setwd(data_dir)

###################
#### Functions ####
plotBox <- function(eve_df, value){
  if(any(grepl("Undetermined", eve_df$id))){
    eve_df <- eve_df[-grep("Undetermined", eve_df$id),]
  }
  eve_df$group <- gsub("-.*", "", eve_df$group.day)
  eve_df$day <- gsub("^.*-", "", eve_df$group.day)
  eve_df$day <- factor(eve_df$day, levels=as.character(sort(as.integer(unique(eve_df$day)))))
  ggp <- ggplot(data=eve_df, aes_string(x='day', y=value, fill='group')) +
    geom_boxplot(position=position_dodge(0.8)) +
    theme_minimal() +
    ylim(0,max(eve_df[,value], na.rm=TRUE)*1.25)
  tukey <- TukeyHSD(aov(data=eve_df, as.formula(paste0(value, "~group.day"))))$group.day
  tukey <- tukey[order(tukey[,'p adj']),]
  return(list("gg"=ggp, "stat"=tukey))
}

formatDat <- function(df, rm_samples=TRUE){
  day_map <- setNames(c(0,5,10,16), c('PreOp', 'Day5', 'Day10', 'Day16'))
  iddf <- as.data.frame(do.call(rbind, strsplit(df$sample.ID, split="-")))
  colnames(df)[1] <- 'id'
  colnames(iddf) <- c('study.id', 'day.raw', 'Mouse')
  iddf$day <- as.integer(day_map[as.character(iddf[,2])])
  iddf$group.day <- with(iddf, paste0(study.id, "-", day))
  
  df_upd <- as.data.frame(cbind(df, iddf))
  if(rm_samples){
    sample_filter <- c('Undetermined', 'HFD-Day16-1')
    rm_idx <- which(df_upd$id %in% sample_filter)
    
    if(length(rm_idx) > 0) df_upd <- df_upd[-rm_idx,]
  }
  df_upd
}

##############
#### Main ####
## QC Alpha plot
richf <- read.table(file.path(qc_dir, "f_teresa_metagenomics_demux.tsv"), h=T, 
                    stringsAsFactors = FALSE, sep="\t")
richr <- read.table(file.path(qc_dir, "r_teresa_metagenomics_demux.tsv"), h=T, 
                    stringsAsFactors = FALSE, sep="\t")

gg_richf <- plotBox(formatDat(richf), value='forward.sequence.count')
gg_richr <- plotBox(formatDat(richr), value='forward.sequence.count')

pdf(file.path(qc_dir, "richness.pdf"), height = 7)
plot_grid(gg_richf$gg, gg_richr$gg, nrow=2)
dev.off()

## Pielous Evennes
evef <- read.table(file.path(data_dir, "evenness_f.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)
ever <- read.table(file.path(data_dir, "evenness_r.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)

gg_evef <- plotBox(evef, value='pielou_evenness')
gg_ever <- plotBox(ever, value='pielou_evenness')

pdf(file.path(data_dir, "evenness.pdf"), height = 7)
plot_grid(gg_evef$gg, gg_ever$gg, nrow=2)
dev.off()

## Faiths Phylogenetic Diversity
fpdf <- read.table(file.path(data_dir, "faithpd_f.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)
fpdr <- read.table(file.path(data_dir, "faithpd_r.tsv"), header=TRUE, 
                   stringsAsFactors = FALSE)

gg_fpdf <- plotBox(fpdf, value='faith_pd')
gg_fpdr <- plotBox(fpdr, value='faith_pd')

plot_grid(gg_fpdf$gg, gg_fpdr$gg, nrow=2)

