library(reshape2)
library(ggplot2)
library(cowplot)

PDIR <- '~/Projects/mcgaha/teresa_metagenomics/16s/'
setwd(PDIR)

getTaxa <- function(taxas, lvl="species", split_regex="\\."){
  pattern <- switch(lvl,
                    genus="g__",
                    phylla="p__",
                    class="c__",
                    order="o__",
                    family="f__",
                    species="s__")
  taxas <- gsub("([;.][gpcofs]__)uncultured_organism", "\\1", taxas)
  taxas <- gsub("([;.][gpcofs]__)uncultured_bacterium", "\\1", taxas)
  taxas <- gsub("([;.][gpcofs]__)uncultured", "\\1", taxas)
  tids <- sapply(strsplit(taxas, split=split_regex), function(i){
    if(any(grepl(pattern, i))) i[grep(pattern, i)] else NA
  })
  na_idx <- (is.na(tids) | grepl("^[sfocpg]__$", tids))
  tids[which(na_idx)] <- sapply(strsplit(taxas[which(na_idx)], split_regex), function(i){
    useful_id <- min(which(i == "__" | grepl("^[sfocpg]__$", i))) - 1
    paste0("X_", paste0(i[useful_id:length(i)], collapse=";"))
  })
  return(tids)
}

drop_samples <- c("HFD-Day16-1", "Undetermined")
imp_n <- 4
lvl_f     <- 'r_level-7'
lvl_taxa  <- if(gsub("^.*-", "", lvl_f) == '6') 'genus' else 'species'

################################
#### Get Important Features ####
longitudinal <- read.table(file.path('longitudinal', 'data', paste0(lvl_f, ".tsv")), header=TRUE,
                           sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
longitudinal$simple_id <- getTaxa(longitudinal$id, lvl_taxa, split_regex = ";")

# order by importance & select top X
longitudinal      <- longitudinal[order(longitudinal$importance, decreasing = TRUE),]
longitudinal_sel  <- longitudinal[c(1:(imp_n+1)),]
imp_feats         <- longitudinal_sel$simple_id


########################################
#### Preprocess the taxa data frame ####
# Filter out samples
taxa <- read.csv(file.path('taxonomy', 'data', paste0(lvl_f, ".csv")))
if(length(drop_samples)>1) taxa <- taxa[-which(taxa$index %in% drop_samples),]

# Separate the taxa sample metadata from abundance matrix
taxa_meta_idx <- grep("^d__|Unassigned", colnames(taxa), invert = TRUE)
taxa_meta <- taxa[,taxa_meta_idx]
taxa_mat  <- taxa[,-taxa_meta_idx]

# Label the abundance matrix
rownames(taxa_mat) <- taxa_meta$index
simple_ids <- getTaxa(colnames(taxa_mat), lvl_taxa)
if(any(is.na(simple_ids))) simple_ids[is.na(simple_ids)] <- paste0("Unk_", c(1:sum(is.na(simple_ids))))
colnames(taxa_mat) <- simple_ids

# Generate relative abundance and melt
taxa_rel <- apply(taxa_mat, 1, function(i) i/sum(i))
rel_taxa_melt <- melt(taxa_rel)
taxa_melt <- melt(t(taxa_mat))
colnames(taxa_melt) <- c('Taxa', 'SampleID', 'Fraction')

# Order the melted data frame by day, then, group, then mouse
meta <- as.data.frame(do.call(rbind, strsplit(as.character(taxa_melt$SampleID), split = "-")))
colnames(meta) <- c('group', 'day', 'mouse')
meta$day <- factor(meta$day, levels = c('PreOp', 'Day5', 'Day10', 'Day16'))
meta$group <- factor(meta$group, levels=c('HFD', 'ND'))
meta$mouse <- factor(meta$mouse, levels=c('1', '2', '3', '4', '5'))
rel_taxa_melt <- rel_taxa_melt[with(meta, order(day, group, mouse)),]
taxa_melt <- taxa_melt[with(meta, order(day, group, mouse)),]
taxa_melt <- cbind(taxa_melt, meta[with(meta, order(day, group, mouse)),])
taxa_melt$SampleID <- factor(as.character(taxa_melt$SampleID),
                             levels=unique(as.character(taxa_melt$SampleID)))
taxa_melt$Rel <- rel_taxa_melt$value

# Only visualize/color the taxas of interest (From longitudinal volatility)
taxa_melt$col <- rep('.', nrow(taxa_melt))
for(toi in imp_feats){
  taxa_melt$col[which(taxa_melt$Taxa == toi)] <- toi
}

########################
#### Visualizations ####
## GGplot of important taxas stacked barplots
taxa_spl <- split(taxa_melt, taxa_melt$group)
ggbps <- lapply(taxa_spl, function(tgrp){
  ggplot(tgrp, aes(fill=col, y=Fraction, x=SampleID)) + 
    facet_grid(cols=vars(day), space='free', scales = "free") +
    geom_bar(position="fill", stat="identity", show.legend = TRUE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank())
})
cow_barp <- plot_grid(plotlist=ggbps, ncol=1)

ggbps_all <- lapply(taxa_spl, function(tgrp){
  taxa_pattern <- paste0("^", substr(lvl_taxa, 1, 1), "__")
  tgrp <- tgrp[grep(taxa_pattern, as.character(tgrp$Taxa)),]
  bigtaxa <- nchar(as.character(tgrp$Taxa)) > 30
  levels(tgrp$Taxa) <- substr(as.character(tgrp$Taxa), 1, 30)
  grp <- unique(as.character(tgrp$group))
  
  ggplot(tgrp, aes(fill=Taxa, y=Fraction, x=SampleID)) + 
    facet_grid(cols=vars(day), space='free', scales = "free") +
    geom_bar(position="fill", stat="identity", show.legend = TRUE) +
    theme_minimal() +
    ggtitle(grp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank()) +
    scale_fill_manual(values=rep(brewer.pal(12, 'Set3'), 5))
})
pdf(file.path("taxonomy", "data", paste0(lvl_f, ".pdf")), width = 20)
legend <- legend <- get_legend(
  # create some space to the left of the legend
  ggbps_all[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)
prow <- plot_grid(ggbps_all[[1]] + theme(legend.position='none'),
                  ggbps_all[[2]] + theme(legend.position='none'), nrow=2)
plot_grid(prow, legend, rel_widths = c(2, 1))
dev.off()

## GGplot spaghetti plots of relative fractions of important taxas
taxa_melt$group_mouse <- with(taxa_melt, paste0(group, "_", mouse))
taxa_melt$group_day <- with(taxa_melt, paste0(group, "_", day))
taxa_spl <- split(taxa_melt, taxa_melt$col)
taxa_spl <- taxa_spl[-grep("^\\.$", names(taxa_spl))]
ggspaghetti <- lapply(taxa_spl, function(tgrp){
  id <- unique(tgrp$col)
  ggspag <- ggplot(tgrp, aes(x=day, y=Rel, group=group_mouse, color=group)) +
    geom_line(size=0.2) +
    theme_minimal() + 
    ggtitle(id) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  red_tgrp <- tgrp[-which(duplicated(tgrp$group_day)),]
  red_tgrp$Rel <-  sapply(split(tgrp, list(tgrp$group, tgrp$day)), function(i){ mean(i$Rel) })
  
  ggspag + 
    geom_line(data=red_tgrp, aes(x=day, y=Rel, group=group, color=group), size=1.5)
})
cow_spag <- plot_grid(plotlist=ggspaghetti, ncol=1)

## GGplot horizontal barplot of feature importances
longitudinal$simple_id <- getTaxa(longitudinal$id, lvl = lvl_taxa, split_regex = ';')
colnames(longitudinal) <- gsub(" ", "_", colnames(longitudinal))
longitudinal <- longitudinal[order(longitudinal$importance, decreasing = FALSE),]
longitudinal$simple_id <- factor(as.character(longitudinal$simple_id),
                                 levels=unique(as.character(longitudinal$simple_id)))

ggstat <- list()
# Feature Importance
ggstat[['importance']] <- ggplot(data=longitudinal, aes(x=simple_id, y=importance)) +
  geom_bar(stat='identity') + 
  coord_flip() + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = rel(0.7)),
        axis.title.y = element_blank())

# Cumulative avg increase and decrease
cumavg <- melt(longitudinal[,c('simple_id', paste0('Cumulative_Avg_', c('Decrease', 'Increase')))])
colnames(cumavg)[3] <- c('Cumulative_Avg_Change')
ggstat[['cumavg']] <- ggplot(data=cumavg, aes(x=simple_id, y=Cumulative_Avg_Change)) + 
  geom_bar(stat='identity') + 
  coord_flip() + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = rel(0.7)),
        axis.title.y = element_blank())

ggstat[['variance']] <- ggplot(data=longitudinal, aes(x=simple_id, y=Global_Variance)) +
  geom_bar(stat='identity') + 
  coord_flip() + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = rel(0.7)),
        axis.title.y = element_blank())
cow_stats <- plot_grid(plotlist=ggstat, ncol=1)





pdf(file.path("longitudinal", "data", paste0(lvl_f, ".pdf")), height = 17, width=12)
cow_spag_stats <- plot_grid(cow_spag, cow_stats, ncol=2)
plot_grid(cow_barp, cow_spag_stats, nrow=2, rel_heights = c(1,2))
dev.off()

pdf(file.path("longitudinal", "data", paste0(lvl_f, "_importance.pdf")), height = 17, width=6)
plot_grid(cow_stats)
dev.off()
