library(reshape2)
library(ggplot2)
library(vsn)
library(Maaslin2)
library(DESeq2)
library(cowplot)
library(vegan)
library(ape)
library(ggrepel)
library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
library(pheatmap)

PDIR <- '/Users/rquevedo/Projects/mcgaha/PASS/shotgun/data/human'
outdir <- file.path(PDIR, "results")
setwd(PDIR)
dir.create(outdir, showWarnings = FALSE)

###################
#### Functions ####
# Takes a metaphlan taxonomic matrix, rarefy's it, and fills in NA
preprocessIt <- function(mat, cnt_data=TRUE){
  mat[is.na(mat)] <- 0
  rare_mat <- mat
  if(cnt_data) rare_mat <- t(rrarefy(x = t(mat),sample = min(colSums(mat))))
  mat <- cbind(data.frame("ensgene"=unlist(getTaxa(rownames(mat), 
                                                   lvl=tolower(taxon_level), 
                                                   split_regex = "\\|"))), 
               rare_mat)
  return(mat)
}

# Susbet the metaphlan matrix for a subset of samples
getSamples <- function(mat, samples, return_ids=FALSE, metad=NULL){
  sample_v <- unlist(samples)
  sample_idx <- sapply(sample_v, function(i) match(i, colnames(mat)))
  if(!is.null(metad)){
    names(sample_idx) <- setNames(metad[,2], metad[,1])[sample_v]
  }
  if(any(is.na(sample_idx))) {
    na_idx <- which(is.na(sample_idx))
    sample_idx <- sample_idx[-na_idx]
    sample_v <- sample_v[-na_idx]
  }
  mat <- mat[,sample_idx]
  colnames(mat) <- names(sample_idx)
  if(return_ids){
    return(list("mat"=mat, "ids"=sample_v))
  } else {
    return(mat)
  }
}

# Extract the top X rowSums rows from matrix i
getTopX <- function(i, topx){
  # Get the most abundant topx taxons
  cum_pct <- rowSums(i)
  topx_ids <- head(sort(cum_pct, decreasing = TRUE), topx)
  topx_idx <- rownames(i) %in% names(topx_ids)
  
  # Compress the non-most abundant taxas into a singular value
  if(any(!topx_idx)){
    colSums(i[which(!topx_idx),])
    i <- rbind(i[which(topx_idx),], (colSums(i[which(!topx_idx),])))
    rownames(i)[topx+1] <- gsub("([kpcofgs])__[a-zA-Z0-9_]+$", "\\1__Other", rownames(i)[topx])
  } else {
    i <- i[which(topx_idx),]
  }
  return(i)
}

# Melt the metaphlan data structure for a subset of samples and subset to a top representative
meltMetaphlan <- function(meta_spl, samples, topx=10, metad=NULL){
  melt_meta <- lapply(meta_spl, function(i){
    # Select only the top X species/genus/etc...
    orig_colids <- colnames(i)
    i <- getSamples(i, samples, metad=metad)
    if(any(is.na(colnames(i)))){
      na_idx <- which(is.na(colnames(i)))
      colnames(i)[na_idx] <- orig_colids[na_idx]
    }
    i <- getTopX(i, topx)
    
    # Format the melted data structure
    mi <- melt(t(i))
    colnames(mi) <- c('Sample', 'ID', 'Percentage')
    mi$Rank <- rank_legend[gsub("^.*([kpcofgs])__.*", "\\1", as.character(mi$ID))]
    mi$ID <- gsub("^.*[kpcofgs]__", "", mi$ID)
    
    return(mi)
  })
  return(melt_meta)
}

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

#############################################
#### Metaphlan - genus/species breakdown ####
rank_legend <- c('k'='Domain',
                 'p'='Phylum',
                 'c'='Class',
                 'o'='Order',
                 'f'='Family',
                 'g'='Genus',
                 's'='Species')
# metadata <- file.path("ref", "hpb_stool_metadata.csv"
metadata <- file.path("ref", "PASS_metadata.csv")
metadata <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[which(nchar(metadata$PASS_ID) > 0),]

samples <- list("HPB"=metadata$PASS_ID)
topx <- 15 # -1  # top taxas to return
taxa_lvl <- 7 # taxa level from 1-7, with 7=species and 6=genus
wins <- 0.999 # winsorization threshold for taxa abundance
min_abundance <- 1

# Read in the metaphlan raw count for healthy HMP project
hmp_cnt <- read.table(file.path("HMP", "metaphlan", "metaphlan_taxonomic_counts.tsv"), sep="\t", 
                            header=TRUE, comment.char = '',
                            stringsAsFactors = FALSE, check.names = FALSE)
rownames(hmp_cnt) <- hmp_cnt[,1]
hmp_cnt <- hmp_cnt[,-c(1:2)]
hmp_meta_levels <- sapply(strsplit(rownames(hmp_cnt), split="\\|"), length)
hmp_cnt_spl <- split(hmp_cnt, hmp_meta_levels)
names(hmp_cnt_spl) <- rank_legend

# Read in the metaphlan raw count data and preformat
metaphlan_cnt <- read.table(file.path("metaphlan", "metaphlan_taxonomic_counts.tsv"), sep="\t", 
                        header=TRUE, comment.char = '',
                        stringsAsFactors = FALSE, check.names = FALSE)
rownames(metaphlan_cnt) <- metaphlan_cnt[,1]
metaphlan_cnt <- metaphlan_cnt[,-c(1:2)]
meta_levels <- sapply(strsplit(rownames(metaphlan_cnt), split="\\|"), length)
metacnt_spl <- split(metaphlan_cnt, meta_levels)
names(metacnt_spl) <- rank_legend

# Read in the metaphlan relative-abundance data and preformat
metaphlan <- read.table(file.path("metaphlan", "metaphlan_taxonomic_profiles.tsv"), sep="\t", 
                        header=TRUE, comment.char = '',
                        stringsAsFactors = FALSE, check.names = FALSE)
rownames(metaphlan) <- metaphlan[,1]
metaphlan <- metaphlan[,-1]
colnames(metaphlan) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(metaphlan))
meta_levels <- sapply(strsplit(rownames(metaphlan), split="\\|"), length)
meta_spl <- split(metaphlan, meta_levels)
names(meta_spl) <- rank_legend

# Read in the metaphlan relative-abundance for healthy HMP project
hmp_metaphlan <- read.table(file.path("HMP", "metaphlan", "merged_metaphlan.relative.tsv"), sep="\t", 
                        header=TRUE, comment.char = '',
                        stringsAsFactors = FALSE, check.names = FALSE)
rownames(hmp_metaphlan) <- hmp_metaphlan[,1]
hmp_metaphlan <- hmp_metaphlan[,-c(1,2)]
colnames(hmp_metaphlan) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(hmp_metaphlan))
hmp_meta_levels <- sapply(strsplit(rownames(hmp_metaphlan), split="\\|"), length)
hmp_meta_spl <- split(hmp_metaphlan, hmp_meta_levels)
names(hmp_meta_spl) <- rank_legend

#########################################
#### DESeq2 between Healthy and PASS ####
library(DESeq2)
library(ggplot2)
library(vegan)

# Extract only chemo data from metadata
submeta <- metadata[,c('PASS_ID', 'Chemo_Naive', 'Sex')]
submeta$PASS_ID <- submeta$PASS_ID
submeta$grp <- 'PASS'
hmpmeta <- data.frame("PASS_ID"=colnames(hmp_meta_spl[[taxon_level]]),
                      "Chemo_Naive"="Yes", "Sex"="U", 'grp'='HMP')
subhmp_meta <- as.data.frame(rbind(submeta, hmpmeta))

# Preprocess the taxon count data
taxon_level <- 'Species'
DA <- lapply(setNames(c('Genus', 'Species'), c('Genus', 'Species')), function(taxon_level){
  ## Function to aggregate PASS and HMP:
  aggregateHmpPass <- function(hmp, pass, meta, rm.pass.colnames=T){
    # Aggregate HMP and PASS taxas
    hmp_pass_tmp <- merge(hmp[[taxon_level]], pass[[taxon_level]], by=0, all=T)
    rownames(hmp_pass_tmp) <- hmp_pass_tmp$Row.names
    hmp_pass <- preprocessIt(hmp_pass_tmp[,-1], cnt_data = F)
    rownames(hmp_pass) <- hmp_pass[,1]
    rm_idx <- which(rowSums(hmp_pass[,-1]==0) == ncol(hmp_pass)-1)
    if(length(rm_idx)>0) hmp_pass <- hmp_pass[-rm_idx,]
    if(rm.pass.colnames) colnames(hmp_pass) <- gsub("^PASS_", "", colnames(hmp_pass))
    
    # ensure same samples in metadata and count matrix
    meta_cons <- meta[which(meta$PASS_ID %in% colnames(hmp_pass)),]
    hmp_pass <- hmp_pass[,which(colnames(hmp_pass) %in% c('ensgene', meta$PASS_ID))]
    rownames(meta_cons) <- meta_cons$PASS_ID
    # rownames(meta_cons) <- rownames(hmp_pass) <- NULL
    ids <- rownames(meta_cons)
    
    # PCA directly on relative abundances
    cnts_rel <- hmp_pass[,-1]
    cnts_rel <- as.data.frame(t(cnts_rel))
    
    list("counts"=hmp_pass[,c('ensgene', ids)], "meta"=meta_cons[ids,])
  }
  
  ### Obtaining normalized and absolute counts
  ## Load in the relative counts for PCA
  hmp_pass_rel <- aggregateHmpPass(hmp=hmp_meta_spl, pass=meta_spl, meta=subhmp_meta)
  # PCA directly on relative abundances
  cnts_rel <- hmp_pass_rel$counts[,-1]
  cnts_rel <- as.data.frame(t(cnts_rel))
  
  ## Load in the absolute counts for Differential Expression  analysis
  hmp_pass_abs <- aggregateHmpPass(hmp=hmp_cnt_spl, pass=metacnt_spl, meta=subhmp_meta)
  
  ## Differential analysis
  dds0 <- DESeqDataSetFromMatrix(countData=hmp_pass_abs$counts[,-1]+1, 
                                 colData=hmp_pass_abs$meta, 
                                 design=~grp)#, tidy = TRUE)
  dds <- DESeq(dds0)
  reslfc <- lfcShrink(dds, coef=resultsNames(dds)[2], type = 'apeglm')
  # obtain normalized counts
  cnts_ntd <- normTransform(dds)
  cnts_vst <- vst(dds, nsub=ceiling(nrow(dds)*0.5), blind=T)
  cnts_rlog <- rlog(dds, blind=T)
  
  ## Compare variance stabilizing transformations:
  pdf(file.path(outdir, "metaphlan", "hmp", paste0(taxon_level, "hmp_mean-sd-norm.pdf")))
  plot_grid(plotlist = list(meanSdPlot(assay(cnts_ntd))$gg + ggtitle("ntd"),
                            meanSdPlot(assay(cnts_vst))$gg + ggtitle("vst"),
                            meanSdPlot(assay(cnts_rlog))$gg + ggtitle("rlog"),
                            meanSdPlot(as.matrix(cnts_rel))$gg + ggtitle("relative-counts")),
            ncol=2)
  dev.off()
  
  
  ### PCA on the relative abundances
  max_pc <- 9
  cnts <- cnts_rel
  tcnts <- as.data.frame(if(class(cnts) == 'DESeq2') t(assay(cnts)) else cnts)
  pca <- prcomp(tcnts, scale=F)
  var_df <- data.frame("PC"=paste0("PC", c(1:length(pca$sdev))), 
                       "sd"=((pca$sdev^2) / sum(pca$sdev^2)))
  
  pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))],
                               "sample"=rownames(tcnts),
                               "grp"=subhmp_meta_cons[rownames(tcnts),]$grp,
                               "chemonaive"=subhmp_meta_cons[rownames(tcnts),]$Chemo_Naive))
  pca_x$grp[is.na(pca_x$grp)] <- 'PASS'
  pca_x$chemonaive[is.na(pca_x$chemonaive)] <- 'Yes'
  for(id in paste0("PC", c(1:max_pc))){
    pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
  }
  pca_x$Name <- rownames(pca_x)
  
  pcs <- paste0("PC", c(1:max_pc))
  pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})
  ggps <- lapply(pcs_l, function(pc){
    ggplot(data=pca_x, aes_string(x=pc[1], y=pc[2], shape='chemonaive', 
                                  color='grp', fill='chemonaive')) +
      geom_point(alpha=0.5) + 
      scale_shape_manual(values=c(17,16)) +
      theme_classic()
  })
  var_df$PC <- factor(var_df$PC, levels=var_df$PC)
  ggscree <- ggplot(data=var_df[1:40,], aes(x=PC, y=sd)) +
    geom_point(alpha=0.9) +
    theme_classic() + 
    ylim(0,1) + ylab("%Var explained") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(file.path(outdir, "metaphlan", "hmp", paste0(taxon_level, "pca_rel-abundance.pdf")), 
      width = 14, height = 14)
  plot_grid(plotlist = ggps, ncol=3)
  dev.off()
  pdf(file.path(outdir, "metaphlan", "hmp", paste0(taxon_level, "pca_scree.pdf")))
  print(ggscree)
  dev.off()
  
  ### Hierarchical clustering to measure sample similarity
  cnts_filt <- cnts[,which(apply(cnts, 2, sd) > 0)]
  cnts_var <- cnts_filt[,names(head(sort(apply(cnts_filt, 2, sd) ^ 2, decreasing = T), 50))]
  cnts_cor <- cor(cnts_var)    ## cor() is a base R function
  cols.cor <- cor(cnts_var, use = "pairwise.complete.obs", method = "pearson")
  # Pairwise correlation between rows (genes)
  rows.cor <- cor(t(cnts_filt), use = "pairwise.complete.obs", method = "pearson")
  
  # Plot the heatmap
  heat.colors <- brewer.pal(6, "Blues")
  pdf(file.path(outdir, "metaphlan", "hmp", paste0(taxon_level, "_hclust-cor.pdf")))
  pheatmap(as.matrix(cnts_var), scale = "none", 
           clustering_distance_cols = as.dist(1 - cols.cor),
           color = heat.colors, border_color=NA, fontsize = 4, 
           fontsize_row = 4, height=20)
  pheatmap(cols.cor, color = heat.colors, border_color=NA, fontsize = 4, 
           fontsize_row = 4, height=20)
  pheatmap(rows.cor, color = heat.colors, border_color=NA, fontsize = 4, 
           fontsize_row = 4, height=20)
  dev.off()
  
  ## Looking at differential taxas according to Maaslin2
  dir.create(file.path(outdir, "metaphlan", "hmp", "maaslin2"), showWarnings = F)
  maaslin_fit <- Maaslin2(input_data=as.data.frame(t(hmp_pass_rel$counts[,-1])),
                          input_metadata=hmp_pass_rel$meta,
                          output=file.path(outdir, "metaphlan", "hmp", "maaslin2", taxon_level),
                          fixed_effects=c("grp"))
  assoc <- maaslin_fit$results
  assoc$log10q <- -1*log10(assoc$qval)
  assoc$feature <- factor(assoc$feature, levels=rev(assoc$feature))
  assoc$dir <- assoc$coef < 0
  pdf(file.path(outdir, "metaphlan", "hmp", paste0(taxon_level, "_maaslin2.pdf")))
  ggplot(data = assoc[which(assoc$qval < 0.01),], aes(x=feature, y=coef, fill=dir)) + 
    geom_bar(stat='identity') + 
    theme_classic() + 
    coord_flip() + 
    ylim(-5.5,5.5) + xlab("") +
    theme(legend.position = "none")
  dev.off()
# 
#   # Format results for ggplots
#   res <- reslrt
#   res <- reslfc
#   res <- res[order(res$padj),]
#   res$taxa <- factor(rownames(res), levels=rev(rownames(res)[order(res$pvalue)]))
#   res$sig <- res$pvalue < 0.05
#   res$log2FoldChange[res$log2FoldChange < -10] <- -10
#   res$log2FoldChange[res$log2FoldChange > 10] <- 10
# 
#   dat <- ggplot(data=head(as.data.frame(res), 20), aes(x=taxa, y=log2FoldChange, fill=sig)) +
#     geom_bar(stat="identity", width = 0.75)  +
#     geom_point(data=head(as.data.frame(res), 20), aes(x=taxa, y=log2FoldChange, size=padj)) +
#     ylim(c(-10, 10)) +
#     coord_flip() +
#     theme_minimal()
# 
#   return(list("gg"=dat, "res"=res, "cnt"=meta_cnt))
})
  
  
##################################
#### DESeq2 differential taxa ####
library(DESeq2)
library(ggplot2)
library(vegan)

# Extract only chemo data from metadata
submeta <- metadata[,c('PASS_ID', 'Chemo_Naive', 'Sex')]
submeta$PASS_ID <- paste0("PASS_", submeta$PASS_ID)

# Preprocess the taxon count data
taxon_level <- 'Genus'
differential <- 'treatment' # ' treatment', '
DA <- lapply(setNames(c('Genus', 'Species'), c('Genus', 'Species')), function(taxon_level){
  meta_cnt <- preprocessIt(metacnt_spl[[taxon_level]])
  
  # ensure same samples in metadata and count matrix
  submeta_cons <- submeta[which(submeta$PASS_ID %in% colnames(metaphlan_cnt)),]
  meta_cnt <- meta_cnt[,which(colnames(meta_cnt) %in% c('ensgene', submeta$PASS_ID))]
  rownames(submeta_cons) <- rownames(meta_cnt) <- NULL
  
  # Differential analysis
  dds <- DESeqDataSetFromMatrix(countData=meta_cnt, 
                                colData=submeta_cons, 
                                design=~Chemo_Naive, tidy = TRUE)
  # dds <- DESeq(dds, test='Wald', fitType='parametric')
  # dds <- DESeq(dds, test='Wald')
  dds <- DESeq(dds, test='LRT', reduced=~1)
  res <- results(dds, cooksCutoff = T)
  
  # Format results for ggplots
  res <- res[order(res$pvalue),]
  res$taxa <- factor(rownames(res), levels=rev(rownames(res)[order(res$pvalue)]))
  res$sig <- res$pvalue < 0.05
  res$log2FoldChange[res$log2FoldChange < -10] <- -10
  res$log2FoldChange[res$log2FoldChange > 10] <- 10
  
  dat <- ggplot(data=head(as.data.frame(res), 20), aes(x=taxa, y=log2FoldChange, fill=sig)) +
    geom_bar(stat="identity", width = 0.75)  + 
    geom_point(data=head(as.data.frame(res), 20), aes(x=taxa, y=log2FoldChange, size=padj)) +
    ylim(c(-10, 10)) +
    coord_flip() + 
    theme_minimal()
  
  return(list("gg"=dat, "res"=res, "cnt"=meta_cnt))
})

# Plot the results of DESeq2
pdf(file.path(outdir, "metaphlan", "differential_abundance.pdf"), width=7)
lapply(DA, function(i) i$gg)
dev.off()

# Plot the raw counts for the top DESeq2 differentially abundant taxa
spl_submeta <- split(submeta, submeta$Chemo_Naive)

taxa_chemo <- lapply(DA, function(da){
  meta_cnt <- da$cnt
  da_x <- as.data.frame(da$res)
  da_x <- da_x[order(da_x$pvalue),]
  sig_taxas <- as.character(na.omit(da_x[da_x$pvalue < 0.05,]$taxa))
  taxa_melts <- lapply(sig_taxas, function(taxa){
    taxa_cnts <- lapply(spl_submeta, function(i){
      cnt <- meta_cnt[match(taxa,meta_cnt$ensgene), 
                      which(colnames(meta_cnt) %in% i$PASS_ID)]
      data.frame("Taxa"=taxa, "Chemo_Naive"=unique(i$Chemo_Naive), 
                 "Sample"=colnames(cnt), "Cnt"=unlist(cnt))
    })
    do.call(rbind, taxa_cnts)
  })
  taxa_melt <- do.call(rbind, taxa_melts)
  
  chemo_map <- setNames(c('Treated', 'Untreated'), c('No', 'Yes'))
  taxa_melt$Chemo_Naive <- factor(chemo_map[taxa_melt$Chemo_Naive], c('Treated', 'Untreated'))
  taxa_melt$Taxa <- factor(as.character(taxa_melt$Taxa), levels=(sig_taxas))
  ggbox <- ggplot(data=taxa_melt, aes(x=Chemo_Naive, y=Cnt, fill=Chemo_Naive)) +
    geom_boxplot(show.legend = T) + 
    facet_grid(rows=vars(Taxa), scales='free', space='free', switch='y') +
    coord_flip() + 
    ylim(c(0,max(taxa_melt$Cnt))) +
    theme_minimal() +
    theme(axis.title.y =element_blank(),
          axis.text.y =element_blank(),
          strip.text.y.left = element_text(angle = 0)) + 
    scale_fill_brewer(palette="Dark2")
  list("gg"=ggbox, "df"=taxa_melts)
})

pdf(file.path(outdir, "metaphlan", "da_plots.pdf"), width=7)
taxa_chemo$Genus$gg
taxa_chemo$Species$gg
dev.off()

####################################
#### PCA of the Taxas abundance ####
gg_pca <- lapply(rank_legend[-1], function(rank){
  meta_l <- getSamples(meta_spl[[rank]], samples, return_ids = TRUE)
  meta_species <- meta_l$mat
  pca <- prcomp(t(as.matrix(meta_species)))
  percent_var <- round(pca$sdev^2 / sum( pca$sdev^2 ),4)  # Calculate variance for each PC
  pc_cutoff <- min(which(cumsum(percent_var) > 0.95))     # Identify 95% variance cutoff
  
  # Generate the PCA plots using DESeq2 plotting function
  dat <- cbind(pca$x[,c('PC1', 'PC2')], 
               data.frame("group"=gsub("[0-9]+$", "", rownames(pca$x)),
                          "name"=meta_l$ids[rownames(pca$x)]))
  gg_pca <- ggplot(dat, aes(x= PC1, y= PC2, color=group, label=name))+
    geom_point() +
    geom_text(aes(label=name),hjust=0, vjust=0) +
    xlim(-100, 100) +
    ylim(-100, 100) +
    ylab(paste0("PC2 (", percent_var[2]*100, "%)")) +
    xlab(paste0("PC1 (", percent_var[1]*100, "%)")) +
    ggtitle(rank) +
    theme_minimal()
  gg_pca
})

pdf(file.path(outdir, "genus_species_pca.pdf"), width = 15, height=15)
plot_grid(plotlist=gg_pca)
dev.off()

###########################
#### Diversity Metrics ####
# Calculate the Bray-Curtis distance matrix
metad <- metadata[,c('PASS_ID', 'UID')]
metad$PASS_ID <- gsub("_.*", "", metad$PASS_ID)
meta_mat_lbld <- t(metaphlan[,metad$PASS_ID])
rownames(meta_mat_lbld) <- setNames(metad[,2], metad[,1])[rownames(meta_mat_lbld)]
meta_mat <- meta_mat_lbld[,which(sapply(strsplit(colnames(meta_mat_lbld), split="\\|"), length) == taxa_lvl)]
bray_mat <- meta_mat[,which(apply(meta_mat, 2, max) >= min_abundance)]
# shannon <- diversity(meta_mat, index = 'shannon')
bray <- (vegdist(bray_mat, method="bray"))

# Calculate the PCoA of the bray-distance
bray_pcoa <- pcoa(bray)
biplot(bray_pcoa, plot.axes = c(2,3))
var_explained <- bray_pcoa$values$Eigenvalues + abs(min( bray_pcoa$values$Eigenvalues))
var_explained <- var_explained / sum(var_explained)

# Plot the PcoA of bray-distance
bray_pcoa_df <- as.data.frame(bray_pcoa$vectors)
bray_pcoa_df$sample <- rownames(bray_pcoa_df)
bray_pcoa_df <- merge(x=bray_pcoa_df, y=metadata, by.x='sample', by.y='UID')

pdf(file.path(outdir, paste0("pcoa_bray.pdf")))
ggplot(data=bray_pcoa_df, aes(x=Axis.1, y=Axis.2, label=sample)) + 
  geom_point(aes(shape=Chemo_Naive, color=Disease_Status)) + 
  geom_text_repel(aes(label = Study_ID), max.overlaps = 50, size = 2.5) + 
  geom_line(aes(group=Study_ID)) +
  xlab(paste0("PCoA1 (", round(var_explained[1],3)*100, "%)")) +
  ylab(paste0("PCoA2 (", round(var_explained[2],3)*100, "%)"))
dev.off()

#############################
#### Clustering of taxas ####
# Establish color functions for metadata
my_palette <- colorRampPalette(c("black", "blue", "cyan", "red", "yellow"))(n = 299)
status_col_fun <- setNames(brewer.pal(length(unique(metadata$Disease_Status)), 'Set1'),
                           unique(metadata$Disease_Status))
chemo_col_fun <- setNames(brewer.pal(length(unique(metadata$Chemo_Naive)), 'PuOr'),
                           unique(metadata$Chemo_Naive))
sex_col_fun <- setNames(brewer.pal(length(unique(metadata$Sex)), 'RdGy'),
                          unique(metadata$Sex))[1:2]
age_col_fun = colorRamp2(c(40, 50, 60, 70, 80, 90), brewer.pal(6, 'Blues')) 

# Identify the top taxa to create clustering plots from
taxa_id <- switch(as.character(taxa_lvl),
                  '1'='kingdom',
                  '2'='phylum',
                  '3'='class',
                  '4'='order',
                  '5'='family',
                  '6'='genus',
                  '7'='species')
taxa_id <- getTaxa(taxas=colnames(meta_mat), lvl = taxa_id, split_regex = "\\|")
taxa_mat <- meta_mat
storage.mode(taxa_mat) <- 'numeric'
colnames(taxa_mat) <- taxa_id
taxa_mat[taxa_mat > quantile(as.numeric(taxa_mat), wins)] <- quantile(as.numeric(taxa_mat), wins)
keep_taxa <- which(apply(taxa_mat, 2, max) >= min_abundance)
taxa_mat <- if(length(keep_taxa) > 0) taxa_mat[,keep_taxa] else taxa_mat
top_taxa <- order(apply(taxa_mat, 2, mean), decreasing = TRUE)

# Generate the heatmap for absolute abundance based on top X most abundant taxa
subset_taxa <- if(topx > -1) taxa_mat[,top_taxa][,c(1:top_taxa)] else taxa_mat[,top_taxa]
ht_abs <- Heatmap(t(subset_taxa), col=my_palette,
              top_annotation=HeatmapAnnotation(status=anno_simple(metadata$Disease_Status, col=status_col_fun, height=unit(0.3, 'cm')),
                                               chemo=anno_simple(metadata$Chemo_Naive, col=chemo_col_fun, height=unit(0.3, 'cm')),
                                               sex=anno_simple(metadata$Sex, col=sex_col_fun, height=unit(0.3, 'cm')),
                                               age=anno_simple(as.numeric(metadata$Age), col=age_col_fun,  height=unit(0.3, 'cm'))),
              column_km=5, row_km=5)

# Generate the heatmap for sample-normalized abundance based on top X most abundant taxa
rel_taxa_mat <- t(apply(taxa_mat,1,function(i) i/sum(i)))
subset_taxa <- if(topx > -1) rel_taxa_mat[,top_taxa][,c(1:top_taxa)] else rel_taxa_mat[,top_taxa]
ht_rel <- Heatmap(t(subset_taxa), col=my_palette,
        top_annotation=HeatmapAnnotation(status=anno_simple(metadata$Disease_Status, col=status_col_fun, height=unit(0.3, 'cm')),
                                         chemo=anno_simple(metadata$Chemo_Naive, col=chemo_col_fun, height=unit(0.3, 'cm')),
                                         sex=anno_simple(metadata$Sex, col=sex_col_fun, height=unit(0.3, 'cm')),
                                         age=anno_simple(as.numeric(metadata$Age), col=age_col_fun,  height=unit(0.3, 'cm'))),
        column_km=5, row_km=5)

# Create the legend diagrams for the different metadata values
stt_leg = Legend(title = "Disease status", legend_gp = gpar(fill=status_col_fun),
                    at = seq_along(status_col_fun), labels = names(status_col_fun))
che_leg = Legend(title = "Chemo naive", legend_gp = gpar(fill=chemo_col_fun),
                    at = seq_along(chemo_col_fun), labels = names(chemo_col_fun))
sex_leg = Legend(title = "Sex", legend_gp = gpar(fill=sex_col_fun),
                    at = seq_along(sex_col_fun), labels = names(sex_col_fun))
age_leg = Legend(title = "Age", col_fun = age_col_fun,  at=seq(40,90, by=10),
                 labels = seq(40,90, by=10))

pdf(file.path(outdir, paste0("taxa", taxa_lvl, "_clusters.pdf")), 
    height = 20, width = 10)
draw(ht_abs, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
draw(ht_rel, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
dev.off()
          
##########################################################
#### Stacked barplot of organisms within each phyllum ####
## Breakdown the composition of organisms across all the mice

melt_meta <- meltMetaphlan(meta_spl, samples, topx = topx, metad=metadata[,c('PASS_ID', 'Specimen_Label')])
gps <- lapply(names(melt_meta), function(r){
  mm <- melt_meta[[r]]
  x <- ggplot(mm, aes(fill=ID, y=Percentage, x=Sample)) + 
    #facet_grid(vars(Rank), scales = "free") +
    geom_bar(position="stack", stat="identity") +
    ggtitle(r) +
    xlab("") +
    theme_minimal() +
    scale_fill_brewer(palette="Set3") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  x
})

ggpie <- lapply(names(melt_meta), function(r){
  mm <- melt_meta[[r]]
  avg_melt <- do.call(rbind, lapply(split(mm, mm$ID), function(i){
    j <- i[1,,drop=F]
    j$Percentage <- mean(i$Percentage)
    j
  }))
  other_idx <- match("Other", rownames(avg_melt))
  avg_melt['Other',]$Percentage <- sum(avg_melt$Percentage) - sum(avg_melt[-other_idx,]$Percentage)
  
  gp <- ggplot(avg_melt, aes(x="", y=Percentage, fill=ID)) +
    geom_bar(width = 1, stat = "identity") +
    theme_minimal() +
    ggtitle(r) +
    scale_fill_brewer(palette="Set3")
  
  gpie <- gp + coord_polar("y", start=0)
  list("box"=gp, "pie"=gpie)
})

pdf(file.path(outdir, "taxa_barplots.pdf"))
lapply(gps, function(i) i)
dev.off()

pdf(file.path(outdir, "avg_taxa_barplots.pdf"))
lapply(ggpie, function(i) i$box)
lapply(ggpie, function(i) i$pie)
dev.off()

###############################################################
#### HMP: Stacked barplot of organisms within each phyllum ####
# Read in the metaphlan relative-abundance data and preformat
hmp <- read.table(file.path("HMP", "hmp1-II_metaphlan2-mtd-qcd.pcl"), sep="\t", 
                  header=TRUE, comment.char = '',
                  stringsAsFactors = FALSE, check.names = FALSE)
rownames(hmp) <- hmp[,1]
hmp <- hmp[,-1]

hmp_meta <- as.data.frame(t(hmp[1:8,]))
hmp <- hmp[-c(1:8),]
stool_idx <- which(hmp_meta$STSite %in% 'Stool')
hmp_stool <- as.matrix(hmp[,stool_idx])
storage.mode(hmp_stool) <- 'numeric'
hmp_stool <- as.data.frame(hmp_stool)

meta_levels <- sapply(strsplit(rownames(hmp_stool), split="\\|"), length)
hmp_spl <- split(hmp_stool, meta_levels)
names(hmp_spl) <- rank_legend

melt_hmp <- meltMetaphlan(hmp_spl[1:7], list("HMP"=colnames(hmp_stool)), topx = topx, 
                           metad=hmp_meta[,c('SNPRNT', 'SNPRNT')])
gps <- lapply(names(melt_hmp), function(r){
  mm <- melt_hmp[[r]]
  mm$Sample <- as.character(mm$Sample)
  x <- ggplot(mm, aes(fill=ID, y=Percentage, x=Sample)) + 
    #facet_grid(vars(Rank), scales = "free") +
    geom_bar(position="stack", stat="identity") +
    ggtitle(r) +
    xlab("") +
    theme_minimal() +
    scale_fill_brewer(palette="Set3") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  x
})

ggpie <- lapply(names(melt_hmp), function(r){
  mm <- melt_hmp[[r]]
  avg_melt <- do.call(rbind, lapply(split(mm, mm$ID), function(i){
    j <- i[1,,drop=F]
    j$Percentage <- mean(i$Percentage)
    j
  }))
  other_idx <- match("Other", rownames(avg_melt))
  avg_melt['Other',]$Percentage <- sum(avg_melt$Percentage) - sum(avg_melt[-other_idx,]$Percentage)
  
  gp <- ggplot(avg_melt, aes(x="", y=Percentage, fill=ID)) +
    geom_bar(width = 1, stat = "identity") +
    theme_minimal() +
    ggtitle(r) +
    scale_fill_brewer(palette="Set3")
  
  gpie <- gp + coord_polar("y", start=0)
  list("box"=gp, "pie"=gpie)
})

pdf(file.path(outdir, "HMP.taxa_barplots.pdf"), width = 25)
lapply(gps, function(i) i)
dev.off()

pdf(file.path(outdir, "HMP.avg_taxa_barplots.pdf"))
lapply(ggpie, function(i) i$box)
lapply(ggpie, function(i) i$pie)
dev.off()
