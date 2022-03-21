library(reshape2)
library(ggplot2)
library(DESeq2)
library(cowplot)
library(vegan)
library(ape)
library(ggrepel)
library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
library(brendaDb)
library(circlize)

brenda.filepath <- DownloadBrenda()
brenda <- ReadBrenda(brenda.filepath)
brenda_spl <- split(brenda, brenda$ID)
brenda_rn <- sapply(brenda_spl, function(i) {
  x <- gsub("^RN\t", "", i[which(i$field=='RECOMMENDED_NAME'),]$description) %>%
    gsub("\n", "", .)
  if(length(x) > 0 ) x else NA
})

PDIR <- '/Users/rquevedo/Projects/mcgaha/PASS/shotgun/data/human'
outdir <- file.path(PDIR, "results", "humann")
setwd(PDIR)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###################
#### Functions ####
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
  cum_pct <- rowSums(i)
  topx_ids <- head(sort(cum_pct, decreasing = TRUE), topx)
  topx_idx <- rownames(i) %in% names(topx_ids)
  if(any(!topx_idx)){
    colSums(i[which(!topx_idx),])
    i <- rbind(i[which(topx_idx),], (colSums(i[which(!topx_idx),])))
    rownames(i)[topx+1] <- gsub("_[a-zA-Z0-9]+$", "_Other", rownames(i)[topx])
  }
  return(i)
}

getFunction <- function(funcs, meta=NULL){
  func_df   <- as.data.frame(do.call(rbind, strsplit(funcs, split="\\|")))
  func_meta <- merge(func_df, meta, by.x=c('V1', 'V2'), by.y=c('ec', 'taxa'), all.x=TRUE)
  func_idx  <- which(func_meta$V1 == func_meta$V2)
  
  id <- apply(func_meta, 1, function(i) {
    species <- i['species']
    simplifySpecies(species)
  })
  
  id <- paste0(id, ": ", func_meta$rn)
  id[func_idx] <- func_meta$rn[func_idx]
  return(list("ids"=id, "func_idx"=func_idx))
}

# Converts species from "s__Akkermansia_muciniphila"  to "s_Akker._mucin."
simplifySpecies <- function(species){
  if(!is.na(species)){
    reduced_species <- sapply(strsplit(gsub("_+", "_", species), split="_")[[1]], function(j) {
      max_j <- if(nchar(j)>5) 5 else nchar(j)
      paste0(substr(j, 1, max_j), if(nchar(j) > 5) ".")
    })
    paste(reduced_species, collapse="_")
  } else {
    NA
  }
}

getAbundanceId <- function(id){
  id=rownames(pathabund)
  id_df <- as.data.frame(do.call(rbind, strsplit(id, split=":")))
  id_df2 <- as.data.frame(do.call(rbind, strsplit(id_df$V2, split="\\|")))
  id_df3 <- cbind(id_df[,1], id_df2)
  colnames(id_df3) <- c('biocycID', 'pathway', 'taxa')
  id_df3$taxa[with(id_df3, which(pathway == taxa))] <- NA
  
  id_df3$species <- sapply(strsplit(id_df3$taxa, split="\\."), function(i) {
    simplifySpecies(i[2])
  })
  id_df3$uid <- with(id_df3, paste0(biocycID, "|", species))
  return(id_df3)
}


###############################
#### 1.a) HUMANN - EC Breakdown ####
metadata <- "hpb_stool_metadata.csv"
metadata <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[which(nchar(metadata$PASS_ID) > 0),]

samples <- list("HPB"=metadata$PASS_ID)
topx <- -1  # top Enzymatic functions to use for hclusting
top_barx <- 20 # top enzymatic functions to highlight in stkaced barplot
min_abundance <- 0.01

# Read in the humann enzymatic (EC) data and preformat
humann <- read.table(file.path("humann", "ecs_relab.tsv"), sep="\t", 
                     header=TRUE, comment.char = '',
                     stringsAsFactors = FALSE, check.names = FALSE)
rownames(humann) <- humann[,1]
humann <- humann[,-1]
colnames(humann) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(humann))

# format the functional meta levels
meta_levels <- do.call(rbind, strsplit(rownames(humann), split="\\|"))
meta_levels <- as.data.frame(meta_levels)
meta_levels$rn <- gsub('\\s', "_", brenda_rn[meta_levels$ec])
meta_levels$species <- gsub("^.*\\.", "", meta_levels$taxa)
meta_levels$genus <- gsub("\\..*", "", meta_levels$taxa)

########################################
#### 1.b) PCA of the Enzymatic Abundance ####
meta_l <- getSamples(humann, samples, return_ids = TRUE)
meta_func <- meta_l$mat
pca <- prcomp(t(as.matrix(meta_func)))
percent_var <- round(pca$sdev^2 / sum( pca$sdev^2 ),4)  # Calculate variance for each PC
pc_cutoff <- min(which(cumsum(percent_var) > 0.95))     # Identify 95% variance cutoff

# Generate the PCA plots using DESeq2 plotting function
dat <- cbind(pca$x[,c('PC1', 'PC2')], 
             data.frame("group"=gsub("[0-9]+$", "", rownames(pca$x)),
                        "name"=meta_l$ids[rownames(pca$x)]))
gg_pca <- ggplot(dat, aes(x= PC1, y= PC2, color=group, label=name))+
  geom_point() +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlim(-0.05, 0.05) +
  ylim(-0.05, 0.05) +
  ylab(paste0("PC2 (", percent_var[2]*100, "%)")) +
  xlab(paste0("PC1 (", percent_var[1]*100, "%)")) +
  theme_minimal()
gg_pca


pdf(file.path(outdir, "ecs_pca.pdf"))
plot_grid(gg_pca)
dev.off()

###########################
#### 1.c) Diversity Metrics ####
# Calculate the Bray-Curtis distance matrix
metad <- metadata[,c('PASS_ID', 'UID')]
metad$PASS_ID <- gsub("_.*", "", metad$PASS_ID)
meta_mat_lbld <- t(humann[,metad$PASS_ID])
rownames(meta_mat_lbld) <- setNames(metad[,2], metad[,1])[rownames(meta_mat_lbld)]
bray_mat <- meta_mat_lbld
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


###########################################
#### 1.d) Clustering of Enzymatic Functions ####
# Establish color functions for metadata
my_palette <- colorRampPalette(c("black", "blue", "cyan", "red", "yellow"))(n = 299)
status_col_fun <- setNames(brewer.pal(length(unique(metadata$Disease_Status)), 'Set1'),
                           unique(metadata$Disease_Status))
chemo_col_fun <- setNames(brewer.pal(length(unique(metadata$Chemo_Naive)), 'PuOr'),
                          unique(metadata$Chemo_Naive))
sex_col_fun <- setNames(brewer.pal(length(unique(metadata$Sex)), 'RdGy'),
                        unique(metadata$Sex))[1:2]
age_col_fun = colorRamp2(c(40, 50, 60, 70, 80, 90), brewer.pal(6, 'Blues')) 

# Identify the top functions to create clustering plots from
func_id <- getFunction(funcs=colnames(meta_mat_lbld), meta = meta_levels)
func_mat <- meta_mat_lbld
storage.mode(func_mat) <- 'numeric'
colnames(func_mat) <- func_id$ids
func_mat_core <- func_mat[,func_id$func_idx]
keep_idx <- which(colSums(func_mat_core)>=min_abundance)
func_mat_core <- func_mat_core[,keep_idx]
top_func <- order(apply(func_mat_core, 2, mean), decreasing = TRUE)

# Consensus clustering
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(func_mat_core,maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title='samples', clusterAlg="hc",distance="pearson",seed=1234,
                               plot="png")
tresults = ConsensusClusterPlus(t(func_mat_core),maxK=20,reps=50,pItem=0.8,pFeature=1,
                               title='functions', clusterAlg="hc",distance="pearson",seed=1234,
                               plot="png")

# Generate the heatmap for absolute abundance based on top X most abundant taxa
subset_func <- if(topx > -1) func_mat_core[,top_func][,c(1:topx)] else func_mat_core[,top_func]
ht_abs <- Heatmap(t(subset_func), col=my_palette,
                  top_annotation=HeatmapAnnotation(status=anno_simple(metadata$Disease_Status, col=status_col_fun, height=unit(0.3, 'cm')),
                                                   chemo=anno_simple(metadata$Chemo_Naive, col=chemo_col_fun, height=unit(0.3, 'cm')),
                                                   sex=anno_simple(metadata$Sex, col=sex_col_fun, height=unit(0.3, 'cm')),
                                                   age=anno_simple(as.numeric(metadata$Age), col=age_col_fun,  height=unit(0.3, 'cm'))),
                  column_km=5, row_km=10)

# Create the legend diagrams for the different metadata values
stt_leg = Legend(title = "Disease status", legend_gp = gpar(fill=status_col_fun),
                 at = seq_along(status_col_fun), labels = names(status_col_fun))
che_leg = Legend(title = "Chemo naive", legend_gp = gpar(fill=chemo_col_fun),
                 at = seq_along(chemo_col_fun), labels = names(chemo_col_fun))
sex_leg = Legend(title = "Sex", legend_gp = gpar(fill=sex_col_fun),
                 at = seq_along(sex_col_fun), labels = names(sex_col_fun))
age_leg = Legend(title = "Age", col_fun = age_col_fun,  at=seq(40,90, by=10),
                 labels = seq(40,90, by=10))

pdf(file.path(outdir, "func_clusters.pdf"),  height = 100, width = 10)
draw(ht_abs, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
# draw(ht_rel, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
dev.off()

##############################################################
#### 1.e) Stacked barplot of organisms within each enzyme ####
## Breakdown the composition of organisms across all the mice
melt_func <- melt(func_mat_core)
colnames(melt_func) <- c('Sample', 'Func', 'Fraction')
melt_func$col <- rep(".", nrow(melt_func))
top_ids <- names(head(sort(colSums(func_mat_core), decreasing = TRUE), top_barx))
top_idx <- which(melt_func$Func %in% top_ids)
melt_func[top_idx,'col'] <- as.character(melt_func[top_idx,'Func'])
melt_func$col <- factor(melt_func$col, levels=rev(c(top_ids, ".")))

ggp <- ggplot(melt_func, aes(fill=col, y=Fraction, x=Sample)) + 
  #facet_grid(vars(Rank), scales = "free") +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=rep(brewer.pal(12, "Set2"), 5), drop=FALSE)

pdf(file.path(outdir, "func_barplots.pdf"), width = 15)
ggp
dev.off()


##########################################
#### 2.a) HUMANN - Pathway Abundance  ####
# Read in the HUMANN pathway abundance data and preformat
pathabund <- read.table(file.path("humann", "pathabundance_relab2.tsv"), sep="\t", 
                        header=TRUE, comment.char = '',
                        stringsAsFactors = FALSE, check.names = FALSE)
rownames(pathabund) <- pathabund[,1]
pathabund <- pathabund[,-1]
colnames(pathabund) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(pathabund))

pathabund_lbld <- t(pathabund[,metad$PASS_ID])
rownames(pathabund_lbld) <- setNames(metad[,2], metad[,1])[rownames(pathabund_lbld)]

##############################################
#### 2.b) Clustering of pathway abundance ####
# Establish color functions for metadata
my_palette <- colorRampPalette(c("black", "blue", "cyan", "red", "yellow"))(n = 299)
status_col_fun <- setNames(brewer.pal(length(unique(metadata$Disease_Status)), 'Set1'),
                           unique(metadata$Disease_Status))
chemo_col_fun <- setNames(brewer.pal(length(unique(metadata$Chemo_Naive)), 'PuOr'),
                          unique(metadata$Chemo_Naive))
sex_col_fun <- setNames(brewer.pal(length(unique(metadata$Sex)), 'RdGy'),
                        unique(metadata$Sex))[1:2]
age_col_fun = colorRamp2(c(40, 50, 60, 70, 80, 90), brewer.pal(6, 'Blues')) 

# Identify the top functions to create clustering plots from
pathway_id <- getAbundanceId(id=rownames(pathabund_lbld))
pathabund_mat <- as.matrix(pathabund_lbld)
storage.mode(pathabund_mat) <- 'numeric'
pathabund_core <- as.data.frame(pathabund_mat[,which(is.na(pathway_id$taxa))])
keep_idx <- which(colSums(pathabund_core)>=min_abundance)
pathabund_core <- pathabund_core[,keep_idx]
colnames(pathabund_core) <- pathway_id$pathway[which(is.na(pathway_id$taxa))][keep_idx]
top_path <- order(apply(pathabund_core, 2, mean), decreasing = TRUE)

# Consensus clustering
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(pathabund_core),maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title='samples', clusterAlg="hc",distance="pearson",seed=1234,
                               plot="png")
tresults = ConsensusClusterPlus(t(pathabund_core),maxK=10,reps=50,pItem=0.8,pFeature=1,
                                title='functions', clusterAlg="hc",distance="pearson",seed=1234,
                                plot="png")

# Generate the heatmap for absolute abundance based on top X most abundant taxa
subset_func <- if(topx > -1) pathabund_core[,top_path][,c(1:topx)] else pathabund_core[,top_path]
ht_abs <- Heatmap(t(subset_func), col=my_palette,
                  top_annotation=HeatmapAnnotation(status=anno_simple(metadata$Disease_Status, col=status_col_fun, height=unit(0.3, 'cm')),
                                                   chemo=anno_simple(metadata$Chemo_Naive, col=chemo_col_fun, height=unit(0.3, 'cm')),
                                                   sex=anno_simple(metadata$Sex, col=sex_col_fun, height=unit(0.3, 'cm')),
                                                   age=anno_simple(as.numeric(metadata$Age), col=age_col_fun,  height=unit(0.3, 'cm'))),
                  column_km=4, row_km=6)

# Create the legend diagrams for the different metadata values
stt_leg = Legend(title = "Disease status", legend_gp = gpar(fill=status_col_fun),
                 at = seq_along(status_col_fun), labels = names(status_col_fun))
che_leg = Legend(title = "Chemo naive", legend_gp = gpar(fill=chemo_col_fun),
                 at = seq_along(chemo_col_fun), labels = names(chemo_col_fun))
sex_leg = Legend(title = "Sex", legend_gp = gpar(fill=sex_col_fun),
                 at = seq_along(sex_col_fun), labels = names(sex_col_fun))
age_leg = Legend(title = "Age", col_fun = age_col_fun,  at=seq(40,90, by=10),
                 labels = seq(40,90, by=10))

pdf(file.path(outdir, "pathway_clusters.pdf"),  height = 50, width = 10)
draw(ht_abs, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
# draw(ht_rel, annotation_legend_list = list(stt_leg, che_leg, sex_leg, age_leg))
dev.off()

##############################################
#### 2.b) Clustering of pathway abundance ####
melt_func <- melt(t(pathabund_core))
colnames(melt_func) <- c('Pathway', 'Sample', 'Fraction')
melt_func$col <- rep(".", nrow(melt_func))
top_ids <- names(head(sort(colSums(pathabund_core), decreasing = TRUE), top_barx))
top_idx <- which(melt_func$Pathway %in% top_ids)
melt_func[top_idx,'col'] <- as.character(melt_func[top_idx,'Pathway'])
melt_func$col <- factor(melt_func$col, levels=rev(c(top_ids, ".")))

ggp <- ggplot(melt_func, aes(fill=col, y=Fraction, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=rep(brewer.pal(8, "Set2"), 5), drop=FALSE)

pdf(file.path(outdir, "pathway_barplots.pdf"), width = 15)
ggp
dev.off()
