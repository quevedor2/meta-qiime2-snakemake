library(reshape2)
library(ggplot2)
library(DESeq2)


PDIR <- '/Users/rquevedo/Projects/mcgaha/PASS/shotgun/data/mouse'
outdir <- '~/Projects/mcgaha/PASS/shotgun/results/functional/'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(PDIR)

###################
#### Functions ####
#### Functions ####
# Susbet the metaphlan matrix for a subset of samples
getSamples <- function(mat, samples, return_ids=FALSE){
  sample_v <- unlist(samples)
  sample_idx <- sapply(sample_v, function(i) match(i, colnames(mat)))
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
getTopX <- function(i, topx, pat_took_a_turn="_[a-zA-Z0-9]+$"){
  cum_pct <- rowSums(as.matrix(i))
  topx_ids <- head(sort(cum_pct, decreasing = TRUE), topx)
  topx_idx <- rownames(i) %in% names(topx_ids)
  if(any(!topx_idx)){
    # colSums(i[which(!topx_idx),])
    i <- rbind(i[which(topx_idx),], (colSums(i[which(!topx_idx),])))
    rownames(i)[topx+1] <- gsub(pat_took_a_turn, "_Other", rownames(i)[topx+1])
  }
  return(i)
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
controls <- list("L.reuteri"="S33",
                 "L.johnsonii"="S36",
                 "L.intestinalis"="S35",
                 "L.murinus"="S34")
mice <- list("BeforeAbx"=paste0("S", c(1:6)),
             "AfterAbx"=paste0("S", c(7:12)),
             "GvgFecal"=paste0("S", c(28:32)),
             "GvgLiLj"=paste0("S", c(18:22)),
             "GvgLmLr"=paste0("S", c(13:17)),
             "GvgPBS"=paste0("S", c(23:27)))
topx <- 7


# Read in the metaphlan data and preformat
# pathabundance_relab2.tsv is the result of a sed command to remove problematic ' characters
pathabundance <- read.table("pathabundance_relab2.tsv", sep="\t", 
                            header=TRUE, comment.char = '',
                            stringsAsFactors = FALSE, check.names = FALSE)
rownames(pathabundance) <- pathabundance[,1]
pathabundance <- pathabundance[,-1]
colnames(pathabundance) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(pathabundance))

##########################################################################
#### Heatmap (and clustering) of top X pathways in bacterial controls ####
# Annotated pathway dataframe
path_l <- getSamples(pathabundance, controls, return_ids = TRUE)

known_pathways <- sapply(gsub("L.", "", names(controls)), function(bacterial_id){
  grep(paste0("Lactobacillus_", bacterial_id), rownames(path_l$mat))
})
path_top <- path_l$mat[unlist(known_pathways),]

# Make a heatmap of the top X functional pathways, separated by groups
path_melt <- melt(t(path_top))
path_melt$functional_id <- gsub("^.*s__", "", path_melt$Var2)
colnames(path_melt) <- c('Sample', 'Pathway', 'Pct', 'functional_id')

gg_hm <- ggplot(path_melt, aes(x=Sample, y=Pathway, fill= Pct)) + 
  facet_grid(rows=vars(functional_id), scales = "free") +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(file.path(outdir, "bacterial_heatmap.pdf"), width = 11, height=14)
gg_hm
dev.off()

###################################################
#### PCA of functional pathways in mice groups ####
path_meta <- do.call(rbind, strsplit(rownames(pathabundance), split="[:|]"))
path_l <- getSamples(pathabundance, mice, return_ids = TRUE)
pat_abund <- path_l$mat
pca <- prcomp(t(as.matrix(pat_abund)))
percent_var <- round(pca$sdev^2 / sum( pca$sdev^2 ),4)  # Calculate variance for each PC
pc_cutoff <- min(which(cumsum(percent_var) > 0.95))     # Identify 95% variance cutoff
screedat <- data.frame("PC"=paste0("PC", c(1:length(percent_var))), 
                       "Variance"=percent_var)
screedat$PC <- factor(as.character(screedat$PC), paste0("PC", c(1:length(percent_var))))

# Generate the PCA plots using DESeq2 plotting function
dat <- cbind(pca$x[,c('PC1', 'PC2')], 
             data.frame("group"=gsub("[0-9]+$", "", rownames(pca$x)),
                        "name"=path_l$ids[rownames(pca$x)]))
gg_pca <- ggplot(dat, aes(x= PC1, y= PC2, color=group, label=name))+
  geom_point() +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlim(-0.15, 0.15) +
  ylim(-0.15, 0.15) +
  ggtitle("PCA-functional pathways") +
  theme_minimal()
gg_scree <- ggplot(screedat, aes(x=PC, y=Variance)) +
  geom_bar(stat='identity') +
  ylim(0, 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(file.path(outdir, "mice_pca.pdf"))
gg_pca
gg_scree
dev.off()

###################################################################
#### Heatmap (and clustering) of top X pathways in mice groups ####
# Annotated pathway dataframe
path_meta <- do.call(rbind, strsplit(rownames(pathabundance), split="[:|]"))
path_l <- getSamples(pathabundance, mice, return_ids = TRUE)
path_top <- getTopX(path_l$mat, topx=100, pat_took_a_turn = ".*")

# Hierarchical clustering of the top functional pathways
hcl_path <- hclust(dist(path_top))

# Make a heatmap of the top X functional pathways, separated by groups
grps <- c("BeforeAbx","AfterAbx","GvgPBS","GvgFecal","GvgLiLj","GvgLmLr") 
path_melt <- melt(t(path_top))
path_melt$Grp <- gsub("[0-9]*$", "", path_melt$Var1)
path_melt$Var1 <- path_l$ids[path_melt$Var1]
colnames(path_melt) <- c('Sample', 'Pathway', 'Pct', 'Grp')
path_melt$Grp <- factor(as.character(path_melt$Grp), grps)

if(use_hcl){
  path_melt$Pathway <- factor(as.character(path_melt$Pathway), rownames(path_top)[hcl_path$order])
}else{
  path_melt$Pathway <- factor(as.character(path_melt$Pathway), rownames(path_top))
}
if(rm_others) path_melt <- path_melt[-grep("_Other", path_melt$Pathway),]

gg_hm <- ggplot(path_melt, aes(Sample, Pathway, fill= Pct)) + 
  facet_grid(cols=vars(Grp), scales = "free") +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu")
pdf(file.path(outdir, "mice_heatmap.pdf"), width = 12, height=14)
gg_hm
dev.off()

##########################################################
#### Find Differentially expressed pathawys per group ####
# Annotated pathway dataframe
path_meta <- do.call(rbind, strsplit(rownames(pathabundance), split="[:|]"))
path_l <- getSamples(pathabundance, mice, return_ids = TRUE)
path_mat <- path_l$mat
qvals <- lapply(grps, function(grp){
  grp1 <- path_mat[,grep(grp, colnames(path_mat))]
  grp2 <- path_mat[,grep(grp, colnames(path_mat), invert = TRUE)]
  pvals <- sapply(setNames(c(1:nrow(grp1)), rownames(grp1)), function(row_idx){
    wilcox.test(unlist(grp1[row_idx,]), unlist(grp2[row_idx,]), alternative = 'greater')$p.val
  })
  qvals <- p.adjust(pvals, method='fdr')
  qvals[which(qvals < 0.05)]
})
names(qvals) <- grps
path_top <- path_mat[unlist(sapply(qvals, names)),]
path_grp <- setNames(gsub("[0-9]*$", "", names(unlist(sapply(qvals, names)))),
                     unlist(sapply(qvals, names)))

# Make a heatmap of the top X functional pathways, separated by groups
grps <- c("BeforeAbx","AfterAbx","GvgPBS","GvgFecal","GvgLiLj","GvgLmLr") 
path_melt <- melt(t(path_top))
path_melt$Grp <- gsub("[0-9]*$", "", path_melt$Var1)
path_melt$Var1 <- path_l$ids[path_melt$Var1]
colnames(path_melt) <- c('Sample', 'Pathway', 'Pct', 'Grp')
path_melt$Grp <- factor(as.character(path_melt$Grp), grps)
path_melt$Pathway <- factor(as.character(path_melt$Pathway), unique(unlist(sapply(qvals, names))))
path_melt$sigGrp <- path_grp[as.character(path_melt$Pathway)]

gg_hm <- ggplot(path_melt, aes(Sample, Pathway, fill= Pct)) + 
  facet_grid(cols=vars(Grp), rows=vars(sigGrp), scales = "free") +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu")
pdf(file.path(outdir, "mice_heatmap-wilcox.pdf"), width = 16, height=14)
gg_hm
dev.off()

###############################################################
#### Find Differentially expressed pathawys between groups ####
path_l <- getSamples(pathabundance, mice, return_ids = TRUE)
path_mat <- path_l$mat

grp_comparisons <- list("LB-inode"=c('GvgLmLr', 'GvgLiLj'),
                        "Fecal"=c('GvgFecal', 'GvgPBS'),
                        "AbxLR"=c('GvgLmLr', 'AfterAbx'),
                        "AbxLI"=c('GvgLiLj', 'AfterAbx'))
gg_hms <- lapply(grp_comparisons, function(grpc){
  grp1 <- path_mat[,grep(grpc[1], colnames(path_mat))]
  grp2 <- path_mat[,grep(grpc[2], colnames(path_mat))]
  pvals <- sapply(setNames(c(1:nrow(grp1)), rownames(grp1)), function(row_idx){
    t.test(unlist(grp1[row_idx,]), unlist(grp2[row_idx,]), 
                alternative = 'greater', na.rm=TRUE)$p.val
  })
  qvals <- p.adjust(pvals, method='fdr')
  sigp <- pvals[which(pvals < 0.05)]
  sigq <- qvals[which(qvals < 0.25)]
  grp_path <- as.data.frame(cbind(grp1, grp2))
  
  path_top <- grp_path[names(sigq),]
  # Make a heatmap of the top X functional pathways, separated by groups
  grps <- c("BeforeAbx","AfterAbx","GvgPBS","GvgFecal","GvgLiLj","GvgLmLr") 
  path_melt <- melt(t(path_top))
  path_melt$Grp <- gsub("[0-9]*$", "", path_melt$Var1)
  path_melt$Var1 <- path_l$ids[path_melt$Var1]
  colnames(path_melt) <- c('Sample', 'Pathway', 'Pct', 'Grp')

  path_melt$Grp <- factor(as.character(path_melt$Grp), grps)
  
  gg_hm <- ggplot(path_melt, aes(Sample, Pathway, fill= Pct)) + 
    facet_grid(cols=vars(Grp), scales = "free") +
    geom_tile() +
    ggtitle(paste(grpc, collapse="_"))
    scale_fill_distiller(palette = "RdPu")
  gg_hm
})
pdf(file.path(outdir, "mice_heatmap-between.pdf"), width = 14, height=17)
lapply(gg_hms, function(i) i)
dev.off()

