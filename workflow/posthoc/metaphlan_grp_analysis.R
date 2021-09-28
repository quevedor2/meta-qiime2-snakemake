library(reshape2)
library(ggplot2)
library(DESeq2)


PDIR <- '/Users/rquevedo/Projects/mcgaha/PASS/shotgun/data'
outdir <- '~/Projects/mcgaha/PASS/shotgun/results/species/'
setwd(PDIR)

###################
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

# Melt the metaphlan data structure for a subset of samples and subset to a top representative
meltMetaphlan <- function(meta_spl, samples, topx=10){
  melt_meta <- lapply(meta_spl, function(i){
    # Select only the top X species/genus/etc...
    i <- getSamples(i, samples)
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


#############################################
#### Metaphlan - genus/species breakdown ####
rank_legend <- c('k'='Domain',
                 'p'='Phylum',
                 'c'='Class',
                 'o'='Order',
                 'f'='Family',
                 'g'='Genus',
                 's'='Species')
controls <- list("L.intestinalis"="S35",
                 "L.johnsonii"="S36",
                 "L.reuteri"="S33",
                 "L.murinus"="S34")
mice <- list("BeforeAbx"=paste0("S", c(1:6)),
             "AfterAbx"=paste0("S", c(7:12)),
             "GvgFecal"=paste0("S", c(28:32)),
             "GvgLiLj"=paste0("S", c(18:22)),
             "GvgLmLr"=paste0("S", c(13:17)),
             "GvgPBS"=paste0("S", c(23:27)))
topx <- 7

# Read in the metaphlan data and preformat
metaphlan <- read.table("metaphlan_taxonomic_profiles.tsv", sep="\t", 
                        header=TRUE, comment.char = '',
                        stringsAsFactors = FALSE, check.names = FALSE)
rownames(metaphlan) <- metaphlan[,1]
metaphlan <- metaphlan[,-1]
colnames(metaphlan) <- gsub("PASS_(S[0-9]+)_.*", "\\1", colnames(metaphlan))
meta_levels <- sapply(strsplit(rownames(metaphlan), split="\\|"), length)
meta_spl <- split(metaphlan, meta_levels)
names(meta_spl) <- rank_legend

############################
#### PCA of the Species ####
gg_pca <- lapply(rank_legend[-1], function(rank){
  meta_l <- getSamples(meta_spl[[rank]], mice, return_ids = TRUE)
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
    ggtitle(rank) +
    theme_minimal()
  gg_pca
})

pdf(file.path(outdir, "genus_species_pca.pdf"))
lapply(gg_pca, function(i) i)
dev.off()

#####################################
#### Comparing changes in groups ####
ref_grp <- 'GvgPBS' #  
ref_grps <- c('BeforeAbx', 'AfterAbx', 'GvgPBS')
grps <- c("BeforeAbx","AfterAbx","GvgPBS","GvgFecal","GvgLiLj","GvgLmLr") # unique(gsub("[0-9]+$", "", colnames(meta_species)))

## Cycle through each rank (e.g. Phylla, Family, Order, Genus, Species)
gg_bps <- lapply(rank_legend[-1], function(rank){
  ## Set the Reference group to compare against (e.g. BeforeAbx, AfterAbx)
  comp_ref_pcts <- lapply(setNames(ref_grps,ref_grps), function(ref_grp){
    # Split the metaphlan species table based on the rank of the phylla
    meta_species <- getSamples(meta_spl[[rank]], mice)
    rownames(meta_species) <- gsub("^.*[kpcofgs]__", "", rownames(meta_species))
    meta_species <- meta_species[-which(rowSums(meta_species)==0),]
    meta_l <- lapply(setNames(grps, grps), function(grp){
      meta_species[,grep(grp, colnames(meta_species))]
    })
    
    # Calculate the t-statistic difference between the comparing group and the reference group
    ref_grp_idx <- match(ref_grp, names(meta_l))
    delta_pct <- lapply(meta_l[-ref_grp_idx], function(i){
      df <- data.frame("delta"=sapply(c(1:nrow(i)), function(row_i) t.test(i[row_i,], meta_l[[ref_grp_idx]][row_i,])$statistic),
                       # rowMeans(i) - rowMeans(meta_l[[ref_grp_idx]]),
                       "Grp"=unique(gsub("[0-9]*$", "", colnames(i))),
                       "Ref"=ref_grp,
                       "Organism"=rownames(i),
                       "Rank"=rank)
      if(any(na.omit(df$delta < -6))) df[which(df$delta < -6),]$delta <- -6
      if(any(na.omit(df$delta > 6))) df[which(df$delta > 6),]$delta <- 6
      return(df)
    })
    delta_pct <- as.data.frame(do.call(rbind,delta_pct))
    
    # Zero out non-comparing groups
    if(ref_grp=='GvgPBS'){
      zero_grp <- c('BeforeAbx', 'AfterAbx')
    } else if(ref_grp=='AfterAbx'){
      zero_grp <- 'BeforeAbx'
    } else {
      zero_grp <- NULL
    }
    if(!is.null(zero_grp)){
      delta_pct[which(delta_pct$Grp %in% zero_grp),]$delta <- 0
    }
    
    return(delta_pct)
  })
  
  delta_pct <- as.data.frame(do.call(rbind, comp_ref_pcts))
  sum_delta <- sapply(split(delta_pct, delta_pct$Organism), function(i) mean((i$delta), na.rm = TRUE))
  sum_delta[is.nan(sum_delta)] <- 0
  delta_pct$Organism <- factor(delta_pct$Organism, names(sort(sum_delta, decreasing = FALSE)))
  delta_pct$Grp <- factor(as.character(delta_pct$Grp), grps)
  delta_pct$Ref <- factor(as.character(delta_pct$Ref), ref_grps)
  
  gg_bps <- ggplot(delta_pct, aes(x=Organism, y=delta, fill=Ref)) +
    facet_grid(cols=vars(Grp), scales = "free") +
    geom_bar(stat='identity', position = position_dodge(width=0.5)) +
    ylim(-6, 6) +
    ggtitle(rank) +
    theme_minimal() +
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(gg_bps)
})

pdf(file.path(outdir, "delta_comp-ref-mice.pdf"), width = 10)
lapply(gg_bps, function(i) i)
dev.off()

##########################################################
#### Stacked barplot of organisms within each phyllum ####
## Breakdown the composition of organisms across all the mice
melt_meta <- meltMetaphlan(meta_spl, mice, topx = topx)
gps <- lapply(names(melt_meta), function(r){
  mm <- melt_meta[[r]]
  x <- ggplot(mm, aes(fill=ID, y=Percentage, x=Sample)) + 
    #facet_grid(vars(Rank), scales = "free") +
    geom_bar(position="stack", stat="identity") +
    ggtitle(r) +
    xlab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  x
})

pdf(file.path(outdir, "mice_breakdown.pdf"))
lapply(gps, function(i) i)
dev.off()



## Breakdown the composition of organisms across all the bacterial controls
melt_meta <- meltMetaphlan(meta_spl, controls, topx = topx)
gps <- lapply(names(melt_meta), function(r){
  mm <- melt_meta[[r]]
  x <- ggplot(mm, aes(fill=ID, y=Percentage, x=Sample)) + 
    #facet_grid(vars(Rank), scales = "free") +
    geom_bar(position="stack", stat="identity") +
    ggtitle(r) +
    xlab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  x
})

pdf(file.path(outdir, "bacterial_breakdown.pdf"))
lapply(gps, function(i) i)
dev.off()

