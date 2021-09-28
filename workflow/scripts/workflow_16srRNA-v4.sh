# https://docs.qiime2.org/2021.4/tutorials/fmt/

# Download files
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake/data
wget 'https://data.qiime2.org/2021.4/tutorials/fmt/fmt-tutorial-demux-1-10p.qza'
wget 'https://data.qiime2.org/2021.4/tutorials/fmt/fmt-tutorial-demux-2-10p.qza'
wget 'https://data.qiime2.org/2021.4/tutorials/fmt/sample_metadata.tsv'

# Sequence QC
## Simple visualization of the demultiplexed data
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake
conda activate qiime2-2021.4
mkdir -p test/seqQC

for ID in fmt-tutorial-demux-1-10p fmt-tutorial-demux-2-10p; do
  qiime demux summarize \
    --i-data data/16S/${ID}.qza \
    --o-visualization test/seqQC/${ID}.qzv
done

## Based on visualization, trim at pleft=13 and ptrunc=150
mkdir -p test/seqQC/DADA2/
for ID in fmt-tutorial-demux-1-10p fmt-tutorial-demux-2-10p; do
  qiime dada2 denoise-single \
    --p-trim-left 13 \
    --p-trunc-len 150 \
    --i-demultiplexed-seqs data/16S/${ID}.qza \
    --o-representative-sequences test/seqQC/DADA2/${ID}.rep_seq.qza \
    --o-table test/seqQC/DADA2/${ID}.table.qza \
    --o-denoising-stats test/seqQC/DADA2/${ID}.stats.qza
done

## Visualizing denoising stats
for ID in fmt-tutorial-demux-1-10p fmt-tutorial-demux-2-10p; do
  qiime metadata tabulate \
    --m-input-file test/seqQC/DADA2/${ID}.stats.qza \
    --o-visualization test/seqQC/DADA2/${ID}.stats.qzv
done

# Merging the denoised data
## Creates a binary .biom file of features
qiime feature-table merge \
  --i-tables test/seqQC/DADA2/fmt-tutorial-demux-1-10p.table.qza \
  --i-tables test/seqQC/DADA2/fmt-tutorial-demux-2-10p.table.qza \
  --o-merged-table test/seqQC/DADA2/table.qza

qiime feature-table summarize \
  --i-table test/seqQC/DADA2/table.qza \
  --o-visualization test/seqQC/DADA2/table.qzv \
  --m-sample-metadata-file data/16S/sample_metadata.tsv

## Create merged fasta file
qiime feature-table merge-seqs \
  --i-data test/seqQC/DADA2/fmt-tutorial-demux-1-10p.rep_seq.qza \
  --i-data test/seqQC/DADA2/fmt-tutorial-demux-2-10p.rep_seq.qza \
  --o-merged-data test/seqQC/DADA2/rep-seqs.qza

qiime feature-table tabulate-seqs \
  --i-data test/seqQC/DADA2/rep-seqs.qza \
  --o-visualization test/seqQC/DADA2/rep-seqs.qzv


# Phylogenetic Diversity Analysis
mkdir -p test/diversity/phylogenetic/
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences test/seqQC/DADA2/rep-seqs.qza \
  --o-alignment test/phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment test/phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree test/phylogeny/unrooted-tree.qza \
  --o-rooted-tree test/phylogeny/rooted-tree.qza

# Alpha and Beta Diversity Analysis
mkdir -p test/diversity/alpha_beta
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny test/phylogeny/rooted-tree.qza \
  --i-table test/seqQC/DADA2/table.qza \
  --p-sampling-depth 500 \
  --m-metadata-file data/16S/sample_metadata.tsv \
  --output-dir test/diversity/alpha_beta/core-metrics-results

## Visualize Alpha Diversity compared to Groups
qiime diversity alpha-group-significance \
  --i-alpha-diversity test/diversity/alpha_beta/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file data/16S/sample_metadata.tsv \
  --o-visualization test/diversity/alpha_beta/core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity test/diversity/alpha_beta/core-metrics-results/evenness_vector.qza \
  --m-metadata-file data/16S/sample_metadata.tsv \
  --o-visualization test/diversity/alpha_beta/core-metrics-results/evenness-group-significance.qzv

## Visualize Beta Diversity comparisons
qiime diversity beta-group-significance \
  --i-distance-matrix test/diversity/alpha_beta/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file data/16S/sample_metadata.tsv \
  --m-metadata-column sample-type \
  --o-visualization test/diversity/alpha_beta/core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

# Taxonomic Analysis
mkdir -p test/taxonomy
## Use a pre-trained Naive Bayes classifier to classify the taxa of each sequence
classifier='silva-138-99-nb-classifier.qza'
qiime feature-classifier classify-sklearn \
  --i-classifier data/16S/classifiers/${classifier} \
  --i-reads test/seqQC/DADA2/rep-seqs.qza \
  --o-classification test/taxonomy/taxonomy.qza

qiime metadata tabulate \
--m-input-file test/taxonomy/taxonomy.qza \
--o-visualization test/taxonomy/taxonomy.qzv

qiime taxa barplot \
  --i-table test/seqQC/DADA2/table.qza \
  --i-taxonomy test/taxonomy/taxonomy.qza \
  --m-metadata-file data/16S/sample_metadata.tsv \
  --o-visualization test/taxonomy/taxa-bar-plots.qzv
