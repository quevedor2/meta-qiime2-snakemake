# https://docs.qiime2.org/2021.4/tutorials/fmt/

# Download files
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake/data/shotgun
module load sratoolkit/2.11.0
prefetch SRR2726671
prefetch SRR2726672
fasterq-dump --split-files SRR2726671/SRR2726671.sra
fasterq-dump --split-files SRR2726672/SRR2726672.sra

# Setup QIIME env
## SHOGUN
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake
conda activate qiime2-2021.4
pip install https://github.com/knights-lab/SHOGUN/archive/master.zip
pip install https://github.com/qiime2/q2-shogun/archive/master.zip
qiime dev refresh-cache

## Biobakery3
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake/data
salloc --partition=all -c 1 -t 8:0:0 --mem 20
conda create --name biobakery3
conda activate biobakery3
conda install -c biobakery biobakery_workflows
conda install tbb=2020.2
conda install networkx=1.11
conda install diamond=0.9.36
#The following packages will be REMOVED:
#  boost-1.68.0-py37h8619c78_1001
#The following packages will be UPDATED:
#  blast                                    2.9.0-h20b68b9_1 --> 2.11.0-pl526he19e7b1_0
#  boost-cpp                            1.68.0-h11c811c_1000 --> 1.70.0-ha2d47e9_1
#  diamond                                 0.9.24-ha888412_1 --> 0.9.34-h56fc30b_0
metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 # Clade markers and database
#run_command(["humann_databases","--download","chocophlan","full",humann_install_folder])
#run_command(["humann_databases","--download","uniref","uniref90_diamond",humann_install_folder])
db_dir='/cluster/home/quever/biobakery_workflows_databases'
humann_install_folder="${db_dir}/humann"
knead_install_folder="${db_dir}/kneaddata_db_human_genome"
humann_databases --download chocophlan full ${humann_install_folder}
humann_databases --download uniref uniref90_diamond ${humann_install_folder}
kneaddata_database --download human_genome bowtie2 ${knead_install_folder}
#biobakery_workflows_databases --install wmgx

wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/HD32R1_subsample.fastq.gz'
wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/HD42R4_subsample.fastq.gz'
wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/HD48R4_subsample.fastq.gz'
wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/LD96R2_subsample.fastq.gz'
wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/LV16R4_subsample.fastq.gz'
wget 'https://github.com/biobakery/biobakery_workflows/raw/master/examples/tutorial/input/LV20R4_subsample.fastq.gz'
 
biobakery_workflows wmgx \
--input input \
--output output_data \
--qc-options "--trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2
 --index mpa_v30_CHOCOPhlAn_201901
 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--contaminate-databases "/cluster/projects/pughlab/projects/cancer_cell_lines/tmp/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--dry-run

 metaphlan /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_data/kneaddata/main/HD32R1_subsample.fastq \
 --input_type fastq \
 --output_file /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_data/metaphlan/main/HD32R1_subsample_taxonomic_profile.tsv \
 --samout /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_data/metaphlan/main/HD32R1_subsample_bowtie2.sam \
 --nproc 1 \
 --no_map \
 --tmp_dir /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_data/metaphlan/main \
 --index mpa_v30_CHOCOPhlAn_201901 \
 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases

 
biobakery_workflows wmgx \
--input shotgun \
--output shotgun_out \
--bypass-strain-profiling \
--pair-identifier _1 \
--input-extension fastq

# Input FASTQ into a qza artifact input
conda activate qiime2-2021.4
cd /mnt/work1/users/home2/quever/git/meta-qiime2-snakemake/data/shotgun
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path sample_metadata.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  

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
