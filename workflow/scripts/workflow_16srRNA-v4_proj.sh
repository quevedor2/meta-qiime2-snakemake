
PROJDIR='/cluster/projects/mcgahalab/data/mcgahalab/PASS/16s'
cd ${PROJDIR}

#############
#### NEW ####
# Manifest and Import
PROJDIR='/cluster/projects/mcgahalab/data/mcgahalab/teresa_metagenomics/16s'
ID='teresa_metagenomics'
PROJDIR='/cluster/projects/mcgahalab/data/mcgahalab/PASS/16s'
ID='PASS'
cd ${PROJDIR}
conda activate ~/mcgahalab/envs/qiime2


## Assembling PE manifest for PHRED33v2
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" \
> data/artifacts/pe-manifest
echo -e "sample-id\tabsolute-filepath" \
> data/artifacts/se-manifest-f
echo -e "sample-id\tabsolute-filepath" \
> data/artifacts/se-manifest-r

for sid in $(cat data/samples.txt); do
  oid=$(echo ${sid} | cut -d"," -f1)
  nid=$(echo ${sid} | cut -d"," -f2)

  echo ${oid}
  r1=$(ls -d data/fastq/* | grep ${oid}"_R1")
  r2=$(ls -d data/fastq/* | grep ${oid}"_R2")
  r1path=$(readlink -f ${r1})
  r2path=$(readlink -f ${r2})

  # Format the names of the fastq samples
  echo -e "${nid}\t${r1path}\t${r2path}" >> data/artifacts/pe-manifest
  echo -e "${nid}\t${r1path}" >> data/artifacts/se-manifest-f
  echo -e "${nid}\t${r2path}" >> data/artifacts/se-manifest-r
done

## Import the manifest files
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path data/artifacts/pe-manifest \
--output-path data/artifacts/${ID}_demux.pe.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path data/artifacts/se-manifest-f \
--output-path data/artifacts/${ID}_demux.se_f.qza \
--input-format SingleEndFastqManifestPhred33V2

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path data/artifacts/se-manifest-r \
--output-path data/artifacts/${ID}_demux.se_r.qza \
--input-format SingleEndFastqManifestPhred33V2



# Sequence QC
## Simple visualization of the demultiplexed data
qcDIR="outputs/pe/seqQC"
mkdir -p ${qcDIR}
qiime demux summarize \
  --i-data data/artifacts/${ID}_demux.pe.qza \
  --o-visualization ${qcDIR}/${ID}_demux.qzv

for direction in f r; do
  qcDIR="outputs/se/${direction}/seqQC"
  mkdir -p ${qcDIR}
  qiime demux summarize \
    --i-data data/artifacts/${ID}_demux.se_${direction}.qza \
    --o-visualization ${qcDIR}/${ID}_demux.qzv
done

# Feature Table Creation
## Based on visualization, trim at pleft=13 and ptrunc=150
dadaDIR="outputs/pe/features/DADA2"
mkdir -p ${dadaDIR}
qiime dada2 denoise-paired \
  --p-trim-left-f 18 \
  --p-trim-left-r 12 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --i-demultiplexed-seqs data/artifacts/${ID}_demux.pe.qza \
  --o-representative-sequences ${dadaDIR}/${ID}.rep_seq.qza \
  --o-table ${dadaDIR}/${ID}.table.qza \
  --o-denoising-stats ${dadaDIR}/${ID}.stats.qza

for direction in f r; do
  dadaDIR="outputs/se/${direction}/features/DADA2"
  mkdir -p ${dadaDIR}

  if [[ ${direction} == 'f' ]]; then
    left=18
  else
    left=12
  fi

  qiime dada2 denoise-single \
    --p-trim-left ${left} \
    --p-trunc-len 150 \
    --i-demultiplexed-seqs data/artifacts/${ID}_demux.se_${direction}.qza \
    --o-representative-sequences ${dadaDIR}/${ID}.rep_seq.qza \
    --o-table ${dadaDIR}/${ID}.table.qza \
    --o-denoising-stats ${dadaDIR}/${ID}.stats.qza
done


## Visualizing denoising stats
dadaDIR="outputs/pe/features/DADA2"
mkdir -p ${dadaDIR}/views
### summarize features
qiime feature-table summarize \
  --i-table ${dadaDIR}/${ID}.table.qza \
  --o-visualization ${dadaDIR}/views/${ID}.table.qzv \
  --m-sample-metadata-file data/metadata/sample_metadata.tsv

### tabulate features
qiime metadata tabulate \
  --m-input-file ${dadaDIR}/${ID}.stats.qza \
  --o-visualization ${dadaDIR}/views/${ID}.stats.qzv

### tabule feature sequences
qiime feature-table tabulate-seqs \
  --i-data ${dadaDIR}/${ID}.rep_seq.qza \
  --o-visualization ${dadaDIR}/views/${ID}.rep-seqs.qzv

for direction in f r; do
  dadaDIR="outputs/se/${direction}/features/DADA2"
  mkdir -p ${dadaDIR}/views
  qiime feature-table summarize \
    --i-table ${dadaDIR}/${ID}.table.qza \
    --o-visualization ${dadaDIR}/views/${ID}.table.qzv \
    --m-sample-metadata-file data/metadata/sample_metadata.tsv

  qiime metadata tabulate \
    --m-input-file ${dadaDIR}/${ID}.stats.qza \
    --o-visualization ${dadaDIR}/views/${ID}.stats.qzv

  qiime feature-table tabulate-seqs \
    --i-data ${dadaDIR}/${ID}.rep_seq.qza \
    --o-visualization ${dadaDIR}/views/${ID}.rep-seqs.qzv
done

# Phylogenetic Diversity Analysis
phyloDIR='outputs/diversity/phylogenetic'
mkdir -p ${phyloDIR}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${dadaDIR}/${ID}.rep_seq.qza \
  --o-alignment ${phyloDIR}/${ID}.aligned-rep-seqs.qza \
  --o-masked-alignment ${phyloDIR}/${ID}.masked-aligned-rep-seqs.qza \
  --o-tree ${phyloDIR}/${ID}.unrooted-tree.qza \
  --o-rooted-tree ${phyloDIR}/${ID}.rooted-tree.qza

for direction in f r; do
  dadaDIR="outputs/se/${direction}/features/DADA2"
  phyloDIR="outputs/se/${direction}/diversity/phylogenetic"
  mkdir -p ${phyloDIR}
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences ${dadaDIR}/${ID}.rep_seq.qza \
    --o-alignment ${phyloDIR}/${ID}.aligned-rep-seqs.qza \
    --o-masked-alignment ${phyloDIR}/${ID}.masked-aligned-rep-seqs.qza \
    --o-tree ${phyloDIR}/${ID}.unrooted-tree.qza \
    --o-rooted-tree ${phyloDIR}/${ID}.rooted-tree.qza
done

scp ${XFER2}/nc-dotplot_anno.pdf .
scp ${XFER2}/nc-dotplot_clusters.pdf .
scp ${XFER2}/nc-dotplot_clusters.pdf .
scp ${XFER2}/nc-dimplot_clusters.pdf .
scp ${XFER2}/nc-featplots.pdf .


# Alpha and Beta Diversity Analysis
abDIR='outputs/diversity/alpha_beta'
mkdir -p ${abDIR}
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${phyloDIR}/${ID}.rooted-tree.qza \
  --i-table ${dadaDIR}/${ID}.table.qza \
  --p-sampling-depth 5500 \
  --m-metadata-file data/metadata/sample_metadata.tsv \
  --output-dir ${abDIR}/${ID}.core-metrics-results

for direction in f r; do
  abDIR="outputs/se/${direction}/diversity/alpha_beta"
  dadaDIR="outputs/se/${direction}/features/DADA2"
  phyloDIR="outputs/se/${direction}/diversity/phylogenetic"
  rm -r ${abDIR}
  mkdir -p ${abDIR}/views
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny ${phyloDIR}/${ID}.rooted-tree.qza \
    --i-table ${dadaDIR}/${ID}.table.qza \
    --p-sampling-depth 5500 \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --output-dir ${abDIR}/${ID}.core-metrics-results

  ## Visualize Alpha Diversity compared to Groups
  qiime diversity alpha-group-significance \
    --i-alpha-diversity ${abDIR}/${ID}.core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --o-visualization ${abDIR}/${ID}.core-metrics-results/faith-pd-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity ${abDIR}/${ID}.core-metrics-results/evenness_vector.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --o-visualization ${abDIR}/${ID}.core-metrics-results/evenness-group-significance.qzv

  ## Visualize Beta Diversity comparisons
  qiime diversity beta-group-significance \
    --i-distance-matrix ${abDIR}/${ID}.core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --m-metadata-column group-day \
    --o-visualization ${abDIR}/${ID}.core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
    --p-pairwise
done



# Taxonomic Analysis
for direction in f r; do
  dadaDIR="outputs/se/${direction}/features/DADA2"
  taxaDIR="outputs/se/${direction}/taxonomy"
  mkdir -p ${taxaDIR}/views
  ## Use a pre-trained Naive Bayes classifier to classify the taxa of each sequence
  classifier_dir='/cluster/projects/mcgahalab/ref/metagenomics/16s_classifiers'
  classifier='silva-138-99-nb-weighted-classifier.qza'

  qiime feature-classifier classify-sklearn \
    --i-classifier ${classifier_dir}/${classifier} \
    --i-reads ${dadaDIR}/${ID}.rep_seq.qza \
    --o-classification ${taxaDIR}/${ID}.taxonomy.qza

  qiime metadata tabulate \
  --m-input-file ${taxaDIR}/${ID}.taxonomy.qza \
  --o-visualization ${taxaDIR}/views/${ID}.taxonomy.qzv

  qiime taxa collapse \
    --i-table ${dadaDIR}/${ID}.table.qza \
    --i-taxonomy ${taxaDIR}/${ID}.taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${taxaDIR}/${ID}.lvl6-collapse.qza

  qiime taxa collapse \
    --i-table ${dadaDIR}/${ID}.table.qza \
    --i-taxonomy ${taxaDIR}/${ID}.taxonomy.qza \
    --p-level 7 \
    --o-collapsed-table ${taxaDIR}/${ID}.lvl7-collapse.qza

  qiime taxa barplot \
    --i-table ${dadaDIR}/${ID}.table.qza \
    --i-taxonomy ${taxaDIR}/${ID}.taxonomy.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --o-visualization ${taxaDIR}/views/${ID}.taxa-bar-plots.qzv
done

# Longitudinal Analysis
## https://docs.qiime2.org/2021.8/tutorials/longitudinal/?highlight=longitudinal
for direction in f r; do
  taxaDIR="outputs/se/${direction}/taxonomy"
  dadaDIR="outputs/se/${direction}/features/DADA2"
  longitudinalDIR="outputs/se/${direction}/longitudinal"
  rmdir ${longitudinalDIR}/lvl6_genus ${longitudinalDIR}/lvl7_species

  qiime longitudinal feature-volatility \
    --i-table ${taxaDIR}/${ID}.lvl6-collapse.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --p-state-column Day \
    --p-individual-id-column group-mouse \
    --p-n-estimators 10 \
    --p-random-state 17 \
    --output-dir ${longitudinalDIR}/lvl6_genus

  qiime longitudinal feature-volatility \
    --i-table ${taxaDIR}/${ID}.lvl7-collapse.qza \
    --m-metadata-file data/metadata/sample_metadata.tsv \
    --p-state-column Day \
    --p-individual-id-column group-mouse \
    --p-n-estimators 10 \
    --p-random-state 17 \
    --output-dir ${longitudinalDIR}/lvl7_species
done
