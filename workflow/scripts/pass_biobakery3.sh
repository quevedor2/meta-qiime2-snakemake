conda activate biobakery3

PROJDIR='/cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun'
cd ${PROJDIR}


biobakery_workflows wmgx \
--input data \
--output results/biobakery3 \
--qc-options "--trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2
              --index mpa_v30_CHOCOPhlAn_201901
              --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--contaminate-databases "/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--dry-run

biobakery_workflows wmgx \
--input data \
--output results/biobakery3 \
--qc-options "--trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2" \
--taxonomic-profiling-options "-t rel_ab_w_read_stats --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--functional-profiling-options "--protein-database /cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/humann/uniref" \
--contaminate-databases "/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--grid-jobs 6 --threads 8 \
--grid-scratch '/cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/scratch'

biobakery_workflows wmgx \
--input input \
--output output_bkup \
--qc-options "--trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2" \
--taxonomic-profiling-options "-t rel_ab_w_read_stats --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--functional-profiling-options "--protein-database /cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/humann/uniref" \
--contaminate-databases "/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--local-jobs 2 --threads 2

PDIR='/cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun'
KNEAD="${PDIR}/results/biobakery3/kneaddata/main"
DATABASE='/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome'
cd ${PDIR}

for id in $(cat samples.txt); do
  echo ${id}
  id='PASS_S10'
  kneaddata --input ${PDIR}/data/${id}.R1.fastq.gz \
  --input ${PDIR}/data/${id}.R2.fastq.gz \
  -db ${DATABASE} \
  --trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2 \
  --output ${KNEAD}
done

PDIR='/cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun'
KNEAD="${PDIR}/results/biobakery3/kneaddata/main"
cd ${PDIR}
for id in $(cat samples.txt); do
  echo ${id}
  metaphlan ${KNEAD}/${id}_paired_1.fastq,${KNEAD}/${id}_paired_2.fastq \
  -t rel_ab_w_read_stats \
  -o ${PDIR}/results/metaphlan/${id}.s1.tsv \
  --input_type fastq \
  --index mpa_v30_CHOCOPhlAn_201901 \
  --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases \
  --biom ${PDIR}/results/metaphlan/${id}.biom \
  --bowtie2out ${PDIR}/results/metaphlan/${id}.s1.bowtie2.bz2 \
  --nproc 4
done

conda activate biobakery3viz
module load texlive/2019
biobakery_workflows wmgx_vis \
--input /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_bkup \
--project-name Demo_Test \
--output /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_bkup/viz

conda activate biobakery3viz
module load texlive/2019
biobakery_workflows wmgx_vis \
--input /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3 \
--project-name PASS \
--output /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/viz

conda activate biobakery3viz
module load texlive/2019
biobakery_workflows stats \
--input /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3 \
--input-metadata $INPUT_METADATA \
--output /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/stats \
--project-name PASS





/cluster/home/quever/miniconda3/envs/biobakery3viz/share/tlpkg
/cluster/home/quever/miniconda3/envs/biobakery3viz/share/texmf-dist/scripts/texlive
/cluster/home/quever/miniconda3/envs/biobakery3viz/lib/site_perl/5.26.2/x86_64-linux-thread-multi
/cluster/home/quever/miniconda3/envs/biobakery3viz/lib/site_perl/5.26.2
/cluster/home/quever/miniconda3/envs/biobakery3viz/lib/5.26.2/x86_64-linux-thread-multi
/cluster/home/quever/miniconda3/envs/biobakery3viz/lib/5.26.2


humann \
--input /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/kneaddata/main/PASS_S23.fastq \
--output /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/humann/main \
--o-log /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/humann/main/PASS_S23.log \
--threads 3 \
--taxonomic-profile /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/metaphlan/main/PASS_S23_taxonomic_profile.tsv  \
--protein-database /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/biobakery_workflows_databases/humann/uniref


Aligning to reference database: uniref90_201901b_full.dmnd

CRITICAL ERROR: Error executing:
/cluster/home/quever/miniconda3/envs/biobakery3/bin/diamond blastx \
--query /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/humann/main/PASS_S23_humann_temp/PASS_S23_bowtie2_unaligned.fa \
--evalue 1.0 \
--threads 3 \
--top 1 \
--outfmt 6 \
--db /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/biobakery_workflows_databases/humann/uniref/uniref90_201901b_full \
--out /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/humann/main/PASS_S23_humann_temp/tmpog5o7htr/diamond_m8_9ikd_0np \
--tmpdir /cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/biobakery3/humann/main/PASS_S23_humann_temp/tmpog5o7htr

Error message returned from diamond :
diamond v0.9.36.137 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

#CPU threads: 3
grep -E "s__|taxonomy" metaphlan_taxonomic_profiles.tsv | sed 's/^.*s__//g' | sed -e 's/# taxon omy/body_site/g' > merged_abundance_table_species.txt
hclust2.py -i merged_abundance_table_species.txt -o abundance_heatmap_species.png \
--f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 \
-l --flabel_size 1 --slabel_size 1 --max_flabel_len 100 --max_slabel_len 100 \
--minv 0.1 --dpi 500

ln -s metaphlan_taxonomic_profiles.tsv  merged_abundance_table_reformatted.txt
export2graphlan.py --skip_rows 1 -i merged_abundance_table_reformatted.txt \
--tree merged_abundance.tree.txt --annotation merged_abundance.annot.txt \
--most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
