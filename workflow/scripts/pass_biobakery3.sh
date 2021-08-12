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
--taxonomic-profiling-options "--index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--functional-profiling-options "--protein-database /cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/humann/uniref" \
--contaminate-databases "/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--grid-jobs 6 --threads 8 \
--grid-scratch '/cluster/projects/mcgahalab/data/mcgahalab/PASS/shotgun/results/scratch'

biobakery_workflows wmgx \
--input input \
--output output_bkup \
--qc-options "--trimmomatic /cluster/home/quever/miniconda3/envs/biobakery3/share/trimmomatic-0.39-2" \
--taxonomic-profiling-options "--index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /cluster/home/quever/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases" \
--functional-profiling-options "--protein-database /cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/humann/uniref" \
--contaminate-databases "/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/kneaddata_db_human_genome" \
--bypass-strain-profiling \
--local-jobs 2 --threads 2



--input data \
--output results/biobakery3 \


module load texlive/2019
biobakery_workflows wmgx_vis \
--input /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_bkup \
--project-name Demo_Test \
--output /cluster/projects/pughlab/projects/cancer_cell_lines/tmp/tutorial/output_bkup/viz

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
