rule metaphlan:
    input:
        fq1="results/kneaddata/main/{sample}.R1_kneaddata_paired_1.fastq.gz",
        fq2="results/kneaddata/main/{sample}.R1_kneaddata_paired_2.fastq.gz",
    output:
        table="results/metaphlan/main/{sample}.s1.tsv",
        biom="results/metaphlan/main/{sample}.biom",
        bowtie2="results/metaphlan/{sample}.s1.bowtie2.bz2",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
        database=directory(config['database']['bowtie2db']),
        cores=config['resources']['metaphlan_nrpoc'],
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/metaphlan/{sample}.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        metaphlan {input.fq1},{input.fq2} \
        -o {output.table} \
        -t rel_ab_w_read_stats \
        --input_type fastq \
        --index mpa_v30_CHOCOPhlAn_201901 \
        --bowtie2db {params.database} \
        --biom {output.biom} \
        --bowtie2out {output.bowtie2} \
        --nproc {params.cores}

        gzip results/kneaddata/main/{params.id}.*fastq;
        """

rule merge_metaphlan:
    input:
        get_all_metaphlan_tables,
    output:
        table="results/metaphlan/merged/merged_metaphlan.tsv",
    params:
        metaphlan_dir="results/metaphlan/main",
        regex="*s1.tsv$",
        conda=config['env']['conda_shell'],
        env=directory(config['env']['r41']),
        database=directory(config['database']['bowtie2db']),
        cores=config['resources']['metaphlan_nrpoc'],
    log:
        "logs/metaphlan/merge.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        Rscript mergeCountsMetaphlan.R \
        --dir {params.metaphlan_dir} \
        --pattern {params.regex} \
        --out {output}
        """
