rule kneaddata_pe:
    input:
        fq1=get_fq1,
        fq2=get_fq2,
    output:
        paired1="results/kneaddata/main/{sample}_{unit}.R1_kneaddata_paired_1.fastq.gz",
        paired2="results/kneaddata/main/{sample}_{unit}.R2_kneaddata_paired_1.fastq.gz",
    params:
        id="{sample}_{unit}"
        trimmomatic=config['tools']['trimmomatic'],
        database=config['database']['kneaddata_human'],
    conda:
        "../envs/biobakery3_core.yaml",
    log:
        "logs/kneaddata/{sample}_{unit}.log",
    shell:
        """
        kneaddata \
        --input {input.fq1} \
        --input {input.fq2} \
        -db {params.database} \
        --trimmomatic {params.trimmomatic} \
        --output {output};

        gzip results/kneaddata/main/{params.id}.*fastq;
        """
