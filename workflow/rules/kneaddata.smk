rule kneaddata_pe:
    input:
        fq1=get_fq1,
        fq2=get_fq2,
    output:
        paired1="results/kneaddata/main/{sample}.knead1.fastq.gz",
        paired2="results/kneaddata/main/{sample}.knead2.fastq.gz",
    params:
        id="{sample}",
        outdir=directory("results/kneaddata/main/"),
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
        trimmomatic=directory(config['tools']['trimmomatic']),
        database=directory(config['database']['kneaddata']),
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/kneaddata/{sample}.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};
        module load pigz/2.6
        
        kneaddata \
        --input {input.fq1} \
        --input {input.fq2} \
        --bypass-trf \
        -db {params.database} \
        --trimmomatic {params.trimmomatic} \
        --output {params.outdir};
        
        echo "Gzipping...";
        pigz -p6 results/kneaddata/main/{params.id}.*fastq;
        echo "Gzipped and done."
        """
