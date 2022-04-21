rule merge_fq:
    input:
        fq1="results/kneaddata/main/{sample}.R1_kneaddata_paired_1.fastq.gz",
        fq2="results/kneaddata/main/{sample}.R1_kneaddata_paired_2.fastq.gz",
    output:
        "results/humann/main/{sample}.merged.fastq.gz",
    params:
        id="{sample}",
    log:
        "logs/humann/{sample}.merge_fq.log",
    shell:
        """
        zcat {input.fq1} {input.fq2} | gzip > {output}
        """

rule humann:
    input:
        metaphlan_table="results/metaphlan/main/{sample}.s1.tsv",
        fq="results/humann/main/{sample}.merged.fastq.gz",
    output:
        genefamilies="results/humann/main/{sample}.merged_genefamilies.tsv",
        pathabundance="results/humann/main/{sample}.merged_pathabundance.tsv",
        pathcoverage="results/humann/main/{sample}.merged_pathcoverage.tsv",
        outdir=temp(directory("results/humann/main/{sample}.merged_humann_temp")),
        log="results/humann/main/{sample}.log",
    params:
        outdir="results/humann/main/",
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
        nucleotide_db=directory(config['database']['nucleotide']),
        protein_db=directory(config['database']['protein']),
        cores=config['resources']['humann_proc'],
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/humann/{sample}.humann.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        humann \
        --threads {params.cores} \
        --input {input.fq} \
        --input-format fastq.gz \
        --output {params.outdir} \
        --nucleotide-database {params.nucleotide_db} \
        --protein-database {params.protein_db} \
        --taxonomic-profile {input.metaphlan_table} \
        --o-log {output.log}
        """

rule regroup_humann:
    input:
        "results/humann/main/{sample}.merged_genefamilies.tsv",
    output:
        "results/humann/main/regrouped/{sample}.ecs.tsv",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/humann/{sample}.regroup.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        humann_regroup_table \
        --input {input} \
        --output {output} \
        --group uniref90_level4ec
        """

rule merge_human_raw:
    input:
        ecs=expand("results/humann/main/regrouped/{sample}.ecs.tsv", sample=samples.index),
        genefamilies=expand("results/humann/main/{sample}.merged_genefamilies.tsv", sample=samples.index),
        pathabundance=expand("results/humann/main/{sample}.merged_pathabundance.tsv", sample=samples.index),
    output:
        ecs="results/humann/merged/ecs.tsv",
        genefamilies="results/humann/merged/genefamilies.tsv",
        pathabundance="results/humann/merged/pathabundance.tsv",
    params:
        ecs="results/humann/main/regrouped/",
        genefamilies="results/humann/main/",
        pathabundance="results/humann/main/",
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/humann/join_table_raw.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        humann_join_tables \
        --input {params.ecs} \
        --output {output.ecs} \
        --file_name ecs;

        humann_join_tables \
        --input {params.genefamilies} \
        --output {output.genefamilies} \
        --file_name genefamilies;

        humann_join_tables \
        --input {params.pathabundance} \
        --output {output.pathabundance} \
        --file_name pathabundance;
        """

rule humann_normalize:
    input:
        genefamilies="results/humann/main/{sample}.merged_genefamilies.tsv",
        pathabundance="results/humann/main/{sample}.merged_pathabundance.tsv",
        ecs="results/humann/main/regrouped/{sample}.ecs.tsv",
    output:
        ecs="results/humann/relab/{sample}.relab.ecs.tsv",
        genefamilies="results/humann/relab/{sample}.relab.genefamilies.tsv",
        pathabundance="results/humann/relab/{sample}.relab.pathabundance.tsv",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/humann/{sample}_renorm.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        humann_renorm_table \
        --input {input.ecs} \
        --output {output.ecs} \
        --units relab \
        --special n;

        humann_renorm_table \
        --input {input.genefamilies} \
        --output {output.genefamilies} \
        --units relab \
        --special n;

        humann_renorm_table \
        --input {input.pathabundance} \
        --output {output.pathabundance} \
        --units relab \
        --special n;
        """


rule merge_human_relab:
    input:
        ecs=expand("results/humann/relab/{sample}.relab.ecs.tsv", sample=samples.index),
        genefamilies=expand("results/humann/relab/{sample}.relab.genefamilies.tsv", sample=samples.index),
        pathabundance=expand("results/humann/relab/{sample}.relab.pathabundance.tsv", sample=samples.index),
    output:
        ecs="results/humann/merged/ecs.relab.tsv",
        genefamilies="results/humann/merged/genefamilies.relab.tsv",
        pathabundance="results/humann/merged/pathabundance.relab.tsv",
    params:
        dir="results/humann/relab/",
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/humann/join_table_relab.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};

        humann_join_tables \
        --input {params.dir} \
        --output {output.ecs} \
        --file_name ecs;

        humann_join_tables \
        --input {params.dir} \
        --output {output.genefamilies} \
        --file_name genefamilies;

        humann_join_tables \
        --input {params.dir} \
        --output {output.pathabundance} \
        --file_name pathabundance;
        """
