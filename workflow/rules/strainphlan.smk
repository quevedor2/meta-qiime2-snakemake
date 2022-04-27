rule generate_markers:
    input:
        sam="results/metaphlan/main/sams/{sample}.sam",
    output:
        clade="results/strainphlan/consensus_markers/{sample}.pkl",
    params:
        dir="results/strainphlan/consensus_markers",
        cores=config['resources']['strainphlan_nproc'],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
    log:
        "logs/strainphlan/generate_markers_{sample}.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        sample2markers.py \
        --input {input.sam} \
        --input_format sam \
        --output_dir {params.dir} \
        --nprocs {params.cores}
        """

rule strainphlan_print_clades:
    input:
        expand("results/strainphlan/consensus_markers/{sample}.pkl", sample=samples.index),
    output:
        clade="results/strainphlan/consensus_markers/clade_list.txt",
    params:
        dir="results/strainphlan/consensus_markers",
        database=config['database']['metaphlan_pkl'],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
    log:
        "logs/strainphlan/print_clades.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        strainphlan \
        --database {params.database} \
        --samples {params.dir}/*.pkl \
        --output_dir {params.dir} \
        --print_clades_only > {output.clade}
        """

rule order_clade:
    input:
        abundance="results/metaphlan/merged/merged_metaphlan.relab.tsv",
        clade="results/strainphlan/consensus_markers/clade_list.txt",
    output:
        "results/strainphlan/consensus_markers/clades_list_order_by_average_abundance.txt",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
    log:
        "logs/strainphlan/order_clade.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        python scripts/order_clade_list.py \
        --abundance {input.abundance} \
        --clade {input.clade} \
        --output {output}
        """

rule strainphlan:
    input:
        pkls=expand("results/strainphlan/consensus_markers/{sample}.pkl", sample=samples.index),
        clades="results/strainphlan/consensus_markers/clades_list_order_by_average_abundance.txt",
    output:
        "results/strainphlan/strainphlan.out",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['biobakery3_core']),
        database=config['database']['metaphlan_pkl'],
        pkldir="results/strainphlan/consensus_markers",
        outdir="results/strainphlan",
        clade="",
        cores=config['resources']['strainphlan_nproc'],
#    conda:
#        "../envs/biobakery3_core.yaml",
    log:
        "logs/strainphlan/strainphlan.log",
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        for i in $(cat {input.clades}); do
            echo ${i} >> .tmp;
            strainphlan \
            --samples {params.pkldir}/*.pkl \
            --database {params.database} \
            --output_dir {params.outdir} \
            --clade ${i} \
            --nprocs {params.cores}
        done;
        
        if [ "$(wc -l < {input.clades})" -eq "$(wc -l < .tmp)" ]; then
            mv tmp {output};
        fi
        """
