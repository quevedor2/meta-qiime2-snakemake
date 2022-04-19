from snakemake.utils import validate
import pandas as pd
import re

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
# Config file
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Samples: List of samples and conditions
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

# List of sample+unit information (e.g. paths, builds, etc.)
units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### test space #####
#u = units.loc[ ('net-037', '1'), ["fq1", "fq2"] ].dropna()
#print(u)
#print([ f"{u.fq1}", f"{u.fq2}" ])
#print("|".join(samples.index))
#print("|".join(units["unit"]))

##### wildcard constraints #####
wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])

'''
##### setting env paths #####
rule get_RlibPath:
    output:
        "results/ref/libpath"
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript -e \"cat(.libPaths(), '\n')\" > {output}"
'''

####### helpers ###########
def get_fq1(wildcards):
    """Get raw FASTQ files from unit sheet."""
    #u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
    u = units.loc[ (wildcards.sample, '1'), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq1}" ]

def get_fq2(wildcards):
    """Get raw FASTQ files from unit sheet."""
    #u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
    u = units.loc[ (wildcards.sample, '1'), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq2}" ]
