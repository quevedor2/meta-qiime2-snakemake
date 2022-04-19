condaprefix=$(readlink -f .snakemake/conda)

snakemake \
--use-conda \
--use-singularity \
--conda-create-envs-only \
--conda-frontend conda \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:///cluster/home/quever/downloads/snakemake-wrappers/' \
--cores 4
