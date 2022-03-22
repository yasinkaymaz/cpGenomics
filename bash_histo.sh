snakemake -s ~/Dropbox/codes/cpGenomics/Denovo.snakefile \
    --configfile ~/Dropbox/codes/cpGenomics/Denovo.config.yml \
    --cores 4 \
    --latency-wait 20 \
    --use-conda -n

snakemake -s ~/Dropbox/codes/cpGenomics/AssemblyCheck.snakefile \
--configfile ~/Dropbox/codes/cpGenomics/AssemblyCheck.config.yml \
--cores 4 \
--latency-wait 20 \
--use-conda -n

snakemake -s ~/Dropbox/codes/cpGenomics/multiGenome.snakefile \
--configfile ~/Dropbox/codes/cpGenomics/multiGenome.config.yml \
--cores 4 \
--latency-wait 2 \
--use-conda -n
