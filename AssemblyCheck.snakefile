
dockerTag = "latest" #FIXME tagged versions

rule all:
    input:
        expand("Reads2Assembly/{species}.Path{path}.bwtindex.1.bt2", path=['1', '2'], species=config["species"]),
        expand("Reads2Assembly/{species}.Path{path}.readsback.sorted.bam", path=['1', '2'], species=config["species"])

rule Bowtie2Index:
    input: "{species}.plastome/embplant_pt.K115.complete.graph1.{path}.path_sequence.fasta"
    output: "Reads2Assembly/{species}.Path{path}.bwtindex.1.bt2"
    params:
        prefix="Reads2Assembly/{species}.Path{path}.bwtindex"
    shell:
        """
        bowtie2-build {input} {params.prefix}
        """

rule Bowtie2map:
    input:
        sample=["trimmed/{species}_R1.fastq", "trimmed/{species}_R2.fastq"]
    output: pipe("Reads2Assembly/{species}.Path{path}.readsback.sam")
    params:
        prefix="Reads2Assembly/{species}.Path{path}.bwtindex"
    message: "now aligning {wildcards.species} Path{wildcards.path} to create {output}"
    shell:
        """
        bowtie2 -x {params.prefix} -1 {input.sample[0]} -2 {input.sample[1]} -S {output}
        """

rule Bamprocess:
    input: "Reads2Assembly/{species}.Path{path}.readsback.sam"
    output: "Reads2Assembly/{species}.Path{path}.readsback.sorted.bam",
            "Reads2Assembly/{species}.Path{path}.readsback.Coverage.stats.txt"
    shell:
        "bash ~/Dropbox/codes/CGTools/bamprocess.sh {input} "
