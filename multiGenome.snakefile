

rule all:
    input: expand("msa/{file}.aligned.fasta", file=config["multifasta"]),
            expand("msa/{file}_pct_dissimilarity.pdf", file=config["multifasta"]),
            expand("msa/PCvars_{file}/{file}.Genes.Ns.Ss.txt", file=config["multifasta"]),
            expand("msa/NCvars_{file}/{file}.NC.Genes.SNVs.txt", file=config["multifasta"])

rule MSA:
    input: "msa/{file}.fasta"
    output: "msa/{file}.aligned.fasta"
    shell:
        """
        mafft {input} > {output}
        """

rule snvFind:
    input:
        msaf= "msa/{file}.aligned.fasta"
    output:
        varpos= temp("msa/{file}_Mismatch_positions.bed"),
        tmpaln= temp("msa/{file}.tmp.file.aln"),
        smvar= "msa/{file}_Mismatch_positions_smoothed.bed"
    params:
        repf = config["RepeatRegions"],
        RefName = config["ReferenceSpecies"],
        #species = str(config["species"][0])+"_cpDNA"
        species = lambda wildcards: f"{config['multifasta'][wildcards.file]}"
    script:
        "scripts/SimilarityPlotter.py"

rule simPlot:
    input: "msa/{file}_Mismatch_positions_smoothed.bed"
    output: "msa/{file}_pct_dissimilarity.pdf"
    shell:
        "Rscript ~/Dropbox/codes/cpGenomics/scripts/plot_similarity.R {input} {output} "

rule PCgeneVars:
    input:
        msaf= "msa/{file}.aligned.fasta"
    output:
        "msa/PCvars_{file}/{file}.Genes.Ns.Ss.txt"
    params:
        #species = str(config["species"][0])+"_cpDNA",
        RefName = config["ReferenceSpecies"],
        geneBedf = config["RefGeneBed"],
        outdir = "msa/PCvars_{file}/",
        species = lambda wildcards: f"{config['multifasta'][wildcards.file]}"
    shell:
        """
        bash ~/Dropbox/codes/CGTools/dNdS_for_CDS.sh {input.msaf} {params.geneBedf} {params.RefName} {params.outdir} &> {params.outdir}/NS.log;
        mv msa/PCvars_{wildcards.file}/Genes.Ns.Ss.txt msa/PCvars_{wildcards.file}/{wildcards.file}.Genes.Ns.Ss.txt
        """

rule NCgeneVars:
    input:
        msaf= "msa/{file}.aligned.fasta"
    output:
        "msa/NCvars_{file}/{file}.NC.Genes.SNVs.txt"
    params:
        #species = str(config["species"][0])+"_cpDNA",
        RefName = config["ReferenceSpecies"],
        geneBedf = config["RefNCGeneBed"],
        outdir = "msa/NCvars_{file}/",
        species = lambda wildcards: f"{config['multifasta'][wildcards.file]}"
    shell:
        """
        bash ~/Dropbox/codes/CGTools/SynChanges_for_NCgenes.sh {input.msaf} {params.geneBedf} {params.RefName} {params.outdir} &> {params.outdir}/Syn.log;
        mv msa/NCvars_{wildcards.file}/NC.Genes.SNVs.txt msa/NCvars_{wildcards.file}/{wildcards.file}.NC.Genes.SNVs.txt
        """
