
dockerTag = "latest" #FIXME tagged versions



rule all:
    input:
        expand("FastQC/{sample}_{num}_fastqc.zip", sample=config["samples"], num=['1', '2']),
        expand("{sample}.plastome/embplant_pt.K115.complete.graph1.1.path_sequence.fasta", sample=config["samples"]),
        #expand("Path{num}.seq.png", num=['1', '2']),
        #expand("Reads2Assembly/Path{num}.bowtie2_idx.1.bt2", num=['1', '2']),
        #"Reads2Assembly/Lens_lamottei.Path1.readsback.sorted.bam"
        #    expand("Reads2Assembly/{sample}.Path{num}.readsback.sorted.bam", num=['1', '2'], sample=config['samples'])
        #expand("{sample}.plastome/Path.{num}.sequence.png", sample=config["samples"], num=['1', '2'])
        #"utdir/fastq_multiqc.html",
        #expand("gatk/{sample}_SnpEff.vcf.gz", sample=config["samples"]),
        #expand("gatk/{sample}_snpEff_summary.html", sample=config["samples"])

rule FastQC:
    input:
        fq=lambda wildcards: f"{config['samples'][wildcards.sample]}_{wildcards.num}.fq.gz"
    output:
        html="FastQC/{sample}_{num}_fastqc.html",
        zip="FastQC/{sample}_{num}_fastqc.zip"
    shell:
        "fastqc "
        "{input.fq} "
        "--outdir FastQC"

rule AdpTrimFastp:
    input:
        sample=["fastq_files/{sample}_1.fq.gz", "fastq_files/{sample}_2.fq.gz"]
    output:
        trimmed=["trimmed/{sample}_R1.fastq", "trimmed/{sample}_R2.fastq"],
        unpaired1=temp("trimmed/{sample}.u1.fastq"),
        unpaired2=temp("trimmed/{sample}.u2.fastq"),
        merged=temp("trimmed/{sample}.merged.fastq"),
        failed=temp("trimmed/{sample}.failed.fastq"),
        html="report/{sample}.html",
        json="report/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    wrapper:
        "v0.86.0/bio/fastp"

rule GetOrganelle:
    input:
        fq1="trimmed/{sample}_R1.fastq",
        fq2="trimmed/{sample}_R2.fastq"
    output:
        Assembly="{sample}.plastome/embplant_pt.K115.complete.graph1.1.path_sequence.fasta"
    params:
        ref=config['GenomeFa']
    shell:
        "get_organelle_from_reads.py -1 {input.fq1} -2 {input.fq2} "
        "-t 4 -o {wildcards.sample}.plastome "
        "-F embplant_pt -R 15 -w 102 --reduce-reads-for-coverage inf --max-reads inf "
        "-s {params.ref} "
        "--overwrite "

# rule Nucmer:
#     input:
#         infa=expand("{sample}.plastome/embplant_pt.K115.complete.graph1.{num}.path_sequence.fasta", sample=config["samples"], num=['1', '2'])
#     output:
#         #expand("{sample}.plastome/Path.{num}.sequence.png", sample=config["samples"], num=['1', '2'])
#         expand("{sample}.plastome/Path{num}.seq.png", num=['1', '2'], sample=config["samples"])
#     params:
#         ref=config['GenomeFa'],
#         prefix=expand("{sample}.plastome/{sample}.plastome", sample=config["samples"])
#     shell:
#         """
#         for i in 1 2;
#             do
#             nucmer --mum {params.prefix}/embplant_pt.K115.complete.graph1.$i.path_sequence.fasta {params.ref} -p Path$i.seq && \
#             delta-filter -l 1000 -q Path$i.seq.delta > Path$i.seq.filtered.delta && \
#             mummerplot -png -f Path$i.seq.filtered.delta -p Path$i.seq;
#         done
#         """
#
#
# rule Bowtie2Index:
#     input:
#         infa=expand("{sample}.plastome/embplant_pt.K115.scaffolds.graph1.{num}.path_sequence.fasta", sample=config["samples"], num=['1', '2'])
#     output:
#         expand("Reads2Assembly/Path{num}.bowtie2_idx.1.bt2", num=['1', '2'])
#     params:
#         prefix=expand("{sample}.plastome", sample=config["samples"])
#     shell:
#         """
#         for i in 1 2;
#             do
#             bowtie2-build {params.prefix}/embplant_pt.K115.scaffolds.graph1.$i.path_sequence.fasta Reads2Assembly/Path$i.bowtie2_idx ;
#         done
#         """
#
#
# rule Bowtie2map:
#     input:
#         sample=["trimmed/{sample}_R1.fastq", "trimmed/{sample}_R2.fastq"]
#         #fq1=lambda wildcards: f"{config['samples'][wildcards.sample]}_1.fq.gz",
#         #fq2=lambda wildcards: f"{config['samples'][wildcards.sample]}_2.fq.gz"
#     output:
#         bams=["Reads2Assembly/{sample}.Path1.readsback.sorted.bam",
#             "Reads2Assembly/{sample}.Path2.readsback.sorted.bam"]
#     message: "now aligning {input} to create {output}"
#     params:
#         sampleprefix=config['samples']
#     shell:
#         """
#         for i in 1 2;
#             do
#             bowtie2 -x Reads2Assembly/Path$i.bowtie2_idx -1 {input.sample[0]} -2 {input.sample[1]} -S {params.sampleprefix}.Path$i.readsback.sam;
#             bash ~/Dropbox/codes/CGTools/bamprocess.sh {params.sampleprefix}.Path$i.readsback.sam;
#         done
#         """

# bowtie2 -x Lens_ervo_idx -1 s7_1.fq.gz -2 s7_2.fq.gz -S Lens_ervo_readsback.sam
#
# bash ~/Dropbox/codes/CGTools/bamprocess.sh Lens_ervo_readsback.sam
#
# #https://www.biostars.org/p/5165/
# samtools depth  Lens_ervo_readsback.sorted.bam |  awk '{sum+=$3} END { print "Average = ", sum/NR}'
#
#
#
# #SNVs in coding genes
# bash ~/Dropbox/codes/cpLensErvo/dNdS_for_CDS.sh \
# ../Lens_cul-ervo_MafftAligned.fa \
# ~/Dropbox/codes/cpLensErvo/data/Lens_culinaris_CDS_genes.bed Lens_culinaris &> NS.log
#
#
# for i in `awk '($2 != 0)' Genes.Ns.Ss.txt |\
# grep -v Gene|\
# cut -f1`;
# do
#   grep -w $i ~/Dropbox/codes/cpLensErvo/data/Lens_culinaris_GFF3.gff3 |\
#   awk '{if($3 == "CDS") print}'|\
#   cut -f9-|\
#   uniq |\
#   cut -d";" -f8 |\
#   sed 's/product=//g' ;
# done
#
# #SNVs in non-coding genes
# bash ~/Dropbox/codes/cpLensErvo/SynChanges_for_NCgenes.sh \
# ../Lens_cul-ervo_MafftAligned.fa \
# ~/Dropbox/codes/cpLensErvo/data/Lens_culinaris_NC_genes.bed Lens_culinaris &> NS.log
