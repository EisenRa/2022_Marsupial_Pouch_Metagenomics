################################################################################
################################################################################
################################################################################
# Snakefile for quality controlling (trimming/mapping to host) reads and assembly
# Raphael Eisenhofer 10/2021
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
          for fn in glob(f"2_Reads/0_Untrimmed/*_1.fastq.gz")]

print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz", sample=SAMPLE),
        expand("2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz", sample=SAMPLE),
        expand("3_Outputs/2_Assemblies/{sample}/{sample}_contigs.fasta", sample=SAMPLE)
################################################################################
### Download reads (do this closer to submission)


################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i = "2_Reads/0_Untrimmed/{sample}_1.fastq.gz",
        r2i = "2_Reads/0_Untrimmed/{sample}_2.fastq.gz"
    output:
        r1o = "2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2o = "2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz",
        fastp_html = "2_Reads/2_fastp_results/{sample}.html",
        fastp_json = "2_Reads/2_fastp_results/{sample}.json"
    conda:
        "1_QC_Assembly.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_fastp.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 40 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --adapter_sequence_r2 CTGTCTCTTATACACATCT \
        &> {log}
        """
################################################################################
### Download host reference genomes and create a bowtie2 index:
# rule download_refs_and_index:
#     output:
#         catted_ref = "1_References/CattedRefs.fna.gz",
#         bt2_index = "1_References/CattedRefs.fna.rev.2.bt2l",
#         temp("GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.fna.gz"),
#         temp("GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz")
#     conda:
#         "1_QC_Assembly.yaml"
#     threads:
#         40
#     message:
#         "Downloading host reference genomes and building a Bowtie2 index"
#     shell:
#         """
#         # Download refs, combine to competitively map:
#         wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/497/805/GCF_900497805.2_bare-nosed_wombat_genome_assembly/GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.fna.gz \
#         wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz \
#         cat 1_References/Lasiorhinus_frefftii_mito.fna.gz GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.fna.gz GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz > {output.catted_ref} \
#
#         # Build bowtie2 index:
#         bowtie2-build \
#             --large-index \
#             --threads {threads} \
#             {output.catted_ref} {output.catted_ref}
#         """
################################################################################
### Index host genomes:
rule index_ref:
    input:
        "1_References/CattedRefs.fna.gz"
    output:
        bt2_index = "1_References/CattedRefs.fna.gz.rev.2.bt2l"
    conda:
        "1_QC_Assembly.yaml"
    threads:
        40
    log:
        "3_Outputs/0_Logs/host_genome_indexing.log"
    message:
        "Indexing host genomes with Bowtie2"
    shell:
        """
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {input} \
        &> {log}
        """
################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i = "2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2i = "2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz",
        catted_ref = "1_References/CattedRefs.fna.gz",
        bt2_index = "1_References/CattedRefs.fna.gz.rev.2.bt2l"
    output:
        all_bam = temp("3_Outputs/1_BAMs/{sample}.bam"),
        host_bam = "3_Outputs/1_BAMs/{sample}_host.bam",
        non_host_r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        non_host_r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz"
    conda:
        "1_QC_Assembly.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_mapping.log"
    message:
        "Mapping {wildcards.sample} reads to host genomes"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
        --time \
        --threads {threads} \
        -x {input.catted_ref} \
        -1 {input.r1i} \
        -2 {input.r2i} \
        | samtools view -@ {threads} -o {output.all_bam} - &&

        # Split extract non-host reads
        samtools view -b -f4 -@ {threads} {output.all_bam} \
        | samtools fastq -@ {threads} -c 6 -1 {output.non_host_r1} -2 {output.non_host_r2} - &&

        # Send host reads to BAM
        samtools view -b -F4 -@ {threads} {output.all_bam} \
        | samtools sort -@ {threads} -o {output.host_bam} - \
        &> {log}
        """
################################################################################
### Perform assembly on individual samples using metaspades:
rule assembly:
    input:
        non_host_r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        non_host_r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz"
    output:
        "3_Outputs/2_Assemblies/{sample}/{sample}_contigs.fasta"
    params:
        workdir = "3_Outputs/2_Assemblies/{sample}"
    conda:
        "1_QC_Assembly.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly.log"
    message:
        "Assembling {wildcards.sample} using metaspades"
    shell:
        """
        metaspades.py \
            -t {threads} \
            -k 21,33,55,77,99 \
            -1 {input.non_host_r1} -2 {input.non_host_r2} \
            -o {params.workdir} \
        &> {log} &&
        mv 3_Outputs/3_Assemblies/{wildcards.sample}/scaffolds.fasta {output}
        """
