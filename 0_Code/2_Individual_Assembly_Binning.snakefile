################################################################################
################################################################################
################################################################################
# Snakefile for individual assembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

ASSEMBLER = 'megahit'

GROUP = [os.path.basename(dir)
         for dir in glob(f"2_Reads/3_Host_removed/*")]

SAMPLE = [os.path.relpath(fn, "2_Reads/3_Host_removed").replace("_non_host_1.fastq.gz", "")
            for group in GROUP
            for fn in glob(f"2_Reads/3_Host_removed/{group}/*_1.fastq.gz")]

print("Detected these sample groups:")
print(GROUP)
print("Detected the following samples:")
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
#        expand("3_Outputs/3_Assembly_Mapping/BAMs/{sample}_coverM.txt", sample=SAMPLE)
        expand("3_Outputs/2_Assemblies/{sample}_QUAST/report.html", sample=SAMPLE)

################################################################################
### Perform assembly on each sample
rule Assembly:
    input:
        r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz",
    output:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
    params:
        workdir = "3_Outputs/2_Assemblies/{sample}",
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly.log"
    message:
        "Assembling {wildcards.sample} using {params.assembler}"
    shell:
        """
        # Set up assembler variable from config file
        export assembler={config[assembler]}

        if [ "$assembler" == "metaspades" ]
        then
        # Run metaspades
            metaspades.py \
                -t {threads} \
                -k 21,33,55,77,99 \
                --only-assembler \
                -1 {input.r1} -2 {input.r2} \
                -o {params.workdir}
                2> {log}

        # Remove contigs shorter than 1,500 bp
            reformat.sh
                in={params.workdir}/scaffolds.fasta \
                out={output.assembly} \
                minlength=1500

        else
        # Run megahit
            metaspades.py \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.workdir}
                2> {log}

        # Move the Coassembly to final destination
            mv {params.workdir}/scaffolds.fasta {output.assembly}

        # Reformat headers
            sed -i 's/ /-/g' {output.assembly}

        # Move the Coassembly to final destination
            mv {params.workdir}/scaffolds.fasta {output.assembly}
        fi
        """
################################################################################
### Create QUAST reports of assemblies
rule QUAST:
    input:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    output:
        report = "3_Outputs/2_Assemblies/{sample}_QUAST/report.html",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    message:
        "Running QUAST on {wildcards.sample} assembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.assembly}
        """
################################################################################
### Index each sample's assembly
rule assembly_index:
    input:
        report = "3_Outputs/2_Assemblies/{sample}_QUAST/report.html"
    output:
        bt2_index = "3_Outputs/2_Assemblies/{sample}_contigs.fasta.rev.2.bt2l",
    params:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.log"
    message:
        "Indexing {wildcards.sample} assembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.assembly} {params.assembly} \
        &> {log}
        """
################################################################################
### Map a sample's reads to it corresponding assembly
rule assembly_mapping:
    input:
        bt2_index = "3_Outputs/2_Assemblies/{sample}_contigs.fasta.rev.2.bt2l"
    output:
        mapped_bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    params:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        r1 = "2_Reads/3_Host_removed/{sample}_1.fastq.gz",
        r2 = "2_Reads/3_Host_removed/{sample}_2.fastq.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.log"
    message:
        "Mapping {wildcards.sample} to its assembly using Bowtie2"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.assembly} \
            -1 {params.r1} \
            -2 {params.r2} \
        | samtools view -@ {threads} -o {output.mapped_bam} -
        """
################################################################################
### Bin each sample's contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        concoct = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/concoct_bins",
        maxbin2 = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/metabat2_bins",
    params:
        outdir = "3_Outputs/3_Assembly_Mapping/Binning/{sample}",
        assembly = "3_Outputs/3_Assemblies/{sample}_contigs.fasta",
        basename = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}",
        memory = "180"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_binning.log"
    message:
        "Binning {wildcards.sample} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_1.fastq)
        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_2.fastq)

        #Symlink BAMs for metaWRAP
        ln -s {input} {params.outdir}/work_files/$(basename {input}); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            -a {params.assembly} \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq
        """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        concoct = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/concoct_bins",
        maxbin2 = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Assembly_Mapping/Binning/{sample}/metabat2_bins",
    output:
        stats = "3_Outputs/3_Assembly_Mapping/Refined_Bins/{sample}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/3_Assembly_Mapping/Refined_Bins/{sample}_metawrap_70_10_bins.contigs"
    params:
        outdir = "3_Outputs/3_Assembly_Mapping/Refined_Bins/{sample}",
        memory = "180",
        sample = "{sample}"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.log"
    message:
        "Refining {wildcards.sample} bins with MetaWRAP's bin refinement module"
    shell:
        """
        metawrap bin_refinement \
            -m {params.memory} \
            -t {threads} \
            -o {params.outdir} \
            -A {input.concoct} \
            -B {input.maxbin2} \
            -C {input.metabat2} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/bin_{params.sample}/g' {output.stats}
        sed -i'' 's/bin/bin_{params.sample}/g' {output.contigmap}
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        stats = "3_Outputs/3_Assembly_Mapping/Refined_Bins/{sample}_metawrap_70_10_bins.stats",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        mapped_bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        "3_Outputs/3_Assembly_Mapping/BAMs/{sample}_coverM.txt"
    params:
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_coverM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_coverM.log"
    message:
        "Calculating assembly mapping rate for {wildcards.sample} with CoverM"
    shell:
        """
        coverm genome \
            -b {input.mapped_bam} \
            --genome-fasta-files {input.assembly} \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
