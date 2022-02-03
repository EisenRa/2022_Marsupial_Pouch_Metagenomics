################################################################################
################################################################################
################################################################################
# Snakefile for coCoassembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

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
        expand("3_Outputs/2_Coassemblies/{group}_QUAST", group=GROUP),
        expand("3_Outputs/3_Coassembly_Mapping/BAMs/{group}_coverM.txt", group=GROUP)

################################################################################
### Perform Coassemblies on each sample group
rule Coassembly:
    input:
        reads = "2_Reads/3_Host_removed/{group}"
    output:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
    params:
        workdir = "3_Outputs/2_Coassemblies/{group}",
        r1_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_1.fastq.gz"),
        r2_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_2.fastq.gz"),
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly.log"
    message:
        "Coassembling {wildcards.group} using {params.assembler}"
    shell:
        """
        # Set up assembler variable from config file
        export assembler={config[assembler]}

        if [ "$assembler" == "metaspades" ]
        then

        # Concatenate reads from the same group for Coassembly
        cat {input.reads}/*_1.fastq.gz > {params.r1_cat}
        cat {input.reads}/*_2.fastq.gz > {params.r2_cat}

        # Run metaspades
            metaspades.py \
                -t {threads} \
                -k 21,33,55,77,99 \
                -1 {params.r1_cat} -2 {params.r2_cat} \
                -o {params.workdir}
                2> {log}

        # Remove contigs shorter than 1,500 bp
            reformat.sh
                in={params.workdir}/scaffolds.fasta \
                out={output.Coassembly} \
                minlength=1500

        else

        # Set up input reads variable for megahit
        R1=$(for i in 2_Reads/3_Host_removed/{wildcards.group}/*_1.fastq.gz; do echo $i | tr '\n' ,; done)
        R2=$(for i in 2_Reads/3_Host_removed/{wildcards.group}/*_2.fastq.gz; do echo $i | tr '\n' ,; done)

        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 $R1 -2 $R2 \
                -f \
                -o {params.workdir}
                2> {log}

        # Move the Coassembly to final destination
            mv {params.workdir}/final.contigs.fa {output.Coassembly}

        # Reformat headers
            sed -i 's/ /-/g' {output.Coassembly}

        fi
        """
################################################################################
### Create QUAST reports of coassemblies
rule QUAST:
    input:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        report = directory("3_Outputs/2_Coassemblies/{group}_QUAST"),
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        20
    message:
        "Running -QUAST on {wildcards.group} coassembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.Coassembly}

        # Rename QUAST files
        for i in {output.report}/*;
            do mv $i {output.report}/{wildcards.group}_$(basename $i);
                done
        """
################################################################################
### Map reads to the coassemblies
rule Coassembly_index:
    input:
        report = "3_Outputs/2_Coassemblies/{group}_QUAST/"
    output:
        bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l",
    params:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_indexing.log"
    message:
        "Indexing {wildcards.group} coassembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.Coassembly} {params.Coassembly} \
        &> {log}
        """
################################################################################
### Map reads to the coassemblies
rule Coassembly_mapping:
    input:
        bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l"
    output:
        mapped_bam = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}"
    params:
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        read_dir = "2_Reads/3_Host_removed/{group}"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_mapping.log"
    message:
        "Mapping {wildcards.group} samples to coassembly using Bowtie2"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        for fq1 in {params.read_dir}/*_1.fastq.gz; do \
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.assembly} \
            -1 $fq1 \
            -2 ${{fq1/_1.fastq.gz/_2.fastq.gz}} \
        | samtools view -@ {threads} -o {output.mapped_bam}/${{fq1/_1.fastq.gz/.bam}} -; done
        """
################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Coassembly_Mapping/BAMs/{group}"
    output:
        concoct = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/metabat2_bins",
    params:
        outdir = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/Binning",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        memory = "16"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_binning.log"
    message:
        "Binning {wildcards.group} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input}/*.bam; do ln -s $bam {params.outdir}/work_files/$(basename $bam); done

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
        concoct = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/metabat2_bins",
    output:
        stats = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.contigs"
    params:
        outdir = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}",
        memory = "16",
        group = "{group}"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Refining {wildcards.group} bins with MetaWRAP's bin refinement module"
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
        sed -i'' '2,$s/bin/bin_{params.group}/g' {output.stats}
        sed -i'' 's/bin/bin_{params.group}/g' {output.contigmap}
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats"
    output:
        "3_Outputs/3_Coassembly_Mapping/BAMs/{group}_coverM.txt"
    params:
        mapped_bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        memory = "16",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Calculating coassembly mapping rate for {wildcards.group} with CoverM"
    shell:
        """
        coverm genome \
            -b {params.mapped_bams}/*.bam \
            --genome-fasta-files {params.assembly} \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
