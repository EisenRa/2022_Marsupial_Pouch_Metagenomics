################################################################################
################################################################################
################################################################################
# Snakefile for MAG dereplication, taxonomic and functional annotation.
# Raphael Eisenhofer 02/2022
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

GROUP = [os.path.basename(dir)
         for dir in glob(f"3_Outputs/5_Refined_Bins/dRep_groups/*")]

MAGS = [os.path.relpath(fn, "3_Outputs/5_Refined_Bins/").replace(".fa.gz", "")
            for group in GROUP
            for fn in glob(f"3_Outputs/5_Refined_Bins/{group}/bins/*.fa.gz")]

print("Detected these sample groups:")
print(GROUP)
print("Detected this many MAGs:")
len(MAGS)
################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("3_Outputs/6_CoverM/{group}_assembly_coverM.txt", group=GROUP)
################################################################################
### Dereplicate refined bins using dRep
rule dereplication:
    input:
        bins = "3_Outputs/5_Refined_Bins/dRep_groups/{group}"
    output:
        drep = "3_Outputs/7_Dereplication/{group}/figures/{group}_Primary_clustering_dendrogram.pdf,
    params:
        ANI = expand("{ANI}", ANI=config['ANI']),
        workdir = "3_Outputs/7_Dereplication/{group}"
    conda:
        "3_dRep.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_dRep.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_dRep.log"
    message:
        "Dereplicating bins for {wildcards.group} that are > {params.ANI} percent indentical"
    shell:
        """
        # Parse/collate metawrap stats files for compatibility with dRep genomeinfo:
        for i in {input.bins}/

        # Run dRep
            dRep dereplicate \
                -p {threads} \
                -comp 70 \
                -sa {params.ANI} \
                -g {input.bins}/bins/*.fa.gz \
                --genomeInfo {input.bins}/genome_info.csv
                2> {log}

        # Rename and compress output
        for bin in {params.workdir}/dereplicated_genomes/*.fa; do
            mv $bin {params.workdir}/dereplicated_genomes/$(basename {group}_"$bin");
                done
        for i in {params.workdir}/figures/*; do
            mv $i {params.workdir}/figures/$(basename {group}_"$i");
                done

        pigz -p {threads} {paras.workdir}/dereplicated_genomes/*.fa
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
        directory("3_Outputs/3_Coassembly_Mapping/BAMs/{group}/Complete")
    params:
        outdir = directory("3_Outputs/3_Coassembly_Mapping/BAMs/{group}"),
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        read_dir = "2_Reads/3_Host_removed/{group}"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
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
        | samtools sort -@ {threads} -o {params.outdir}/$(basename ${{fq1/_1.fastq.gz/.bam}}); done

        #Create output file for snakemake
        mkdir -p {output}
        """
################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Coassembly_Mapping/BAMs/{group}/Complete"
    output:
        "3_Outputs/4_Binning/{group}/Done.txt"
    params:
        concoct = "3_Outputs/4_Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/4_Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/4_Binning/{group}/metabat2_bins",
        outdir = "3_Outputs/4_Binning/{group}",
        bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        memory = "180"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_binning.log"
    message:
        "Binning {wildcards.group} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir -p {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {params.bams}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {params.bams}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {params.bams}/*.bam; do ln -s `pwd`/$bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            -a {params.assembly} \
            -l 1500 \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq

        # Create dummy file for refinement input
        echo "Binning complete" > {output}
        """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        "3_Outputs/4_Binning/{group}/Done.txt"
    output:
        stats = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.contigs"
    params:
        concoct = "3_Outputs/4_Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/4_Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/4_Binning/{group}/metabat2_bins",
        outdir = "3_Outputs/5_Refined_Bins/{group}",
        memory = "180",
        group = "{group}"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        40
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
            -A {params.concoct} \
            -B {params.maxbin2} \
            -C {params.metabat2} \
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
        "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats"
    output:
        "3_Outputs/6_CoverM/{group}_assembly_coverM.txt"
    params:
        mapped_bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        binning_files = "3_Outputs/4_Binning/{group}",
        refinement_files = "3_Outputs/5_Refined_Bins/{group}",
        memory = "180",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
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

        # Clean up metaWRAP temp files
        rm -rf {params.binning_files}/work_files
        rm -f {params.binning_files}/*/*.fa
        rm -rf {params.refinement_files}/work_files
        pigz -p {threads} {params.refinement_files}/concoct_bins/*
        pigz -p {threads} {params.refinement_files}/metabat2_bins/*
        pigz -p {threads} {params.refinement_files}/maxbin2_bins/*

        """
