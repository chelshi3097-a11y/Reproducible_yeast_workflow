import pandas as pd
import os

# Load sample metadata
samples = pd.read_csv("metadata/samples.csv")
SAMPLES = samples["Run_Accession"].tolist()
sample_links = dict(zip(samples["Run_Accession"], samples["Download_link"]))

# Ensure shell uses bash
shell.executable("/bin/bash")

# -----------------------------
# Rule: all
# -----------------------------
rule all:
    input:
        expand("data/raw_fastq/{run}.fastq.gz", run=SAMPLES),
        expand("data/processed_fastq/{run}_fastqc.html", run=SAMPLES),
        "data/processed_fastq/multiqc_report.html",
        "data/reference/STAR_index",
        expand("data/star_aligned/{run}.Aligned.sortedByCoord.out.bam", run=SAMPLES),
        "data/featurecounts/featureCounts_gene_counts.txt",
        "data/deseq2_results/deseq2_results.csv",
        "data/deseq2_results/MA_plot.png",
        "data/deseq2_results/Volcano_plot.png",
        "data/deseq2_results/PCA_plot.png",
        "data/deseq2_results/Top20_DE_genes_heatmap.png",
        "data/deseq2_results/Sample_distance_heatmap.png"

# -----------------------------
# Rule: download raw FASTQ
# -----------------------------
rule download_fastq:
    output:
        "data/raw_fastq/{run}.fastq.gz"
    params:
        link=lambda wildcards: sample_links[wildcards.run]
    shell:
        """
        mkdir -p data/raw_fastq
        wget -c {params.link} -O {output}
        """

# -----------------------------
# Rule: FastQC
# -----------------------------
rule fastqc:
    input:
        "data/raw_fastq/{run}.fastq.gz"
    output:
        "data/processed_fastq/{run}_fastqc.html"
    threads: 2
    shell:
        """
        mkdir -p data/processed_fastq
        fastqc -t {threads} {input} -o data/processed_fastq
        """

# -----------------------------
# Rule: MultiQC
# -----------------------------
rule multiqc:
    input:
        expand("data/processed_fastq/{run}_fastqc.html", run=SAMPLES)
    output:
        "data/processed_fastq/multiqc_report.html"
    shell:
        """
        mkdir -p data/processed_fastq
        multiqc data/processed_fastq -o data/processed_fastq
        """

# -----------------------------
# Rule: download genome FASTA
# -----------------------------
rule download_genome:
    output:
        "data/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
    shell:
        """
        mkdir -p data/reference
        wget -c ftp://ftp.ensemblgenomes.org/pub/fungi/release-56/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O {output}.gz
        gunzip -c {output}.gz > {output}
        """

# -----------------------------
# Rule: download GTF
# -----------------------------
rule download_gtf:
    output:
        "data/reference/Saccharomyces_cerevisiae.R64-1-1.56.gtf"
    shell:
        """
        mkdir -p data/reference
        wget -c ftp://ftp.ensemblgenomes.org/pub/fungi/release-56/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.56.gtf.gz -O {output}.gz
        gunzip -c {output}.gz > {output}
        """

# -----------------------------
# Rule: STAR genome index with dynamic overhang
# -----------------------------
rule star_index:
    input:
        genome="data/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
        gtf="data/reference/Saccharomyces_cerevisiae.R64-1-1.56.gtf",
        fastq=lambda wildcards: f"data/raw_fastq/{SAMPLES[0]}.fastq.gz"  # use first sample for read length
    output:
        directory("data/reference/STAR_index")
    threads: 4
    run:
        import gzip

        # Calculate read length from first FASTQ
        fq_file = input.fastq
        with gzip.open(fq_file, "rt") as fh:
            next(fh)  # skip header
            seq = next(fh).strip()
            read_length = len(seq)

        sjdbOverhang = read_length - 1
        print(f"Detected read length = {read_length}, using sjdbOverhang = {sjdbOverhang}")

        # Make output directory if it doesn't exist
        os.makedirs(output[0], exist_ok=True)

        # Run STAR
        shell(f"""
            STAR --runThreadN {threads} \
                 --runMode genomeGenerate \
                 --genomeDir {output[0]} \
                 --genomeFastaFiles {input.genome} \
                 --sjdbGTFfile {input.gtf} \
                 --sjdbOverhang {sjdbOverhang}
        """)


# -----------------------------
# Rule: STAR alignment
# -----------------------------
rule star_align:
    input:
        fastq="data/raw_fastq/{run}.fastq.gz",
        index="data/reference/STAR_index"
    output:
        bam="data/star_aligned/{run}.Aligned.sortedByCoord.out.bam"
    threads: 4
    shell:
        """
        mkdir -p data/star_aligned
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.fastq} \
             --readFilesCommand zcat \
             --outFileNamePrefix data/star_aligned/{wildcards.run}. \
             --outSAMtype BAM SortedByCoordinate
        """

# -----------------------------
# Rule: featureCounts
# -----------------------------
rule featurecounts:
    input:
        bam=lambda wildcards: expand("data/star_aligned/{run}.Aligned.sortedByCoord.out.bam", run=SAMPLES),
        gtf="data/reference/Saccharomyces_cerevisiae.R64-1-1.56.gtf"
    output:
        counts="data/featurecounts/featureCounts_gene_counts.txt"
    threads: 4
    shell:
        """
        mkdir -p data/featurecounts
        featureCounts -T {threads} -s 1 \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam}
        """

# -----------------------------
# Rule: DESeq2 analysis
# -----------------------------
rule deseq2:
    input:
        counts="data/featurecounts/featureCounts_gene_counts.txt",
        metadata="metadata/samples.csv"
    output:
        results="data/deseq2_results/deseq2_results.csv",
        ma="data/deseq2_results/MA_plot.png",
        volcano="data/deseq2_results/Volcano_plot.png",
        pca="data/deseq2_results/PCA_plot.png",
        heatmap="data/deseq2_results/Top20_DE_genes_heatmap.png",
        sample_dist="data/deseq2_results/Sample_distance_heatmap.png"
    log:
        "logs/deseq2.log"
    shell:
        """
        Rscript --vanilla scripts/deseq2_yeast.R {input.counts} {input.metadata} \
        {output.results} {output.ma} {output.volcano} {output.pca} \
        {output.heatmap} {output.sample_dist} > {log} 2>&1
        """

