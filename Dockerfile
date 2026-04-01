FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# -------------------------
# System dependencies + Python for Snakemake
# -------------------------
RUN apt-get update && apt-get install -y \
    ca-certificates \
    curl \
    wget \
    git \
    build-essential \
    software-properties-common \
    dirmngr \
    gnupg \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3 \
    python3-dev \
    python3-pip \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# -------------------------
# Install R (CRAN)
# -------------------------
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" && \
    apt-get update && \
    apt-get install -y r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# -------------------------
# Install Bioconductor + RNA-seq R packages
# -------------------------
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(ask=FALSE, update=FALSE)" && \
    Rscript -e "BiocManager::install(c('DESeq2','edgeR','tximport','GenomicFeatures','GenomicRanges','SummarizedExperiment','AnnotationDbi','limma'))"

# -------------------------
# Install additional CRAN packages for DESeq2 plots
# -------------------------
RUN Rscript -e "install.packages(c('pheatmap','ggrepel'), repos='https://cloud.r-project.org')"

# -------------------------
# Upgrade pip tooling
# -------------------------
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel

# -------------------------
# Install Snakemake + Python deps
# -------------------------
RUN pip3 install --no-cache-dir \
    snakemake==7.32.4 \
    pulp==2.7.0 \
    pandas \
    multiqc

# -------------------------
# Install RNA-seq CLI tools
# -------------------------
RUN apt-get update && apt-get install -y \
    fastqc \
    subread \
    rna-star \
    && rm -rf /var/lib/apt/lists/*

# -------------------------
# Set working directory
# -------------------------
WORKDIR /work

CMD ["/bin/bash"]
