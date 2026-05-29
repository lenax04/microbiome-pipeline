FROM ubuntu:22.04

LABEL maintainer="Dawid X <dawidx1233>"
LABEL description="Dockerfile for MicroSnake: reproducible microbiome pipeline"

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    r-base-core \
    wget \
    curl \
    git \
    fastqc \
    fastp \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    snakemake \
    multiqc \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    jinja2

# Install R packages
RUN R -e "install.packages(c('ggplot2', 'vegan'), repos='https://cloud.r-project.org/')"

WORKDIR /app
COPY . /app

CMD ["snakemake", "--cores", "2"]
