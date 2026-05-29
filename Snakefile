# MicroSnake: A reproducible Snakemake workflow for 16S rRNA amplicon and shotgun metagenomics analysis
# Author: Dawid X (dawidx1233)
# License: MIT

import os
import pandas as pd

# Load configuration
configfile: "config/config.yaml"

# Load samples
samples_df = pd.read_csv(config["samples_tsv"], sep="\t").set_index("sample_id", drop=False)
SAMPLES = list(samples_df.index)

# Helper function to get fastq files
def get_fastq(wildcards, group):
    return samples_df.loc[wildcards.sample, f"fastq_{group}"]

# Define target files based on run mode
def get_targets():
    targets = []
    # Quality Control targets (always run)
    targets.append("results/qc/multiqc_report.html")
    
    # 16S rRNA Amplicon targets
    if config["run_modes"]["amplicon"]:
        targets.append("results/amplicon/dada2/asv_table.tsv")
        targets.append("results/amplicon/dada2/representative_seqs.fasta")
        targets.append("results/amplicon/diversity/alpha_diversity.tsv")
        targets.append("results/amplicon/diversity/pcoa_bray_curtis.png")
        targets.append("results/amplicon/diversity/taxa_barplot.png")
        targets.append("results/amplicon/differential_abundance/deseq2_results.tsv")
        
    # Shotgun Metagenomics targets
    if config["run_modes"]["shotgun"]:
        targets.append("results/shotgun/metaphlan/merged_taxonomic_profile.tsv")
        targets.append("results/shotgun/humann/merged_pathway_abundance.tsv")
        targets.append("results/shotgun/diversity/shotgun_taxa_barplot.png")
        
    return targets

rule all:
    input:
        get_targets()

# Include sub-workflows/rules
include: "rules/qc.smk"
include: "rules/amplicon.smk"
include: "rules/shotgun.smk"
include: "rules/diversity.smk"
include: "rules/report.smk"
