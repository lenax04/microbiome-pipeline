# MicroSnake Pipeline Architecture

## Directed Acyclic Graph (DAG)

The following Mermaid diagram illustrates the complete Snakemake DAG for a typical MicroSnake run with 4 samples:

```mermaid
graph TD
    A1[Sample_1 R1/R2 FASTQ] --> B1(fastp: Sample_1)
    A2[Sample_2 R1/R2 FASTQ] --> B2(fastp: Sample_2)
    A3[Sample_3 R1/R2 FASTQ] --> B3(fastp: Sample_3)
    A4[Sample_4 R1/R2 FASTQ] --> B4(fastp: Sample_4)
    
    B1 --> C1(FastQC: Sample_1)
    B2 --> C2(FastQC: Sample_2)
    B3 --> C3(FastQC: Sample_3)
    B4 --> C4(FastQC: Sample_4)
    
    C1 --> D(MultiQC)
    C2 --> D
    C3 --> D
    C4 --> D
    
    B1 --> E(DADA2 ASV calling)
    B2 --> E
    B3 --> E
    B4 --> E
    
    E --> F(Taxonomic Classification: SILVA)
    F --> G(Alpha Diversity)
    F --> H(Beta Diversity PCoA)
    F --> I(Differential Abundance: DESeq2)
    
    B1 --> J1(MetaPhlAn4: Sample_1)
    B2 --> J2(MetaPhlAn4: Sample_2)
    B3 --> J3(MetaPhlAn4: Sample_3)
    B4 --> J4(MetaPhlAn4: Sample_4)
    
    J1 --> K(Merge MetaPhlAn4 Profiles)
    J2 --> K
    J3 --> K
    J4 --> K
    
    J1 --> L1(HUMAnN3: Sample_1)
    J2 --> L2(HUMAnN3: Sample_2)
    J3 --> L3(HUMAnN3: Sample_3)
    J4 --> L4(HUMAnN3: Sample_4)
    
    L1 --> M(Merge HUMAnN3 Pathways)
    L2 --> M
    L3 --> M
    L4 --> M
    
    D --> N(HTML Report)
    G --> N
    H --> N
    I --> N
    K --> N
    M --> N
```

## Rule Descriptions

| Rule | Input | Output | Tool |
| --- | --- | --- | --- |
| `fastp` | Raw FASTQ (R1, R2) | Trimmed FASTQ, QC JSON/HTML | fastp 0.23.4 |
| `fastqc` | Trimmed FASTQ | FastQC HTML/ZIP reports | FastQC 0.12.1 |
| `multiqc` | fastp JSON, FastQC ZIP | Aggregated HTML report | MultiQC 1.21 |
| `dada2` | Trimmed FASTQ (all samples) | ASV table, representative sequences | DADA2 1.30.0 (R) |
| `diversity_amplicon` | ASV table, metadata | Alpha diversity TSV, PCoA PNG, barplot PNG | vegan, ggplot2 (R) |
| `diff_abundance` | ASV table, metadata | DESeq2 results TSV | DESeq2 (R) |
| `metaphlan` | Trimmed FASTQ (per sample) | Taxonomic profile TXT | MetaPhlAn4 4.1.0 |
| `merge_metaphlan` | Per-sample profiles | Merged profile TSV | MetaPhlAn4 |
| `humann` | Trimmed FASTQ, MetaPhlAn profile | Pathway abundance TSV | HUMAnN3 3.8 |
| `merge_humann` | Per-sample pathway files | Merged pathway TSV | HUMAnN3 |
| `generate_html_report` | All results | Standalone HTML report | Python/Jinja2 |
