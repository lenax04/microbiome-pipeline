# 16S rRNA Amplicon Rules

rule dada2:
    input:
        trimmed_r1 = expand("results/qc/trimmed/{sample}_1.trimmed.fastq.gz", sample=SAMPLES),
        trimmed_r2 = expand("results/qc/trimmed/{sample}_2.trimmed.fastq.gz", sample=SAMPLES),
        samples_tsv = config["samples_tsv"]
    output:
        asv_table = "results/amplicon/dada2/asv_table.tsv",
        rep_seqs = "results/amplicon/dada2/representative_seqs.fasta",
        stats = "results/amplicon/dada2/dada2_stats.tsv"
    threads: config["amplicon"]["threads"]
    conda:
        "../envs/dada2.yaml"
    log:
        "results/logs/dada2.log"
    shell:
        """
        Rscript scripts/dada2_pipeline.R \
            --samples {input.samples_tsv} \
            --input_dir results/qc/trimmed \
            --output_asv {output.asv_table} \
            --output_fasta {output.rep_seqs} \
            --output_stats {output.stats} \
            --truncLen1 {config[amplicon][dada2][truncLen_R1]} \
            --truncLen2 {config[amplicon][dada2][truncLen_R2]} \
            --maxEE1 {config[amplicon][dada2][maxEE_R1]} \
            --maxEE2 {config[amplicon][dada2][maxEE_R2]} \
            --threads {threads} > {log} 2>&1
        """
