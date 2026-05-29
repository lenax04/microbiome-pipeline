# Diversity and Statistical Analysis Rules

rule diversity_amplicon:
    input:
        asv_table = "results/amplicon/dada2/asv_table.tsv",
        samples_tsv = config["samples_tsv"]
    output:
        alpha = "results/amplicon/diversity/alpha_diversity.tsv",
        pcoa = "results/amplicon/diversity/pcoa_bray_curtis.png",
        barplot = "results/amplicon/diversity/taxa_barplot.png"
    conda:
        "../envs/diversity.yaml"
    log:
        "results/logs/diversity_amplicon.log"
    shell:
        """
        Rscript scripts/diversity_analysis.R \
            --asv {input.asv_table} \
            --metadata {input.samples_tsv} \
            --output_alpha {output.alpha} \
            --output_pcoa {output.pcoa} \
            --output_barplot {output.barplot} \
            --taxa_level {config[amplicon][diversity][taxa_level]} > {log} 2>&1
        """

rule diff_abundance:
    input:
        asv_table = "results/amplicon/dada2/asv_table.tsv",
        samples_tsv = config["samples_tsv"]
    output:
        results = "results/amplicon/differential_abundance/deseq2_results.tsv"
    conda:
        "../envs/diversity.yaml"
    log:
        "results/logs/diff_abundance.log"
    shell:
        """
        Rscript scripts/differential_abundance.R \
            --asv {input.asv_table} \
            --metadata {input.samples_tsv} \
            --output {output.results} \
            --formula "{config[amplicon][differential_abundance][formula]}" \
            --p_cutoff {config[amplicon][differential_abundance][p_cutoff]} \
            --log2fc_cutoff {config[amplicon][differential_abundance][log2fc_cutoff]} > {log} 2>&1
        """

rule diversity_shotgun:
    input:
        merged_profile = "results/shotgun/metaphlan/merged_taxonomic_profile.tsv",
        samples_tsv = config["samples_tsv"]
    output:
        barplot = "results/shotgun/diversity/shotgun_taxa_barplot.png"
    conda:
        "../envs/diversity.yaml"
    log:
        "results/logs/diversity_shotgun.log"
    shell:
        """
        Rscript scripts/diversity_analysis.R \
            --asv {input.merged_profile} \
            --metadata {input.samples_tsv} \
            --output_barplot {output.barplot} \
            --is_shotgun TRUE > {log} 2>&1
        """
