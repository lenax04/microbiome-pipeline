# Report Generation Rules

rule generate_html_report:
    input:
        multiqc = "results/qc/multiqc_report.html",
        alpha = "results/amplicon/diversity/alpha_diversity.tsv",
        pcoa = "results/amplicon/diversity/pcoa_bray_curtis.png",
        barplot = "results/amplicon/diversity/taxa_barplot.png",
        deseq = "results/amplicon/differential_abundance/deseq2_results.tsv"
    output:
        report = "results/report.html"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python3 scripts/generate_report.py \
            --multiqc {input.multiqc} \
            --alpha {input.alpha} \
            --pcoa {input.pcoa} \
            --barplot {input.barplot} \
            --deseq {input.deseq} \
            --output {output.report}
        """
