# Quality Control Rules

rule fastp:
    input:
        r1 = lambda wildcards: get_fastq(wildcards, "1"),
        r2 = lambda wildcards: get_fastq(wildcards, "2")
    output:
        r1_clean = "results/qc/trimmed/{sample}_1.trimmed.fastq.gz",
        r2_clean = "results/qc/trimmed/{sample}_2.trimmed.fastq.gz",
        html = "results/qc/fastp/{sample}.html",
        json = "results/qc/fastp/{sample}.json"
    threads: config["qc"]["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -h {output.html} -j {output.json} \
              --qualified_quality_phred {config[qc][fastp][qualified_quality_phred]} \
              --unqualified_percent_limit {config[qc][fastp][unqualified_percent_limit]} \
              --length_required {config[qc][fastp][length_required]} \
              --thread {threads}
        """

rule fastqc:
    input:
        r1 = "results/qc/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = "results/qc/trimmed/{sample}_2.trimmed.fastq.gz"
    output:
        zip1 = "results/qc/fastqc/{sample}_1.trimmed_fastqc.zip",
        zip2 = "results/qc/fastqc/{sample}_2.trimmed_fastqc.zip",
        html1 = "results/qc/fastqc/{sample}_1.trimmed_fastqc.html",
        html2 = "results/qc/fastqc/{sample}_2.trimmed_fastqc.html"
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -o results/qc/fastqc/ {input.r1} {input.r2}
        """

rule multiqc:
    input:
        fastp_jsons = expand("results/qc/fastp/{sample}.json", sample=SAMPLES),
        fastqc_zips = expand("results/qc/fastqc/{sample}_1.trimmed_fastqc.zip", sample=SAMPLES)
    output:
        report = "results/qc/multiqc_report.html"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc results/qc/ -o results/qc/ -n multiqc_report.html --force
        """
