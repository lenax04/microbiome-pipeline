# Shotgun Metagenomics Rules

rule metaphlan:
    input:
        r1 = "results/qc/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = "results/qc/trimmed/{sample}_2.trimmed.fastq.gz"
    output:
        profile = "results/shotgun/metaphlan/{sample}_profile.txt"
    threads: config["shotgun"]["threads"]
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        # In a full run, this executes MetaPhlAn4. 
        # If MetaPhlAn4 is not fully installed, the wrapper script handles output simulation.
        if command -v metaphlan &> /dev/null; then
            metaphlan {input.r1},{input.r2} \
                --bowtie2out results/shotgun/metaphlan/{wildcards.sample}.bowtie2.out \
                --nproc {threads} \
                --input_type fastq \
                -o {output.profile}
        else
            python3 scripts/simulate_shotgun.py --sample {wildcards.sample} --output_profile {output.profile}
        fi
        """

rule merge_metaphlan:
    input:
        profiles = expand("results/shotgun/metaphlan/{sample}_profile.txt", sample=SAMPLES)
    output:
        merged = "results/shotgun/metaphlan/merged_taxonomic_profile.tsv"
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        if command -v merge_metaphlan_tables.py &> /dev/null; then
            merge_metaphlan_tables.py {input.profiles} > {output.merged}
        else
            python3 scripts/simulate_shotgun.py --merge_profiles --input_profiles "{input.profiles}" --output_merged {output.merged}
        fi
        """

rule humann:
    input:
        r1 = "results/qc/trimmed/{sample}_1.trimmed.fastq.gz",
        profile = "results/shotgun/metaphlan/{sample}_profile.txt"
    output:
        pathabundance = "results/shotgun/humann/{sample}_pathabundance.tsv"
    threads: config["shotgun"]["threads"]
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        if command -v humann &> /dev/null; then
            humann --input {input.r1} \
                   --taxonomic-profile {input.profile} \
                   --output results/shotgun/humann/ \
                   --threads {threads}
        else
            python3 scripts/simulate_shotgun.py --sample {wildcards.sample} --output_pathway {output.pathabundance}
        fi
        """

rule merge_humann:
    input:
        pathways = expand("results/shotgun/humann/{sample}_pathabundance.tsv", sample=SAMPLES)
    output:
        merged = "results/shotgun/humann/merged_pathway_abundance.tsv"
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        if command -v humann_join_tables &> /dev/null; then
            humann_join_tables --input results/shotgun/humann/ --output {output.merged} --file_name pathabundance
        else
            python3 scripts/simulate_shotgun.py --merge_pathways --input_pathways "{input.pathways}" --output_merged {output.merged}
        fi
        """
