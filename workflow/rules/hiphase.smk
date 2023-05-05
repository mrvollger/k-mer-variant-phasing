rule clean_bam:
    input:
        bam=HIFI_BAM,
    output:
        bam=pipe("temp/{sm}/{sm}.hiphase.bam"),
    conda:
        CONDA
    resources:
        mem_mb=32 * 1024,
    threads: 32
    shell:
        """
        reset-bam-read-groups.py \
            -t {threads} -r {wildcards.sm} \
            -i {input.bam} \
            -o {output.bam} 
        """

rule hiphase:
    input:
        vcf=VCF,
        bam=rules.clean_bam.output.bam,
        ref=REFERENCE,
    output:
        vcf="results/{sm}/{sm}.hiphase.vcf.gz",
        bam="results/{sm}/{sm}.hiphase.bam",
        summary="results/{sm}/{sm}.hiphase.summary.tsv",
        stats="results/{sm}/{sm}.hiphase.stats.tsv",
        blocks="results/{sm}/{sm}.hiphase.blocks.tsv",
    conda:
        CONDA
    resources:
        mem_mb=32 * 1024,
    threads: 32
    shell:
        """
        hiphase -t {threads} \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --reference {input.ref} \
            --output-bam {output.bam} \
            --output-vcf {output.vcf} \
            --summary-file {output.summary} \
            --stats-file {output.stats} \
            --blocks-file {output.blocks} 
        """
