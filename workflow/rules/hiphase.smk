rule clean_bam:
    input:
        bam=HIFI_BAM,
    output:
        bam=temp("temp/{sm}/{sm}.hiphase.bam"),
        bai=temp("temp/{sm}/{sm}.hiphase.bam.bai"),
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 16
    params:
        script=workflow.source_path("../scripts/reset-bam-read-groups.py"),
    shell:
        """
        python {params.script} -v \
            -t {threads} \
            -r {wildcards.sm} \
            -s {wildcards.sm} \
            -i {input.bam} \
            -o {output.bam} 
        samtools index -@ {threads} {output.bam}
        """


rule hiphase:
    input:
        vcf=VCF,
        bam=get_hifi_bam,
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
