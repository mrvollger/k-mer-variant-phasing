rule hiphase:
    input:
        vcf=VCF,
        bam=HIFI_BAM,
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
        mem_mb=64 * 1024,
    threads: 64
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
