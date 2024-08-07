rule clean_bam:
    input:
        bam=get_hifi_bam,
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


rule clean_vcf:
    input:
        vcf=get_input_vcf,
    output:
        vcf=temp("temp/{sm}/{sm}.hiphase.vcf.gz"),
        tbi=temp("temp/{sm}/{sm}.hiphase.vcf.gz.tbi"),
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 16
    shell:
        """
        bcftools reheader \
            --threads {threads} \
            -s <(echo {wildcards.sm}) \
            -o {output.vcf} {input.vcf} 
        tabix {output.vcf}
        """


# Currently, DeepVariant, pbsv, and TRGT are the three supported input types.
rule hiphase:
    input:
        gvcf=get_variant_file,
        tbi=get_variant_tbi,
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=REFERENCE,
        pbsv_vcf=rules.pbsv_index.output.vcf,
    output:
        gvcf="results/{sm}/hiphase/{sm}.gvcf.gz",
        pbsv_vcf="results/{sm}/hiphase/{sm}.pbsv.vcf.gz",
        summary="results/{sm}/hiphase/summary.tsv",
        stats="results/{sm}/hiphase/stats.tsv",
        blocks="results/{sm}/hiphase/blocks.tsv",
        haptag="results/{sm}/hiphase/read-level-phasing.tsv",
        bam="results/{sm}/hiphase/{sm}.bam" if NO_PARENTAL else [],
    conda:
        CONDA
    resources:
        mem_mb=96 * 1024,
    threads: 16
    params:
        bam="--output-bam" if NO_PARENTAL else "",
    benchmark:
        "benchmark/{sm}/hiphase/bench.txt"
    shell:
        """
        hiphase -t {threads} \
            --ignore-read-groups \
            --bam {input.bam} \
            --reference {input.ref} \
            --vcf {input.gvcf} \
            --output-vcf {output.gvcf} \
            --vcf {input.pbsv_vcf} \
            --output-vcf {output.pbsv_vcf} \
            {params.bam} {output.bam} \
            --haplotag-file {output.haptag} \
            --summary-file {output.summary} \
            --stats-file {output.stats} \
            --blocks-file {output.blocks} 
        """


rule hiphase_vcf:
    input:
        gvcf=rules.hiphase.output.gvcf,
    output:
        vcf="results/{sm}/hiphase/{sm}.vcf.gz",
        index="results/{sm}/hiphase/{sm}.vcf.gz.csi",
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 4
    shell:
        """
        bcftools view \
            --threads {threads} \
            -e'FILTER="."' -AA \
            --write-index \
            -o {output.vcf} \
            {input.gvcf}
        """
