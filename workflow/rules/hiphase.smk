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


rule clean_vcf:
    input:
        vcf=VCF,
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


rule hiphase:
    input:
        vcf=get_vcf,
        tbi=get_tbi,
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=REFERENCE,
    output:
        vcf="results/{sm}/hiphase/{sm}.vcf.gz",
        #bam="results/{sm}/hiphase/{sm}.bam",
        summary="results/{sm}/hiphase/summary.tsv",
        stats="results/{sm}/hiphase/stats.tsv",
        blocks="results/{sm}/hiphase/blocks.tsv",
        haptag="results/{sm}/hiphase/read-level-phasing.tsv",
    conda:
        CONDA
    resources:
        mem_mb=32 * 1024,
    threads: 32
    shell:
        """
        hiphase -t {threads} \
            --ignore-read-groups \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --reference {input.ref} \
            --output-vcf {output.vcf} \
            --haplotag-file {output.haptag} \
            --summary-file {output.summary} \
            --stats-file {output.stats} \
            --blocks-file {output.blocks} 
        """
        #--output-bam {output.bam} \
        
'''
rule hiphase_read_list:
    input:
        bam=rules.hiphase.output.bam,
    output:
        tsv="results/{sm}/hiphase/read-phase-blocks.tsv.gz",
    conda:
        CONDA
    threads: 8
    params:
        script=workflow.source_path("../scripts/get-hap-and-phaseblock.py"),
    shell:
        """
        python {params.script} \
            -t {threads} \
            -i {input.bam} \
            | bgzip -@ {threads} \
            > {output.tsv}
        """
'''