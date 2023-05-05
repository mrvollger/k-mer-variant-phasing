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

rule hiphase_read_lists:
    input:
        bam=rules.hiphase.output.bam
    output:
        "results/{sm}/hiphase/{tag}.hiphase.reads.tsv"
    conda:
        CONDA
    threads: 8 
    shell: 
        """
        ( \
            printf "read\tphase_block\thap\n"; \
            samtools view -@ {threads} -d HP:{wildcards.tag}  {input.bam} \
                | cut -f1,12- \
                | sed -E 's/^(\w+).*PS:i:(\w+)\tHP:i:(\w+).*/\1\t\2\t\3/' \
        ) | bgzip -@ {threads} > {output}
        """