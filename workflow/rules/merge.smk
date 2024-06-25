rule merge_kmer_and_variant_phasing:
    input:
        kmer=rules.kmer_read_list.output.tsv,
        variant=rules.hiphase.output.haptag,
    output:
        tsv="results/{sm}/combined_phasing.tsv.gz",
    conda:
        CONDA
    threads: 4
    resources:
        mem_mb=64 * 1024,
    params:
        script=workflow.source_path("../scripts/merge-kmer-and-variant-phasing.py"),
    benchmark:
        "benchmark/{sm}/merge_kmer_and_variant_phasing/bench.txt"
    log:
        "results/{sm}/merge-kmer-and-variant-phasing.log",
    shell:
        """
        python {params.script} \
            -v -t {threads} \
            {input.kmer} {input.variant} \
            | bgzip -@ {threads} \
            > {output.tsv} 2> {log}
        """


rule haplotagged_bam:
    input:
        bam=get_hifi_bam,
        tsv=rules.merge_kmer_and_variant_phasing.output.tsv,
    output:
        bam="results/{sm}/{sm}.haplotagged.bam",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=32 * 1024,
    params:
        script=workflow.source_path("../scripts/merged-haplotag.py"),
    benchmark:
        "benchmark/{sm}/haplotagged_bam/{sm}.bench.txt"
    shell:
        """
        python {params.script} \
            --sample-name {wildcards.sm} \
            -v -t {threads} \
            -i {input.bam} -r {input.tsv} \
            -o {output.bam}
        """


rule haplotagged_bai:
    input:
        bam=rules.haplotagged_bam.output.bam,
    output:
        bai=f"{rules.haplotagged_bam.output.bam}.bai",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """


rule haplotagged_vcf:
    input:
        bam=rules.haplotagged_bam.output.bam,
        bai=rules.haplotagged_bai.output.bai,
        vcf=rules.hiphase_vcf.output.vcf,
        ref=get_ref,
    output:
        vcf="results/{sm}/{sm}.haplotagged.vcf.gz",
    conda:
        CONDA
    threads: 4
    resources:
        mem_mb=32 * 1024,
        runtime=8 * 60,
    benchmark:
        "benchmark/{sm}/haplotagged_vcf/{sm}.bench.txt"
    shell:
        """
        whatshap haplotagphase \
            -r {input.ref} \
            -o {output.vcf} \
            {input.vcf} {input.bam}
        """


rule haplotagged_vcf_index:
    input:
        vcf=rules.haplotagged_vcf.output.vcf,
    output:
        tbi=f"{rules.haplotagged_vcf.output.vcf}.tbi",
    conda:
        CONDA
    threads: 1
    resources:
        mem_mb=12 * 1024,
    shell:
        """
        tabix {input.vcf}
        """


rule haplotagged_gvcf:
    input:
        vcf=rules.haplotagged_vcf.output.vcf,
        gvcf=rules.hiphase_gvcf.output.gvcf,
    output:
        gvcf="results/{sm}/{sm}.haplotagged.gvcf.gz",
        tbi="results/{sm}/{sm}.haplotagged.gvcf.gz.tbi",
    conda:
        CONDA
    threads: 1
    resources:
        mem_mb=32 * 1024,
    params:
        script=workflow.source_path("../scripts/combine-phased-vcf-and-gvcf.py"),
    shell:
        """
        python {params.script} {input.vcf} {input.gvcf} {output.gvcf}
        tabix {output.gvcf}
        """
