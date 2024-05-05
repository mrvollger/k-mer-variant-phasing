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


rule haplotaged_bam:
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
        "benchmark/{sm}/haplotaged_bam/{sm}.bench.txt"
    shell:
        """
        python {params.script} \
            --sample-name {wildcards.sm} \
            -v -t {threads} \
            -i {input.bam} -r {input.tsv} \
            -o {output.bam}
        """


rule haplotaged_bai:
    input:
        bam=rules.haplotaged_bam.output.bam,
    output:
        bai=f"{rules.haplotaged_bam.output.bam}.bai",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """



rule haplotaged_vcf:
    input:
        bam=rules.haplotaged_bam.output.bam,
        vcf=rules.hiphase.output.vcf,
        ref=get_ref,
    output:
        vcf="results/{sm}/{sm}.haplotagged.vcf.gz",
    conda:
        CONDA
    threads: 4
    resources:
        mem_mb=48 * 1024,
    benchmark:
        "benchmark/{sm}/haplotaged_vcf/{sm}.bench.txt"
    shell:
        """
        whatshap haplotagphase \
            -r {input.ref} \
            -o {output.vcf} \
            {input.vcf} {input.bam}
        """