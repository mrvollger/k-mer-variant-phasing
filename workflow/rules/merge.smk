rule merge_kmer_and_variant_phasing:
    input:
        kmer=rules.canu_read_list.output.tsv,
        variant=rules.hiphase_read_list.output.tsv,
    output:
        tsv="results/{sm}/combined_phasing.tsv.gz",
    conda:
        CONDA
    threads: 4
    resources:
        mem_mb=64*1024,
    params:
        script=workflow.source_path("../scripts/merge-kmer-and-variant-phasing.py"),
    shell:
        """
        python {params.script} \
            -v -t {threads} \
            {input.kmer} {input.variant} \
            | bgzip -@ {threads} \
            > {output.tsv}
        """


rule haplotaged_bam:
    input:
        bam=HIFI_BAM,
        tsv=rules.merge_kmer_and_variant_phasing.output.tsv,
    output:
        bam="results/{sm}/{sm}.haplotagged.bam",
        bai="results/{sm}/{sm}.haplotagged.bam.bai",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=32*1024,
    params:
        script=workflow.source_path("../scripts/merged-haplotag.py"),
    shell:
        """
        python {params.script} \
            -v -t {threads} \
            -i {input.bam} -r {input.tsv} \
            -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """
