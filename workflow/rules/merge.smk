
rule merge_kmer_and_variant_phasing:
    input:
        kmer=rules.canu_read_list.output.tsv,
        variant=rules.hiphase_read_list.output.tsv,
    output:
        tsv="results/{sm}/combined_phasing.tsv.gz",
    conda:
        CONDA
    threads: 8
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