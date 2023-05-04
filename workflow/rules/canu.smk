rule hifi_fasta:
    input:
        HIFI_BAM,
    output:
        fasta="temp/{sm}/hifi.fa.gz",
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 16
    shell:
        """
        samtools fasta -@ {threads} {input} | bgzip -@ {threads} > {output.fasta}
        """
