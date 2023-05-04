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

rule canu_phase:
    input:
        fasta=rules.hifi_fasta.output.fasta,
    output:
        pat=temp("temp/{sm}/canu_phase/haplotype/haplotype-pat.fasta.gz"),
        mat=temp("temp/{sm}/canu_phase/haplotype/haplotype-mat.fasta.gz"),
        unk=temp("temp/{sm}/canu_phase/haplotype/haplotype-unk.fasta.gz"),
        odir=temp(directory("temp/{sm}/canu_phase")),
    conda:
        CONDA
    resources:
        mem_mb=64 * 1024,
    threads: 64
    shell:
        """
        canu -haplotype \
            maxThreads={threads} \
            useGrid=false \
            -p asm -d {output.dir} \
            -genomeSize={params.genomeSize} \
            -haplotypemat {input.mat} \
            -haplotypepat {input.pat} \
            -pacbio-raw {input.fasta}
        """


