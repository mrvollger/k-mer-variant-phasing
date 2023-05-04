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

rule collect_reads:
    input:
        reads=lambda wc: config[wc.hap],
    output:
        fasta="{temp}/{pre}.{hap}.fasta",
    conda:
        "envs/env.yaml"
    threads: 16
    shell:
        """
        if [[ {input.reads} =~ .*\.(fofn) ]]; then
            cat $(cat {input.reads}) | seqtk seq -l 80 -A | bgzip -@ 8 > {output.fasta}
        elif [[ {input.reads} =~ .*\.(bam|sam|cram) ]]; then
            samtools fasta -@ {threads} {input.reads} | bgzip -@ {threads} > {output.fasta}
        fi
        """

def get_reads(hap):
    if hap == "mat":
        reads = MAT_DATA
    elif hap == "pat":
        reads = PAT_DATA
    else:
        raise ValueError("hap must be mat or pat")

    if reads.endswith(".fofn") or reads.endswith("bam") or reads.endswith("sam") or reads.endswith("cram"):
        return expand(rules.collect_reads.output.fasta, sm=sm, hap=hap)
    else:
        return reads
    

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


