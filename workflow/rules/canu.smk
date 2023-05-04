rule hifi_fasta:
    input:
        bam=HIFI_BAM,
    output:
        fasta=temp("temp/{sm}/hifi.fa.gz"),
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
        reads=get_reads,
    output:
        fasta=temp("temp/{sm}/{hap}.fa.gz"),
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


rule canu_phase:
    input:
        pat=get_pat,
        mat=get_mat,
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