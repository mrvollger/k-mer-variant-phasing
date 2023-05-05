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
        CONDA
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
        unk=temp("temp/{sm}/canu_phase/haplotype/haplotype-unknown.fasta.gz"),
        outdir=temp(directory("temp/{sm}/canu_phase")),
    conda:
        CONDA
    resources:
        mem_mb=64 * 1024,
    params:
        genomeSize="3.1G",
    threads: 64
    shell:
        """
        canu -haplotype \
            maxThreads={threads} \
            useGrid=false \
            -p asm -d {output.outdir} \
            -genomeSize={params.genomeSize} \
            -haplotypemat {input.mat} \
            -haplotypepat {input.pat} \
            -pacbio-raw {input.fasta}
        """

rule canu_read_list:
    input:
        pat=rules.canu_phase.output.pat,
        mat=rules.canu_phase.output.mat,
        unk=rules.canu_phase.output.unk,
    output:
        txt="results/{sm}/canu/read-haplotypes.tsv",
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 8
    shell:
        """
        printf "read\\tkmer_hap\\n" > {output.txt}

        bgzip -cd@8 {input.pat} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tpat/' \
            >> {output.txt}

        bgzip -cd@8 {input.mat} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tmat/' \
            >> {output.txt}

        bgzip -cd@8 {input.unk} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tunk/' \
            >> {output.txt}
        """