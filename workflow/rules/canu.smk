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
    conda:
        CONDA
    resources:
        mem_mb=24 * 1024,
    params:
        genomeSize="3.1G",
        meryl_gb="24G",
        meryl_threads=8,
        grid=" --account=stergachislab --partition=compute --nodes=1 --export=all --parsable "
    threads: 8
    shell:
        """
        OUTDIR="temp/{wildcards.sm}/canu_phase"
        canu -haplotype \
            merylMemory={params.meryl_gb} \
            merylThreads={params.meryl_threads} \
            useGrid=true \
            gridOptions="{params.grid}" \
            -p asm -d ${{OUTDIR}} \
            -genomeSize={params.genomeSize} \
            -haplotypemat {input.mat} \
            -haplotypepat {input.pat} \
            -pacbio-raw {input.fasta}
        rm -rf ${{OUTDIR}}/canu-logs
        rm -rf ${{OUTDIR}}/canu-scripts
        rm -rf ${{OUTDIR}}/haplotype/*-kmers
        rm -rf ${{OUTDIR}}/asm*
        """
        #maxThreads={threads} \
        #useGrid=false \
        #maxMemory={params.meryl_gb} 


rule canu_read_list:
    input:
        pat=rules.canu_phase.output.pat,
        mat=rules.canu_phase.output.mat,
        unk=rules.canu_phase.output.unk,
    output:
        tsv="results/{sm}/canu/read-haplotypes.tsv",
    conda:
        CONDA
    resources:
        mem_mb=8 * 1024,
    threads: 8
    shell:
        """
        printf "read\\tkmer_hap\\n" > {output.tsv}

        bgzip -cd@8 {input.pat} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tpat/' \
            >> {output.tsv}

        bgzip -cd@8 {input.mat} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tmat/' \
            >> {output.tsv}

        bgzip -cd@8 {input.unk} \
            | grep '^>' | cut -f 1 | sed 's/^>//' | sed 's/$/\\tunk/' \
            >> {output.tsv}
        """
