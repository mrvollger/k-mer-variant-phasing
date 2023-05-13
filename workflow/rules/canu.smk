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
        mem_mb=256 * 1024,
    params:
        genomeSize="3.1G",
        meryl_gb="256G",
        meryl_threads=32,
        grid=" --account=stergachislab --partition=compute --nodes=1 --export=all --parsable ",
    threads: 48
    shell:
        """
        OUTDIR="temp/{wildcards.sm}/canu_phase"
        canu -haplotype \
            maxThreads={threads} \
            useGrid=false \
            maxMemory={params.meryl_gb} \
            merylMemory={params.meryl_gb} \
            merylThreads={threads} \
            -p asm -d ${{OUTDIR}} \
            -genomeSize={params.genomeSize} \
            -haplotypemat {input.mat} \
            -haplotypepat {input.pat} \
            -pacbio-raw {input.fasta}
        # rm -rf ${{OUTDIR}}/canu-logs ${{OUTDIR}}/canu-scripts ${{OUTDIR}}/haplotype/*-kmers ${{OUTDIR}}/asm*
        """
        #useGrid=true \
        #gridOptions="{params.grid}" \


rule make_meryl:
    input:
        reads=get_meryl_input,
    output:
        meryl=directory("temp/{sm}/kmer_phase/{individual}/"),
        done=temp("temp/{sm}/kmer_phase/{individual}/done.txt"),
    conda:
        CONDA
    resources:
        mem_mb=70 * 1024,
    threads: 32
    shell:
        """
        which meryl
        rm -rf {output.meryl}
        meryl threads={threads} memory=64 k=21 count {input.reads} output {output.meryl}
        echo "done" > {output.done}
        """


rule hapmers:
    input:
        mat="temp/{sm}/kmer_phase/mat/",
        pat="temp/{sm}/kmer_phase/pat/",
        pro="temp/{sm}/kmer_phase/pro/",
    output:
        mat=directory("temp/{sm}/kmer_phase/pat.hapmer.meryl/"),
        pat=directory("temp/{sm}/kmer_phase/mat.hapmer.meryl/"),
        done=temp("temp/{sm}/kmer_phase/done.txt"),
    conda:
        CONDA
    resources:
        mem_mb=64 * 1024,
    threads: 32
    shell:
        """
        rm -rf {output.mat} {output.pat}
        MAT=$(realpath {input.mat})
        PAT=$(realpath {input.pat})
        PRO=$(realpath {input.pro})
        pushd {input.mat} && cd ..
        $MERQURY/trio/hapmers.sh $MAT $PAT $PRO
        popd
        echo "done" > {output.done}
        """


rule split_haplotype:
    input:
        pat=rules.hapmers.output.pat,
        mat=rules.hapmers.output.mat,
        fasta=rules.hifi_fasta.output.fasta,
    output:
        pat=temp("temp/{sm}/kmer_phase/haplotype/haplotype-pat.fasta.gz"),
        mat=temp("temp/{sm}/kmer_phase/haplotype/haplotype-mat.fasta.gz"),
        unk=temp("temp/{sm}/kmer_phase/haplotype/haplotype-unknown.fasta.gz"),
    conda:
        CONDA
    resources:
        mem_mb=32 * 1024,
    params:
        min_rl=1000,
    threads: 32
    shell:
        """
        splitHaplotype \
            -cl {params.min_rl} \
            -memory {resources.mem_mb} \
            -threads {threads} \
            -H {input.mat} 1 {output.mat} \
            -H {input.pat} 1 {output.pat} \
            -A {output.unk} \
            -R {input.fasta}
        """


rule kmer_read_list:
    input:
        pat=rules.split_haplotype.output.pat,
        mat=rules.split_haplotype.output.mat,
        unk=rules.split_haplotype.output.unk,
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
