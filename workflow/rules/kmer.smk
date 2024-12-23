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
        if [[ {input.reads} =~ .*\\.(fofn) ]]; then
            cat $(cat {input.reads}) | seqtk seq -l 80 -A | bgzip -@ 8 > {output.fasta}
        elif [[ {input.reads} =~ .*\\.(bam|sam|cram) ]]; then
            samtools fasta -@ {threads} {input.reads} | bgzip -@ {threads} > {output.fasta}
        fi
        """


rule meryl:
    input:
        reads=get_meryl_input,
    output:
        meryl=directory("temp/{sm}/kmer_phase/{sm}.{individual}.meryl/"),
        done=temp("temp/{sm}/kmer_phase/{sm}.{individual}.done.txt"),
    conda:
        CONDA
    resources:
        mem_mb=148 * 1024,
    params:
        k_mer_size=K_MER_SIZE,
    threads: K_MER_THREADS
    benchmark:
        "benchmark/{sm}/meryl/{individual}.txt"
    shell:
        """
        which meryl
        rm -rf {output.meryl}
        meryl threads={threads} memory=124 k={params.k_mer_size} count {input.reads} output {output.meryl}
        echo "done" > {output.done}
        """


rule hapmers:
    input:
        mat="temp/{sm}/kmer_phase/{sm}.mat.meryl/",
        pat="temp/{sm}/kmer_phase/{sm}.pat.meryl/",
        pro="temp/{sm}/kmer_phase/{sm}.pro.meryl/",
    output:
        pat=directory("temp/{sm}/kmer_phase/hapmers/{sm}.pat.hapmer.meryl/"),
        mat=directory("temp/{sm}/kmer_phase/hapmers/{sm}.mat.hapmer.meryl/"),
        done=temp("temp/{sm}/kmer_phase/hapmers/done.txt"),
    conda:
        CONDA
    resources:
        mem_mb=64 * 1024,
    threads: K_MER_THREADS
    benchmark:
        "benchmark/{sm}/hapmers/bench.txt"
    shell:
        """
        # setup 
        RUNDIR=$(dirname {output.done})
        MAT=$(realpath {input.mat})/
        PAT=$(realpath {input.pat})/
        PRO=$(realpath {input.pro})/

        # clean
        rm -rf $RUNDIR
        mkdir -p $RUNDIR

        # hapmers
        pushd $RUNDIR
        $MERQURY/trio/hapmers.sh $MAT $PAT $PRO
        popd
        
        # finish
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
        mem_mb=64 * 1024,
    params:
        min_rl=1000,
    threads: K_MER_THREADS
    benchmark:
        "benchmark/{sm}/split_haplotype/bench.txt"
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
        tsv="results/{sm}/kmer/read-haplotypes.tsv",
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
