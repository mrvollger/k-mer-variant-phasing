rule pbindex:
    input:
        bam=ancient(HIFI_BAM),
    output:
        pbi=f"{HIFI_BAM}.pbi",
    conda:
        CONDA
    threads: 8
    shell:
        """
        pbindex -j {threads} {input.bam} 
        """


rule zmws:
    input:
        bam=HIFI_BAM,
        pbi=f"{HIFI_BAM}.pbi",
    output:
        txt=temp("temp/{sm}/zmws/zmws.txt"),
    conda:
        CONDA
    threads: 8
    shell:
        """
        zmwfilter -j {threads} --show-all {input.bam} > {output.txt}
        """


rule split_zmws:
    input:
        txt=rules.zmws.output.txt,
    output:
        txt=temp("temp/{sm}/zmws/{scatteritem}.txt"),
    params:
        split_zmws=workflow.source_path("../scripts/split_zmws.py"),
    threads: 1
    conda:
        CONDA
    shell:
        """
        python {params.split_zmws} \
            --scatteritem {wildcards.scatteritem} \
            {input.txt} -o {output.txt}
        """


rule pbmm2_chunk:
    input:
        bam=HIFI_BAM,
        pbi=f"{HIFI_BAM}.pbi",
        ref=get_ref,
        fai=get_fai,
        zmws=rules.split_zmws.output.txt,
    output:
        tmp=temp("temp/{sm}/align/zmw.{scatteritem}.bam"),
        bam=temp("temp/{sm}/align/{scatteritem}.bam"),
        bai=temp("temp/{sm}/align/{scatteritem}.bam.bai"),
    conda:
        CONDA
    resources:
        disk_mb=32 * 1024,
        runtime=120,
        mem_mb=32 * 1024,
    threads: 8
    shell:
        """
        zmwfilter -j {threads} --include {input.zmws} {input.bam} {output.tmp}
        pbmm2 align \
            -j {threads} \
            --preset CCS --sort \
            --sort-memory 1G \
            --log-level INFO \
            --strip \
            --sample "{wildcards.sm}" \
            --unmapped \
            {input.ref} {output.tmp} {output.bam} 
        """


rule pbmm2_merge:
    input:
        bams=expand(
            rules.pbmm2_chunk.output.bam,
            scatteritem=custom_scatteritems(),
            allow_missing=True,
        ),
    output:
        #bam="results/{sm}/pbmm2/{sm}.cram",
        #bai="results/{sm}/pbmm2/{sm}.cram.crai",
        bam="results/{sm}/pbmm2/{sm}.bam",
        bai="results/{sm}/pbmm2/{sm}.bam.bai",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=64 * 1024,
    shell:
        """
        #pbmerge -j {threads} {input.bams} -o {output.bam}
        samtools merge -@ {threads} -cp --write-index -o {output.bam}##idx##{output.bai} {input.bams}
        """


# unused cram
"""
samtools merge \
    -@ {threads} \
    -O CRAM --output-fmt-option embed_ref=1 \
    --reference {input.ref} \
    -cp --write-index -o {output.bam} {input.bams}
"""


rule pbmm2_index:
    input:
        bam=rules.pbmm2_merge.output.bam,
    output:
        pbi="results/{sm}/pbmm2/{sm}.bam.pbi",
    conda:
        CONDA
    threads: 16
    resources:
        mem_mb=32 * 1024,
    shell:
        """
    pbindex -j {threads} {input.bam}
        """


rule pbmm2:
    input:
        expand(rules.pbmm2_merge.output, sm=SAMPLE),
        expand(rules.pbmm2_index.output, sm=SAMPLE),
