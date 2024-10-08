rule discover_pbsv:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
    output:
        svsig="results/{sm}/pbsv/{sm}.pbsv.svsig.gz",
    threads: 16
    resources:
        mem_mb=64 * 1024,
        runtime=16 * 60,
    conda:
        CONDA
    shell:
        """
        pbsv discover \
            --ccs \
            --sample {wildcards.sm} \
            {input.bam} \
            {output.svsig}
        """


rule run_pbsv:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
        svsig=rules.discover_pbsv.output.svsig,
    output:
        vcf=temp("temp/{sm}/pbsv/{sm}.pbsv.vcf"),
    threads: 16
    resources:
        mem_mb=64 * 1024,
        runtime=16 * 60,
    conda:
        CONDA
    shell:
        """
        pbsv call {input.ref} \
            --ccs {input.svsig} {output.vcf}
        """


rule pbsv_index:
    input:
        vcf=rules.run_pbsv.output.vcf,
    output:
        vcf="results/{sm}/pbsv/{sm}.pbsv.vcf.gz",
        tbi="results/{sm}/pbsv/{sm}.pbsv.vcf.gz.tbi",
    conda:
        CONDA
    threads: 4
    shell:
        """
        bgzip -@ {threads} -c {input.vcf} > {output.vcf}
        tabix {output.vcf}
        """


rule run_sniffles:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
    output:
        vcf="results/{sm}/sniffles/{sm}.sniffles.vcf",
    threads: 16
    resources:
        mem_mb=64 * 1024,
    conda:
        CONDA
    shell:
        """
        sniffles --input {input.bam} \
            --sample-id {wildcards.sm} \
            --vcf {output.vcf} \
            --reference {input.ref}
        """


def get_sv_caller_outputs(wc):
    if not SV_CALLERS:
        return []
    rtn = expand(rules.pbsv_index.output, sm=SAMPLE)
    rtn += expand(rules.run_sniffles.output, sm=SAMPLE)
    return rtn
