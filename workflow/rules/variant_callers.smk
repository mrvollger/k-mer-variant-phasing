rule run_pbsv:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
        sm=SAMPLE,
    output:
        svsig="results/{sm}/pbsv/{sm}.pbsv.svsig.gz",
        vcf="results/{sm}/pbsv/{sm}.pbsv.vcf",
    threads: 16
    resources:
        mem_mb=64 * 1024,
    conda:
        CONDA
    shell:
        """
        pbsv discover \
            {input.bam} \
            {output.svsig}
            
        pbsv call {input.ref} \
        --ccs {output.svsig} {output.vcf}
            
        """


rule run_sniffles:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
        sm=SAMPLE,
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
        --vcf {output.vcf} \
        --reference {input.ref}
            
        """


def get_sv_caller_outputs():
    rtn = expand(rules.run_pbsv.output, sm=SAMPLE)
    rtn += expand(rules.run_sniffles.output, sm=SAMPLE)
    print(f"SV_CALLERS: {SV_CALLERS}\nSV_CALLER_OUTPUTS: {rtn}")
    if not SV_CALLERS:
        return []
    return rtn
