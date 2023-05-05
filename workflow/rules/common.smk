def get_cleaned_reads(sm, hap):
    if hap == "mat":
        reads = MAT_DATA
    elif hap == "pat":
        reads = PAT_DATA
    else:
        raise ValueError("hap must be mat or pat")

    if (
        reads.endswith(".fofn")
        or reads.endswith("bam")
        or reads.endswith("sam")
        or reads.endswith("cram")
    ):
        return expand(rules.collect_reads.output.fasta, sm=sm, hap=hap)
    else:
        return reads


def get_mat(wc):
    return get_cleaned_reads(wc.sm, "mat")


def get_pat(wc):
    return get_cleaned_reads(wc.sm, "pat")


def get_reads(wc):
    if wc.hap == "mat":
        return MAT_DATA
    return PAT_DATA


def get_hifi_bam(wc):
    if config.get("clean_bam"):
        return expand(rules.clean_bam.output.bam, sm=wc.sm, allow_missing=True)
    return HIFI_BAM

def get_vcf(wc):
    if config.get("clean_vcf"):
        return expand(rules.clean_vcf.output.vcf, sm=wc.sm, allow_missing=True)
    return VCF