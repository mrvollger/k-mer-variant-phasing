def get_cleaned_reads(sm, hap):
    if hap == "mat":
        reads = MAT_DATA
    elif hap == "pat":
        reads = PAT_DATA
    else:
        raise ValueError("hap must be mat or pat")
    # check format
    valid_extensions = [".fa", ".fasta", ".fq", ".fastq"]
    valid_extensions += [x + ".gz" for x in valid_extensions]
    valid_extensions += [".fofn", ".bam", ".sam", ".cram"]
    assert reads.endswith(tuple(valid_extensions)), (
        f"reads must end with one of {valid_extensions}, " f"but {reads} does not"
    )
    # get fasta
    if (
        reads.endswith(".fofn")
        or reads.endswith("bam")
        or reads.endswith("sam")
        or reads.endswith("cram")
    ):
        return expand(rules.collect_reads.output.fasta, sm=sm, hap=hap)
    else:
        return reads


def get_ref(wc):
    return REFERENCE


def get_fai(wc):
    return f"{get_ref(wc)}.fai"

REGIONS=["chr20_10000000-10010000"]
def get_region(wc):
    print(wc)
    return "chr20:10000000-10010000" 

def get_mat(wc):
    return get_cleaned_reads(wc.sm, "mat")


def get_pat(wc):
    return get_cleaned_reads(wc.sm, "pat")


def get_reads(wc):
    reads = PAT_DATA
    if wc.hap == "mat":
        reads = MAT_DATA
    return reads


def get_hifi_bam(wc):
    # if config.get("clean_bam"):
    #    return expand(rules.clean_bam.output.bam, sm=wc.sm, allow_missing=True)[0]
    return HIFI_BAM


def get_hifi_bai(wc):
    return f"{get_hifi_bam(wc)}.bai"


def get_vcf(wc):
    if config.get("clean_vcf"):
        return expand(rules.clean_vcf.output.vcf, sm=wc.sm, allow_missing=True)[0]
    return VCF


def get_tbi(wc):
    return f"{get_vcf(wc)}.tbi"


def get_meryl_input(wc):
    if wc.individual == "pro":
        return rules.hifi_fasta.output.fasta
    return get_cleaned_reads(wc.sm, wc.individual)
