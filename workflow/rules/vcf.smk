BIN_VERSION="1.5.0"

rule deepvariant:
    input:
        bam=get_hifi_bai,
        ref=REFERENCE,
    output:
        vcf=temp("temp/{sm}/deepvariant/{sm}.deepvariant.vcf.gz"),
        gvcf=temp("temp/{sm}/deepvariant/{sm}.deepvariant.gvcf.gz"),
    threads: 8
    resources:
        mem_mb=64 * 1024,
    singularity:
        f"docker://google/deepvariant:{BIN_VERSION}"
        #f"docker://google/deepvariant:{BIN_VERSION}-gpu"
    params:
        model_type="PACBIO", # WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
        bin_version=BIN_VERSION,
    shell:
        """
        singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
            docker://google/deepvariant:"{params.bin_version}" \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref={input.ref} \
            --reads={input.bam} \
            --regions "chr20:10,000,000-10,010,000" \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads}
        """

#--intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir" \ **Optional.
# singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \
# docker://google/deepvariant:"${BIN_VERSION}-gpu" \