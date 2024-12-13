rule deepvariant_chunk:
    input:
        bam=get_hifi_bam,
        bai=get_hifi_bai,
        ref=get_ref,
        fai=get_fai,
    output:
        # no longer in 1.8.0
        #html=temp("temp/{sm}/deepvariant/{rgn}/{sm}.deepvariant.visual_report.html"),
        vcf=temp("temp/{sm}/deepvariant/{rgn}/{sm}.deepvariant.vcf.gz"),
        vcf_tbi=temp("temp/{sm}/deepvariant/{rgn}/{sm}.deepvariant.vcf.gz.tbi"),
        gvcf=temp("temp/{sm}/deepvariant/{rgn}/{sm}.deepvariant.gvcf.gz"),
        gvcf_tbi=temp("temp/{sm}/deepvariant/{rgn}/{sm}.deepvariant.gvcf.gz.tbi"),
    threads: 8
    resources:
        mem_mb=32 * 1024,
    singularity:
        f"docker://google/deepvariant:{BIN_VERSION}"
    params:
        model_type="PACBIO",  # WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA
        bin_version=BIN_VERSION,
        rgn=get_region,
    shell:
        """
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref={input.ref} \
            --reads={input.bam} \
            --regions {params.rgn} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads} \
            --call_variants_extra_args="allow_empty_examples=true"
        """


# --intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir" \ **Optional.
# f"docker://google/deepvariant:{BIN_VERSION}-gpu"
# singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \


rule deepvariant_merge:
    input:
        vcfs=expand(
            rules.deepvariant_chunk.output.vcf, rgn=make_regions(), allow_missing=True
        ),
        gvcfs=expand(
            rules.deepvariant_chunk.output.gvcf, rgn=make_regions(), allow_missing=True
        ),
    output:
        vcf="results/{sm}/deepvariant/{sm}.deepvariant.vcf.gz",
        vcf_tbi="results/{sm}/deepvariant/{sm}.deepvariant.vcf.gz.tbi",
        gvcf="results/{sm}/deepvariant/{sm}.deepvariant.gvcf.gz",
        gvcf_tbi="results/{sm}/deepvariant/{sm}.deepvariant.gvcf.gz.tbi",
        tmp=temp("temp/{sm}/deepvariant/{sm}.deepvariant.vcf.gz"),
        vcfs=temp("temp/{sm}/deepvariant/vcfs.txt"),
    threads: 16
    resources:
        mem_mb=64 * 1024,
    conda:
        CONDA
    shell:
        """
        > {output.vcfs}
        for vcf in {input.vcfs}; do
            sample=$(bcftools head $vcf | tail -n 1 | cut -f 10)
            if [[ $sample == "default" ]]; then
                echo "skipping, no variants in $vcf"
            else
                echo $vcf >> {output.vcfs}
            fi
        done

        bcftools concat $(cat {output.vcfs}) -o {output.tmp}
        bcftools reheader --threads {threads} \
            -s <(echo {wildcards.sm}) {output.tmp} \
            -o {output.vcf} 
        bcftools index -t {output.vcf} 


        > {output.vcfs}
        for vcf in {input.gvcfs}; do
            sample=$(bcftools head $vcf | tail -n 1 | cut -f 10)
            if [[ $sample == "default" ]]; then
                echo "skipping, no variants in $vcf"
            else
                echo $vcf >> {output.vcfs}
            fi
        done
        bcftools concat $(cat {output.vcfs}) -o {output.tmp}
        bcftools reheader --threads {threads} \
            -s <(echo {wildcards.sm}) {output.tmp} \
            -o {output.gvcf} 
        bcftools index -t {output.gvcf} 
       """


rule deepvariant:
    input:
        expand(rules.deepvariant_merge.output, sm=SAMPLE),
