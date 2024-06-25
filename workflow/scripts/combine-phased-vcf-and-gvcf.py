#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
from cyvcf2 import VCF, Writer
from tqdm import tqdm


def add_phasing_format_to_header(vcf):
    vcf.add_format_to_header(
        dict(ID="PS", Number=1, Type="Integer", Description="Phase set identifier")
    )
    vcf.add_format_to_header(
        dict(ID="PF", Number=1, Type="String", Description="Phasing flag")
    )


def get_tag(rec):
    return (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])


def run(vcf: Path, gvcf: Path, outfile: Path, break_n=None):
    phased_vcf_data = {}
    vcf = VCF(vcf)
    logging.info("Reading and storing phased vcf")
    for idx, rec in enumerate(tqdm(vcf)):
        if idx == break_n:
            break
        # check if phased
        if rec.genotypes[0][2]:
            # error when there are multiple genotypes
            if len(rec.genotypes) != 1:
                sys.stderr.write("\n")
                logging.error(f"Multiple genotypes in phased vcf: {rec}")
                # and now fail
                sys.exit(1)
            # save the record to insert into gvcf
            phased_vcf_data[get_tag(rec)] = rec

    logging.info("Reading and modifying gvcf")
    gvcf = VCF(gvcf)
    add_phasing_format_to_header(gvcf)
    o_gvcf = Writer(outfile, gvcf, mode="wz")
    o_gvcf.add_to_header(
        f"{sys.argv[0]}##Command={}".format(" ".join(sys.argv)),
    )
    change_count = 0
    for idx, rec in enumerate(tqdm(gvcf)):
        tag = get_tag(rec)
        if tag in phased_vcf_data:
            # error when there are multiple genotypes
            if len(rec.genotypes) != 1:
                sys.stderr.write("\n")
                logging.error(f"Multiple genotypes in phased vcf: {rec}")
                # and now fail
                sys.exit(1)
            vcf_rec = phased_vcf_data[tag]
            new_gt = vcf_rec.genotypes
            rec.genotypes = new_gt
            change_count += 1
        o_gvcf.write_record(rec)

    logging.info(
        f"Total number of phased records in phased vcf: {len(phased_vcf_data)}"
    )
    logging.info(f"Total number of records in gvcf: {idx+1}")
    logging.info(f"Total number of gvcf records changed: {change_count}")


def main(
    vcf: Path,
    gvcf: Path,
    outfile: Optional[Path] = None,
    *,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param vcf: phased vcf
    :param gvcf: unphased gvcf
    :param outfile: phased gvcf
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    if outfile is None:
        outfile = sys.stdout
    # Set logging level
    logging.basicConfig(
        level=logging.DEBUG if verbose > 0 else logging.INFO,
    )
    run(vcf, gvcf, outfile)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
