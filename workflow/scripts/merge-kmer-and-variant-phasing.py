#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Docs here
"""

__author__ = "Mitchell R. Vollger"
__credits__ = ["Mitchell R. Vollger"]
__maintainer__ = "Mitchell R. Vollger"
__email__ = "mrvollger@gmail.com"
__status__ = "Development"

MAX_ROWS = 100

import logging
import argparse
import sys
import pandas as pd
import numpy as np

pd.set_option("display.width", 150)
pd.set_option("display.max_columns", 150)

READ_COL = "read"
PATERNAL = "pat"
MATERNAL = "mat"
UNKNOWN = "unk"


def read_kmer(file):
    df = pd.read_csv(file, sep="\t")
    df.drop_duplicates(inplace=True)
    logging.info(f"Read {len(df):,} sequences from {file}")
    logging.info(f"K-mer based read counts:\n{df.kmer_hap.value_counts()}")
    return df


def read_variant(file):
    df = pd.read_csv(file, sep="\t")
    df.drop_duplicates(inplace=True)
    logging.info(f"Read {len(df):,} sequences from {file}")
    logging.info(
        f"Read {len(df.phase_block.unique()):,} phase blocks from variant based phasing"
    )
    logging.info(f"Variant based read counts:\n{df.variant_hap.value_counts()}")
    return df


def parse():
    """Change all reads (and header) in a bam file to have one read group (RG)"""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "kmer",
        help="Input files with paternal maternal and unknown haplotypes for each read",
    )
    parser.add_argument(
        "variant",
        help="Input file",
    )
    parser.add_argument(
        "-k",
        "--prioritize-kmer",
        help="keep phased kmer reads even if they disagree with variant calls",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--strict",
        help="Set any reads that disagree in variant and kmer phasing to unknown",
        action="store_true",
    )
    parser.add_argument(
        "-m",
        "--max-frac-disagree",
        help="Maximum fraction of reads in a phaseblock that can disagree with k-mers assignments before the whole phaseblock is set to unknown",
        default=0.05,
        type=float,
    )
    parser.add_argument("-o", "--output", help="Output file", default=sys.stdout)
    parser.add_argument(
        "-t", "--threads", help="Number of threads to use", type=int, default=8
    )
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def assign_per_hap(group_df, max_frac_disagree):
    cur_variant_hap = group_df.variant_hap.iloc[0]
    # if variant reads are unphased use the kmer phasing
    if cur_variant_hap == UNKNOWN or cur_variant_hap.isna():
        group_df.merged_hap = group_df.kemr_hap
        group_df.fraction_disagreement = 0.0
        return group_df

    # check which kmer haplotype is more common in the phaseblock
    kmer_counts = group_df.kmer_hap.value_counts()
    paternal = kmer_counts.get(PATERNAL, 0)
    maternal = kmer_counts.get(MATERNAL, 0)

    winner = UNKNOWN
    looser = UNKNOWN
    if paternal > maternal:
        winner = PATERNAL
        looser = MATERNAL
    elif maternal > 0:
        winner = MATERNAL
        looser = PATERNAL
    # set all to unphased if too many disagree
    fraction_disagreement = kmer_counts.get(looser, 0) / group_df.shape[0]
    if fraction_disagreement > max_frac_disagree:
        # logging.info("Too many reads disagree. Setting all to unphased.")
        winner = UNKNOWN
    # set the merged haplotype to the winner between paternal and maternal
    group_df.merged_hap = winner
    group_df.fraction_disagreement = fraction_disagreement
    return group_df


def assign_per_phase_block(group_df, max_frac_disagree):
    cur_phase_block = group_df.phase_block.iloc[0]
    # if the current phaseblock is unknown use the kmer phasing
    if cur_phase_block == UNKNOWN:
        group_df.merged_hap = group_df.kemr_hap
        group_df.fraction_disagreement = 0.0
        return group_df

    # now go through H1, H2 in the phaseblock
    group_df = group_df.groupby("variant_hap", group_keys=False).apply(
        lambda row: assign_per_hap(row, max_frac_disagree)
    )

    is_h1 = group_df.variant_hap == 1
    is_h2 = group_df.variant_hap == 2
    h1_assignment = group_df.merged_hap[is_h1].iloc[0]
    h2_assignment = group_df.merged_hap[is_h2].iloc[0]
    # mat < pat < unk
    assignments = sorted([h1_assignment, h2_assignment])
    if assignments[0] != UNKNOWN and assignments[1] == UNKNOWN:
        logging.info(
            f"Within a phase block one haplotype was assigned maternal or paternal but the other was unassigned.\n{group_df}"
        )
        # TODO sync up the unknown haplotype

    return group_df


def main():
    args = parse()
    kmer_df = read_kmer(args.kmer)
    variant_df = read_variant(args.variant)
    merged_df = kmer_df.merge(variant_df, on=READ_COL, how="outer")
    merged_df.fillna(UNKNOWN, inplace=True)
    n_reads = merged_df.shape[0]

    merged_df["merged_hap"] = UNKNOWN
    merged_df["fraction_disagreement"] = 0.0
    merged_df = (
        merged_df[merged_df.variant_hap != UNKNOWN]
        .groupby(["phase_block", "variant_hap"], group_keys=False)
        .apply(lambda row: assign_per_phase_block(row, args.max_frac_disagree))
    )

    # say the number of disagreements
    can_disagree = (
        (merged_df.variant_hap != UNKNOWN)
        & (merged_df.kmer_hap != UNKNOWN)
        & (merged_df.merge_hap != UNKNOWN)
    )
    disagreements = (merged_df.merged_hap != merged_df.kmer_hap) & can_disagree
    logging.info(
        f"{disagreements.sum()/n_reads:.2%} of reads disagree in variant versus kmer based phasing."
    )

    # resolve disagreements
    if args.prioritize_kmer:
        logging.info("Prioritizing k-mer based phasing when there is disagreements.")
        merged_df.merged_hap[disagreements] = merged_df.kmer_hap[disagreements]
    if args.strict:
        logging.info(
            "Strict mode. Setting reads to unphased if there is a disagreement."
        )
        merged_df.merged_hap[disagreements] = UNKNOWN

    # say the number of reassigned reads
    can_be_improved = (merged_df.variant_hap != UNKNOWN) | (
        merged_df.kmer_hap != UNKNOWN
    )
    improved = merged_df[can_be_improved].merged_hap != UNKNOWN
    logging.info(f"{improved.sum()/n_reads:.2%} of previously unphased reads were phased.")

    # log the phasing rates
    z = (merged_df.variant_hap != UNKNOWN).sum()
    logging.info(f"Variant based phasing rate: {z/n_reads:.2%}")
    z = (merged_df.kmer_hap != UNKNOWN).sum()
    logging.info(f"K-mer phasing rate: {z/n_reads:.2%}")
    z = (merged_df.merged_hap != UNKNOWN).sum()
    logging.info(f"Merged phasing rate: {z/n_reads:.2%}")

    # write the output
    merged_df["hap"] = merged_df.merged_hap
    merged_df[
        [
            READ_COL,
            "hap",
            "kmer_hap",
            "variant_hap",
            "phase_block",
            "fraction_disagreement",
        ]
    ].to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
