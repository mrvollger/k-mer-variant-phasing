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
pd.set_option('display.width', 150)
pd.set_option('display.max_columns', 150)

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
        "-m",
        "--max-frac-disagree",
        help="Maximum fraction of reads in a phaseblock that can disagree with k-mers assignments before the whole phaseblock is set to unknown",
        default=0.2,
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


def assign_based_on_kmer(group_df, max_frac_disagree=0.2):
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
    if kmer_counts.get(looser, 0) / group_df.shape[0] > max_frac_disagree:
        logging.info(
            "Too many reads disagree within phase block. Setting all to unphased."
        )
        winner = UNKNOWN
    # set the merged haplotype to the winner between paternal and maternal
    group_df.merged_hap = winner
    return group_df


def main():
    args = parse()
    kmer_df = read_kmer(args.kmer)
    variant_df = read_variant(args.variant)
    merged_df = kmer_df.merge(variant_df, on=READ_COL, how="left")

    merged_df["merged_hap"] = UNKNOWN
    merged_df = (
        merged_df[merged_df.variant_hap != UNKNOWN]
        .groupby(["phase_block", "variant_hap"], group_keys=False)
        .apply(lambda row: assign_based_on_kmer(row, args.max_frac_disagree))
    )
    # count disagreements
    has_kmer = merged_df.kmer_hap != UNKNOWN
    disagreements = merged_df[has_kmer].kmer_hap != merged_df[has_kmer].merged_hap
    logging.info(
        f"{disagreements.sum()/merged_df[has_kmer].shape[0]:.2%} of variant based calls that disagree with the phased kmer calls."
    )
    if args.prioritize_kmer:
        logging.info("Prioritizing kmer phasing.")

    # say # of reassignments
    phasing_added = (merged_df[~has_kmer].merged_hap != UNKNOWN).sum()
    logging.info(
        f"{phasing_added:,} unassigned sequences were assigned a haplotype based on variant based phasing"
    )

    # assign final haplotype
    merged_df["hap"] = merged_df.kmer_hap
    can_reassign_kmer = (
        merged_df.merged_hap != UNKNOWN
    )  # & (merged_df.kmer_hap == UNKNOWN)
    merged_df.loc[can_reassign_kmer, "hap"] = merged_df.merged_hap[can_reassign_kmer]

    # make final outputs
    out = kmer_df.merge(
        merged_df[[READ_COL, "merged_hap", "hap", "variant_hap", "phase_block"]],
        on=READ_COL,
        how="left",
    )
    # add any reads in the variant df missing from the kmer df
    out = out.merge(variant_df[[READ_COL]], on=READ_COL, how="outer")
    logging.info(f"\n{out}\n")

    # set NA haps to the kmer value (if it exists)
    out.loc[out.hap.isna(), "hap"] = out.kmer_hap[out.hap.isna()]
    # set all NAs to unknown
    out.fillna(UNKNOWN, inplace=True)

    # drop ambiguous reads from phasing
    if not args.prioritize_kmer and not args.max_frac_disagree:
        known = (out.kmer_hap.isin([MATERNAL, PATERNAL])) & (
            out.merged_hap.isin([MATERNAL, PATERNAL])
        )
        disagreements = known & (out.kmer_hap != out.merged_hap)
        out.loc[disagreements, "hap"] = UNKNOWN
        logging.info(
            f"{disagreements.sum()/out.shape[0]:.2%} of total reads changed to unknown due to variant and kmer disagreements"
        )

    z = (out.variant_hap != UNKNOWN).sum()
    logging.info(f"Variant based phasing rate: {z/len(out):.2%}")
    z = (out.kmer_hap != UNKNOWN).sum()
    logging.info(f"K-mer phasing rate: {z/len(out):.2%}")
    z = (out.hap != UNKNOWN).sum()
    logging.info(f"Merged phasing rate: {z/len(out):.2%}")

    out[[READ_COL, "hap", "kmer_hap", "variant_hap", "phase_block"]].to_csv(
        args.output, index=False, sep="\t"
    )


if __name__ == "__main__":
    main()
