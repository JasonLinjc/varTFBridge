#!/usr/bin/env python
"""
Filter snpRes_hg38 files to create credible set common variant files.

Keeps common variants with PIP > threshold (configurable).
Annotates each variant with LCS credible-set-level PEP (PEP_cs) and
credible set ID (CS_id) by matching SNP IDs to exploded LCS files.

Usage:
    python scripts/comvar_filter_credible_set.py \
        --snpres-dir data/GWFM_erythroids/snpRes_hg38 \
        --lcs-dir data/GWFM_erythroids/lcs \
        --out-dir data/GWFM_erythroids/credible_set_snpRes
"""

import argparse
import glob
import os

import pandas as pd
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser(description="Filter snpRes to credible set variants")
    p.add_argument("--snpres-dir", required=True, help="Directory containing lifted-over .snpRes files")
    p.add_argument("--lcs-dir", default="data/GWFM_erythroids/lcs",
                   help="Directory containing .lcs files (default: data/GWFM_erythroids/lcs)")
    p.add_argument("--out-dir", default="data/GWFM_erythroids/credible_set_snpRes",
                   help="Output directory (default: data/GWFM_erythroids/credible_set_snpRes)")
    p.add_argument("--pip-threshold", type=float, default=0.1, help="Minimum PIP (default: 0.1)")
    return p.parse_args()


def load_lcs_annotations(lcs_file):
    """Load LCS file, explode SNPs, return SNP -> (PEP_cs, CS_id) mapping.

    When a SNP appears in multiple credible sets, keep the highest PEP.
    """
    lcs = pd.read_csv(lcs_file, sep=r"\s+")
    lcs["SNP"] = lcs["SNP"].str.split(",")
    lcs = lcs.explode("SNP")
    # Keep highest PEP per SNP
    lcs = lcs.sort_values("PEP", ascending=False).drop_duplicates(subset="SNP", keep="first")
    return lcs.set_index("SNP")[["PEP", "CS"]]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    snpres_files = sorted(glob.glob(os.path.join(args.snpres_dir, "*.snpRes")))
    if not snpres_files:
        print("ERROR: No .snpRes files found in", args.snpres_dir)
        return

    print(f"Found {len(snpres_files)} snpRes files")
    print(f"LCS directory: {args.lcs_dir}")
    print(f"Filtering: PIP > {args.pip_threshold}")
    print()

    summary = []
    for snpres_file in tqdm(snpres_files, desc="Filtering", unit="file", ncols=60):
        trait = os.path.basename(snpres_file).replace(".snpRes", "")
        df = pd.read_csv(snpres_file, sep="\t")

        # Load LCS PEP_cs for this trait
        lcs_file = os.path.join(args.lcs_dir, f"{trait}.lcs")
        if os.path.exists(lcs_file):
            lcs_annot = load_lcs_annotations(lcs_file)
            df["PEP_cs"] = df["Name"].map(lcs_annot["PEP"])
            df["CS_id"] = df["Name"].map(lcs_annot["CS"])
        else:
            tqdm.write(f"  Warning: {lcs_file} not found, PEP_cs/CS_id will be empty")
            df["PEP_cs"] = float("nan")
            df["CS_id"] = float("nan")

        filtered = df[df["PIP"] > args.pip_threshold].copy()

        out_path = os.path.join(args.out_dir, f"{trait}_credible_set_hg38.csv")
        filtered.to_csv(out_path, index=False)
        n_with_cs = filtered["PEP_cs"].notna().sum()
        summary.append((trait, len(df), len(filtered), n_with_cs))

    print()
    print(f"{'Trait':<12} {'Total':>12} {'Filtered':>10} {'w/ PEP_cs':>10}")
    print("-" * 46)
    for trait, total, n_filtered, n_cs in summary:
        print(f"{trait:<12} {total:>12,} {n_filtered:>10,} {n_cs:>10,}")

    print()
    print(f"Output: {args.out_dir}/")
    print("Done!")


if __name__ == "__main__":
    main()
