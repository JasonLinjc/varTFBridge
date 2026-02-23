#!/usr/bin/env python
"""
VAR2TFBS Step 2: Predict common variant effects on TF binding sites.

Trait-agnostic: takes the merged BED file from Step 1 (containing all unique
common variants overlapping FOODIE footprints across traits) and predicts
variant effects on TF binding using FIMO motif scanning against JASPAR PWMs.

For each variant, extracts ref/alt sequences from the reference genome
(footprint ± ext_bp) and classifies TF binding changes as Create, Disrupt,
Increase, Decrease, or Unchange.

Input:
    - Merged BED from Step 1: GWFM_variants_in_{cell}.merged.hg38.bed
      (columns: Chromosome, Start, End, SNP, footprint_region)
    - Allele source: per-trait overlap CSV directory or single CSV with SNP, A1, A2

Usage:
    python scripts/comvar_var2tfbs.py \
        --input-bed results/comvar_footprint_overlap_credible/GWFM_variants_in_K562.merged.hg38.bed \
        --allele-src results/comvar_footprint_overlap_credible/K562.merged.hg38 \
        --ref-genome data/reference/hg38.fa \
        --jaspar-meme data/JASPAR_MEME/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
        --out-dir results/comvar_var2tfbs_results
"""

import argparse
import os

import numpy as np
import pandas as pd
import pyfaidx
from kipoiseq import Interval
from memelite import fimo
from memelite.io import read_meme
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Pretty printing helpers
# ---------------------------------------------------------------------------

def print_header(text):
    width = 62
    print()
    print("=" * width)
    print(f"  {text}")
    print("=" * width)


def print_step(step_num, total, text):
    print(f"\n[Step {step_num}/{total}] {text}")
    print("-" * 50)


def print_info(text):
    print(f"  -> {text}")


def print_success(text):
    print(f"  [OK] {text}")


# ---------------------------------------------------------------------------
# Reference genome helper
# ---------------------------------------------------------------------------

class FastaStringExtractor:
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval) -> str:
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(
            interval.chrom,
            max(interval.start, 0),
            min(interval.end, chromosome_length),
        )
        sequence = str(
            self.fasta.get_seq(
                trimmed_interval.chrom,
                trimmed_interval.start + 1,
                trimmed_interval.stop,
            ).seq
        ).upper()
        pad_upstream = "N" * max(-interval.start, 0)
        pad_downstream = "N" * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def parse_footprint_region(region_str):
    """Parse 'chr1:100-200' into (chrom, start, end)."""
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def load_allele_info(allele_src):
    """Load A1/A2 alleles from a CSV file or directory of CSVs.

    If allele_src is a directory, reads all CSVs and merges to get the
    union of SNP allele info (needed when no single trait covers all variants).
    """
    import glob as _glob

    if os.path.isdir(allele_src):
        csv_files = sorted(_glob.glob(os.path.join(allele_src, "*.csv")))
        dfs = []
        for f in csv_files:
            dfs.append(pd.read_csv(f, usecols=["SNP", "A1", "A2"]))
        df = pd.concat(dfs, ignore_index=True)
    else:
        df = pd.read_csv(allele_src, usecols=["SNP", "A1", "A2"])
    df = df.drop_duplicates(subset="SNP", keep="first")
    return df.set_index("SNP")[["A1", "A2"]].to_dict("index")


def load_merged_bed(bed_path):
    """Load merged BED file into a DataFrame."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["Chromosome", "Start", "End", "SNP", "footprint_region"],
    )
    return df


def extract_sequences(df, seq_extractor, ext_bp):
    """Add Ref_seq and Alt_seq columns to the variant DataFrame.

    For each variant, the extracted window is footprint ± ext_bp.
    Determines which of A1/A2 is the reference allele by checking against
    the reference genome sequence at the variant position.
    """
    ref_seqs, alt_seqs, seq_names = [], [], []
    n_a1_ref, n_a2_ref, n_mismatch = 0, 0, 0
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Extracting sequences", ncols=60):
        fp_chrom, fp_start, fp_end = parse_footprint_region(row["footprint_region"])
        var_pos = int(row["Start"])  # 0-based variant position (hg38)

        # Check which allele matches the reference genome
        a1, a2 = row["A1"], row["A2"]
        ref_base = seq_extractor.extract(
            Interval(row["Chromosome"], var_pos, var_pos + len(a2))
        )
        if ref_base == a2:
            ref_var, alt_var = a2, a1
            n_a2_ref += 1
        elif ref_base == a1:
            ref_var, alt_var = a1, a2
            n_a1_ref += 1
        else:
            # Check if A1 matches (for indels, compare with A1 length)
            ref_base_a1 = seq_extractor.extract(
                Interval(row["Chromosome"], var_pos, var_pos + len(a1))
            )
            if ref_base_a1 == a1:
                ref_var, alt_var = a1, a2
                n_a1_ref += 1
            else:
                # Default to A2=ref if neither matches exactly
                ref_var, alt_var = a2, a1
                n_mismatch += 1

        if len(ref_var) > 1:
            pad = len(ref_var)
        else:
            pad = 0
        window_start = fp_start - ext_bp - pad
        window_end = fp_end + ext_bp + pad

        var_pos_in_window = var_pos - window_start

        ref_seq = seq_extractor.extract(Interval(fp_chrom, window_start, window_end))
        alt_seq = (
            ref_seq[:var_pos_in_window]
            + alt_var
            + ref_seq[var_pos_in_window + len(ref_var) :]
        )

        # Sequence name: rsID|chrom:pos(0-based):ref:alt|footprint_region
        seq_name = (
            f"{row['SNP']}|{row['Chromosome']}:{var_pos}:{ref_var}:{alt_var}"
            f"|{row['footprint_region']}"
        )

        ref_seqs.append(ref_seq)
        alt_seqs.append(alt_seq)
        seq_names.append(seq_name)

    print(f"  Allele assignment: A2=ref {n_a2_ref:,}, A1=ref {n_a1_ref:,}, mismatch {n_mismatch:,}")

    df = df.copy()
    df["Ref_seq"] = ref_seqs
    df["Alt_seq"] = alt_seqs
    df["seq_name"] = seq_names
    return df


def write_fasta(df, out_path, allele="ref"):
    """Write ref or alt FASTA from the DataFrame."""
    col = "Ref_seq" if allele == "ref" else "Alt_seq"
    dedup = df.drop_duplicates(subset=["seq_name"]).reset_index(drop=True)
    with open(out_path, "w") as fh:
        for _, row in dedup.iterrows():
            fh.write(f">{row['seq_name']}\n{row[col]}\n")
    return out_path


# ---------------------------------------------------------------------------
# FIMO motif scanning
# ---------------------------------------------------------------------------

def run_fimo(motifs, fasta_path, ext_bp, threshold=0.0001):
    """Run FIMO and post-process hits.

    Returns a DataFrame of hits where the variant falls within the motif.
    """
    hits = fimo(motifs, fasta_path, threshold=threshold)
    dfs = [h for h in hits if not h.empty]
    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)
    # motif_name format from memelite: "MA0001.1 ARNT"
    df["motif_id"] = df["motif_name"].str.split(" ").str[0] + df["strand"]
    df["motif_name"] = df["motif_name"].str.split(" ").str[1] + df["strand"]
    df["motif_info"] = (
        df["motif_id"] + "-" + df["motif_name"]
        + "_" + df["start"].astype(str)
        + "_" + df["end"].astype(str)
    )
    df["var-motif"] = df["sequence_name"] + "-" + df["motif_info"]

    # Parse sequence_name: rsID|chrom:pos:ref:alt|footprint_region
    df["rsID"] = df["sequence_name"].str.split("|").str[0]
    df["var_locus"] = df["sequence_name"].str.split("|").str[1]
    df["foodie_id"] = df["sequence_name"].str.split("|").str[2]
    df["pos_var"] = df["var_locus"].str.split(":").str[1].astype(int)

    fp_parts = df["foodie_id"].str.split(":")
    df["chrom_fp"] = fp_parts.str[0]
    fp_coords = fp_parts.str[1].str.split("-")
    df["start_fp"] = fp_coords.str[0].astype(int)
    df["end_fp"] = fp_coords.str[1].astype(int)

    # Map motif coordinates back to genomic coordinates
    df["start_motif"] = df["start"] - ext_bp + df["start_fp"]
    df["end_motif"] = df["end"] - ext_bp + df["start_fp"]

    # Keep only hits where the variant falls within the motif
    df["var_in_motif"] = (df["pos_var"] >= df["start_motif"]) & (
        df["pos_var"] < df["end_motif"]
    )
    df = df[df["var_in_motif"]].reset_index(drop=True)
    return df


def classify_tf_changes(ref_alt_df):
    """Classify TF binding changes based on ref/alt FIMO p-values.

    Categories: Create, Disrupt, Increase, Decrease, Unchange, NoTFBS.
    """
    ref_p = ref_alt_df["p-value_ref"].fillna(1)
    alt_p = ref_alt_df["p-value_alt"].fillna(1)

    conditions = [
        (ref_p != 1) & (alt_p != 1) & (ref_p > alt_p),   # Increase
        (ref_p != 1) & (alt_p != 1) & (ref_p < alt_p),   # Decrease
        (ref_p != 1) & (alt_p != 1) & (ref_p == alt_p),  # Unchange
        (ref_p == 1) & (alt_p != 1),                       # Create
        (ref_p != 1) & (alt_p == 1),                       # Disrupt
    ]
    choices = ["Increase", "Decrease", "Unchange", "Create", "Disrupt"]
    ref_alt_df["TF_change"] = np.select(conditions, choices, default="NoTFBS")
    return ref_alt_df


def deduplicate_tf_hits(ref_alt_df):
    """Corrected TFBS: keep only the best p-value per TF per variant.

    For each rsID, takes ref and alt hits separately, keeps only the best
    (lowest) p-value per TF, then re-merges on TF with outer join and
    re-classifies TF binding changes.
    """
    results = []
    for rsid, var_tfbs in ref_alt_df.groupby("rsID"):
        foodie_id = var_tfbs["foodie_id"].iloc[0]

        # Ref: keep best p-value per TF
        ref_cols = ["TF", "motif_name_ref", "motif_id_ref", "p-value_ref", "start_ref", "end_ref"]
        var_ref = var_tfbs[ref_cols].dropna(subset=["p-value_ref"])
        if not var_ref.empty:
            var_ref = (
                var_ref.sort_values(by=["TF", "p-value_ref"])
                .drop_duplicates(subset=["TF"], keep="first")
                .reset_index(drop=True)
            )

        # Alt: keep best p-value per TF
        alt_cols = ["TF", "motif_name_alt", "motif_id_alt", "p-value_alt", "start_alt", "end_alt"]
        var_alt = var_tfbs[alt_cols].dropna(subset=["p-value_alt"])
        if not var_alt.empty:
            var_alt = (
                var_alt.sort_values(by=["TF", "p-value_alt"])
                .drop_duplicates(subset=["TF"], keep="first")
                .reset_index(drop=True)
            )

        # Re-merge on TF
        if var_ref.empty and var_alt.empty:
            continue
        elif var_ref.empty:
            merged = var_alt.copy()
        elif var_alt.empty:
            merged = var_ref.copy()
        else:
            merged = var_ref.merge(var_alt, on="TF", how="outer").reset_index(drop=True)

        merged["rsID"] = rsid
        merged["foodie_id"] = foodie_id
        results.append(merged)

    if not results:
        return pd.DataFrame()

    deduped = pd.concat(results, ignore_index=True)
    deduped = classify_tf_changes(deduped)
    return deduped


def var2tfbs(motifs, ref_fasta, alt_fasta, ext_bp, fimo_threshold=0.0001):
    """Run FIMO on ref and alt FASTA, merge, and classify TF changes."""
    print("Running FIMO on ref sequences ...")
    ref_hits = run_fimo(motifs, ref_fasta, ext_bp, threshold=fimo_threshold)
    print("Running FIMO on alt sequences ...")
    alt_hits = run_fimo(motifs, alt_fasta, ext_bp, threshold=fimo_threshold)

    if ref_hits.empty and alt_hits.empty:
        return pd.DataFrame()

    if ref_hits.empty:
        ref_alt = alt_hits.rename(columns=lambda c: c + "_alt" if c not in ["var-motif"] else c)
        ref_alt["TF_change"] = "Create"
        ref_alt["sequence_name"] = ref_alt.get("sequence_name_alt", "")
        return ref_alt

    if alt_hits.empty:
        ref_alt = ref_hits.rename(columns=lambda c: c + "_ref" if c not in ["var-motif"] else c)
        ref_alt["TF_change"] = "Disrupt"
        ref_alt["sequence_name"] = ref_alt.get("sequence_name_ref", "")
        return ref_alt

    print("Merging ref/alt hits ...")
    ref_alt = ref_hits.merge(
        alt_hits, suffixes=["_ref", "_alt"], on=["var-motif"], how="outer"
    ).reset_index(drop=True)

    ref_alt["sequence_name"] = ref_alt["sequence_name_ref"].fillna(
        ref_alt["sequence_name_alt"]
    )

    # Parse variant/footprint info from sequence_name
    ref_alt["rsID"] = ref_alt["sequence_name"].str.split("|").str[0]
    ref_alt["variant_id"] = ref_alt["sequence_name"].str.split("|").str[1]
    ref_alt["foodie_id"] = ref_alt["sequence_name"].str.split("|").str[2]
    ref_alt["TF"] = (
        ref_alt["motif_name_ref"]
        .fillna(ref_alt["motif_name_alt"])
        .str[:-1]
        .str.upper()
    )

    ref_alt = classify_tf_changes(ref_alt)

    # Corrected TFBS: deduplicate to keep best p-value per TF per variant
    print("Deduplicating: keeping best p-value per TF per variant ...")
    ref_alt_corrected = deduplicate_tf_hits(ref_alt)

    return ref_alt_corrected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="VAR2TFBS Step 2: predict common variant effects on TF binding (trait-agnostic)"
    )
    p.add_argument(
        "--input-bed", required=True,
        help="Merged BED file from Step 1 "
             "(e.g. comvar_footprint_overlap_credible/GWFM_variants_in_K562.merged.hg38.bed)",
    )
    p.add_argument(
        "--allele-src", required=True,
        help="Directory of per-trait overlap CSVs or single CSV with SNP, A1, A2 columns "
             "(e.g. comvar_footprint_overlap_credible/K562.merged.hg38)",
    )
    p.add_argument("--ref-genome", required=True, help="Path to hg38.fa reference genome")
    p.add_argument("--jaspar-meme", required=True, help="Path to JASPAR MEME motif file")
    p.add_argument(
        "--out-dir", default="./results/comvar_var2tfbs_results",
        help="Output directory (default: ./results/comvar_var2tfbs_results)",
    )
    p.add_argument("--ext-bp", type=int, default=30, help="Sequence extension in bp (default: 30)")
    p.add_argument("--fimo-threshold", type=float, default=0.0001, help="FIMO p-value threshold (default: 0.0001)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    fasta_dir = os.path.join(args.out_dir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    # Derive cell name from BED filename (e.g. GWFM_variants_in_K562.merged.hg38.bed -> K562)
    bed_basename = os.path.basename(args.input_bed)
    cell = bed_basename.replace("GWFM_variants_in_", "").split(".")[0]

    total_steps = 3

    print_header("Common Variant VAR2TFBS Analysis")
    print_info(f"Cell:            {cell}")
    print_info(f"Input BED:       {args.input_bed}")
    print_info(f"Allele source:   {args.allele_src}")
    print_info(f"Reference genome: {args.ref_genome}")
    print_info(f"JASPAR motifs:   {args.jaspar_meme}")
    print_info(f"Seq extension:   {args.ext_bp} bp")
    print_info(f"FIMO threshold:  {args.fimo_threshold}")

    # Step 1: Load variants and allele info
    print_step(1, total_steps, "Loading variants and allele info")
    df = load_merged_bed(args.input_bed)
    print_info(f"{len(df):,} variants from merged BED")

    allele_dict = load_allele_info(args.allele_src)
    df["A1"] = df["SNP"].map(lambda x: allele_dict.get(x, {}).get("A1", ""))
    df["A2"] = df["SNP"].map(lambda x: allele_dict.get(x, {}).get("A2", ""))
    n_missing = (df["A1"] == "").sum()
    if n_missing > 0:
        print_info(f"WARNING: {n_missing:,} variants missing allele info, skipping")
        df = df[df["A1"] != ""].reset_index(drop=True)
    print_success(f"{len(df):,} variants with allele info")

    # Step 2: Extract ref/alt sequences
    print_step(2, total_steps, "Extracting ref/alt sequences")
    print_info("Loading reference genome ...")
    seq_extractor = FastaStringExtractor(args.ref_genome)
    print_info("Loading JASPAR motifs ...")
    motifs = read_meme(args.jaspar_meme)

    df = extract_sequences(df, seq_extractor, args.ext_bp)
    df = df.drop_duplicates(subset=["seq_name"]).reset_index(drop=True)
    seq_extractor.close()

    ref_fasta = os.path.join(fasta_dir, f"{cell}_ref_seqExt{args.ext_bp}bp.fa")
    alt_fasta = os.path.join(fasta_dir, f"{cell}_alt_seqExt{args.ext_bp}bp.fa")
    write_fasta(df, ref_fasta, allele="ref")
    write_fasta(df, alt_fasta, allele="alt")
    print_success(f"Wrote {ref_fasta}")
    print_success(f"Wrote {alt_fasta}")

    # Step 3: FIMO motif scanning + TF change classification
    print_step(3, total_steps, "FIMO motif scanning + TF change classification")
    result = var2tfbs(motifs, ref_fasta, alt_fasta, args.ext_bp, args.fimo_threshold)

    if result.empty:
        print_info("No motif hits found.")
        return

    result["cell"] = cell
    out_path = os.path.join(args.out_dir, f"{cell}_var2tfbs.csv")
    result.to_csv(out_path, index=False)

    print()
    print("=" * 62)
    print(f"  Results: {len(result):,} variant-TF pairs across {result['rsID'].nunique()} variants")
    print(f"    Create:    {(result['TF_change'] == 'Create').sum():,}")
    print(f"    Disrupt:   {(result['TF_change'] == 'Disrupt').sum():,}")
    print(f"    Increase:  {(result['TF_change'] == 'Increase').sum():,}")
    print(f"    Decrease:  {(result['TF_change'] == 'Decrease').sum():,}")
    print(f"    Unchange:  {(result['TF_change'] == 'Unchange').sum():,}")
    print(f"  Output: {out_path}")
    print("=" * 62)
    print()


if __name__ == "__main__":
    main()
