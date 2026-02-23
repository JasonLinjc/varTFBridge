#!/usr/bin/env python
"""
Rare variant VAR2TFBS: Identify driver rare variants from burden test
leave-one-out analysis and predict their effects on TF binding sites.

Pipeline:
1. Load burden test results per trait and find significant footprints (p < 5e-8)
2. For each significant footprint, use leave-one-out (LOO) results to find the
   driver variant — the variant whose removal causes the largest increase in
   burden test p-value (i.e., highest LOO p-value)
3. Filter out footprints where no variant has MAC > min_carrier threshold
4. Run FIMO motif scanning on driver variants to predict TF binding effects

Usage:
    python scripts/rarevar_var2tfbs.py \
        --burden-dir data/burdentest_erythroids \
        --loo-file data/leaveoneout_results/K562.leave_one_out.all_traits.20251120.csv \
        --ref-genome data/reference/hg38.fa \
        --jaspar-meme data/JASPAR_MEME/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
        --out-dir rarevar_var2tfbs_results
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyfaidx
from kipoiseq import Interval
from memelite import fimo as run_fimo_engine
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
# Reference genome helper (shared with comvar_var2tfbs.py)
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
# Driver variant identification
# ---------------------------------------------------------------------------

def load_significant_footprints(burden_dir, sig_threshold=None):
    """Load burden test results and return significant footprints per trait.

    If sig_threshold is None, uses Bonferroni correction: 0.05 / N_footprints.
    """
    burden_files = sorted(glob.glob(os.path.join(burden_dir, "*.xlsx")))
    if not burden_files:
        raise FileNotFoundError(f"No .xlsx files found in {burden_dir}")

    # Determine threshold from first file if not specified
    first_df = pd.read_excel(burden_files[0])
    n_footprints = len(first_df)
    if sig_threshold is None:
        sig_threshold = 0.05 / n_footprints
    print_info(f"{n_footprints:,} footprints tested, significance threshold: p < {sig_threshold:.2e}")

    sig_records = []
    for f in burden_files:
        # Trait name: e.g. HC.K562_FOODIE_fps_rare_var_result.2025-05-29.xlsx -> HC
        trait = os.path.basename(f).split(".")[0]
        df = pd.read_excel(f)
        sig = df[df["p_regenie"] < sig_threshold]
        for _, row in sig.iterrows():
            sig_records.append({
                "trait": trait,
                "footprint": row["ID"],
                "burden_p": row["p_regenie"],
                "MAC": row["MAC"],
            })

    sig_df = pd.DataFrame(sig_records)
    print_info(f"{len(sig_df)} significant trait-footprint pairs across {sig_df['trait'].nunique()} traits")
    print_success(f"{sig_df['footprint'].nunique()} unique footprints")
    return sig_df


def find_driver_variants(sig_df, loo_df, min_carrier=30):
    """For each significant footprint × trait, find the driver variant from LOO.

    The driver is the variant whose removal causes the largest increase in
    burden test p-value (highest LOO p-value → most signal loss when removed).

    Filters:
    - Remove footprints where no variant has MAC > min_carrier
    - Keep only the driver variant per footprint (deduplicated across traits)
    """
    drivers = []
    removed_no_carrier = 0

    for _, sig_row in sig_df.iterrows():
        trait = sig_row["trait"]
        footprint = sig_row["footprint"]

        loo = loo_df[(loo_df["trait"] == trait) & (loo_df["foodie"] == footprint)]
        if loo.empty:
            continue

        # Filter: at least one variant with MAC > min_carrier
        if not (loo["MAC_rare"] > min_carrier).any():
            removed_no_carrier += 1
            continue

        # Driver = variant whose removal causes highest LOO p-value
        driver_idx = loo["p_regenie"].idxmax()
        driver = loo.loc[driver_idx]

        drivers.append({
            "trait": trait,
            "footprint": footprint,
            "burden_p": sig_row["burden_p"],
            "driver_variant": driver["ID_rare"],
            "driver_loo_p": driver["p_regenie"],
            "driver_MAC": driver["MAC_rare"],
            "driver_ALLELE0": driver["ALLELE0_rare"],
            "driver_ALLELE1": driver["ALLELE1_rare"],
        })

    if removed_no_carrier > 0:
        print_info(f"Removed {removed_no_carrier} footprint-trait pairs (no variant with MAC > {min_carrier})")

    driver_df = pd.DataFrame(drivers)
    print_info(f"{len(driver_df)} trait-footprint pairs with driver variants")
    print_success(f"{driver_df['driver_variant'].nunique()} unique driver variants across {driver_df['footprint'].nunique()} footprints")
    return driver_df


def prepare_driver_bed(driver_df):
    """Create a BED-like DataFrame for driver variants (deduplicated).

    Parses DRAGEN variant IDs (1-based) into 0-based BED coordinates.
    Returns DataFrame with columns: Chromosome, Start, End, SNP, footprint_region, A1, A2
    """
    records = []
    seen = set()
    for _, row in driver_df.iterrows():
        variant_id = row["driver_variant"]
        footprint = row["footprint"]
        key = (variant_id, footprint)
        if key in seen:
            continue
        seen.add(key)

        # Parse DRAGEN ID: DRAGEN:chr16:88810611:G:A
        parts = variant_id.split(":")
        chrom = parts[1]
        pos_1based = int(parts[2])
        ref_allele = parts[3]
        alt_allele = parts[4]

        # Convert to 0-based BED coordinates
        pos_0based = pos_1based - 1

        records.append({
            "Chromosome": chrom,
            "Start": pos_0based,
            "End": pos_0based + len(ref_allele),
            "SNP": variant_id,
            "footprint_region": footprint,
            "A1": alt_allele,   # A1 = alt allele (REGENIE convention)
            "A2": ref_allele,   # A2 = ref allele
        })

    bed_df = pd.DataFrame(records)
    print_success(f"{len(bed_df)} unique variant-footprint pairs for FIMO analysis")
    return bed_df


# ---------------------------------------------------------------------------
# Sequence extraction (adapted from comvar_var2tfbs.py)
# ---------------------------------------------------------------------------

def parse_footprint_region(region_str):
    """Parse 'chr1:100-200' into (chrom, start, end)."""
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def extract_sequences(df, seq_extractor, ext_bp):
    """Extract ref/alt sequences for each variant.

    For each variant, the extracted window is footprint ± ext_bp.
    Determines which of A1/A2 is the reference allele by checking against
    the reference genome sequence at the variant position.
    """
    ref_seqs, alt_seqs, seq_names = [], [], []
    n_a1_ref, n_a2_ref, n_mismatch = 0, 0, 0

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Extracting sequences", ncols=60):
        fp_chrom, fp_start, fp_end = parse_footprint_region(row["footprint_region"])
        var_pos = int(row["Start"])  # 0-based variant position

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
            ref_base_a1 = seq_extractor.extract(
                Interval(row["Chromosome"], var_pos, var_pos + len(a1))
            )
            if ref_base_a1 == a1:
                ref_var, alt_var = a1, a2
                n_a1_ref += 1
            else:
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

        # Sequence name: variantID|chrom:pos(0-based):ref:alt|footprint_region
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
# FIMO motif scanning (shared with comvar_var2tfbs.py)
# ---------------------------------------------------------------------------

def run_fimo(motifs, fasta_path, ext_bp, threshold=0.0001):
    """Run FIMO and post-process hits.

    Returns a DataFrame of hits where the variant falls within the motif.
    """
    hits = run_fimo_engine(motifs, fasta_path, threshold=threshold)
    dfs = [h for h in hits if not h.empty]
    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)
    df["motif_id"] = df["motif_name"].str.split(" ").str[0] + df["strand"]
    df["motif_name"] = df["motif_name"].str.split(" ").str[1] + df["strand"]
    df["motif_info"] = (
        df["motif_id"] + "-" + df["motif_name"]
        + "_" + df["start"].astype(str)
        + "_" + df["end"].astype(str)
    )
    df["var-motif"] = df["sequence_name"] + "-" + df["motif_info"]

    # Parse sequence_name: variantID|chrom:pos:ref:alt|footprint_region
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
    """Classify TF binding changes based on ref/alt FIMO p-values."""
    ref_p = ref_alt_df["p-value_ref"].fillna(1)
    alt_p = ref_alt_df["p-value_alt"].fillna(1)

    conditions = [
        (ref_p != 1) & (alt_p != 1) & (ref_p > alt_p),
        (ref_p != 1) & (alt_p != 1) & (ref_p < alt_p),
        (ref_p != 1) & (alt_p != 1) & (ref_p == alt_p),
        (ref_p == 1) & (alt_p != 1),
        (ref_p != 1) & (alt_p == 1),
    ]
    choices = ["Increase", "Decrease", "Unchange", "Create", "Disrupt"]
    ref_alt_df["TF_change"] = np.select(conditions, choices, default="NoTFBS")
    return ref_alt_df


def deduplicate_tf_hits(ref_alt_df):
    """Keep only the best p-value per TF per variant."""
    results = []
    for rsid, var_tfbs in ref_alt_df.groupby("rsID"):
        foodie_id = var_tfbs["foodie_id"].iloc[0]

        ref_cols = ["TF", "motif_name_ref", "motif_id_ref", "p-value_ref", "start_ref", "end_ref"]
        var_ref = var_tfbs[ref_cols].dropna(subset=["p-value_ref"])
        if not var_ref.empty:
            var_ref = (
                var_ref.sort_values(by=["TF", "p-value_ref"])
                .drop_duplicates(subset=["TF"], keep="first")
                .reset_index(drop=True)
            )

        alt_cols = ["TF", "motif_name_alt", "motif_id_alt", "p-value_alt", "start_alt", "end_alt"]
        var_alt = var_tfbs[alt_cols].dropna(subset=["p-value_alt"])
        if not var_alt.empty:
            var_alt = (
                var_alt.sort_values(by=["TF", "p-value_alt"])
                .drop_duplicates(subset=["TF"], keep="first")
                .reset_index(drop=True)
            )

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

    print("Deduplicating: keeping best p-value per TF per variant ...")
    ref_alt_corrected = deduplicate_tf_hits(ref_alt)

    return ref_alt_corrected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Rare variant VAR2TFBS: identify driver variants and predict TF binding effects"
    )
    p.add_argument(
        "--burden-dir", required=True,
        help="Directory of burden test result Excel files (one per trait)",
    )
    p.add_argument(
        "--loo-file", required=True,
        help="Leave-one-out results CSV (all traits combined)",
    )
    p.add_argument("--ref-genome", required=True, help="Path to hg38.fa reference genome")
    p.add_argument("--jaspar-meme", required=True, help="Path to JASPAR MEME motif file")
    p.add_argument(
        "--out-dir", default="./rarevar_var2tfbs_results",
        help="Output directory (default: ./rarevar_var2tfbs_results)",
    )
    p.add_argument("--ext-bp", type=int, default=30, help="Sequence extension in bp (default: 30)")
    p.add_argument("--fimo-threshold", type=float, default=0.0001, help="FIMO p-value threshold (default: 0.0001)")
    p.add_argument("--sig-threshold", type=float, default=None, help="Burden test significance threshold (default: Bonferroni 0.05/N_footprints)")
    p.add_argument("--min-carrier", type=int, default=30, help="Minimum MAC for at least one variant in footprint (default: 30)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    fasta_dir = os.path.join(args.out_dir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    total_steps = 4

    print_header("Rare Variant VAR2TFBS Analysis")
    print_info(f"Burden test dir: {args.burden_dir}")
    print_info(f"LOO file:        {args.loo_file}")
    print_info(f"Reference genome: {args.ref_genome}")
    print_info(f"JASPAR motifs:   {args.jaspar_meme}")
    print_info(f"Significance:    p < {args.sig_threshold}")
    print_info(f"Min carrier:     MAC > {args.min_carrier}")
    print_info(f"Seq extension:   {args.ext_bp} bp")
    print_info(f"FIMO threshold:  {args.fimo_threshold}")

    # Step 1: Load significant footprints from burden tests
    print_step(1, total_steps, "Loading significant footprints from burden tests")
    sig_df = load_significant_footprints(args.burden_dir, args.sig_threshold)

    # Step 2: Identify driver variants from LOO analysis
    print_step(2, total_steps, "Identifying driver variants from leave-one-out analysis")
    loo_df = pd.read_csv(args.loo_file)
    print_info(f"LOO file: {len(loo_df):,} rows, {loo_df['trait'].nunique()} traits")
    driver_df = find_driver_variants(sig_df, loo_df, min_carrier=args.min_carrier)

    driver_out = os.path.join(args.out_dir, "driver_variants_summary.csv")
    driver_df.to_csv(driver_out, index=False)
    print_success(f"Saved driver summary: {driver_out}")

    # Step 3: Extract ref/alt sequences
    print_step(3, total_steps, "Extracting ref/alt sequences")
    bed_df = prepare_driver_bed(driver_df)

    print_info("Loading reference genome ...")
    seq_extractor = FastaStringExtractor(args.ref_genome)
    print_info("Loading JASPAR motifs ...")
    motifs = read_meme(args.jaspar_meme)

    bed_df = extract_sequences(bed_df, seq_extractor, args.ext_bp)
    bed_df = bed_df.drop_duplicates(subset=["seq_name"]).reset_index(drop=True)
    seq_extractor.close()

    ref_fasta = os.path.join(fasta_dir, "K562_rarevar_ref_seqExt{ext}bp.fa".format(ext=args.ext_bp))
    alt_fasta = os.path.join(fasta_dir, "K562_rarevar_alt_seqExt{ext}bp.fa".format(ext=args.ext_bp))
    write_fasta(bed_df, ref_fasta, allele="ref")
    write_fasta(bed_df, alt_fasta, allele="alt")
    print_success(f"Wrote {ref_fasta}")
    print_success(f"Wrote {alt_fasta}")

    # Step 4: FIMO motif scanning + TF change classification
    print_step(4, total_steps, "FIMO motif scanning + TF change classification")
    result = var2tfbs(motifs, ref_fasta, alt_fasta, args.ext_bp, args.fimo_threshold)

    if result.empty:
        print_info("No motif hits found.")
        return

    result["cell"] = "K562"
    out_path = os.path.join(args.out_dir, "K562_rarevar_var2tfbs.csv")
    result.to_csv(out_path, index=False)

    print()
    print("=" * 62)
    print(f"  Results: {len(result):,} variant-TF pairs across {result['rsID'].nunique()} driver variants")
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
