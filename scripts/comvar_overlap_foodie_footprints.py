#!/usr/bin/env python3
"""
Overlap GWFM common variants with FOODIE footprint BED files.

For each footprint BED file, outputs per-trait CSVs of overlapping common
variants with PIP, PEP, PEP_cs, and CS_id annotations. Accepts both
credible set CSV and snpRes (tab-delimited) formats with auto-detection.

Usage:
    # Using credible set files:
    python scripts/comvar_overlap_foodie_footprints.py \
        --snp-dir data/GWFM_erythroids/credible_set_snpRes \
        --footprint-dir data/FOODIE_footprints \
        --out-dir results/comvar_footprint_overlap

    # Using full snpRes_hg38 files:
    python scripts/comvar_overlap_foodie_footprints.py \
        --snp-dir data/GWFM_erythroids/snpRes_hg38 \
        --snp-suffix .snpRes \
        --footprint-dir data/FOODIE_footprints \
        --out-dir results/comvar_footprint_overlap_snpRes

Requirements:
    - bedtools (must be on PATH)
    - pandas
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import time

import pandas as pd

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


# ── Pretty printing helpers ──────────────────────────────────────────────────

def print_header(text):
    width = 62
    print()
    print(f"{'=' * width}")
    print(f"  {text}")
    print(f"{'=' * width}")


def print_step(step_num, total, text):
    print(f"\n[Step {step_num}/{total}] {text}")
    print(f"{'-' * 50}")


def print_success(text):
    print(f"  [OK] {text}")


def print_info(text):
    print(f"  -> {text}")


def progress_bar(iterable, desc="Processing", total=None):
    """Wrap an iterable with tqdm if available, otherwise print progress."""
    if HAS_TQDM:
        return tqdm(iterable, desc=f"  {desc}", total=total, bar_format="  {l_bar}{bar:30}{r_bar}")
    else:
        items = list(iterable)
        n = len(items)
        for i, item in enumerate(items):
            pct = (i + 1) / n * 100
            bar_len = 30
            filled = int(bar_len * (i + 1) / n)
            bar = "#" * filled + "-" * (bar_len - filled)
            print(f"\r  {desc}: [{bar}] {pct:5.1f}% ({i+1}/{n})", end="", flush=True)
            yield item
        print()


# ── Core functions ───────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="VAR2TFBS Step 1: Overlap GWFM common variants with FOODIE footprint BED files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using credible set files (filtered by PIP/PEP):
  python scripts/comvar_overlap_foodie_footprints.py \\
      --snp-dir data/GWFM_erythroids/credible_set_snpRes \\
      --footprint-dir data/FOODIE_footprints

  # Using full snpRes_hg38 files (all ~13M common variants):
  python scripts/comvar_overlap_foodie_footprints.py \\
      --snp-dir data/GWFM_erythroids/snpRes_hg38 \\
      --snp-suffix .snpRes \\
      --footprint-dir data/FOODIE_footprints \\
      --out-dir results/comvar_footprint_overlap_snpRes
        """,
    )
    parser.add_argument(
        "--snp-dir",
        required=True,
        help="Directory containing GWFM variant files (CSV or tab-delimited snpRes).",
    )
    parser.add_argument(
        "--footprint-dir",
        required=True,
        help="Directory containing FOODIE footprint BED files.",
    )
    parser.add_argument(
        "--out-dir",
        default="./results/comvar_footprint_overlap",
        help="Output directory for common variant overlap results (default: ./results/comvar_footprint_overlap).",
    )
    parser.add_argument(
        "--pip-threshold",
        type=float,
        default=0,
        help="Minimum max-PIP (across all traits) to include in output (default: 0).",
    )
    parser.add_argument(
        "--snp-suffix",
        default="_credible_set_hg38.csv",
        help="Suffix to strip from SNP filenames to get trait name (default: _credible_set_hg38.csv; use .snpRes for snpRes files).",
    )
    parser.add_argument(
        "--lcs-dir",
        default=None,
        help="Directory containing .lcs files for PEP_cs annotation (default: auto-detect from snp-dir).",
    )
    parser.add_argument(
        "--bed-pattern",
        default="*.bed",
        help="Glob pattern for footprint BED files (default: *.bed).",
    )
    return parser.parse_args()


def check_bedtools():
    """Verify bedtools is available."""
    try:
        result = subprocess.run(
            ["bedtools", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        version = result.stdout.strip()
        print_success(f"bedtools found: {version}")
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("Error: bedtools not found. Please install bedtools and ensure it is on PATH.", file=sys.stderr)
        sys.exit(1)


REQUIRED_SNP_COLUMNS = {"Chromosome", "Start", "End", "SNP", "PIP", "PEP", "A1", "A2", "freq"}

# snpRes column names -> standard names used by this script
COLUMN_ALIASES = {
    "Name": "SNP",
    "Chromosome_hg38": "Chromosome",
    "Start_hg38": "Start",
    "End_hg38": "End",
    "A1Frq": "freq",
}


def _detect_sep(filepath):
    """Auto-detect delimiter (tab or comma) by reading the first line."""
    with open(filepath) as f:
        first_line = f.readline()
    return "\t" if "\t" in first_line else ","


def validate_snp_file(snp_file):
    """Check that a SNP CSV file has all required columns (after alias mapping)."""
    try:
        header = pd.read_csv(snp_file, nrows=0, sep=_detect_sep(snp_file))
    except Exception as e:
        print(f"Error: Cannot read {snp_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Apply column aliases to check
    mapped_cols = {COLUMN_ALIASES.get(c, c) for c in header.columns}
    missing = REQUIRED_SNP_COLUMNS - mapped_cols
    if missing:
        print(f"Error: {os.path.basename(snp_file)} is missing required columns: {', '.join(sorted(missing))}", file=sys.stderr)
        print(f"  Required columns: {', '.join(sorted(REQUIRED_SNP_COLUMNS))}", file=sys.stderr)
        print(f"  Found columns:    {', '.join(header.columns[:10])}{'...' if len(header.columns) > 10 else ''}", file=sys.stderr)
        sys.exit(1)


def validate_bed_file(bed_file):
    """Check that a BED file has at least 3 tab-separated columns."""
    with open(bed_file) as f:
        first_line = f.readline().strip()
    if not first_line:
        print(f"Error: {os.path.basename(bed_file)} is empty.", file=sys.stderr)
        sys.exit(1)
    fields = first_line.split("\t")
    if len(fields) < 3:
        print(f"Error: {os.path.basename(bed_file)} does not appear to be a valid BED file.", file=sys.stderr)
        print(f"  Expected tab-separated with >=3 columns (chr, start, end), got {len(fields)} field(s).", file=sys.stderr)
        print(f"  First line: {first_line[:100]}", file=sys.stderr)
        sys.exit(1)


def discover_files(snp_dir, snp_suffix, footprint_dir, bed_pattern):
    """Find SNP trait files and footprint BED files."""
    snp_files = sorted(glob.glob(os.path.join(snp_dir, f"*{snp_suffix}")))
    if not snp_files:
        print(f"Error: No SNP files matching *{snp_suffix} in {snp_dir}", file=sys.stderr)
        print(f"  Check that --snp-dir and --snp-suffix are correct.", file=sys.stderr)
        sys.exit(1)

    bed_files = sorted(glob.glob(os.path.join(footprint_dir, bed_pattern)))
    if not bed_files:
        print(f"Error: No BED files matching {bed_pattern} in {footprint_dir}", file=sys.stderr)
        print(f"  Check that --footprint-dir and --bed-pattern are correct.", file=sys.stderr)
        sys.exit(1)

    # Validate input file formats
    print_info("Validating input file formats...")
    validate_snp_file(snp_files[0])
    for bed_file in bed_files:
        validate_bed_file(bed_file)
    print_success("All input files passed format checks.")

    traits = []
    for f in snp_files:
        trait = os.path.basename(f).replace(snp_suffix, "")
        traits.append(trait)

    print_info(f"SNP directory:       {snp_dir}")
    print_info(f"Footprint directory: {footprint_dir}")
    print_info(f"Trait files ({len(snp_files)}):    {', '.join(traits)}")
    print_info(f"Footprint BEDs ({len(bed_files)}): {', '.join(os.path.basename(b) for b in bed_files)}")
    return snp_files, traits, bed_files


def make_variant_bed(snp_files, tmpdir):
    """Merge all trait CSVs into a deduplicated 4-column BED (chr, start, end, SNP).

    When using credible set files, each trait may have a different subset of
    variants (filtered by PIP threshold), so we must take the union across all
    trait files to avoid missing variants that only appear in some traits.
    """
    bed_path = os.path.join(tmpdir, "variants.bed")
    print_info(f"Reading variants from {len(snp_files)} trait file(s)...")
    frames = []
    for snp_file in snp_files:
        df = pd.read_csv(snp_file, sep=_detect_sep(snp_file))
        df = df.rename(columns=COLUMN_ALIASES)
        df = df[["Chromosome", "Start", "End", "SNP"]].dropna(subset=["Chromosome", "Start", "End"])
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["SNP"])
    merged["Start"] = merged["Start"].astype(int)
    merged["End"] = merged["End"].astype(int)
    merged = merged.sort_values(["Chromosome", "Start"])
    merged.to_csv(
        bed_path, sep="\t", header=False, index=False
    )
    print_success(f"{len(merged):,} unique variants loaded into BED format.")
    return bed_path


def run_bedtools_intersect(variant_bed, footprint_bed, output_path):
    """Run bedtools intersect and return overlapping variant SNP IDs."""
    result = subprocess.run(
        ["bedtools", "intersect", "-a", variant_bed, "-b", footprint_bed, "-wa", "-wb"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip()
        if stderr and "Error" in stderr:
            print(f"  Warning: bedtools reported: {stderr[:200]}", file=sys.stderr)

    with open(output_path, "w") as f:
        f.write(result.stdout)

    overlap_records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) >= 8:
            snp_id = fields[3]
            fp_name = fields[7] if len(fields) > 7 else f"{fields[4]}:{fields[5]}-{fields[6]}"
            overlap_records.append((snp_id, fp_name))

    if not overlap_records:
        return set(), {}

    snp_ids = set(r[0] for r in overlap_records)

    fp_map = {}
    for snp_id, fp_name in overlap_records:
        fp_map.setdefault(snp_id, []).append(fp_name)
    fp_map = {k: ";".join(sorted(set(v))) for k, v in fp_map.items()}

    return snp_ids, fp_map


def _load_lcs_annotations(lcs_dir, trait):
    """Load LCS PEP_cs and CS_id mappings for a trait. Returns DataFrame with PEP, CS columns."""
    lcs_file = os.path.join(lcs_dir, f"{trait}.lcs")
    if not os.path.exists(lcs_file):
        return pd.DataFrame(columns=["PEP", "CS"], dtype=float)
    lcs = pd.read_csv(lcs_file, sep=r"\s+")
    lcs["SNP"] = lcs["SNP"].str.split(",")
    lcs = lcs.explode("SNP")
    lcs = lcs.sort_values("PEP", ascending=False).drop_duplicates(subset="SNP", keep="first")
    return lcs.set_index("SNP")[["PEP", "CS"]]


def load_trait_pips(snp_files, traits, snp_ids, lcs_dir=None):
    """Load PIP, PEP, and PEP_cs values for overlapping SNPs from all trait files."""
    pip_data = {}
    pep_data = {}
    pep_cs_data = {}
    cs_id_data = {}
    base_info = None
    n = len(traits)

    for i, (snp_file, trait) in enumerate(zip(snp_files, traits)):
        pct = (i + 1) / n * 100
        filled = int(30 * (i + 1) / n)
        bar = "#" * filled + "-" * (30 - filled)
        print(f"\r  Loading trait PIPs: [{bar}] {pct:5.1f}% ({i+1}/{n}) {trait:<10}", end="", flush=True)

        df = pd.read_csv(snp_file, sep=_detect_sep(snp_file))
        df = df.rename(columns=COLUMN_ALIASES)
        df = df[["SNP", "PIP", "PEP", "Chromosome", "Start", "End", "A1", "A2", "freq"]]
        # Deduplicate SNPs (same SNP can appear in multiple credible sets);
        # keep the entry from the credible set with the highest PEP.
        df = df.sort_values("PEP", ascending=False).drop_duplicates(subset="SNP", keep="first")
        subset = df[df["SNP"].isin(snp_ids)].set_index("SNP")
        pip_data[trait] = subset["PIP"]
        pep_data[trait] = subset["PEP"]

        # Load PEP_cs and CS_id from LCS
        if lcs_dir:
            lcs_annot = _load_lcs_annotations(lcs_dir, trait)
            pep_cs_data[trait] = lcs_annot["PEP"].reindex(subset.index)
            cs_id_data[trait] = lcs_annot["CS"].reindex(subset.index)
        else:
            pep_cs_data[trait] = pd.Series(dtype=float)
            cs_id_data[trait] = pd.Series(dtype=float)

        # Accumulate base coordinate info from all traits (some variants
        # only appear in certain traits' credible sets)
        new_info = subset[["Chromosome", "Start", "End", "A1", "A2", "freq"]].copy()
        new_info["Start"] = new_info["Start"].astype(int)
        new_info["End"] = new_info["End"].astype(int)
        if base_info is None:
            base_info = new_info
        else:
            # Add variants not yet in base_info
            missing = new_info.loc[~new_info.index.isin(base_info.index)]
            if len(missing) > 0:
                base_info = pd.concat([base_info, missing])

    print()
    return base_info, pip_data, pep_data, pep_cs_data, cs_id_data


def build_result_table(base_info, pip_data, traits, fp_map, pip_threshold):
    """Assemble the final combined result DataFrame."""
    result = base_info.copy()
    for trait in traits:
        result[f"PIP_{trait}"] = pip_data[trait]

    pip_cols = [f"PIP_{t}" for t in traits]
    result["max_PIP"] = result[pip_cols].max(axis=1)
    result["footprint_region"] = result.index.map(lambda x: fp_map.get(x, ""))
    result = result.sort_values("max_PIP", ascending=False)

    if pip_threshold > 0:
        before = len(result)
        result = result[result["max_PIP"] >= pip_threshold]
        print_info(f"PIP filter >= {pip_threshold}: {before:,} -> {len(result):,} variants")

    return result


def build_trait_table(base_info, pip_data, pep_data, pep_cs_data, cs_id_data, trait, fp_map, pip_threshold):
    """Assemble a per-trait result DataFrame with PIP, PEP, PEP_cs, CS_id."""
    result = base_info.copy()
    result["PIP"] = pip_data[trait]
    result["PEP"] = pep_data[trait]
    if trait in pep_cs_data and len(pep_cs_data[trait]) > 0:
        result["PEP_cs"] = pep_cs_data[trait]
    else:
        result["PEP_cs"] = float("nan")
    if trait in cs_id_data and len(cs_id_data[trait]) > 0:
        result["CS_id"] = cs_id_data[trait]
    else:
        result["CS_id"] = float("nan")
    result["footprint_region"] = result.index.map(lambda x: fp_map.get(x, ""))
    result = result.dropna(subset=["PIP"])
    result = result.sort_values("PIP", ascending=False)

    if pip_threshold > 0:
        result = result[result["PIP"] >= pip_threshold]

    return result


def main():
    start_time = time.time()
    args = parse_args()

    print_header("VAR2TFBS Step 1: Overlap GWFM Common Variants with FOODIE Footprints")

    total_steps = 4

    # Step 1: Validate environment and discover files
    print_step(1, total_steps, "Checking environment and discovering files")
    check_bedtools()
    snp_files, traits, bed_files = discover_files(
        args.snp_dir, args.snp_suffix, args.footprint_dir, args.bed_pattern
    )
    os.makedirs(args.out_dir, exist_ok=True)
    print_info(f"Output directory: {os.path.abspath(args.out_dir)}")
    if args.pip_threshold > 0:
        print_info(f"PIP threshold: >= {args.pip_threshold}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 2: Convert variants to BED and run intersections
        print_step(2, total_steps, "Intersecting variants with footprints")
        variant_bed = make_variant_bed(snp_files, tmpdir)

        all_snp_ids = set()
        footprint_results = {}

        for bed_file in bed_files:
            fp_label = os.path.basename(bed_file).replace(".bed", "")
            print_info(f"Running bedtools intersect with {fp_label}...")
            raw_overlap = os.path.join(tmpdir, f"overlap_{fp_label}.bed")
            snp_ids, fp_map = run_bedtools_intersect(variant_bed, bed_file, raw_overlap)
            print_success(f"{fp_label}: {len(snp_ids):,} overlapping variants found.")
            all_snp_ids |= snp_ids
            footprint_results[fp_label] = (snp_ids, fp_map)

        if not all_snp_ids:
            print("\n  No overlapping variants found across any footprint. Nothing to output.")
            return

        print_info(f"Total unique overlapping variants: {len(all_snp_ids):,}")

        # Step 3: Load PIP/PEP values from all trait files
        print_step(3, total_steps, f"Loading PIP/PEP values from {len(traits)} trait files")
        lcs_dir = args.lcs_dir
        base_info, pip_data, pep_data, pep_cs_data, cs_id_data = load_trait_pips(
            snp_files, traits, all_snp_ids, lcs_dir=lcs_dir
        )

        # Step 4: Write output files (per trait per footprint)
        print_step(4, total_steps, "Writing output files")

        for fp_label, (snp_ids, fp_map) in footprint_results.items():
            if not snp_ids:
                print_info(f"{fp_label}: no overlapping variants, skipping.")
                continue

            print_info(f"Processing {fp_label}...")
            fp_base = base_info.loc[base_info.index.isin(snp_ids)].copy()

            # Per-trait CSV files
            fp_dir = os.path.join(args.out_dir, fp_label)
            os.makedirs(fp_dir, exist_ok=True)

            all_trait_snps = set()
            for trait in traits:
                fp_pips = {trait: pip_data[trait].loc[pip_data[trait].index.isin(snp_ids)]}
                fp_peps = {trait: pep_data[trait].loc[pep_data[trait].index.isin(snp_ids)]}
                fp_pep_cs = {trait: pep_cs_data[trait].loc[pep_cs_data[trait].index.isin(snp_ids)] if trait in pep_cs_data and len(pep_cs_data[trait]) > 0 else pd.Series(dtype=float)}
                fp_cs_id = {trait: cs_id_data[trait].loc[cs_id_data[trait].index.isin(snp_ids)] if trait in cs_id_data and len(cs_id_data[trait]) > 0 else pd.Series(dtype=float)}
                result = build_trait_table(fp_base, fp_pips, fp_peps, fp_pep_cs, fp_cs_id, trait, fp_map, args.pip_threshold)

                if len(result) > 0:
                    out_path = os.path.join(fp_dir, f"{trait}_{fp_label}.csv")
                    result.to_csv(out_path)
                    all_trait_snps |= set(result.index)

            print_success(f"{fp_dir}/  ({len(traits)} traits, {len(all_trait_snps):,} unique variants)")

            # Combined deduplicated BED across all traits, annotated with footprint region
            bed_out = os.path.join(args.out_dir, f"GWFM_variants_in_{fp_label}.bed")
            bed_snps = fp_base.loc[fp_base.index.isin(all_trait_snps)]
            bed_df = bed_snps[["Chromosome", "Start", "End"]].copy()
            bed_df.insert(3, "SNP", bed_snps.index)
            bed_df["footprint_region"] = bed_df["SNP"].map(lambda x: fp_map.get(x, ""))
            bed_df = bed_df.drop_duplicates().sort_values(["Chromosome", "Start"])
            bed_df.to_csv(bed_out, sep="\t", header=False, index=False)
            print_success(f"{bed_out}  ({len(bed_df):,} unique variants)")

        # Summary
        elapsed = time.time() - start_time
        print_header("SUMMARY")
        print(f"  Traits ({len(traits)}):  {', '.join(traits)}")
        print(f"  Total unique variants in any footprint: {len(all_snp_ids):,}")
        for fp_label, (snp_ids, _) in footprint_results.items():
            print(f"    {fp_label}: {len(snp_ids):,} variants")
        if len(footprint_results) > 1:
            shared = set.intersection(*(s for s, _ in footprint_results.values()))
            print(f"    Shared across all: {len(shared):,} variants")
        if args.pip_threshold > 0:
            print(f"  PIP threshold: >= {args.pip_threshold}")
        print(f"  Output: {os.path.abspath(args.out_dir)}")
        print(f"  Completed in {elapsed:.1f}s")
        print()


if __name__ == "__main__":
    main()
