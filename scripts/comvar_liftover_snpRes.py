#!/usr/bin/env python
"""
Liftover GWFM common variant snpRes files from hg19 to hg38.

Adds Chromosome_hg38, Start_hg38, End_hg38 columns to each snpRes file.
Since all snpRes files share the same common variant positions, liftover is
performed once using the first file, then applied to all.

Usage:
    python scripts/comvar_liftover_snpRes.py \
        --snpres-dir data/GWFM_erythroids/snpRes \
        --chain data/reference/hg19ToHg38.over.chain.gz \
        --out-dir data/GWFM_erythroids/snpRes_hg38
"""

import argparse
import glob
import os
import subprocess
import tempfile
from functools import partial
from multiprocessing import Pool, cpu_count

import pandas as pd
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser(description="Liftover snpRes files from hg19 to hg38")
    p.add_argument("--snpres-dir", required=True, help="Directory containing .snpRes files")
    p.add_argument("--chain", required=True, help="Path to hg19ToHg38.over.chain.gz")
    p.add_argument("--out-dir", required=True, help="Output directory for lifted-over files")
    p.add_argument("--workers", type=int, default=0,
                   help="Number of parallel workers (default: number of CPU cores)")
    return p.parse_args()


def process_one_trait(snpres_file, out_dir, mapped_dict):
    """Process a single snpRes file: read, add hg38 columns, write."""
    trait = os.path.basename(snpres_file).replace(".snpRes", "")
    out_path = os.path.join(out_dir, f"{trait}.snpRes")

    df = pd.read_csv(snpres_file, sep=r"\s+")
    df["Chromosome_hg38"] = df["Name"].map(lambda x: mapped_dict.get(x, {}).get("Chromosome_hg38", ""))
    df["Start_hg38"] = df["Name"].map(lambda x: mapped_dict.get(x, {}).get("Start_hg38", ""))
    df["End_hg38"] = df["Name"].map(lambda x: mapped_dict.get(x, {}).get("End_hg38", ""))

    df.to_csv(out_path, sep="\t", index=False)
    n_mapped = df["Chromosome_hg38"].ne("").sum()
    return trait, n_mapped, len(df)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    snpres_files = sorted(glob.glob(os.path.join(args.snpres_dir, "*.snpRes")))
    if not snpres_files:
        print("ERROR: No .snpRes files found in", args.snpres_dir)
        return

    n_workers = args.workers if args.workers > 0 else cpu_count()
    print(f"Found {len(snpres_files)} snpRes files, using {n_workers} workers")

    # --- Step 1: Read positions from first file to build liftover input ---
    ref_file = snpres_files[0]
    ref_name = os.path.basename(ref_file)
    print(f"\n[1/3] Reading positions from {ref_name} ...")
    ref_df = pd.read_csv(ref_file, sep=r"\s+", usecols=["Name", "Chrom", "Position"])
    print(f"  {len(ref_df):,} variants")

    # Create BED: chr, start (0-based), end (1-based), name
    bed_df = pd.DataFrame()
    bed_df["chrom"] = "chr" + ref_df["Chrom"].astype(str)
    bed_df["start"] = ref_df["Position"].astype(int) - 1  # hg19 Position is 1-based; BED is 0-based
    bed_df["end"] = ref_df["Position"].astype(int)
    bed_df["name"] = ref_df["Name"]

    bed_in = tempfile.mktemp(suffix=".bed")
    bed_df.to_csv(bed_in, sep="\t", header=False, index=False)
    print(f"  Wrote BED with {len(bed_df):,} entries")

    # --- Step 2: Run liftOver ---
    print("\n[2/3] Running liftOver ...")
    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as f_out, \
         tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as f_unmapped:
        bed_out = f_out.name
        bed_unmapped = f_unmapped.name

    cmd = ["liftOver", bed_in, args.chain, bed_out, bed_unmapped]
    subprocess.run(cmd, check=True)

    # Count unmapped
    n_unmapped = 0
    with open(bed_unmapped) as f:
        for line in f:
            if not line.startswith("#"):
                n_unmapped += 1

    # Read mapped results
    mapped = pd.read_csv(
        bed_out, sep="\t", header=None,
        names=["Chromosome_hg38", "Start_hg38", "End_hg38", "Name"],
    )
    print(f"  Mapped:   {len(mapped):,} variants")
    print(f"  Unmapped: {n_unmapped:,} variants")

    # Build lookup: Name -> (Chromosome_hg38, Start_hg38, End_hg38)
    mapped_dict = mapped.set_index("Name")[
        ["Chromosome_hg38", "Start_hg38", "End_hg38"]
    ].to_dict("index")

    # Clean up temp files
    os.unlink(bed_in)
    os.unlink(bed_out)
    os.unlink(bed_unmapped)

    # --- Step 3: Add hg38 columns to each snpRes file (parallel) ---
    print(f"\n[3/3] Writing lifted-over snpRes files to {args.out_dir} ...")
    worker_fn = partial(process_one_trait, out_dir=args.out_dir, mapped_dict=mapped_dict)

    with Pool(n_workers) as pool:
        results = list(tqdm(
            pool.imap_unordered(worker_fn, snpres_files),
            total=len(snpres_files),
            desc="Writing",
            unit="file",
            ncols=60,
        ))

    for trait, n_mapped, n_total in sorted(results):
        print(f"  {trait}: {n_mapped:,}/{n_total:,} mapped")

    print("\nDone!")


if __name__ == "__main__":
    main()
