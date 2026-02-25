#!/usr/bin/env python
"""
Score VAR2TFBS variants using AlphaGenome.

Extracts unique variants from VAR2TFBS output, resolves coordinates/alleles,
and scores each variant using AlphaGenome's variant scoring API with
recommended CenterMaskScorer configurations. Supports both common and rare
variant ID formats, with checkpoint/resume for large batches.

Usage (common variants):
    python scripts/alphag_score_variants.py \
        --var2tfbs results/comvar_var2tfbs_results/K562_var2tfbs.csv \
        --allele-src results/comvar_footprint_overlap_credible/K562.merged.hg38 \
        --ref-genome data/reference/hg38.fa \
        --api-key $ALPHAGENOME_API_KEY \
        --cell K562 \
        --prefix comvar \
        --out-dir results/alphag_scores

Usage (rare variants):
    python scripts/alphag_score_variants.py \
        --var2tfbs results/rarevar_var2tfbs_results/K562_rarevar_var2tfbs.csv \
        --api-key $ALPHAGENOME_API_KEY \
        --cell K562 \
        --prefix rarevar \
        --out-dir results/alphag_scores
"""

import argparse
import glob as _glob
import os
import re
import sys

import pandas as pd
import pyfaidx
from kipoiseq import Interval
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


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
# Constants
# ---------------------------------------------------------------------------

CELL_ONTOLOGY = {
    "K562": "EFO:0002067",
    "GM12878": "EFO:0002784",
}

SEQ_LENGTH_MAP = {
    "16KB": dna_client.SEQUENCE_LENGTH_16KB,
    "100KB": dna_client.SEQUENCE_LENGTH_100KB,
    "500KB": dna_client.SEQUENCE_LENGTH_500KB,
    "1MB": dna_client.SEQUENCE_LENGTH_1MB,
}

# Rare variant ID pattern: PREFIX:chrCHROM:POS:REF:ALT
_RARE_PATTERN = re.compile(
    r"^[A-Za-z]+:(chr(?:X|Y|M|\d+)):(\d+):([ACGTNacgtn]+):([ACGTNacgtn]+)$"
)


# ---------------------------------------------------------------------------
# Reference genome helper (same as comvar_var2tfbs.py)
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
# Variant parsing
# ---------------------------------------------------------------------------

def parse_rare_variant_id(rsid):
    """Parse rare variant ID (DRAGEN:chrCHROM:POS:REF:ALT) -> genome.Variant.

    Position is 1-based, alleles are already ref/alt.
    Returns genome.Variant or None if format doesn't match.
    """
    match = _RARE_PATTERN.match(rsid)
    if match is None:
        return None
    chrom, pos, ref, alt = match.groups()
    return genome.Variant(
        chromosome=chrom,
        position=int(pos),
        reference_bases=ref.upper(),
        alternate_bases=alt.upper(),
    )


def load_allele_info(allele_src):
    """Load variant coordinate, allele, and PIP info from allele source.

    Returns:
        allele_info: dict SNP -> {Chromosome, Start, A1, A2}
        max_pip: dict SNP -> max PIP across all traits (None if PIP not available)
    """
    usecols = ["SNP", "Chromosome", "Start", "A1", "A2"]
    if os.path.isdir(allele_src):
        csv_files = sorted(_glob.glob(os.path.join(allele_src, "*.csv")))
        # Check if PIP column exists
        sample_cols = pd.read_csv(csv_files[0], nrows=0).columns
        has_pip = "PIP" in sample_cols
        if has_pip:
            usecols = usecols + ["PIP"]
        dfs = []
        for f in csv_files:
            dfs.append(pd.read_csv(f, usecols=usecols))
        df = pd.concat(dfs, ignore_index=True)
    else:
        sample_cols = pd.read_csv(allele_src, nrows=0).columns
        has_pip = "PIP" in sample_cols
        if has_pip:
            usecols = usecols + ["PIP"]
        df = pd.read_csv(allele_src, usecols=usecols)

    # Max PIP per variant across all traits
    max_pip = None
    if has_pip:
        max_pip = df.groupby("SNP")["PIP"].max().to_dict()

    df = df.drop_duplicates(subset="SNP", keep="first")
    allele_info = df.set_index("SNP")[["Chromosome", "Start", "A1", "A2"]].to_dict("index")
    return allele_info, max_pip


def resolve_ref_alt(chrom, pos_0based, a1, a2, seq_extractor):
    """Determine ref/alt from A1/A2 by checking reference genome."""
    ref_base = seq_extractor.extract(
        Interval(chrom, pos_0based, pos_0based + len(a2))
    )
    if ref_base == a2:
        return a2, a1
    ref_base_a1 = seq_extractor.extract(
        Interval(chrom, pos_0based, pos_0based + len(a1))
    )
    if ref_base_a1 == a1:
        return a1, a2
    return a2, a1


def extract_unique_variants(var2tfbs_df, allele_info=None, seq_extractor=None):
    """Extract unique variants from VAR2TFBS and resolve to genome.Variant objects.

    Returns:
        resolved: list of (rsid_str, genome.Variant)
        skipped: list of rsid_str that could not be resolved
    """
    unique_rsids = var2tfbs_df["rsID"].unique()
    resolved = []
    skipped = []

    for rsid in tqdm(unique_rsids, desc="Resolving variants", ncols=60):
        # Try rare variant format first
        variant = parse_rare_variant_id(rsid)
        if variant is not None:
            resolved.append((rsid, variant))
            continue

        # Common variant: need allele_info lookup
        if allele_info is None or rsid not in allele_info:
            skipped.append(rsid)
            continue

        info = allele_info[rsid]
        chrom = info["Chromosome"]
        pos_0based = int(info["Start"])
        a1, a2 = info["A1"], info["A2"]

        if seq_extractor is not None:
            ref, alt = resolve_ref_alt(chrom, pos_0based, a1, a2, seq_extractor)
        else:
            ref, alt = a2, a1

        variant = genome.Variant(
            chromosome=chrom,
            position=pos_0based + 1,  # genome.Variant uses 1-based position
            reference_bases=ref,
            alternate_bases=alt,
        )
        resolved.append((rsid, variant))

    return resolved, skipped


# ---------------------------------------------------------------------------
# Scorer selection
# ---------------------------------------------------------------------------

def build_scorers(output_type_names):
    """Build variant scorers from output type names."""
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected = []
    for name in output_type_names:
        key = name.upper()
        if key in all_scorers:
            selected.append(all_scorers[key])
        else:
            raise ValueError(
                f"Unknown output type: {name}. "
                f"Available: {list(all_scorers.keys())}"
            )
    return selected


# ---------------------------------------------------------------------------
# Scoring with checkpoints
# ---------------------------------------------------------------------------

def score_variants_with_checkpoints(
    model, rsid_variant_pairs, scorers, seq_length, batch_size, out_dir, cell, prefix
):
    """Score variants one at a time with batch checkpointing.

    Saves a checkpoint TSV every batch_size variants. On re-run, resumes
    from the last checkpoint.
    """
    checkpoint_dir = os.path.join(out_dir, "checkpoints")
    os.makedirs(checkpoint_dir, exist_ok=True)

    # Load existing checkpoints
    all_dfs = []
    completed_rsids = set()
    cp_pattern = os.path.join(checkpoint_dir, f"{cell}_{prefix}_batch_*.tsv")
    checkpoint_files = sorted(_glob.glob(cp_pattern))
    for cp_file in checkpoint_files:
        cp_df = pd.read_csv(cp_file, sep="\t", low_memory=False)
        all_dfs.append(cp_df)
        completed_rsids.update(cp_df["rsID"].unique())

    if completed_rsids:
        print_info(f"Resuming: {len(completed_rsids)} variants already scored")

    remaining = [
        (rsid, var) for rsid, var in rsid_variant_pairs
        if rsid not in completed_rsids
    ]

    if not remaining:
        print_info("All variants already scored")
        if all_dfs:
            return pd.concat(all_dfs, ignore_index=True)
        return pd.DataFrame()

    # Score in batches
    batch_idx = len(checkpoint_files)
    batch_results = []
    batch_rsids = []

    for i, (rsid, variant) in enumerate(
        tqdm(remaining, desc="Scoring variants", ncols=60)
    ):
        interval = variant.reference_interval.resize(seq_length)

        try:
            result = model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
            batch_results.append(result)
            batch_rsids.append((rsid, variant))
        except Exception as e:
            print_info(f"ERROR scoring {rsid}: {e}")
            continue

        # Save checkpoint at batch boundary
        if len(batch_results) >= batch_size or i == len(remaining) - 1:
            if batch_results:
                tidy_df = variant_scorers.tidy_scores(batch_results)
                if tidy_df is not None and not tidy_df.empty:
                    # Map variant_id back to rsID
                    vid_to_rsid = {}
                    for r, v in batch_rsids:
                        vid_to_rsid[str(v)] = r
                    tidy_df["variant_id"] = tidy_df["variant_id"].astype(str)
                    tidy_df["rsID"] = tidy_df["variant_id"].map(vid_to_rsid)

                    os.makedirs(checkpoint_dir, exist_ok=True)
                    cp_path = os.path.join(
                        checkpoint_dir, f"{cell}_{prefix}_batch_{batch_idx:04d}.tsv"
                    )
                    tidy_df.to_csv(cp_path, sep="\t", index=False)
                    all_dfs.append(tidy_df)
                    batch_idx += 1

                    print_info(
                        f"Checkpoint {batch_idx}: "
                        f"{len(batch_rsids)} variants, "
                        f"{len(tidy_df):,} score rows"
                    )

                batch_results = []
                batch_rsids = []

    if not all_dfs:
        return pd.DataFrame()
    return pd.concat(all_dfs, ignore_index=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Score VAR2TFBS variants using AlphaGenome"
    )
    p.add_argument(
        "--var2tfbs", required=True,
        help="VAR2TFBS output CSV",
    )
    p.add_argument(
        "--api-key", default=None,
        help="AlphaGenome API key (or set ALPHAGENOME_API_KEY env variable)",
    )
    p.add_argument(
        "--cell", default="K562", choices=["K562", "GM12878"],
        help="Cell type (default: K562)",
    )
    p.add_argument(
        "--allele-src", default=None,
        help="Directory of per-trait overlap CSVs or single CSV with "
             "SNP, Chromosome, Start, A1, A2 columns (for common variants)",
    )
    p.add_argument(
        "--ref-genome", default=None,
        help="Path to hg38.fa reference genome (for ref/alt determination)",
    )
    p.add_argument(
        "--out-dir", default="./results/alphag_scores",
        help="Output directory (default: ./results/alphag_scores)",
    )
    p.add_argument(
        "--output-types", nargs="+", default=None,
        help="AlphaGenome output types to score "
             "(default: all recommended scorers). "
             f"Available: {list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys())}",
    )
    p.add_argument(
        "--seq-length", default="1MB",
        choices=list(SEQ_LENGTH_MAP.keys()),
        help="Sequence length for predictions (default: 1MB)",
    )
    p.add_argument(
        "--batch-size", type=int, default=50,
        help="Number of variants per checkpoint (default: 50)",
    )
    p.add_argument(
        "--prefix", default="comvar",
        help="Output filename prefix: comvar or rarevar (default: comvar)",
    )
    p.add_argument(
        "--tf-change-filter", nargs="*", default=None,
        help="Only score variants with these TF_change categories "
             "(e.g., Create Disrupt). Default: all variants.",
    )
    p.add_argument(
        "--pip-threshold", type=float, default=0.1,
        help="Minimum max-PIP across traits to include a variant "
             "(requires --allele-src with PIP column). Default: 0.1",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Resolve API key: --api-key > env var > interactive prompt
    if args.api_key is None:
        args.api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if args.api_key is None:
        args.api_key = input("Enter AlphaGenome API key: ").strip()
    if not args.api_key:
        print("ERROR: No API key provided. Use --api-key, set ALPHAGENOME_API_KEY, or enter interactively.")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    total_steps = 4
    seq_length = SEQ_LENGTH_MAP[args.seq_length]

    print_header("AlphaGenome Variant Scoring")
    print_info(f"VAR2TFBS input:  {args.var2tfbs}")
    print_info(f"Cell type:       {args.cell}")
    print_info(f"Prefix:          {args.prefix}")
    print_info(f"Allele source:   {args.allele_src or 'not provided (rare variant mode)'}")
    print_info(f"Reference genome:{args.ref_genome or 'not provided'}")
    print_info(f"Output types:    {args.output_types or 'all recommended'}")
    print_info(f"Sequence length: {args.seq_length} ({seq_length:,} bp)")
    print_info(f"Batch size:      {args.batch_size}")
    print_info(f"PIP threshold:   {args.pip_threshold}")

    # --- Step 1: Load VAR2TFBS results ---
    print_step(1, total_steps, "Loading VAR2TFBS results")
    var2tfbs_df = pd.read_csv(args.var2tfbs)
    print_info(
        f"{len(var2tfbs_df):,} rows, "
        f"{var2tfbs_df['rsID'].nunique():,} unique variants"
    )

    if args.tf_change_filter:
        var2tfbs_df = var2tfbs_df[
            var2tfbs_df["TF_change"].isin(args.tf_change_filter)
        ].reset_index(drop=True)
        print_info(
            f"After TF_change filter ({args.tf_change_filter}): "
            f"{var2tfbs_df['rsID'].nunique():,} unique variants"
        )

    # --- Step 2: Resolve variant coordinates ---
    print_step(2, total_steps, "Resolving variant coordinates")
    allele_info = None
    max_pip = None
    seq_extractor = None

    if args.allele_src:
        print_info(f"Loading allele info from {args.allele_src}")
        allele_info, max_pip = load_allele_info(args.allele_src)
        print_info(f"  {len(allele_info):,} variants in allele source")

    # Filter by max PIP across traits
    if max_pip is not None and args.pip_threshold > 0:
        passing_rsids = {snp for snp, pip in max_pip.items() if pip >= args.pip_threshold}
        n_before = var2tfbs_df["rsID"].nunique()
        var2tfbs_df = var2tfbs_df[var2tfbs_df["rsID"].isin(passing_rsids)].reset_index(drop=True)
        n_after = var2tfbs_df["rsID"].nunique()
        print_info(
            f"PIP >= {args.pip_threshold} filter: "
            f"{n_after:,} / {n_before:,} variants pass"
        )

    if args.ref_genome:
        print_info(f"Loading reference genome: {args.ref_genome}")
        seq_extractor = FastaStringExtractor(args.ref_genome)

    resolved, skipped = extract_unique_variants(
        var2tfbs_df, allele_info, seq_extractor
    )

    if seq_extractor is not None:
        seq_extractor.close()

    print_success(f"{len(resolved):,} variants resolved")
    if skipped:
        print_info(f"WARNING: {len(skipped):,} variants could not be resolved")
        for s in skipped[:10]:
            print_info(f"  skipped: {s}")
        if len(skipped) > 10:
            print_info(f"  ... and {len(skipped) - 10} more")

    if not resolved:
        print_info("No variants to score. Exiting.")
        return

    # --- Step 3: Score variants with AlphaGenome ---
    print_step(3, total_steps, "Scoring variants with AlphaGenome")
    print_info("Connecting to AlphaGenome API ...")
    model = dna_client.create(args.api_key)

    if args.output_types:
        scorers = build_scorers(args.output_types)
    else:
        scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    print_info(f"Using {len(scorers)} scorer(s)")

    tidy_df = score_variants_with_checkpoints(
        model=model,
        rsid_variant_pairs=resolved,
        scorers=scorers,
        seq_length=seq_length,
        batch_size=args.batch_size,
        out_dir=args.out_dir,
        cell=args.cell,
        prefix=args.prefix,
    )

    if tidy_df.empty:
        print_info("No scores returned. Exiting.")
        return

    # --- Step 4: Save results ---
    print_step(4, total_steps, "Saving results")

    # Full scores (all biosamples)
    out_full = os.path.join(
        args.out_dir, f"{args.cell}_{args.prefix}_alphag_scores.tsv"
    )
    tidy_df.to_csv(out_full, sep="\t", index=False)
    print_info(f"Full scores: {out_full}")
    print_info(f"  {len(tidy_df):,} rows, {tidy_df['rsID'].nunique():,} variants")

    # Cell-type-specific scores
    ontology = CELL_ONTOLOGY[args.cell]
    cell_mask = (
        tidy_df["ontology_curie"].str.contains(ontology, na=False)
        | tidy_df["biosample_name"].str.contains(args.cell, na=False)
    )
    cell_df = tidy_df[cell_mask].reset_index(drop=True)

    out_cell = os.path.join(
        args.out_dir, f"{args.cell}_{args.prefix}_alphag_scores_{args.cell}.tsv"
    )
    cell_df.to_csv(out_cell, sep="\t", index=False)
    print_info(f"{args.cell}-only scores: {out_cell}")
    print_info(f"  {len(cell_df):,} rows, {cell_df['rsID'].nunique():,} variants")

    print()
    print("=" * 62)
    print("  AlphaGenome Scoring Complete")
    print(f"    Variants scored:   {tidy_df['rsID'].nunique():,}")
    n_types = tidy_df["output_type"].nunique() if "output_type" in tidy_df.columns else "N/A"
    print(f"    Output types:      {n_types}")
    print(f"    Total score rows:  {len(tidy_df):,} (full), {len(cell_df):,} ({args.cell})")
    n_tracks = tidy_df["track_name"].nunique() if "track_name" in tidy_df.columns else "N/A"
    print(f"    Unique tracks:     {n_tracks}")
    print(f"    Full output:       {out_full}")
    print(f"    {args.cell} output:      {out_cell}")
    print("=" * 62)
    print()


if __name__ == "__main__":
    main()
