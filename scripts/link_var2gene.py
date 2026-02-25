#!/usr/bin/env python
"""
Link variants to target genes via ABC-FP-Max scores.

For each variant overlapping a FOODIE footprint, finds the footprint's
overlapping ATAC enhancer region in the ABC-FP predictions and assigns
the target gene with the highest ABC-FP score. Annotates with TF binding
changes and cell-type-specific TF RNA expression.

Workflow:
    1. Load VAR2TFBS results (common or rare) — variant-to-footprint + TF info
    2. Intersect FOODIE footprints with ATAC enhancer regions (bedtools)
    3. Load ABC-FP predictions, compute ABC-FP-Max score, link to genes
    4. Annotate with TF RNA expression and output results

Usage:
    python scripts/link_var2gene.py \
        --var2tfbs results/comvar_var2tfbs_results/K562_var2tfbs.csv \
        --footprint-bed data/FOODIE_footprints/K562.merged.hg38.bed \
        --enhancer-bed data/ABC_FP_results/K562_FOODIE_ATAC/Neighborhoods/EnhancerList.bed \
        --abc-predictions data/ABC_FP_results/K562_FOODIE_ATAC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
        --tf-expr data/gene_expr/K562_ENCFF485RIA_gene.tsv \
        --cell K562 \
        --out-dir results/var2gene_results
"""

import argparse
import os
import subprocess

import pandas as pd


def print_header(text):
    print()
    print("=" * 62)
    print(f"  {text}")
    print("=" * 62)


def print_step(step, total, text):
    print(f"\n  [{step}/{total}] {text}")
    print("  " + "-" * 58)


def print_info(text):
    print(f"    {text}")


def load_var2tfbs(path):
    """Load VAR2TFBS results CSV."""
    df = pd.read_csv(path)
    required = ["rsID", "foodie_id"]
    if not all(c in df.columns for c in required):
        raise ValueError(f"VAR2TFBS file must contain columns: {required}")
    return df


def intersect_footprints_enhancers(footprint_bed, enhancer_bed):
    """Use bedtools intersect to map FOODIE footprints to ATAC enhancer regions."""
    cmd = [
        "bedtools", "intersect",
        "-a", footprint_bed,
        "-b", enhancer_bed,
        "-wa", "-wb",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        fp_name = fields[3]
        # Footprint BED has 6 cols, enhancer BED has 4 cols → enhancer name at index 9
        enh_name = fields[-1]
        rows.append({"foodie_id": fp_name, "enhancer_name": enh_name})

    return pd.DataFrame(rows).drop_duplicates()


def load_abc_predictions(path):
    """Load ABC-FP predictions and compute ABC-FP-Max score."""
    df = pd.read_csv(path, sep="\t")
    required = ["name", "TargetGene", "ABC.Score"]
    if not all(c in df.columns for c in required):
        raise ValueError(f"ABC predictions must contain columns: {required}")

    # ABC-FP-Max score: activity_base * hic_contact, normalized per gene
    if "activity_base" in df.columns and "hic_contact" in df.columns:
        df["ABC.Score.FP"] = df["activity_base"] * df["hic_contact"]
        gene_sum = df.groupby("TargetGene")["ABC.Score.FP"].transform("sum")
        df["ABC.Score.FP"] = df["ABC.Score.FP"].where(gene_sum == 0, df["ABC.Score.FP"] / gene_sum)

    return df


def load_tf_expression(path, cell):
    """Load TF RNA expression (ENCODE RNA-seq, max isoform TPM per gene).

    For co-binding TFs (e.g. GATA1::TAL1), returns max TPM of components.
    """
    df = pd.read_csv(path, sep='\t', usecols=['gene_name', 'tpm'])
    df = df.groupby('gene_name')['tpm'].max().reset_index()
    expr_lookup = dict(zip(df['gene_name'].str.upper(), df['tpm']))

    def lookup_tf_expr(tf_name):
        components = tf_name.split('::')
        tpms = [expr_lookup.get(c.upper()) for c in components]
        tpms = [t for t in tpms if t is not None]
        return max(tpms) if tpms else None

    return expr_lookup, lookup_tf_expr


def parse_args():
    p = argparse.ArgumentParser(
        description="Link variants to target genes via ABC-FP-Max scores"
    )
    p.add_argument(
        "--var2tfbs", required=True,
        help="VAR2TFBS output CSV (e.g. results/comvar_var2tfbs_results/K562_var2tfbs.csv)",
    )
    p.add_argument(
        "--footprint-bed", required=True,
        help="FOODIE footprint BED file (e.g. data/FOODIE_footprints/K562.merged.hg38.bed)",
    )
    p.add_argument(
        "--enhancer-bed", required=True,
        help="ABC-FP EnhancerList BED (e.g. data/ABC_FP_results/K562_FOODIE_ATAC/Neighborhoods/EnhancerList.bed)",
    )
    p.add_argument(
        "--abc-predictions", required=True,
        help="ABC-FP predictions TSV (e.g. ABC_FP_results/.../EnhancerPredictionsAllPutative.tsv.gz)",
    )
    p.add_argument(
        "--tf-expr", default=None,
        help="TF RNA expression TSV (e.g. data/gene_expr/K562_ENCFF485RIA_gene.tsv)",
    )
    p.add_argument(
        "--cell", default="K562",
        help="Cell type name for TF expression column lookup (default: K562)",
    )
    p.add_argument(
        "--abc-threshold", type=float, default=0.0,
        help="Minimum ABC score threshold for gene assignment (default: 0, keep all)",
    )
    p.add_argument(
        "--prefix", default="comvar",
        help="Output filename prefix, e.g. comvar or rarevar (default: comvar)",
    )
    p.add_argument(
        "--out-dir", default="./results/var2gene_results",
        help="Output directory (default: ./results/var2gene_results)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    total_steps = 4
    tpm_col = f"TF_{args.cell}_rna_tpm"

    print_header("Variant-to-Gene Linking via ABC-FP-Max")
    print_info(f"VAR2TFBS input:    {args.var2tfbs}")
    print_info(f"Footprint BED:     {args.footprint_bed}")
    print_info(f"Enhancer BED:      {args.enhancer_bed}")
    print_info(f"ABC predictions:   {args.abc_predictions}")
    print_info(f"TF expression:     {args.tf_expr or 'not provided'}")
    print_info(f"Cell type:         {args.cell}")
    print_info(f"ABC threshold:     {args.abc_threshold}")

    # Step 1: Load VAR2TFBS results
    print_step(1, total_steps, "Loading VAR2TFBS results")
    var2tfbs = load_var2tfbs(args.var2tfbs)

    # Extract TF info: rsID, TF, TF_change, foodie_id (excluding Unchange)
    tf_cols = ["rsID", "foodie_id", "TF", "TF_change"]
    tf_cols = [c for c in tf_cols if c in var2tfbs.columns]
    var_tf = var2tfbs[var2tfbs["TF_change"] != "Unchange"][tf_cols].drop_duplicates()
    print_info(f"{len(var_tf):,} variant-TF-footprint entries (TF_change != Unchange)")
    print_info(f"{var_tf['rsID'].nunique():,} unique variants")

    # Unique variant-footprint pairs for gene linking
    var_fp = var_tf[["rsID", "foodie_id"]].drop_duplicates()

    # Annotate with TF expression
    if args.tf_expr:
        print_info(f"Loading TF expression from {args.tf_expr}")
        _, lookup_tf_expr = load_tf_expression(args.tf_expr, args.cell)
        var_tf[tpm_col] = var_tf["TF"].apply(lookup_tf_expr)
        print_info(f"{(var_tf[tpm_col] > 0).sum():,} entries with TF expression > 0")

    # Step 2: Map footprints to ATAC enhancer regions
    print_step(2, total_steps, "Intersecting footprints with ATAC enhancer regions")
    fp_enh = intersect_footprints_enhancers(args.footprint_bed, args.enhancer_bed)
    print_info(f"{len(fp_enh):,} footprint-enhancer overlaps")
    print_info(f"{fp_enh['foodie_id'].nunique():,} footprints mapped to enhancers")

    # Step 3: Load ABC-FP predictions and merge
    print_step(3, total_steps, "Loading ABC-FP predictions and linking to genes")
    abc = load_abc_predictions(args.abc_predictions)
    print_info(f"{len(abc):,} enhancer-gene predictions loaded")

    # Bridge: variant -> footprint -> enhancer -> gene
    var_enh = var_fp.merge(fp_enh, on="foodie_id", how="inner")
    print_info(f"{len(var_enh):,} variant-enhancer pairs after footprint-enhancer join")

    abc_cols = ["name", "class", "TargetGene", "TargetGeneExpression",
                "TargetGeneIsExpressed", "distance", "isSelfPromoter",
                "hic_contact", "activity_base", "ABC.Score", "ABC.Score.FP"]
    abc_cols = [c for c in abc_cols if c in abc.columns]

    var_gene = var_enh.merge(abc[abc_cols], left_on="enhancer_name", right_on="name", how="inner")
    print_info(f"{len(var_gene):,} variant-gene links before filtering")

    if args.abc_threshold > 0:
        var_gene = var_gene[var_gene["ABC.Score"] > args.abc_threshold].reset_index(drop=True)
        print_info(f"{len(var_gene):,} variant-gene links after ABC threshold > {args.abc_threshold}")

    # Step 4: Assign top gene per variant and merge TF info
    print_step(4, total_steps, "Assigning top target gene per variant")

    # Top gene per variant: highest ABC.Score (the ABC model's built-in score)
    var_gene_top = (
        var_gene
        .sort_values(["rsID", "ABC.Score"], ascending=[True, False])
        .drop_duplicates(subset=["rsID"], keep="first")
        .reset_index(drop=True)
    )
    print_info(f"{len(var_gene_top):,} variants with assigned target genes")
    print_info(f"{var_gene_top['TargetGene'].nunique():,} unique target genes")

    # Merge TF info back: variant-TF-gene table
    merge_cols = ["rsID", "enhancer_name", "class", "TargetGene",
                  "distance", "isSelfPromoter", "hic_contact",
                  "activity_base", "ABC.Score"]
    if "ABC.Score.FP" in var_gene_top.columns:
        merge_cols.append("ABC.Score.FP")
    var_tf_gene = var_tf.merge(var_gene_top[merge_cols], on="rsID", how="inner")
    print_info(f"{len(var_tf_gene):,} variant-TF-gene entries in final table")

    # Count variants without gene assignment
    no_gene = var_fp["rsID"].nunique() - var_gene_top["rsID"].nunique()
    if no_gene > 0:
        print_info(f"{no_gene:,} variants could not be linked to any gene")

    # Save outputs — derive cell name from input filename
    basename = os.path.splitext(os.path.basename(args.var2tfbs))[0]
    basename = basename.replace("_rarevar_var2tfbs", "").replace("_var2tfbs", "")

    # ABC-FP-Full: all variant-gene links with TF info
    full_path = os.path.join(args.out_dir, f"{basename}_{args.prefix}_ABC-FP-Full.csv")
    out_cols = ["rsID", "TF", "TF_change"]
    if tpm_col in var_tf_gene.columns:
        out_cols.append(tpm_col)
    out_cols += ["TargetGene", "ABC.Score"]
    if "ABC.Score.FP" in var_tf_gene.columns:
        out_cols.append("ABC.Score.FP")
    out_cols += ["distance", "foodie_id"]
    out_cols = [c for c in out_cols if c in var_tf_gene.columns]
    var_tf_gene[out_cols].to_csv(full_path, index=False)

    # ABC-FP-Max: top gene per variant (highest ABC-FP score)
    max_path = os.path.join(args.out_dir, f"{basename}_{args.prefix}_ABC-FP-Max.csv")
    max_cols = ["rsID", "foodie_id", "enhancer_name", "class", "TargetGene",
                "distance", "isSelfPromoter", "hic_contact", "activity_base",
                "ABC.Score"]
    if "ABC.Score.FP" in var_gene_top.columns:
        max_cols.append("ABC.Score.FP")
    max_cols = [c for c in max_cols if c in var_gene_top.columns]
    var_gene_top[max_cols].to_csv(max_path, index=False)

    print()
    print("=" * 62)
    print(f"  Results")
    print(f"    Variants with target genes: {len(var_gene_top):,}")
    print(f"    Unique target genes:        {var_gene_top['TargetGene'].nunique():,}")
    print(f"    ABC-FP-Full (all links):    {full_path}")
    print(f"    ABC-FP-Max (top gene):      {max_path}")
    print("=" * 62)
    print()


if __name__ == "__main__":
    main()
