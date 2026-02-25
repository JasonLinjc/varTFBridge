#!/usr/bin/env python
"""Merge var2tfbs + snpRes overlap + var2gene + AlphaGenome into a final linkage table.

Each row represents one rsID x TF x trait combination, showing the full chain:
  variant -> TF binding change -> target gene -> trait association

Filters:
  - Affected TF_change only (Create/Disrupt/Increase/Decrease)
  - Credible set variants only (PEP_cs not NaN)

Output: results/K562_comvar2grn.csv
"""

import os
import glob
import argparse
import pandas as pd


def load_snpres_overlap(overlap_dir):
    """Load per-trait snpRes overlap files and filter to credible set variants."""
    frames = []
    for fn in sorted(glob.glob(os.path.join(overlap_dir, '*_K562.merged.hg38.csv'))):
        df = pd.read_csv(fn)
        trait = os.path.basename(fn).split('_')[0]
        df['trait'] = trait
        frames.append(df)
    snpres = pd.concat(frames, ignore_index=True)
    snpres = snpres.rename(columns={'SNP': 'rsID'})
    # Filter to credible set variants
    snpres = snpres[snpres['PEP_cs'].notna()].reset_index(drop=True)
    print(f'  snpRes overlap: {snpres["rsID"].nunique():,} unique variants, '
          f'{len(snpres):,} variant-trait rows (credible set)')
    return snpres


def load_a1effect(credset_dir):
    """Load A1Effect from credible set source files."""
    frames = []
    for fn in sorted(glob.glob(os.path.join(credset_dir, '*_credible_set_hg38.csv'))):
        df = pd.read_csv(fn, usecols=['Name', 'A1Effect'])
        trait = os.path.basename(fn).split('_')[0]
        df['trait'] = trait
        df = df.rename(columns={'Name': 'rsID'})
        frames.append(df)
    a1eff = pd.concat(frames, ignore_index=True).drop_duplicates()
    print(f'  A1Effect: {a1eff["rsID"].nunique():,} unique variants across {a1eff["trait"].nunique()} traits')
    return a1eff


def load_var2tfbs(var2tfbs_path):
    """Load var2tfbs results, filter to affected TF changes."""
    var2tfbs = pd.read_csv(var2tfbs_path)
    var2tfbs = var2tfbs[var2tfbs['TF_change'] != 'Unchange'].reset_index(drop=True)
    keep_cols = ['TF', 'p-value_ref', 'p-value_alt', 'rsID', 'foodie_id', 'TF_change']
    var2tfbs = var2tfbs[keep_cols]
    print(f'  var2tfbs: {var2tfbs["rsID"].nunique():,} variants, '
          f'{len(var2tfbs):,} rsID-TF pairs (affected only)')
    return var2tfbs


def load_tf_expression(tf_expr_path):
    """Load TF expression (ENCODE RNA-seq, max isoform TPM per gene).

    For co-binding TFs (e.g. GATA1::TAL1), returns max TPM of components.
    """
    tf_expr = pd.read_csv(tf_expr_path, sep='\t', usecols=['gene_name', 'tpm'])
    tf_expr = tf_expr.groupby('gene_name')['tpm'].max().reset_index()
    expr_lookup = dict(zip(tf_expr['gene_name'].str.upper(), tf_expr['tpm']))
    return expr_lookup


def load_var2gene(var2gene_path):
    """Load variant-to-gene linkage (ABC-FP-Max)."""
    var2gene = pd.read_csv(var2gene_path)
    var2gene = var2gene.rename(columns={'class': 'enhancer_class'})
    keep_cols = ['rsID', 'foodie_id', 'TargetGene', 'ABC.Score', 'ABC.Score.FP',
                 'distance', 'enhancer_class']
    var2gene = var2gene[keep_cols]
    print(f'  var2gene: {var2gene["rsID"].nunique():,} variants -> '
          f'{var2gene["TargetGene"].nunique():,} genes')
    return var2gene


def load_alphagenome(alphag_path):
    """Load AlphaGenome K562 scores, extract H3K27ac and ATAC (DIFF_LOG2_SUM)."""
    alphag = pd.read_csv(alphag_path, sep='\t')

    # H3K27ac
    h3k = alphag[
        (alphag['output_type'] == 'CHIP_HISTONE') &
        (alphag['histone_mark'] == 'H3K27ac') &
        (alphag['variant_scorer'].str.contains('DIFF_LOG2_SUM', na=False))
    ][['rsID', 'raw_score', 'quantile_score']].rename(columns={
        'raw_score': 'alphag_H3K27ac_score',
        'quantile_score': 'alphag_H3K27ac_quantile'
    })

    # ATAC
    atac = alphag[
        (alphag['output_type'] == 'ATAC') &
        (alphag['variant_scorer'].str.contains('DIFF_LOG2_SUM', na=False))
    ][['rsID', 'raw_score', 'quantile_score']].rename(columns={
        'raw_score': 'alphag_ATAC_score',
        'quantile_score': 'alphag_ATAC_quantile'
    })

    alphag_scores = h3k.merge(atac, on='rsID', how='outer')
    print(f'  AlphaGenome: {alphag_scores["rsID"].nunique():,} variants with K562 scores')
    return alphag_scores


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--project-root', default='.', help='Project root directory')
    args = parser.parse_args()

    root = args.project_root
    results_dir = os.path.join(root, 'results')
    data_dir = os.path.join(root, 'data')

    print('Loading data...')

    # Step 1: Per-trait variant data (credible set only)
    snpres = load_snpres_overlap(
        os.path.join(results_dir, 'comvar_footprint_overlap_snpRes', 'K562.merged.hg38')
    )

    # Step 2: A1Effect from credible set source files
    a1eff = load_a1effect(
        os.path.join(data_dir, 'GWFM_erythroids', 'credible_set_snpRes')
    )
    snpres = snpres.merge(a1eff, on=['rsID', 'trait'], how='left')

    # Step 3: var2tfbs (affected only)
    var2tfbs = load_var2tfbs(
        os.path.join(results_dir, 'comvar_var2tfbs_results', 'K562_var2tfbs.csv')
    )

    # Step 4: Merge variant-trait x TF binding
    print('\nMerging...')
    merged = snpres.merge(
        var2tfbs,
        left_on=['rsID', 'footprint_region'],
        right_on=['rsID', 'foodie_id'],
        how='inner'
    )
    merged = merged.drop(columns=['foodie_id'])
    print(f'  After var-trait x TF merge: {len(merged):,} rows, '
          f'{merged["rsID"].nunique():,} variants')

    # Step 5: TF expression (ENCODE K562 RNA-seq, max isoform TPM)
    expr_lookup = load_tf_expression(
        os.path.join(data_dir, 'gene_expr', 'K562_ENCFF485RIA_gene.tsv')
    )
    # For co-binding TFs (e.g. GATA1::TAL1), use max TPM of components
    def lookup_tf_expr(tf_name):
        components = tf_name.split('::')
        tpms = [expr_lookup.get(c.upper()) for c in components]
        tpms = [t for t in tpms if t is not None]
        return max(tpms) if tpms else None
    merged['TF_K562_rna_tpm'] = merged['TF'].apply(lookup_tf_expr)

    # Step 6: Gene linkage
    var2gene = load_var2gene(
        os.path.join(results_dir, 'var2gene_results', 'K562_comvar_ABC-FP-Max.csv')
    )
    merged = merged.merge(
        var2gene,
        left_on=['rsID', 'footprint_region'],
        right_on=['rsID', 'foodie_id'],
        how='left'
    )
    if 'foodie_id' in merged.columns:
        merged = merged.drop(columns=['foodie_id'])

    # Step 7: AlphaGenome scores
    alphag_path = os.path.join(results_dir, 'alphag_scores', 'K562_comvar_alphag_scores_K562.tsv')
    if os.path.exists(alphag_path):
        alphag_scores = load_alphagenome(alphag_path)
        merged = merged.merge(alphag_scores, on='rsID', how='left')
    else:
        print(f'  Warning: AlphaGenome scores not found at {alphag_path}, skipping')

    # Step 8: Sort and save
    # Natural chromosome sort
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24
    merged['_chrom_sort'] = merged['Chromosome'].map(chrom_order)
    merged = merged.sort_values(
        by=['_chrom_sort', 'Start', 'trait', 'TF']
    ).drop(columns=['_chrom_sort']).reset_index(drop=True)

    # Column order
    col_order = [
        'rsID', 'Chromosome', 'Start', 'End', 'A1', 'A2', 'freq', 'footprint_region',
        'trait', 'PIP', 'PEP', 'PEP_cs', 'CS_id', 'A1Effect',
        'TF', 'TF_change', 'p-value_ref', 'p-value_alt', 'TF_K562_rna_tpm',
        'TargetGene', 'ABC.Score', 'ABC.Score.FP', 'distance', 'enhancer_class',
        'alphag_H3K27ac_score', 'alphag_H3K27ac_quantile',
        'alphag_ATAC_score', 'alphag_ATAC_quantile',
    ]
    col_order = [c for c in col_order if c in merged.columns]
    merged = merged[col_order]

    out_path = os.path.join(results_dir, 'K562_comvar2grn.csv')
    merged.to_csv(out_path, index=False)

    # Summary
    print(f'\n=== Summary ===')
    print(f'Output: {out_path}')
    print(f'Rows: {len(merged):,}')
    print(f'Unique variants: {merged["rsID"].nunique():,}')
    print(f'Unique TFs: {merged["TF"].nunique():,}')
    print(f'Unique genes: {merged["TargetGene"].nunique():,}')
    print(f'Traits: {sorted(merged["trait"].unique())}')
    print(f'\nTF_change distribution:')
    print(merged['TF_change'].value_counts().to_string())
    print(f'\nWith gene assignment: '
          f'{merged["TargetGene"].notna().sum():,}/{len(merged):,} rows')
    print(f'With AlphaGenome H3K27ac: '
          f'{merged.get("alphag_H3K27ac_score", pd.Series(dtype=float)).notna().sum():,}/{len(merged):,} rows')

    # Spot check rs112233623
    rs_check = merged[merged['rsID'] == 'rs112233623']
    if len(rs_check) > 0:
        print(f'\nSpot check rs112233623:')
        print(f'  Traits: {sorted(rs_check["trait"].unique())}')
        print(f'  TFs: {rs_check[["TF", "TF_change"]].drop_duplicates().to_string(index=False)}')
        gene = rs_check['TargetGene'].dropna().unique()
        print(f'  Gene: {gene}')


if __name__ == '__main__':
    main()
