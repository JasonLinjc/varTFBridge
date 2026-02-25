#!/usr/bin/env python
"""Annotate PIP>0.7 common variants with K562 AlphaGenome TF ChIP scores.

For each variant with an affected TF binding change (FIMO), looks up the
AlphaGenome TF ChIP DIFF_LOG2_SUM score for that TF.

For co-binding TFs (e.g., GATA1::TAL1), splits into two rows — one per
component TF — each with its own AlphaGenome score and expression.

Output: results/K562_comvar_pip70_alphag_tf_chip.csv
"""

import os
import glob
import argparse
import pandas as pd


def load_pip70_variants(credset_dir):
    """Load variants with max PIP > 0.7 across all traits."""
    frames = []
    for fn in sorted(glob.glob(os.path.join(credset_dir, '*_credible_set_hg38.csv'))):
        df = pd.read_csv(fn, usecols=['Name', 'PIP'])
        frames.append(df)
    pip_all = pd.concat(frames).groupby('Name')['PIP'].max().reset_index()
    pip_all.columns = ['rsID', 'max_PIP']
    pip70 = pip_all[pip_all['max_PIP'] > 0.7].reset_index(drop=True)
    print(f'  PIP > 0.7 variants: {len(pip70):,}')
    return pip70


def load_alphag_tf_chip(alphag_path):
    """Load AlphaGenome K562 TF ChIP scores (DIFF_LOG2_SUM) as a lookup dict."""
    alphag = pd.read_csv(alphag_path, sep='\t')
    chip_tf = alphag[
        (alphag['output_type'] == 'CHIP_TF') &
        (alphag['variant_scorer'].str.contains('DIFF_LOG2_SUM', na=False))
    ][['rsID', 'transcription_factor', 'raw_score', 'quantile_score']].copy()
    print(f'  AlphaGenome TF ChIP: {chip_tf["rsID"].nunique():,} variants, '
          f'{chip_tf["transcription_factor"].nunique()} TFs')
    return chip_tf


def load_var2tfbs(var2tfbs_path):
    """Load FIMO-predicted TF binding changes (affected only)."""
    var2tfbs = pd.read_csv(var2tfbs_path)
    var2tfbs = var2tfbs[var2tfbs['TF_change'] != 'Unchange'].reset_index(drop=True)
    keep = var2tfbs[['rsID', 'TF', 'TF_change', 'p-value_ref', 'p-value_alt']].drop_duplicates()
    print(f'  var2tfbs (affected): {keep["rsID"].nunique():,} variants, '
          f'{len(keep):,} rsID-TF pairs')
    return keep


def load_tf_expression(tf_expr_path):
    """Load TF expression as a lookup dict: gene_name_upper -> max TPM across isoforms."""
    tf_expr = pd.read_csv(tf_expr_path, sep='\t', usecols=['gene_name', 'tpm'])
    # Keep highest TPM isoform per gene
    tf_expr = tf_expr.groupby('gene_name')['tpm'].max().reset_index()
    tf_expr['gene_upper'] = tf_expr['gene_name'].str.upper()
    expr_lookup = dict(zip(tf_expr['gene_upper'], tf_expr['tpm']))
    print(f'  TF expression (max isoform TPM): {len(expr_lookup):,} genes')
    return expr_lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--project-root', default='.', help='Project root directory')
    args = parser.parse_args()

    root = args.project_root
    results_dir = os.path.join(root, 'results')
    data_dir = os.path.join(root, 'data')

    print('Loading data...')

    # Step 1: PIP > 0.7 variants
    pip70 = load_pip70_variants(
        os.path.join(data_dir, 'GWFM_erythroids', 'credible_set_snpRes')
    )

    # Step 2: AlphaGenome TF ChIP scores
    chip_tf = load_alphag_tf_chip(
        os.path.join(results_dir, 'alphag_scores', 'K562_comvar_alphag_scores_K562.tsv')
    )

    # Step 3: FIMO TF binding changes
    var2tfbs = load_var2tfbs(
        os.path.join(results_dir, 'comvar_var2tfbs_results', 'K562_var2tfbs.csv')
    )

    # Step 3b: TF expression (ENCODE K562 RNA-seq, max isoform TPM per gene)
    expr_lookup = load_tf_expression(
        os.path.join(data_dir, 'gene_expr', 'K562_ENCFF485RIA_gene.tsv')
    )

    # Step 4: Filter var2tfbs to PIP > 0.7 variants
    print('\nMerging...')
    var2tfbs_pip70 = var2tfbs.merge(pip70, on='rsID', how='inner')
    print(f'  var2tfbs x PIP>0.7: {var2tfbs_pip70["rsID"].nunique():,} variants, '
          f'{len(var2tfbs_pip70):,} rsID-TF pairs')

    # Step 5: Split co-binding TFs and look up AlphaGenome scores
    # Build lookup: (rsID, TF_upper) -> (raw_score, quantile_score)
    chip_tf['TF_upper'] = chip_tf['transcription_factor'].str.upper()
    chip_tf_dedup = chip_tf.drop_duplicates(subset=['rsID', 'TF_upper'], keep='first')
    alphag_lookup = chip_tf_dedup.set_index(['rsID', 'TF_upper'])[
        ['raw_score', 'quantile_score']
    ].to_dict('index')

    rows = []
    for _, row in var2tfbs_pip70.iterrows():
        rsid = row['rsID']
        tf_motif = row['TF']
        base = {
            'rsID': rsid,
            'max_PIP': row['max_PIP'],
            'TF_motif': tf_motif,
            'TF_change': row['TF_change'],
            'p-value_ref': row['p-value_ref'],
            'p-value_alt': row['p-value_alt'],
        }

        # Split co-binding TFs into individual rows
        components = tf_motif.split('::')
        for tf_name in components:
            tf_upper = tf_name.upper()
            s = alphag_lookup.get((rsid, tf_upper), {})
            rows.append({
                **base,
                'TF_alphag': tf_name,
                'TF_alphag_K562_rna_tpm': expr_lookup.get(tf_upper),
                'alphag_TF_chip_score': s.get('raw_score'),
                'alphag_TF_chip_quantile': s.get('quantile_score'),
            })

    merged = pd.DataFrame(rows)

    # Sort by rsID then absolute quantile score descending
    merged['_sort_q'] = merged['alphag_TF_chip_quantile'].abs()
    merged = merged.sort_values(
        by=['rsID', '_sort_q'], ascending=[True, False]
    ).drop(columns=['_sort_q']).reset_index(drop=True)

    col_order = [
        'rsID', 'max_PIP', 'TF_motif', 'TF_change', 'p-value_ref', 'p-value_alt',
        'TF_alphag', 'TF_alphag_K562_rna_tpm',
        'alphag_TF_chip_score', 'alphag_TF_chip_quantile',
    ]
    merged = merged[col_order]

    out_path = os.path.join(results_dir, 'K562_comvar_pip70_alphag_tf_chip.csv')
    merged.to_csv(out_path, index=False)

    # Summary
    n_cobind = merged[merged['TF_motif'].str.contains('::', na=False)].shape[0]
    n_with_alphag = merged['alphag_TF_chip_score'].notna().sum()
    print(f'\n=== Summary ===')
    print(f'Output: {out_path}')
    print(f'Rows: {len(merged):,}')
    print(f'Unique variants: {merged["rsID"].nunique():,}')
    print(f'Unique TF motifs: {merged["TF_motif"].nunique():,}')
    print(f'Co-binding TF rows (:: expanded): {n_cobind:,}')
    print(f'With AlphaGenome score: {n_with_alphag:,}/{len(merged):,}')
    print(f'\nTF_change distribution:')
    print(merged['TF_change'].value_counts().to_string())

    # Show co-binding examples
    cobind_rows = merged[merged['TF_motif'].str.contains('::', na=False)].head(10)
    if len(cobind_rows) > 0:
        print(f'\nCo-binding TF examples (two rows per motif):')
        print(cobind_rows[['rsID', 'TF_motif', 'TF_change',
                           'TF_alphag', 'TF_alphag_K562_rna_tpm',
                           'alphag_TF_chip_score', 'alphag_TF_chip_quantile'
                           ]].to_string(index=False))


if __name__ == '__main__':
    main()
