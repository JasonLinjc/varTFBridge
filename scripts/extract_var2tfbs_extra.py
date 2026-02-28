#!/usr/bin/env python
"""Extract TFBS affected by variants that are NOT scored by AlphaGenome TF ChIP.

Reads comvar2grn and AlphaGenome scores, finds TF binding changes where the
TF has no AlphaGenome K562 ChIP-seq track. Filters to max PIP above threshold.

Usage:
  python scripts/extract_var2tfbs_extra.py                  # default PIP > 0.7
  python scripts/extract_var2tfbs_extra.py --pip-threshold 0.5

Output: results/K562_comvar_var2tfbsExtra.csv
"""

import os
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--project-root', default='.', help='Project root directory')
    parser.add_argument('--pip-threshold', type=float, default=0.7,
                        help='Max PIP threshold (default: 0.7)')
    args = parser.parse_args()

    root = args.project_root
    results_dir = os.path.join(root, 'results')

    # Load comvar2grn
    grn = pd.read_csv(os.path.join(results_dir, 'K562_comvar2grn.csv'))
    print(f'comvar2grn: {len(grn):,} rows, {grn["rsID"].nunique()} variants')

    # Load AlphaGenome TF names (K562 ChIP-seq)
    alphag = pd.read_csv(
        os.path.join(results_dir, 'alphag_scores', 'K562_comvar_alphag_scores_K562.tsv'),
        sep='\t'
    )
    alphag_tfs = set(
        alphag[
            (alphag['output_type'] == 'CHIP_TF') &
            (alphag['variant_scorer'].str.contains('DIFF_LOG2_SUM', na=False))
        ]['transcription_factor'].str.upper().unique()
    )
    print(f'AlphaGenome K562 TFs: {len(alphag_tfs)}')

    # Filter to affected TF changes only
    affected = grn[grn['TF_change'] != 'Unchange'].copy()

    # Check if ALL component TFs are missing from AlphaGenome
    # (for co-binding like GATA1::TAL1, only flag if NONE of the components are scored)
    def all_components_missing(tf_name):
        components = tf_name.split('::')
        return all(comp.upper() not in alphag_tfs for comp in components)

    affected = affected[affected['TF'].apply(all_components_missing)].copy()
    print(f'Affected rows with TF not in AlphaGenome: {len(affected):,}')

    # Compute max PIP per rsID and filter by threshold
    pip_thresh = args.pip_threshold
    max_pip = grn.groupby('rsID')['PIP'].max().rename('max_PIP')
    affected = affected.merge(max_pip, on='rsID', how='left')
    affected = affected[affected['max_PIP'] > pip_thresh].reset_index(drop=True)
    print(f'After max_PIP > {pip_thresh} filter: {len(affected):,} rows, '
          f'{affected["rsID"].nunique()} variants')

    # Select and order columns
    col_order = [
        'rsID', 'max_PIP', 'Chromosome', 'Start', 'End', 'A1', 'A2',
        'TF', 'TF_change', 'p-value_ref', 'p-value_alt', 'TF_K562_rna_tpm',
        'TargetGene', 'ABC.Score.FP', 'trait', 'PIP',
    ]
    out = affected[col_order].sort_values(
        by=['max_PIP', 'rsID', 'TF'], ascending=[False, True, True]
    ).reset_index(drop=True)

    out_path = os.path.join(results_dir, 'K562_comvar_var2tfbsExtra.csv')
    out.to_csv(out_path, index=False)

    print(f'\n=== Summary ===')
    print(f'Output: {out_path}')
    print(f'Rows: {len(out):,}')
    print(f'Unique variants: {out["rsID"].nunique():,}')
    print(f'Unique TFs: {out["TF"].nunique():,}')
    print(f'\nTF_change distribution:')
    print(out['TF_change'].value_counts().to_string())
    print(f'\nTop 10 TFs by frequency:')
    print(out['TF'].value_counts().head(10).to_string())


if __name__ == '__main__':
    main()
