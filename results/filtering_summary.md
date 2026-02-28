# varTFBridge Filtering Summary (K562)

## Common Variant Pipeline

| Step | Variants | Genes | Filter |
|------|----------|-------|--------|
| GWFM SNPs (hg38) | 13,065,104 | — | All SNPs across 15 traits |
| Credible set | 60,450 | — | PIP > 0.1 (13 erythroid traits) |
| Footprint overlap | 961 | — | Inside K562 FOODIE footprint |
| VAR2TFBS (FIMO hit) | 847 | — | Variant falls within a TF motif |
| VAR2TFBS (TF changed) | 839 | — | Create / Disrupt / Increase / Decrease |
| Var-to-Gene | 839 | 232 | ABC-FP-Max gene assignment |
| comvar2grn | 224 | 208 | Credible set filter (PEP_cs not NaN) |
| **comvar2grn (final)** | **209** | **195** | **+ ABC.Score > 0.015** |

### Common variant TFBS altered

- **1,154** unique variant–TF pairs across **330** TFs
- TF binding change: Create (457), Disrupt (421), Decrease (141), Increase (135)

### Common variant filtering strategy

1. **PIP > 0.1**: Retain variants in SBayesRC credible sets across 13 erythroid blood traits.
2. **FOODIE footprint overlap**: Variant must physically overlap a K562 TF footprint region (bedtools intersect).
3. **FIMO motif scanning**: Variant position must fall within a JASPAR TF binding motif hit (p < 1e-4).
4. **TF binding change**: At least one TF motif must show a change between ref and alt allele (Create, Disrupt, Increase, or Decrease). Variants with only Unchange hits are excluded.
5. **ABC-FP-Max gene assignment**: Each variant is linked to its top target gene via the Activity-by-Contact model with FOODIE footprint weighting.
6. **Credible set annotation**: Only variant-trait pairs where PEP_cs is not NaN (variant belongs to a local credible set for that trait) are retained.
7. **ABC.Score > 0.015**: Enhancer-gene links with ABC.Score <= 0.015 are filtered out as weak regulatory connections.

## Rare Variant Pipeline

| Step | Count | Genes | Filter |
|------|-------|-------|--------|
| K562 FOODIE footprints | 188,484 | — | All footprints tested in burden analysis |
| Significant footprints | 24 | — | Burden test Bonferroni (p < 2.65e-07) |
| MAC > 30 footprints | 21 | — | At least one LOO variant with MAC > 30 |
| Driver variants (LOO) | 21 | — | Best leave-one-out variant per footprint |
| VAR2TFBS (FIMO hit) | 19 | — | Variant falls within a TF motif |
| Var-to-Gene | 19 | 15 | ABC-FP-Max gene assignment (ABC.Score > 0.015) |
| **rarevar2grn (final)** | **19** | **15** | **Full linkage table** |

### Rare variant TFBS altered

- **143** unique variant–TF pairs across **99** TFs
- TF binding change: Create (60), Disrupt (55), Decrease (20), Increase (8)

### Rare variant filtering strategy

1. **Burden test**: Rare variant burden test across 188,484 K562 FOODIE footprints per trait. Significance determined by Bonferroni correction: -log10(0.05 / 188,484) = 6.58.
2. **MAC > 30**: Footprints must contain at least one rare variant with minor allele count > 30 in the leave-one-out (LOO) analysis, ensuring sufficient statistical power.
3. **Driver variant identification**: For each significant footprint, the LOO variant with the strongest effect (largest drop in burden test significance upon removal) is selected as the driver.
4. **FIMO motif scanning**: Driver variant position must fall within a JASPAR TF binding motif hit. 2 driver variants are lost at this step (no motif at variant position).
5. **ABC-FP-Max gene assignment**: Each driver variant is linked to its top target gene via ABC model (ABC.Score > 0.015). 1 variant loses its gene link at this step (best ABC = 0.012, below threshold), but all 19 variants with FIMO hits are retained because they all have ABC.Score > 0.015.
