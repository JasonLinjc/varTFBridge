# varTFBridge

<p align="center">
  <img src="logo.png" alt="Logo">
</p>

## Genome-wide maps of transcription factor footprints identify noncoding variants rewiring gene regulatory networks

**varTFBridge** is an integrative framework that combines transcription factor (TF) footprinting data with genome-wide association analyses to identify causal noncoding variants and elucidate their regulatory mechanisms in TF-mediated gene regulation.

## Overview

Common-variant genome-wide association studies have identified thousands of noncoding loci associated with human diseases and complex traits; however, interpreting their functional mechanisms remains a major challenge. varTFBridge addresses this by:

- Leveraging high-resolution FOODIE (single-molecule deaminase footprinting) TF footprints
- Integrating both common and rare variant association analyses
- Employing footprint-gene linking models (ABC-FP-Max)
- Utilizing AlphaGenome for variant effect prediction
- Prioritizing causal noncoding variants that rewire gene regulatory networks

## Key Features

### Two-Stage Analysis Framework

**Stage 1: Variant Association**
- **Common variants (MAF ≥ 0.1%)**: GWAS followed by genome-wide fine-mapping (GWFM) with SBayesRC
- **Rare variants (MAF < 0.1%)**: Footprint-based burden tests with leave-one-variant-out analysis

**Stage 2: Functional Dissection**
- **VAR2TFBS**: Predicts how variants affect TF binding affinity using JASPAR-based position weight matrices
- **ABC-FP-Max**: Links variants to target genes through Activity-by-Contact scoring adapted for TF footprints
- **AlphaGenome**: Assesses variant effects across multiple epigenomic layers (histone modifications, TF binding, chromatin accessibility)

## Results Highlights

Using 490,640 UK Biobank whole-genome sequences across 13 erythroid traits:

- K562 FOODIE footprints show ~70-fold heritability enrichment for erythroid traits (comprising <0.5% of the genome)
- Identified **206 common variants** and **18 rare variants** linked to TF binding sites and target genes
- Successfully recapitulated the causal variant **rs112233623**, revealing how disruption of GATA1/TAL1 co-binding alters CCND3 regulation to drive variation in red blood cell count

## Installation

```bash
git clone https://github.com/JasonLinjc/varTFBridge.git
cd varTFBridge
```

## Dependencies

- Python 3.8+
- FIMO (from MEME Suite)
- AlphaGenome (v0.5.1)
- REGENIE (v3.3)
- SBayesRC / GCTB
- Bismark (v0.24.2)
- Trim Galore (v0.6.10)

## Data

- **FOODIE footprints**: K562 and GM12878 cell lines available in this repository
- **TF binding motifs**: JASPAR 2025 core non-redundant vertebrates
- **UK Biobank WGS data**: 490,640 participants (requires authorized access)

## Usage

### FOODIE Footprint Calling
```bash
# See https://github.com/sunneyxie-lab/bulk-foodie-pipeline
```

### VAR2TFBS Analysis
VAR2TFBS takes as input:
- Reference genome
- Variant list (BED format)
- FOODIE footprints (BED format)
- Position Weight Matrices (PWMs)

And outputs variant effects on TF binding sites (creation, disruption, increase, decrease).

### ABC-FP-Max Predictions
Adapted from the [ABC model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) to link TF footprints to target genes using:
- Chromatin accessibility (ATAC-seq)
- Hi-C contact frequency
- Footprint activity scores

## Methods

| Component | Description |
|-----------|-------------|
| **FOODIE** | Single-molecule deaminase footprinting for near-base-resolution TF binding detection |
| **GWFM** | Genome-wide fine-mapping using SBayesRC for posterior inclusion probabilities |
| **VAR2TFBS** | FIMO-based scanning to assess variant effects on TF binding motifs |
| **ABC-FP-Max** | Footprint-to-gene linkage scoring combining activity and chromatin contact |
| **AlphaGenome** | Deep learning model for cell-type-specific variant effect prediction |


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback, please open an issue on GitHub or contact Jiecong Lin.
