# varTFBridge: genome-wide maps of transcription factor footprints identify noncoding variants rewiring gene regulatory networks

<p align="center">
  <img src="workflow.png" alt="Workflow" width="800">
</p>

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
- Identified **248 common variants** and **19 rare variants** linked to TF binding sites and target genes
- Successfully recapitulated the causal variant **rs112233623**, revealing how disruption of GATA1/TAL1 co-binding alters CCND3 regulation to drive variation in red blood cell count

## Installation

```bash
git clone https://github.com/JasonLinjc/varTFBridge.git
cd varTFBridge
```

## Dependencies

**Python packages** (install via `pip install -r requirements.txt`):
- pandas, numpy, pyfaidx, kipoiseq, memelite, tqdm

**External tools**:
- Python 3.8+
- bedtools (must be on PATH)
- FIMO (from MEME Suite)
- AlphaGenome (v0.5.1)
- REGENIE (v3.3)
- SBayesRC / GCTB
- Bismark (v0.24.2)
- Trim Galore (v0.6.10)

## Scripts

| Pipeline               | Script                                       | Description                                                                         |
|------------------------|----------------------------------------------|-------------------------------------------------------------------------------------|
| Data preparation       | `scripts/comvar_liftover_snpRes.py`          | Liftover GWFM `.snpRes` files from hg19 to hg38 coordinates                        |
| Data preparation       | `scripts/comvar_filter_credible_set.py`      | Filter variants by PIP threshold and annotate with LCS credible-set info            |
| Common variant VAR2TFBS | `scripts/comvar_overlap_foodie_footprints.py` | Step 1 — Overlap GWFM common variants with FOODIE footprint BED files               |
| Common variant VAR2TFBS | `scripts/comvar_var2tfbs.py`                  | Step 2 — Predict variant effects on TF binding via FIMO motif scanning              |
| Rare variant VAR2TFBS  | `scripts/rarevar_var2tfbs.py`                | Identify driver rare variants from burden-test LOO and predict TF binding effects   |
| Variant-to-gene with ABC-FP-Max  | `scripts/link_var2gene.py`                   | Link variants to target genes via ABC-FP-Max footprint-gene scores                  |

## Data

- **FOODIE footprints**: K562 and GM12878 cell lines available in this repository
- **TF binding motifs**: JASPAR 2025 core non-redundant vertebrates
- **UK Biobank WGS data**: 490,640 participants (requires authorized access)

## Usage

### FOODIE Footprint Calling

```bash
# See https://github.com/sunneyxie-lab/bulk-foodie-pipeline
```

### GWFM Fine-Mapping Pipeline

Genome-wide fine-mapping (GWFM) uses [GCTB 2.5.4](https://gctbhub.cloud.edu.au/software/gctb/#Genome-wideFine-mappinganalysis) (`--impute-summary` + `--gwfm`) to perform Bayesian fine-mapping on GWAS summary statistics.

**Input pipeline**: REGENIE GWAS (hg38) → rsID mapping → liftover to hg19 → `.ma` format (MHC excluded)

**Output files per trait**:

| File          | Description                    | Key columns                                                                                    |
|---------------|--------------------------------|------------------------------------------------------------------------------------------------|
| **`.snpRes`** | Genome-wide SNP results (hg19) | Index, Name, Chrom, Position, A1, A2, A1Frq, A1Effect, SE, VarExplained, PEP, Pi1–Pi5, PIP   |
| **`.lcs`**    | Local credible sets            | CS, Size, PIP, PGV, PGVenrich, PEP, SNP (comma-separated), ENSGID_hg19, GeneName_hg19        |
| **`.lcsRes`** | Credible set summary           | PIP/PEP thresholds, # sets, avg size, estimated causal variants, variance explained            |

**Blood cell traits analysed** (13 erythroid + 2 others):

| Trait   | Description                                | Trait   | Description                                |
|---------|--------------------------------------------|---------|--------------------------------------------|
| HC      | Hemoglobin concentration                   | MCH     | Mean corpuscular hemoglobin                |
| HP      | Hemoglobin percentage                      | MCHC    | Mean corpuscular hemoglobin concentration  |
| HLDRC   | High light scatter reticulocyte count      | MCV     | Mean corpuscular volume                    |
| HLSRP   | High light scatter reticulocyte percentage | MSCV    | Mean sphered cell volume                   |
| IRF     | Immature reticulocyte fraction             | RBC     | Red blood cell count                       |
| RC      | Reticulocyte count                         | RBCDW   | RBC distribution width                     |
| RP      | Reticulocyte percentage                    |         |                                            |

### snpRes Liftover (hg19 → hg38)

The `.snpRes` files from GCTB use hg19 coordinates. Liftover to hg38 and add `Chromosome_hg38`, `Start_hg38`, `End_hg38` columns:

```bash
python scripts/comvar_liftover_snpRes.py \
    --snpres-dir data/GWFM_erythroids/snpRes \
    --chain data/reference/hg19ToHg38.over.chain.gz \
    --out-dir data/GWFM_erythroids/snpRes_hg38
```

Since all traits share the same variant set (~13M), the liftover is performed once and applied to all files. Supports multiprocessing with `--workers N`.

### Common Variant Credible Set Preparation

Filter snpRes_hg38 common variants by PIP threshold and annotate with LCS credible set information (PEP_cs, CS_id):

```bash
python scripts/comvar_filter_credible_set.py \
    --snpres-dir data/GWFM_erythroids/snpRes_hg38 \
    --lcs-dir data/GWFM_erythroids/lcs \
    --out-dir data/GWFM_erythroids/credible_set_snpRes \
    --pip-threshold 0.1
```

Output: `data/GWFM_erythroids/credible_set_snpRes/{trait}_credible_set_hg38.csv` with columns including PEP_cs (credible-set-level PEP from LCS) and CS_id (credible set ID).

### VAR2TFBS Analysis

#### Step 1: Overlap Common Variants with FOODIE Footprints

Identify GWFM common variants that fall within FOODIE TF footprints. The script accepts both credible set CSV and snpRes (tab-delimited) formats, with auto-detection of delimiter and column name mapping. Outputs per-trait CSVs with PIP, PEP, PEP_cs, and CS_id annotations.

**Using credible set files** (filtered by PIP/PEP thresholds):

```bash
python scripts/comvar_overlap_foodie_footprints.py \
    --snp-dir data/GWFM_erythroids/credible_set_snpRes \
    --footprint-dir data/FOODIE_footprints \
    --lcs-dir data/GWFM_erythroids/lcs \
    --out-dir results/comvar_footprint_overlap
```

**Using full snpRes_hg38 files** (all ~13M common variants per trait):

```bash
python scripts/comvar_overlap_foodie_footprints.py \
    --snp-dir data/GWFM_erythroids/snpRes_hg38 \
    --snp-suffix .snpRes \
    --footprint-dir data/FOODIE_footprints \
    --lcs-dir data/GWFM_erythroids/lcs \
    --out-dir results/comvar_footprint_overlap_snpRes
```

| Option             | Default                            | Description                                         |
|--------------------|------------------------------------|-----------------------------------------------------|
| `--snp-dir`        | (required)                         | Directory of GWFM common variant files (CSV or snpRes) |
| `--footprint-dir`  | (required)                         | Directory of FOODIE footprint BED files             |
| `--out-dir`        | `./results/comvar_footprint_overlap` | Output directory                                  |
| `--pip-threshold`  | `0`                                | Minimum PIP to include per trait                    |
| `--snp-suffix`     | `_credible_set_hg38.csv`          | Suffix to strip for trait names                     |
| `--lcs-dir`        | (optional)                         | Directory of `.lcs` files for PEP_cs/CS_id annotation |

Output: per-trait CSVs (`{out_dir}/{footprint}/{trait}_{footprint}.csv`) with columns SNP, Chromosome, Start, End, A1, A2, freq, PIP, PEP, PEP_cs, CS_id, footprint_region, plus combined BED files per footprint.

#### Step 2: Predict Variant Effects on TF Binding

Trait-agnostic: takes the merged BED from Step 1 (all unique variants across traits) and predicts variant effects on TF binding using FIMO motif scanning against JASPAR PWMs. Classifies TF binding changes as Create, Disrupt, Increase, Decrease, or Unchange.

```bash
python scripts/comvar_var2tfbs.py \
    --input-bed results/comvar_footprint_overlap_credible/GWFM_variants_in_K562.merged.hg38.bed \
    --allele-src results/comvar_footprint_overlap_credible/K562.merged.hg38 \
    --ref-genome data/reference/hg38.fa \
    --jaspar-meme data/JASPAR_MEME/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
    --out-dir results/comvar_var2tfbs_results
```

| Option            | Default                              | Description                                                    |
|-------------------|--------------------------------------|----------------------------------------------------------------|
| `--input-bed`     | (required)                           | Merged BED from Step 1 (e.g. `GWFM_variants_in_K562.merged.hg38.bed`) |
| `--allele-src`    | (required)                           | Directory of per-trait CSVs or single CSV with SNP, A1, A2    |
| `--ref-genome`    | (required)                           | Path to hg38.fa reference genome                              |
| `--jaspar-meme`   | (required)                           | Path to JASPAR MEME motif file                                |
| `--out-dir`       | `./results/comvar_var2tfbs_results`  | Output directory                                              |
| `--ext-bp`        | `30`                                 | Sequence extension in bp around footprint                     |
| `--fimo-threshold` | `0.0001`                            | FIMO p-value threshold                                        |

Output: `{out_dir}/{cell}_var2tfbs.csv` with ref/alt FIMO hits, TF change classification (Create/Disrupt/Increase/Decrease/Unchange), and FASTA files in `{out_dir}/fasta/`.

### Rare Variant VAR2TFBS Analysis

Identifies driver rare variants from footprint-based burden test leave-one-out (LOO) analysis and predicts their effects on TF binding. For each significant footprint (Bonferroni-corrected p < 0.05/N_footprints), the driver variant is the one whose removal causes the largest increase in burden test p-value. Footprints where no variant has more than 30 carriers (MAC > 30) are excluded.

```bash
python scripts/rarevar_var2tfbs.py \
    --burden-dir data/burdentest_erythroids \
    --loo-file data/leaveoneout_results/K562.leave_one_out.all_traits.20251120.csv \
    --ref-genome data/reference/hg38.fa \
    --jaspar-meme data/JASPAR_MEME/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
    --out-dir results/rarevar_var2tfbs_results
```

| Option            | Default                                | Description                                               |
|-------------------|----------------------------------------|-----------------------------------------------------------|
| `--burden-dir`    | (required)                             | Directory of burden test result Excel files (one per trait) |
| `--loo-file`      | (required)                             | Leave-one-out results CSV (all traits combined)           |
| `--ref-genome`    | (required)                             | Path to hg38.fa reference genome                          |
| `--jaspar-meme`   | (required)                             | Path to JASPAR MEME motif file                            |
| `--out-dir`       | `./results/rarevar_var2tfbs_results`   | Output directory                                          |
| `--sig-threshold` | Bonferroni (0.05/N)                    | Burden test significance threshold                        |
| `--min-carrier`   | `30`                                   | Minimum MAC for at least one variant in footprint         |
| `--ext-bp`        | `30`                                   | Sequence extension in bp around footprint                 |
| `--fimo-threshold` | `0.0001`                              | FIMO p-value threshold                                    |

Output: `driver_variants_summary.csv` (driver variants per trait-footprint) and `K562_rarevar_var2tfbs.csv` (TF binding effect predictions).

### ABC-FP-Max Predictions

Adapted from the [ABC model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) to link TF footprints to target genes using chromatin accessibility (ATAC-seq), Hi-C contact frequency, and footprint activity scores. See [ABC-FP/README.md](ABC-FP/README.md) for license details.

**1. Configure biosamples** — edit `ABC-FP/config/config_FOODIE_ATAC.tsv` with paths to your input files:

| Column           | Description                                      |
|------------------|--------------------------------------------------|
| `biosample`      | Sample name (e.g. `K562_FOODIE_ATAC`)            |
| `narrowPeaks`    | FOODIE footprint BED file                        |
| `ATAC`           | ATAC-seq BAM file                                |
| `HiC_file`       | Hi-C contact matrix (`.hic` format)              |
| `HiC_type`       | Hi-C file type (`hic`)                           |
| `HiC_resolution` | Hi-C resolution in bp (e.g. `5000`)              |

**2. Configure parameters** — edit `ABC-FP/config/config.yaml` to set `results_dir` and reference file paths.

**3. Run the pipeline**:

```bash
cd ABC-FP
snakemake -j <cores> --use-conda
```

Output: `{results_dir}/{biosample}/Predictions/EnhancerPredictionsAllPutative.tsv.gz` containing ABC scores for all putative enhancer-gene links.

### Variant-to-Gene Linking

Links variants to target genes by bridging: **variant → footprint → enhancer → gene** using ABC-FP-Max scores. For each variant, the target gene is the one with the highest ABC-FP score. Also annotates with TF binding changes and cell-type-specific TF RNA expression.

```bash
python scripts/link_var2gene.py \
    --var2tfbs results/comvar_var2tfbs_results/K562_var2tfbs.csv \
    --footprint-bed data/FOODIE_footprints/K562.merged.hg38.bed \
    --enhancer-bed ABC_FP_results/K562_FOODIE_ATAC/Neighborhoods/EnhancerList.bed \
    --abc-predictions ABC_FP_results/K562_FOODIE_ATAC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --tf-expr data/gene_expr/TF_K562_GM12878_expression.csv \
    --cell K562 \
    --out-dir results/var2gene_results
```

| Option              | Default                        | Description                                            |
|---------------------|--------------------------------|--------------------------------------------------------|
| `--var2tfbs`        | (required)                     | VAR2TFBS output CSV (common or rare)                   |
| `--footprint-bed`   | (required)                     | FOODIE footprint BED file                              |
| `--enhancer-bed`    | (required)                     | ABC-FP EnhancerList BED from Neighborhoods             |
| `--abc-predictions` | (required)                     | ABC-FP EnhancerPredictionsAllPutative TSV              |
| `--tf-expr`         | (optional)                     | TF RNA expression CSV for cell-type annotation         |
| `--cell`            | `K562`                         | Cell type for TF expression column lookup              |
| `--abc-threshold`   | `0`                            | Minimum ABC score for gene assignment                  |
| `--prefix`          | `comvar`                       | Output filename prefix (`comvar` or `rarevar`)         |
| `--out-dir`         | `./results/var2gene_results`   | Output directory                                       |

Output: `{cell}_{prefix}_ABC-FP-Full.csv` (variant-TF-gene table with rsID, TF, TF_change, TF expression, TargetGene, ABC.Score.FP, distance) and `{cell}_{prefix}_ABC-FP-Max.csv` (one row per variant with top ABC-FP-Max gene).

## Methods

| Component        | Description                                                                                             |
|------------------|---------------------------------------------------------------------------------------------------------|
| **FOODIE**       | Single-molecule deaminase footprinting for near-base-resolution TF binding detection                    |
| **GWFM**         | Genome-wide fine-mapping using SBayesRC producing global (GCS) and local (LCS) credible sets with PIP, PEP, and PGV |
| **VAR2TFBS**     | FIMO-based scanning to assess variant effects on TF binding motifs                                      |
| **ABC-FP-Max**   | Footprint-to-gene linkage scoring combining activity and chromatin contact                              |
| **AlphaGenome**  | Deep learning model for cell-type-specific variant effect prediction                                    |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback, please open an issue on GitHub or contact Jiecong Lin (jieconglin@cpl.ac.cn) and Yajie Zhao (yajiezhao@cpl.ac.cn).

<p align="center">
  <img src="logo.png" alt="Logo" width="400">
</p>
