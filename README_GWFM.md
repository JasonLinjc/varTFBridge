# MCH Fine-mapping Project

Genetic fine-mapping analysis of blood cell traits using GCTB (GCTA-COJO). This project processes GWAS summary statistics, performs variant mapping and coordinate liftover, and runs Bayesian fine-mapping to identify credible sets of causal variants.

## Fine-mapping Usage Summary

**Tool:** GCTB 2.5.4 (`--impute-summary` + `--gwfm`)

**Input:** REGENIE GWAS summary statistics (hg38) → rsID mapping → liftover to hg19 → `.ma` format (MHC excluded)

**Run:** `sbatch code/sbatch_finemapping_all.sh` (or `step1_run_gbmm.sh <gwas_hg19.txt> <ncpu>` for a single trait)

**Parameters:** PIP threshold 0.9, PEP threshold 0.7 (configurable in GCTB)

---

## Output Summary

| File | Description | Key columns / fields |
|------|-------------|----------------------|
| **`.lcs`** | SNP-level credible sets | CS, Size, PIP, PGV, PGVenrich, PEP, SNP, ENSGID_hg19, GeneName_hg19 |
| **`.lcsRes`** | Credible set summary | PIP/PEP thresholds, # 1-SNP vs multi-SNP sets, total SNPs, avg set size, estimated causal variants, power, variance explained |
| **`.gcs`** | Gene-level credible sets | Gene-centric credible sets |
| **`.snpRes`** | Genome-wide SNP results | Index, Name, Chrom, Position, A1, A2, A1Frq, A1Effect, SE, VarExplained, PEP, Pi1–Pi5, PIP, GelmanRubin_R |
| **`.gwas.ma`** | GCTB input | SNP, A1, A2, freq, b, se, p, N |
| **`.imputed.ma`** | Imputed summary stats | Same format as .gwas.ma, LD-imputed |

**Collected folders:** `collected_lcs/` and `collected_snpRes/` hold all trait outputs in one place.

---

## Project Structure

```
.
├── code/                    # Scripts and utilities
│   ├── get_ma.pl           # Convert REGENIE output to GCTB .ma format (excludes MHC)
│   ├── map_to_rs_direct.sh # Map variants to dbSNP rsIDs
│   ├── convert_regenie.sh  # LiftOver hg38 → hg19
│   ├── step1_run_gbmm.sh   # GCTB imputation + fine-mapping
│   └── sbatch_finemapping_all.sh  # SLURM batch pipeline
├── result_jc/              # Fine-mapping results by trait
│   ├── <trait>/            # Per-trait outputs (e.g., HC/, MCH/, RBC/)
│   │   ├── *.lcs          # Credible set results
│   │   ├── *.lcsRes       # Summary statistics
│   │   └── *.gcs          # Gene-level credible sets
│   ├── collected_lcs/      # All .lcs files collected in one folder
│   └── collected_snpRes/  # All .snpRes files collected in one folder
├── LDSC_blood_trait_results_rmMHC_20250528/  # LDSC heritability results
└── heritability.ipynb     # LDSC result parsing and analysis
```

## Blood Cell Traits

| Trait | Description |
|-------|-------------|
| HC | Hemoglobin concentration |
| HP | Hemoglobin percentage |
| HLDRC | High light scatter reticulocyte count |
| HLSRP | High light scatter reticulocyte percentage |
| IRF | Immature reticulocyte fraction |
| LP | Low peroxidase |
| Lymphocyte | Lymphocyte count |
| MCH | Mean corpuscular hemoglobin |
| MCHC | Mean corpuscular hemoglobin concentration |
| MCV | Mean corpuscular volume |
| MSCV | Mean sphered cell volume |
| RBC | Red blood cell count |
| RBCDW | RBC distribution width |
| RC | Reticulocyte count |
| RP | Reticulocyte percentage |

## Output Files

- **`.lcs`** — Credible set results: SNP, PIP (posterior inclusion probability), PGV, gene annotations (ENSGID, GeneName)
- **`.lcsRes`** — Summary: number of credible sets, average size, estimated causal variants, power, variance explained
- **`.gcs`** — Gene-level credible sets
- **`.snpRes`** — Genome-wide SNP-level results: PIP, PEP, effect sizes, variance explained per variant

## Pipeline Overview

1. **Input**: REGENIE GWAS summary statistics (`.site_all_site_all_var_result.*.txt`)
2. **Map to rsID**: `map_to_rs_direct.sh` — annotate variants with dbSNP rsIDs
3. **LiftOver**: `convert_regenie.sh` — hg38 → hg19 for GCTB LD reference
4. **Format**: `get_ma.pl` — convert to GCTB `.ma` format, exclude MHC region
5. **Fine-mapping**: GCTB `--impute-summary` then `--gwfm` — imputation and Bayesian fine-mapping

## Script Reference

### `get_ma.pl`

Converts REGENIE GWAS summary statistics to GCTB `.ma` format. Excludes variants in the MHC region (chr6, from `MHC_region.hg19.bed`).

**Usage:**
```bash
perl get_ma.pl <gwas_file>
```

**Arguments:**
| Argument   | Description |
|------------|-------------|
| `gwas_file` | Tab-delimited REGENIE output with rsID column (from `map_to_rs_direct.sh` + liftover) |

**Required input columns:** `CHROM`, `GENPOS`, `ALLELE1`, `ALLELE0`, `A1FREQ`, `BETA`, `SE`, `p_regenie`, `N`, `rsID`

**Output:** `{trait}.gwas.ma` in current directory (trait = basename of input without extension)

**Output format:** Tab-delimited with columns `SNP`, `A1`, `A2`, `freq`, `b`, `se`, `p`, `N`

**Note:** Requires `MHC_region.hg19.bed` at `/cpl/pipeline/GCTB_finemapping/MHC_region.hg19.bed`

---

### `map_to_rs_direct.sh`

Maps variants to dbSNP rsIDs using CHROM, POS, REF, ALT. Uses `chr2refseq.txt` to convert chromosome names to RefSeq accessions for dbSNP lookup.

**Usage:**
```bash
./map_to_rs_direct.sh <input.txt> <dbsnp.vcf.gz> <output.txt>
```

**Arguments:**
| Argument   | Description |
|------------|-------------|
| `input.txt` | REGENIE GWAS summary file (tab-delimited, hg38 coordinates) |
| `dbsnp.vcf.gz` | dbSNP VCF (e.g. `GCF_000001405.40.gz`) |
| `output.txt` | Output file with rsID column appended |

**Required input columns:** `CHROM`, `GENPOS`, `ALLELE1`, `ALLELE0` (used for REF/ALT lookup)

**Output:** Same as input with extra `rsID` column; variants without dbSNP match get `.`

**Requirements:**
- `chr2refseq.txt` in current directory (chr → RefSeq mapping, e.g. `1  NC_000001.11`)
- `bcftools` in PATH

**Note:** Creates `$TMPDIR/dbsnp_map.txt` (default `/tmp/finemapping/`) for lookup; first run may be slow while building the map.

---

### `convert_regenie.sh`

Lifts GWAS coordinates from hg38 to hg19 using UCSC liftOver. Updates `CHROM` and `GENPOS` in place; other columns unchanged.

**Usage:**
```bash
./convert_regenie.sh <input_file> <chain_file> <output_file>
```

**Arguments:**
| Argument     | Description |
|--------------|-------------|
| `input_file` | Tab-delimited GWAS file (hg38); first row = header |
| `chain_file` | UCSC chain file (e.g. `hg38ToHg19.over.chain.gz`) |
| `output_file`| Output file with hg19 coordinates |

**Required input format:**
- Header row
- Column 1: `CHROM` (1–22, X, Y)
- Column 2: `GENPOS` (hg38 position)
- Column 3: Variant ID (used to track rows through liftover)

**Output:** Same structure with `CHROM` and `GENPOS` updated to hg19; unmapped variants dropped

**Unmapped variants:** Written to `tmp_unmapped.bed` in current directory

**Requirements:** `liftOver` at `/cpl/software/liftover/liftOver`

---

### `step1_run_gbmm.sh`

Runs GCTB imputation and fine-mapping. Step 1: impute summary stats; Step 2: run `--gwfm` (genome-wide fine-mapping).

**Usage:**
```bash
./step1_run_gbmm.sh <gwas_file> <ncpu>
```

**Arguments:**
| Argument   | Description |
|------------|-------------|
| `gwas_file` | REGENIE-style file in hg19 with rsID (output of `convert_regenie.sh`); used as input to `get_ma.pl` |
| `ncpu`     | Number of threads for GCTB |

**Workflow:**
1. Calls `get_ma.pl` to produce `{trait}.gwas.ma`
2. Runs `gctb --impute-summary` → `{trait}.imputed.ma`
3. Runs `gctb --gwfm` with annotation and gene map → `{trait}.lcs`, `{trait}.lcsRes`, `{trait}.gcs`

**Outputs (in current directory):**
| File | Description |
|------|-------------|
| `{trait}.gwas.ma` | GCTB input format |
| `{trait}.imputed.ma` | Imputed summary stats |
| `{trait}.lcs` | SNP-level credible sets |
| `{trait}.lcsRes` | Credible set summary |
| `{trait}.gcs` | Gene-level credible sets |

**Hardcoded paths:**
- LD matrix: `/cpl/databases/GCTB/2.5.4/blk1588`
- Annotation: `/cpl/databases/GCTB/2.5.4/annot_bolt_clean.txt`
- Gene map: `/cpl/pipeline/GCTB_finemapping/gene_map_hg38_hg19.txt`
- GCTB: `/cpl/software/GCTB/2.5.4/gctb`

**Note:** The `--gwfm` call uses trait name `RC`; for other traits you may need to edit this or pass the trait name as a parameter.

---

### `sbatch_finemapping_all.sh`

SLURM batch script that runs the full fine-mapping pipeline for all matching input files.

**Usage:**
```bash
sbatch code/sbatch_finemapping_all.sh
```

**SLURM settings:** 1 node, 32 CPUs, logs to `log/finemapping.err` and `log/finemapping.out`

**Input:** All files matching:
```
{BASE_DIR}/**/*.site_all_site_all_var_result.2025-08-06.txt
```
Default `BASE_DIR`: `/cpl/output/20250909_MCH_finemapping/run/`

**Per-file workflow:**
1. `map_to_rs_direct.sh` → `{prefix}.hg38.txt`
2. Filter to variants with rsID: `awk '$18!="."{print $0}'` → `{prefix}.hg38_rs.txt`
3. `convert_regenie.sh` → `{prefix}.hg19.txt`
4. `step1_run_gbmm.sh` → GCTB outputs

**Output directory:** `/cpl/output/20250909_MCH_finemapping/result/{prefix}/`

**Requirements:**
- Conda env `GATK` (for bcftools)
- Scripts `map_to_rs_direct.sh`, `convert_regenie.sh`, `step1_run_gbmm.sh` in `PATH` or run from directory containing them
- `log/` directory for SLURM logs

---

## Dependencies

- **GCTB** 2.5.4 (GCTA-COJO)
- **bcftools** (dbSNP lookup)
- **liftOver** (UCSC)
- **Perl** (get_ma.pl)
- LD reference: `blk1588` eigen decomposition
- Annotation: `annot_bolt_clean.txt`
- Gene map: `gene_map_hg38_hg19.txt`

## Running the Pipeline

```bash
# Submit full pipeline via SLURM
sbatch code/sbatch_finemapping_all.sh
```

The batch script processes all matching input files under `run/`, runs mapping, liftover, and fine-mapping for each trait, and writes results to `result/`.

## LDSC Heritability

LDSC results (MHC-excluded) are in `LDSC_blood_trait_results_rmMHC_20250528/`. Use `heritability.ipynb` to parse `.log` files for h², intercept, and enrichment ratios.
