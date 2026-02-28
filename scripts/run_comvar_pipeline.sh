#!/usr/bin/env bash
# Run the full common variant analysis pipeline.
#
# Usage:
#   bash scripts/run_comvar_pipeline.sh              # run all steps
#   bash scripts/run_comvar_pipeline.sh --from 3     # resume from step 3
#   bash scripts/run_comvar_pipeline.sh --only 4     # run step 4 only
#   bash scripts/run_comvar_pipeline.sh --skip-alphag # skip AlphaGenome scoring (step 6)
#
# Steps:
#   1. Liftover snpRes from hg19 to hg38
#   2. Filter credible set variants (PIP > 0.1)
#   3. Overlap variants with FOODIE footprints
#   4. Predict TF binding effects (FIMO motif scanning)
#   5. Link variants to target genes (ABC-FP-Max)
#   6. Score variants with AlphaGenome API (optional, requires API key)
#   7. Merge linkage table (comvar2grn)
#   8. Annotate PIP>0.7 variants with AlphaGenome TF ChIP scores
#   9. Extract TFBS changes missing from AlphaGenome
set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
PYTHON="${PYTHON:-/Users/jieconglin/miniconda3/envs/foodie_env/bin/python}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Data paths
SNPRES_DIR="$ROOT/data/GWFM_erythroids/snpRes"
SNPRES_HG38_DIR="$ROOT/data/GWFM_erythroids/snpRes_hg38"
CHAIN="$ROOT/data/reference/hg19ToHg38.over.chain.gz"
LCS_DIR="$ROOT/data/GWFM_erythroids/lcs"
CREDSET_DIR="$ROOT/data/GWFM_erythroids/credible_set_snpRes"

FOOTPRINT_DIR="$ROOT/data/FOODIE_footprints"
FOOTPRINT_BED="$FOOTPRINT_DIR/K562.merged.hg38.bed"

REF_GENOME="$ROOT/data/reference/hg38.fa"
JASPAR_MEME="$ROOT/data/JASPAR_MEME/JASPAR2024_CORE_non-redundant_pfms_meme.txt"

OVERLAP_DIR="$ROOT/results/comvar_footprint_overlap_credible"
VAR2TFBS_DIR="$ROOT/results/comvar_var2tfbs_results"
VAR2GENE_DIR="$ROOT/results/var2gene_results"
ALPHAG_DIR="$ROOT/results/alphag_scores"

ENHANCER_BED="$ROOT/data/ABC_FP_results/K562_FOODIE_ATAC/Neighborhoods/EnhancerList.bed"
ABC_PRED="$ROOT/data/ABC_FP_results/K562_FOODIE_ATAC/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
TF_EXPR="$ROOT/data/gene_expr/K562_ENCFF485RIA_gene.tsv"

CELL="K562"

# ── Parse arguments ──────────────────────────────────────────────────────────
FROM_STEP=1
ONLY_STEP=0
SKIP_ALPHAG=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)    FROM_STEP="$2"; shift 2 ;;
        --only)    ONLY_STEP="$2"; shift 2 ;;
        --skip-alphag) SKIP_ALPHAG=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

should_run() {
    local step=$1
    if [[ $ONLY_STEP -ne 0 ]]; then
        [[ $step -eq $ONLY_STEP ]]
    else
        [[ $step -ge $FROM_STEP ]]
    fi
}

# ── Pipeline ─────────────────────────────────────────────────────────────────

if should_run 1; then
    echo "=== Step 1/9: Liftover snpRes (hg19 → hg38) ==="
    "$PYTHON" "$ROOT/scripts/comvar_liftover_snpRes.py" \
        --snpres-dir "$SNPRES_DIR" \
        --chain "$CHAIN" \
        --out-dir "$SNPRES_HG38_DIR"
    echo ""
fi

if should_run 2; then
    echo "=== Step 2/9: Filter credible set (PIP > 0.1) ==="
    "$PYTHON" "$ROOT/scripts/comvar_filter_credible_set.py" \
        --snpres-dir "$SNPRES_HG38_DIR" \
        --lcs-dir "$LCS_DIR" \
        --out-dir "$CREDSET_DIR" \
        --pip-threshold 0.1
    echo ""
fi

if should_run 3; then
    echo "=== Step 3/9: Overlap variants with FOODIE footprints ==="
    "$PYTHON" "$ROOT/scripts/comvar_overlap_foodie_footprints.py" \
        --snp-dir "$CREDSET_DIR" \
        --footprint-dir "$FOOTPRINT_DIR" \
        --out-dir "$OVERLAP_DIR" \
        --snp-suffix "_credible_set_hg38.csv" \
        --lcs-dir "$LCS_DIR"
    echo ""
fi

if should_run 4; then
    echo "=== Step 4/9: Predict TF binding effects (FIMO) ==="
    "$PYTHON" "$ROOT/scripts/comvar_var2tfbs.py" \
        --input-bed "$OVERLAP_DIR/GWFM_variants_in_${CELL}.merged.hg38.bed" \
        --allele-src "$OVERLAP_DIR/${CELL}.merged.hg38" \
        --ref-genome "$REF_GENOME" \
        --jaspar-meme "$JASPAR_MEME" \
        --out-dir "$VAR2TFBS_DIR"
    echo ""
fi

if should_run 5; then
    echo "=== Step 5/9: Link variants to target genes (ABC-FP-Max) ==="
    "$PYTHON" "$ROOT/scripts/link_var2gene.py" \
        --var2tfbs "$VAR2TFBS_DIR/${CELL}_var2tfbs.csv" \
        --footprint-bed "$FOOTPRINT_BED" \
        --enhancer-bed "$ENHANCER_BED" \
        --abc-predictions "$ABC_PRED" \
        --tf-expr "$TF_EXPR" \
        --cell "$CELL" \
        --prefix comvar \
        --out-dir "$VAR2GENE_DIR"
    echo ""
fi

if should_run 6; then
    if $SKIP_ALPHAG; then
        echo "=== Step 6/9: AlphaGenome scoring — SKIPPED (--skip-alphag) ==="
    elif [[ -z "${ALPHAGENOME_API_KEY:-}" ]]; then
        echo "=== Step 6/9: AlphaGenome scoring — SKIPPED (no ALPHAGENOME_API_KEY) ==="
    else
        echo "=== Step 6/9: Score variants with AlphaGenome ==="
        "$PYTHON" "$ROOT/scripts/alphag_score_variants.py" \
            --var2tfbs "$VAR2TFBS_DIR/${CELL}_var2tfbs.csv" \
            --cell "$CELL" \
            --allele-src "$OVERLAP_DIR/${CELL}.merged.hg38" \
            --ref-genome "$REF_GENOME" \
            --prefix comvar \
            --out-dir "$ALPHAG_DIR"
    fi
    echo ""
fi

if should_run 7; then
    echo "=== Step 7/9: Merge linkage table (comvar2grn) ==="
    "$PYTHON" "$ROOT/scripts/merge_comvar_linkage_table.py" \
        --project-root "$ROOT"
    echo ""
fi

if should_run 8; then
    echo "=== Step 8/9: Annotate PIP>0.7 with AlphaGenome TF ChIP ==="
    "$PYTHON" "$ROOT/scripts/annotate_comvar_alphag_tf_chip.py" \
        --project-root "$ROOT"
    echo ""
fi

if should_run 9; then
    echo "=== Step 9/9: Extract TFBS changes missing from AlphaGenome ==="
    "$PYTHON" "$ROOT/scripts/extract_var2tfbs_extra.py" \
        --project-root "$ROOT"
    echo ""
fi

echo "=== Common variant pipeline complete ==="
