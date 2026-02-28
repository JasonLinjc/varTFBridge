#!/usr/bin/env bash
# Run the full rare variant analysis pipeline.
#
# Usage:
#   bash scripts/run_rarevar_pipeline.sh              # run all steps
#   bash scripts/run_rarevar_pipeline.sh --from 2     # resume from step 2
#   bash scripts/run_rarevar_pipeline.sh --only 1     # run step 1 only
#   bash scripts/run_rarevar_pipeline.sh --skip-alphag # skip AlphaGenome scoring (step 3)
#
# Steps:
#   1. Identify driver variants and predict TF binding effects (LOO + FIMO)
#   2. Link variants to target genes (ABC-FP-Max)
#   3. Score variants with AlphaGenome API (optional, requires API key)
#   4. Merge linkage table (rarevar2grn)
set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
PYTHON="${PYTHON:-/Users/jieconglin/miniconda3/envs/foodie_env/bin/python}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Data paths
BURDEN_DIR="$ROOT/data/burdentest_erythroids"
LOO_FILE="$ROOT/data/leaveoneout_results/K562.leave_one_out.all_traits.20251120.csv"
REF_GENOME="$ROOT/data/reference/hg38.fa"
JASPAR_MEME="$ROOT/data/JASPAR_MEME/JASPAR2024_CORE_non-redundant_pfms_meme.txt"

FOOTPRINT_BED="$ROOT/data/FOODIE_footprints/K562.merged.hg38.bed"
ENHANCER_BED="$ROOT/data/ABC_FP_results/K562_FOODIE_ATAC/Neighborhoods/EnhancerList.bed"
ABC_PRED="$ROOT/data/ABC_FP_results/K562_FOODIE_ATAC/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
TF_EXPR="$ROOT/data/gene_expr/K562_ENCFF485RIA_gene.tsv"

VAR2TFBS_DIR="$ROOT/results/rarevar_var2tfbs_results"
VAR2GENE_DIR="$ROOT/results/var2gene_results"
ALPHAG_DIR="$ROOT/results/alphag_scores"

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
    echo "=== Step 1/4: Identify driver variants + TF binding effects (LOO + FIMO) ==="
    "$PYTHON" "$ROOT/scripts/rarevar_var2tfbs.py" \
        --burden-dir "$BURDEN_DIR" \
        --loo-file "$LOO_FILE" \
        --ref-genome "$REF_GENOME" \
        --jaspar-meme "$JASPAR_MEME" \
        --out-dir "$VAR2TFBS_DIR"
    echo ""
fi

if should_run 2; then
    echo "=== Step 2/4: Link variants to target genes (ABC-FP-Max) ==="
    "$PYTHON" "$ROOT/scripts/link_var2gene.py" \
        --var2tfbs "$VAR2TFBS_DIR/${CELL}_rarevar_var2tfbs.csv" \
        --footprint-bed "$FOOTPRINT_BED" \
        --enhancer-bed "$ENHANCER_BED" \
        --abc-predictions "$ABC_PRED" \
        --tf-expr "$TF_EXPR" \
        --cell "$CELL" \
        --prefix rarevar \
        --out-dir "$VAR2GENE_DIR"
    echo ""
fi

if should_run 3; then
    if $SKIP_ALPHAG; then
        echo "=== Step 3/4: AlphaGenome scoring — SKIPPED (--skip-alphag) ==="
    elif [[ -z "${ALPHAGENOME_API_KEY:-}" ]]; then
        echo "=== Step 3/4: AlphaGenome scoring — SKIPPED (no ALPHAGENOME_API_KEY) ==="
    else
        echo "=== Step 3/4: Score variants with AlphaGenome ==="
        "$PYTHON" "$ROOT/scripts/alphag_score_variants.py" \
            --var2tfbs "$VAR2TFBS_DIR/${CELL}_rarevar_var2tfbs.csv" \
            --cell "$CELL" \
            --prefix rarevar \
            --out-dir "$ALPHAG_DIR"
    fi
    echo ""
fi

if should_run 4; then
    echo "=== Step 4/4: Merge linkage table (rarevar2grn) ==="
    "$PYTHON" "$ROOT/scripts/merge_rarevar_linkage_table.py" \
        --project-root "$ROOT"
    echo ""
fi

echo "=== Rare variant pipeline complete ==="
