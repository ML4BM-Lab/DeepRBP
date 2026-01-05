#!/usr/bin/env bash
set -euo pipefail

#########################
# CONFIG PATHS
#########################

REPO_ROOT="/scratch/jsanchoz/DeepRBP"

# Código R
R_PLOT_DIR="${REPO_ROOT}/src/deeprbp/explainability_module/real_knockdowns/plot_results"

# Inputs
GETBM="${REPO_ROOT}/data/training_module/selected_genes_rbps/getBM.csv"
SCORES_BASE="${REPO_ROOT}/output/results/real_knockdowns/dl_kout_t_stat"
LIMMA_BASE="${REPO_ROOT}/output/results/real_knockdowns/de_limma/logFC1_pvalue0.01"

# Output
OUT_DIR="${REPO_ROOT}/output/results/real_knockdowns/figures_scores"
mkdir -p "${OUT_DIR}"

#########################
# PARAMETERS
#########################

EXPERIMENT_LIST="PRJEB39343,GSE75491,GSE136366"
RBP_INTEREST_LIST="MBNL1,RBM47,TARDBP"

LOGFC_THRESH=1
P_CUT_TYPE="pvalue"
P_CUT_VALUE=0.01

echo "▶️ Post-score analysis"
echo "   Scores: ${SCORES_BASE}"
echo "   Limma:  ${LIMMA_BASE}"
echo "   Output: ${OUT_DIR}"
echo

#########################
# RUN
#########################

Rscript "${R_PLOT_DIR}/generate_realkd_plot.R" \
  --output_dir "${OUT_DIR}" \
  --path_getBM "${GETBM}" \
  --path_results "${SCORES_BASE}" \
  --path_de_limma "${LIMMA_BASE}" \
  --experiment_list "${EXPERIMENT_LIST}" \
  --rbp_interest_list "${RBP_INTEREST_LIST}" \
  --logfc_thresh "${LOGFC_THRESH}" \
  --p_cut_type "${P_CUT_TYPE}" \
  --p_cut_value "${P_CUT_VALUE}"

echo "✅ Post-score analysis finished successfully"
