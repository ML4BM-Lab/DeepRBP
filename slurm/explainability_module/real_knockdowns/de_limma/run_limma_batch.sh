#!/usr/bin/env bash
set -euo pipefail

#########################
# CONFIG PATHS
#########################

# Raíz del repositorio
REPO_ROOT="/scratch/jsanchoz/DeepRBP"

# ----------------------
# Código R (limma core)
# ----------------------
R_DE_LIMMA_DIR="${REPO_ROOT}/src/deeprbp/explainability_module/de_limma"

# ----------------------
# Código R (post-plots)
# ----------------------
R_PLOT_DIR="${REPO_ROOT}/src/deeprbp/explainability_module/real_knockdowns/plot_results"

# ----------------------
# Datos
# ----------------------
BASE_DATA="${REPO_ROOT}/data/explainability_module/real_kds"
GETBM="${REPO_ROOT}/data/training_module/selected_genes_rbps/getBM.csv"

# ----------------------
# Salida
# ----------------------
BASE_OUT="${REPO_ROOT}/output/results/real_knockdowns/de_limma"

#########################
# THRESHOLDS
#########################

LOGFC_THRESH=1
P_CUT_TYPE="pvalue"      # fdr | pvalue
P_CUT_VALUE=0.01

OUT_ROOT="${BASE_OUT}/logFC${LOGFC_THRESH}_${P_CUT_TYPE}${P_CUT_VALUE}"
mkdir -p "${OUT_ROOT}"

echo "▶️ Limma config:"
echo "   |logFC| > ${LOGFC_THRESH}"
echo "   ${P_CUT_TYPE} < ${P_CUT_VALUE}"
echo "   Output → ${OUT_ROOT}"
echo

#########################
# FUNCTION
#########################
run_limma() {
  local dataset="$1"
  local ctrl="$2"
  local kd="$3"
  local legend="$4"

  local out_dir="${OUT_ROOT}/${dataset}"

  echo "▶️ Running limma for ${dataset}"

  Rscript "${R_DE_LIMMA_DIR}/run_voom-limma_kd.R" \
    --path "${BASE_DATA}/${dataset}" \
    --output_dir "${out_dir}" \
    --path_getBM "${GETBM}" \
    --condition_control "${ctrl}" \
    --condition_kd "${kd}" \
    --show_legend "${legend}" \
    --logfc_thresh "${LOGFC_THRESH}" \
    --p_cut_type "${P_CUT_TYPE}" \
    --p_cut_value "${P_CUT_VALUE}"
}

#########################
# DATASETS
#########################

run_limma "PRJEB39343" "HFE145siNeg_1" "HFE145siMBNL1" FALSE
run_limma "GSE75491"   "Control"      "RBM47 KD"      FALSE
run_limma "GSE136366"  "Rescued_tdp43" "tdp43_ko"      TRUE

#########################
# OPTIONAL: CONCATENATE VOLCANOS
#########################

echo "▶️ Concatenating volcano plots"

Rscript "${R_PLOT_DIR}/just_concatenate_plots.R" \
  --base_dir "${OUT_ROOT}" \
  --logfc_thresh "${LOGFC_THRESH}" \
  --p_cut_type "${P_CUT_TYPE}" \
  --p_cut_value "${P_CUT_VALUE}"

echo "✅ Limma batch finished successfully"
