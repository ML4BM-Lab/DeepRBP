#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_realkd_all_datasets
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_realkd_all_datasets.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"
echo

########################################
# ENVIRONMENT PYTHON / CONDA
########################################
module load Miniforge3
conda activate DeepRBP
echo "Python path: $(which python)"
echo "Conda env: $CONDA_DEFAULT_ENV"

############################
# RUTAS BASE
############################
CONFIG_PATH="/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_real_knockdowns.yaml"
BASE_PROC="/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds"
RESULTS_BASE="/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns"

# Entorno python
export PYTHONPATH=/scratch/jsanchoz/DeepRBP
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1

############################
# FUNCIÓN AUXILIAR
############################
run_realkd_dataset () {
  local SCALER_MODE="$1"      # tcga | fit_new
  local OUT_ROOT="$2"         # root para ese modo
  local DATASET_LABEL="$3"    # nombre de carpeta de salida (PRJEB..., GSE..._FUS, etc.)
  local PROCESSED_DIR="$4"    # carpeta processed_xxx

  local OUT="${OUT_ROOT}/${DATASET_LABEL}"
  mkdir -p "${OUT}"

  local TS
  TS="$(date +'%Y%m%d_%H%M%S')"
  local LOG="${OUT}/${TS}_${DATASET_LABEL}_${SCALER_MODE}_realkd.log"

  echo
  echo "========================================"
  echo "Dataset      : ${DATASET_LABEL}"
  echo "Scaler mode  : ${SCALER_MODE}"
  echo "Output dir   : ${OUT}"
  echo "Log file     : ${LOG}"
  echo "========================================"

  run-deeprbp-realkd \
    --config_path "${CONFIG_PATH}" \
    --processed_data_dir "${PROCESSED_DIR}" \
    --output_dir "${OUT}" \
    --scaler_mode "${SCALER_MODE}" 2>&1 | tee -a "${LOG}"

  echo "Exit code: ${PIPESTATUS[0]}" | tee -a "${LOG}"
}

############################
# BUCLE POR MODOS DE SCALER
############################

for SCALER_MODE in tcga fit_new; do
  if [[ "${SCALER_MODE}" == "tcga" ]]; then
    OUT_ROOT="${RESULTS_BASE}/using_tcga_scaler/dl_kout_t_stat"
  else
    OUT_ROOT="${RESULTS_BASE}/refited_scaler/dl_kout_t_stat"
  fi

  echo
  echo "########################################"
  echo ">>> Running all datasets with scaler_mode=${SCALER_MODE}"
  echo "    Output root: ${OUT_ROOT}"
  echo "########################################"
  echo

  # 1) PRJEB39343 (MBNL1 KD)
  run_realkd_dataset \
    "${SCALER_MODE}" \
    "${OUT_ROOT}" \
    "PRJEB39343" \
    "${BASE_PROC}/PRJEB39343/processed_MBNL1"

  # 2) GSE75491 (RBM47 KD)
  run_realkd_dataset \
    "${SCALER_MODE}" \
    "${OUT_ROOT}" \
    "GSE75491" \
    "${BASE_PROC}/GSE75491/processed"

  # 3) GSE136366 (TDP-43 KO)
  run_realkd_dataset \
    "${SCALER_MODE}" \
    "${OUT_ROOT}" \
    "GSE136366" \
    "${BASE_PROC}/GSE136366/processed"

done

echo
echo "########################################"
echo "Job finished at: $(date)"
echo "########################################"