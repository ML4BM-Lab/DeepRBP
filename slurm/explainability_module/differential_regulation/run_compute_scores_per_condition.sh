#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/scores_%x.out

# ============================================================
# Usage:
# sbatch run_compute_scores_per_condition.sh \
#   <DATASET_ID> <CONFIG_PATH> <OUTPUT_DIR> <SELECT_CATEGORY> <CONDITIONS>
#
# Example:
# sbatch run_compute_scores_per_condition.sh \
#   TCGA \
#   src/deeprbp/configs/config_diff_regulation.yaml \
#   /scratch/jsanchoz/DeepRBP/output/diff_reg/TCGA-Liver \
#   Liver_Hepatocellular_Carcinoma \
#   Primary_Tumor,Solid_Tissue_Normal
# ============================================================

set -euo pipefail

# ----------------------------
# Parse arguments
# ----------------------------
if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <DATASET_ID> <CONFIG_PATH> <OUTPUT_DIR> <SELECT_CATEGORY> <CONDITIONS>"
  exit 1
fi

DATASET_ID="$1"
CONFIG_PATH="$2"
OUTPUT_DIR="$3"
SELECT_CATEGORY="$4"
CONDITIONS="$5"

# ----------------------------
# Fixed configuration
# ----------------------------
CKPT_PATH="/scratch/jsanchoz/DeepRBP/pretrained_model/model.ckpt"

echo "########################################"
echo "Dataset        : ${DATASET_ID}"
echo "Config         : ${CONFIG_PATH}"
echo "Output dir     : ${OUTPUT_DIR}"
echo "Category       : ${SELECT_CATEGORY}"
echo "Conditions     : ${CONDITIONS}"
echo "Started at     : $(date)"
echo "Local time     : $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

# ----------------------------
# Environment setup
# ----------------------------
module purge
module load Miniforge3
conda activate DeepRBP

python --version
export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python3.9 -c "import deeprbp; print('DeepRBP import OK')"

# ----------------------------
# Run
# ----------------------------
CMD=(run-deeprbp-differential-regulation-scores
  --config_path "$CONFIG_PATH"
  --model_ckpt_path "$CKPT_PATH"
  --output_dir "$OUTPUT_DIR"
  --select_category "$SELECT_CATEGORY"
  --conditions "$CONDITIONS"
)

echo "Running command:"
printf ' %q' "${CMD[@]}"; echo
"${CMD[@]}"

echo "########################################"
echo "Finished at: $(date)"
echo "########################################"
