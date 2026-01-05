#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_explainer
#SBATCH --cpus-per-task=4
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o output/logs/run_explainer_%j.out
#SBATCH -e output/logs/run_explainer_%j.err
#SBATCH --mail-type=START,END,FAIL
#SBATC﻿H --mail-user=you@yourdomain.com

set -euo pipefail

echo "########################################"
echo "Starting job at: $(date)"
echo "Local time (Europe/Madrid): $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

# ----------------------------
# User config (edit these)
# ----------------------------
# Path to your explainability YAML (single-run mode; no categories needed)
CONFIG_PATH="src/deeprbp/configs/config_model_explain.yaml"

# Pretrained (or your own) model + scaler bundle
CKPT_PATH="pretrained_model/model.ckpt"
SCALER_DIR="pretrained_model"

# Your dataset folder (must contain RBPs_log2p_tpm.csv, trans_log2p_tpm.csv, gn_tpm.csv
# phenotype_metadata.csv is optional and ignored in single-run mode)
DATASET_DIR="data/my_dataset"

# Where results will be saved
OUTPUT_DIR="output/results/explainer_all_samples"

# Toggle hidden-layer attributions (DeepLIFT only)
ANALYZE_HL="${ANALYZE_HL:-false}"

# Optional: save per-sample TxRBP tensor (Tx × RBP × Samples)
# Controlled by YAML key save_per_sample_scores: true/false

# ----------------------------
# Environment setup (adapt to your cluster)
# ----------------------------
module purge
module load Miniforge3

# Activate your DeepRBP environment
# (edit to your env path/name)
source activate /data/jsanchoz/conda-env/DeepRBP

python --version
conda info --envs

# Make sure the package is importable
export PYTHONPATH="$(pwd)/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('DeepRBP import OK')"

# ----------------------------
# Sanity checks
# ----------------------------
mkdir -p output/logs "$OUTPUT_DIR"

for f in "RBPs_log2p_tpm.csv" "trans_log2p_tpm.csv" "gn_tpm.csv"; do
  if [[ ! -f "${DATASET_DIR}/${f}" ]]; then
    echo "[ERROR] Missing required file: ${DATASET_DIR}/${f}"
    exit 1
  fi
done

# ----------------------------
# Run
# ----------------------------
# NOTE: single-run mode (no --select_category)
# Make sure your YAML points test_path_files to $DATASET_DIR
# or override it by editing the YAML before submission.
CMD=(run-deeprbp-explainer
  --config_path "$CONFIG_PATH"
  --model_ckpt_path "$CKPT_PATH"
  --scaler_dir "$SCALER_DIR"
  --output_dir "$OUTPUT_DIR"
)

if [[ "$ANALYZE_HL" == "true" ]]; then
  CMD+=(--analyze_hidden_layer)
  echo "[INFO] analyze_hidden_layer: ENABLED"
else
  echo "[INFO] analyze_hidden_layer: DISABLED"
fi

echo "Running:"
printf ' %q' "${CMD[@]}"; echo

"${CMD[@]}"

echo "########################################"
echo "Finished at: $(date)"
echo "########################################"
