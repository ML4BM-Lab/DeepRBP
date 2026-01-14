#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=deeprbp_diff_reg_scores
#SBATCH --cpus-per-task=4
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o output/logs/diff_reg_scores_%j.out
#SBATCH -e output/logs/diff_reg_scores_%j.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=you@yourdomain.com

set -euo pipefail

echo "########################################"
echo "Starting job at: $(date)"
echo "Local time (Europe/Madrid): $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

# ----------------------------
# User configuration
# ----------------------------
CONFIG_PATH="src/deeprbp/configs/config_diff_regulation.yaml"
CKPT_PATH="/path/to/deeprbp_predictor.ckpt"
OUTPUT_DIR="/scratch/jsanchoz/DeepRBP/output/diff_reg/Liver"

SELECT_CATEGORY="Liver_Hepatocellular_Carcinoma"
CONDITIONS="Primary_Tumor,Solid_Tissue_Normal"

# ----------------------------
# Environment setup
# ----------------------------
module purge
module load Miniforge3

source activate /data/jsanchoz/conda-env/DeepRBP

python --version
conda info --envs

export PYTHONPATH="$(pwd)/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python -c "import deeprbp; print('DeepRBP import OK')"

mkdir -p output/logs "$OUTPUT_DIR"

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

echo "Running:"
printf ' %q' "${CMD[@]}"; echo

"${CMD[@]}"

echo "########################################"
echo "Finished at: $(date)"
echo "########################################"