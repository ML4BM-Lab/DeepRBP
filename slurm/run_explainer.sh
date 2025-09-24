#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_explainer🚀
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_explainer.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

# ---------- Config ----------
CONFIG_PATH="/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml"
CKPT_PATH="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt"
SCALER_DIR="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data"
OUTPUT_DIR="/scratch/jsanchoz/DeepRBP/output/results/explainability/explainer_dl_kout_t_stat"
SELECT_CATEGORIES="Acute_Myeloid_Leukemia,Kidney_Chromophobe,Liver_Hepatocellular_Carcinoma"

ANALYZE_HL="${ANALYZE_HL:-false}" # check if calculate scores also for last hidden layer (by default: false)

module purge
module load Miniforge3
source activate /data/jsanchoz/conda-env/DeepRBP

python --version
conda info --envs

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

CMD=(python3.9 -m deeprbp.explainability_module.main_explainer
  --config_path "$CONFIG_PATH"
  --model_ckpt_path "$CKPT_PATH"
  --scaler_dir "$SCALER_DIR"
  --output_dir "$OUTPUT_DIR"
  --select_category "$SELECT_CATEGORIES"
)

if [[ "$ANALYZE_HL" == "true" ]]; then
  CMD+=(--analyze_hidden_layer)
  echo "[INFO] analyze_hidden_layer: ENABLED"
else
  echo "[INFO] analyze_hidden_layer: DISABLED"
fi

echo "Running:"
printf ' %q' "${CMD[@]}"; echo

# ---------- Execute ----------
"${CMD[@]}"

echo "########################################"
echo "Finished at: $(date)"
echo "########################################"

# python3.9 -m deeprbp.explainability_module.main_explainer \
#   --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml" \
#   --model_ckpt_path "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt" \
#   --scaler_dir "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data" \
#   --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainability/explainer_dl_kout_t_stat" \
#   --select_category "Acute_Myeloid_Leukemia,Kidney_Chromophobe,Liver_Hepatocellular_Carcinoma"
