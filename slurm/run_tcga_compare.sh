#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=tcga_compare
#SBATCH --cpus-per-task=1
#SBATCH --mem=75gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_tcga_compare.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

echo "########################################"
echo " TCGA comparative regulation job started "
echo " Started at: $(date)"
echo " Local time (Hondarribia): $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

# ---------- CONFIG ----------
CONFIG_PATH="/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml"
CKPT_PATH="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt"
SCALER_DIR="/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data"
OUTPUT_DIR="/scratch/jsanchoz/DeepRBP/output/results/explain_tcga_compare/LIHC"
SELECT_CATEGORY="Liver_Hepatocellular_Carcinoma"
NORMAL_CONDITION="Solid_Tissue_Normal"
TUMOR_CONDITION="Primary_Tumor"
GETBM_PATH="/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv"

module purge
module load Miniforge3
conda activate /home/jsanchoz/.conda/envs/DeepRBP

python --version
conda info --envs

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python -c "import deeprbp; print('[OK] deeprbp package found')"

########################################
# STEP 1/3 – compute scores + Wilcoxon
########################################
CMD1=(python3.9 -m deeprbp.explainability_module.tcga_normal_vs_tumor.cli.run_compute_scores_and_tests
  --config_path "$CONFIG_PATH"
  --model_ckpt_path "$CKPT_PATH"
  --scaler_dir "$SCALER_DIR"
  --output_dir "$OUTPUT_DIR"
  --select_category "$SELECT_CATEGORY"
  --normal_condition "$NORMAL_CONDITION"
  --tumor_condition "$TUMOR_CONDITION"
)

echo "[INFO] Running STEP 1 (compute_scores_and_tests):"
printf ' %q' "${CMD1[@]}"; echo
"${CMD1[@]}"

########################################
# STEP 2/3 – Section A: RBP programs
########################################
# (opcional) thresholds si quieres sobreescribir los defaults config de Section A
# TX_FDR_THR=0.01
# RBP_FDR_THR=0.01
# LOG2FC_STRICT_THR=1.0
# LOG2FC_MODERATE_THR=0.5
# PRIOR_PROP_SIG_THR=0.20
# PRIOR_EXTREME_PROP_SIG_MAX=0.05
# PRIOR_TOP_N_GLOBAL=15
# PRIOR_TOP_N_EXTREME=10

# CMD2=(python3.9 -m deeprbp.explainability_module.tcga_normal_vs_tumor.cli.run_normal_vs_tumor_rbp_programs
#   --output_dir "$OUTPUT_DIR"
#   --getBM_path "$GETBM_PATH"
# )
# opcional: --rbps_focus "SLU7,SRSF3,SRSF1"

echo "[INFO] Running STEP 2 (Section A - RBP programs):"
printf ' %q' "${CMD2[@]}"; echo
"${CMD2[@]}"

# (en el futuro, STEP 3/3 – tumor_only)
# CMD3=(python3.9 -m deeprbp.explainability_module.tcga_normal_vs_tumor.cli.run_tumor_only_rbp_targets
#   --output_dir "$OUTPUT_DIR"
#   --getBM_path "$GETBM_PATH"
#   --rbps_focus "SRSF1,SRSF3,SLU7"
#   --min_abs_score 0.0
#   --min_n_genes 0
#   --top_n 50
# )
# echo "[INFO] Running STEP 3 (Section B - tumor-only):"
# printf ' %q' "${CMD3[@]}"; echo
# "${CMD3[@]}"

echo "########################################"
echo " DeepRBP - TCGA LIHC job finished at: $(date)"
echo "########################################"