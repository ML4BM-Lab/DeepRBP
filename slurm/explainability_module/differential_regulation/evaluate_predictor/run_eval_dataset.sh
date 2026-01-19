#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=deeprbp_eval
#SBATCH --gres=gpu:1
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/%u/DeepRBP/output/logs/%x.out

set -euo pipefail

DATASET_ID="$1"
PROCESSED_DIR="$2"
OUT_BASE="$3"

if [[ -z "$DATASET_ID" || -z "$PROCESSED_DIR" || -z "$OUT_BASE" ]]; then
  echo "Usage: $0 <DATASET_ID> <processed_dir> <output_base>"
  exit 1
fi

module purge
module load Miniforge3
conda activate DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

echo "########################################"
echo "DeepRBP evaluation — ${DATASET_ID}"
echo "Started at: $(date)"
echo "########################################"

# Detectar si hay subgrupos
SUBDIRS=($(find "${PROCESSED_DIR}" -mindepth 1 -maxdepth 1 -type d))

if [[ ${#SUBDIRS[@]} -eq 0 ]]; then
  echo "▶ Evaluating ALL samples together"

  GROUP_NAME="ALL"
  GROUP_DIR="${PROCESSED_DIR}"
  OUTDIR="${OUT_BASE}/ALL"
  mkdir -p "$OUTDIR"

  CONFIG="/tmp/config_${DATASET_ID}_ALL_eval.yaml"

  cat > "$CONFIG" <<EOF
test_path_files: "${GROUP_DIR}"

getBM_path: "/scratch/jsanchoz/DeepRBP/data/annotation/getBM_gencode_v23.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"

cuda: True
val_batch_size: 256
plot_results: True
seed: 42

scaler_mode: "tcga"
scaler_dir: "/scratch/jsanchoz/DeepRBP/pretrained_model"
EOF

  python -m deeprbp.training_module.evaluate_predictor \
    --config_path "$CONFIG" \
    --model_checkpoint "/scratch/jsanchoz/DeepRBP/pretrained_model/model.ckpt" \
    --output_dir "$OUTDIR" \
    --num_workers 4 \
    --verbose 1

else
  for GROUP_DIR in "${SUBDIRS[@]}"; do
    GROUP_NAME=$(basename "$GROUP_DIR")
    echo "▶ Evaluating group: ${GROUP_NAME}"

    OUTDIR="${OUT_BASE}/${GROUP_NAME}"
    mkdir -p "$OUTDIR"
    CONFIG="/tmp/config_${DATASET_ID}_${GROUP_NAME}_eval.yaml"

    cat > "$CONFIG" <<EOF
test_path_files: "${GROUP_DIR}"

getBM_path: "/scratch/jsanchoz/DeepRBP/data/annotation/getBM_gencode_v23.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"

cuda: True
val_batch_size: 256
plot_results: True
seed: 42

scaler_mode: "tcga"
scaler_dir: "/scratch/jsanchoz/DeepRBP/pretrained_model"
EOF

    python -m deeprbp.training_module.evaluate_predictor \
      --config_path "$CONFIG" \
      --model_checkpoint "/scratch/jsanchoz/DeepRBP/pretrained_model/model.ckpt" \
      --output_dir "$OUTDIR" \
      --num_workers 4 \
      --verbose 1

    echo "✔ Finished ${GROUP_NAME}"
  done
fi

echo "✅ Evaluation completed for ${DATASET_ID}"
