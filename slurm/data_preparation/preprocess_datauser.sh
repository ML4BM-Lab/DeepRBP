#!/bin/bash
#SBATCH --job-name=preprocess_user_data
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

echo "########################################"
echo " DeepRBP user-data preprocessing job"
echo " Job ID     : ${SLURM_JOB_ID:-N/A}"
echo " Started at : $(date)"
echo "########################################"

# -------------------------
# Environment
# -------------------------
module purge
module load Miniforge3
conda activate DeepRBP

python --version
export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python -c "import deeprbp; print('[OK] deeprbp package found')"

echo ""
echo "[preprocess-user-data] Arguments:"
printf ' %q' "$@"; echo
echo ""

# Quick hint about grouping mode
if printf '%s\n' "$@" | grep -q -- '--group_col'; then
  echo "[preprocess-user-data] Mode: group-wise preprocessing (group_col provided)"
else
  echo "[preprocess-user-data] Mode: single-dataset preprocessing (no group_col)"
fi
echo ""

# -------------------------
# Run preprocessing
# -------------------------
preprocess-user-data "$@"

echo ""
echo "########################################"
echo " Finished at: $(date)"
echo "########################################"
