#!/bin/bash
#SBATCH --job-name=preprocess_user_data
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/preprocess_user_data_%j.out

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

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python -c "import deeprbp; print('[OK] deeprbp package found')"

echo ""
echo "[preprocess-user-data] Command-line arguments:"
echo "  $@"
echo ""

# -------------------------
# Run preprocessing
# -------------------------
preprocess-user-data "$@"

echo ""
echo "########################################"
echo " Finished at: $(date)"
echo "########################################"
