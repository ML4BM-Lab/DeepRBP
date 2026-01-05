#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=xlong
#SBATCH --job-name=run_hyper_optuna_job2
#SBATCH --gres=gpu:1
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_hyper_optuna_job2.out
#SBATCH --mail-type=ALL       
#SBATCH --mail-user=jsanchoz@unav.es

# Timestamp para logs únicos
TIMESTAMP=$(date "+%Y%m%d_%H%M%S")

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "Log timestamp: $TIMESTAMP"
echo "########################################"

module load Miniforge3
source activate /data/jsanchoz/conda-env/DeepRBP
 
# Check Python version
python --version

# Check active conda environment
conda info --envs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

echo "🚀 Lanzando proceso en GPU 0"

LOG_FILE="/scratch/jsanchoz/DeepRBP/output/logs/run_hyper_optuna_gpu0_${TIMESTAMP}_definitive_edition.out"

CUDA_VISIBLE_DEVICES=0 \
python3.9 -m deeprbp.training_module.hyperparameter_optimization.grid_search_optuna \
    --storage_path "/scratch/jsanchoz/DeepRBP/output/results/hyperparameter_optimization_SLURM/optuna.db" \
    --n_trials 200 \
    --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml" \
    --output_dir "/scratch/jsanchoz/DeepRBP/output/results/hyperparameter_optimization_SLURM_gpu0_${TIMESTAMP}_definitive_edition" \
    --num_workers 0 \
    --min_delta 0.001 \
    --patience 30 \
    --gpu_id 0 > "$LOG_FILE" 2>&1

wait   
echo "✅ All process are finished"

