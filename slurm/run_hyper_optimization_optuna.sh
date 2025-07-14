#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=🔬run_hyper_optuna
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --constraint=a100-sxm4 
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_hyper_optuna.out
#SBATCH --mail-type=ALL       
#SBATCH --mail-user=jsanchoz@unav.es

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module load Miniforge3
conda activate /data/jsanchoz/conda-env/DeepRBP
#PYTHON_EXEC="/data/jsanchoz/conda-env/DeepRBP/bin/python" (esto está haciendo que falle)

# Check Python version
python --version

# Check active conda environment
conda info --envs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

# Lanzar un proceso por GPU
for i in {0..3}; do
    echo "➤ Lanzando proceso en GPU $i"
    LOG_FILE="/scratch/jsanchoz/DeepRBP/output/logs/run_hyper_optuna_gpu_${i}.out"
    CUDA_VISIBLE_DEVICES=$i \
    python -m deeprbp.training_module.hyperparameter_optimization.grid_search_optuna \
        --storage_path "/scratch/jsanchoz/DeepRBP/output/results/hyperparameter_optimization_SLURM/optuna.db" \
        --n_trials 250 \
        --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml" \
        --output_dir "/scratch/jsanchoz/DeepRBP/output/results/hyperparameter_optimization_SLURM" \
        --num_workers 0 \
        --min_delta 0.001 \
        --patience 30 \
        --gpu_id $i > "$LOG_FILE" 2>&1 &   
done

wait   
echo "✅ All process are finished"

