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

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

# Run the training script with torchrun
torchrun \
    --nproc_per_node=$SLURM_NTASKS \
    --nnodes=$SLURM_JOB_NUM_NODES \
    --node_rank=$SLURM_NODEID \
    --master_addr=$(hostname) \
    --master_port=$(shuf -i 20000-30000 -n 1) \
    -m deeprbp.training_module.hyperparameter_optimization.grid_search_optuna \
        --config_path '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml' \
        --n_trials 33 \
        --output_dir '/scratch/jsanchoz/DeepRBP/output/results/hyperpameter_optimization_SLURM' \
        --num_workers 0 \
        --min_delta 0.001 \
        --patience 30