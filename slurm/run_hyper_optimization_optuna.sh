#!/bin/bash
#SBATCH --qos=test
#SBATCH --job-name=run_hyper_optuna
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=3gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_hyper_optuna.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.training_module.hyperparameter_optimization.grid_search_optuna \
        --config_path_file '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml' \
        --output_dir '/scratch/jsanchoz/DeepRBP/stuff/' \
        --val_batch_size 32 \
        --n_trials 4


