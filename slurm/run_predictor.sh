#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=🚀run_predictor
#SBATCH --gres=gpu:4 # prueba con 4
#SBATCH --ntasks-per-node=4 # prueba con 4 This needs to match Trainer(devices=...), must be number of gpus
#SBATCH --constraint=a100-sxm4 # --constraint=rtx3090
#SBATCH --mem=90gb
#SBATCH --nodes=1 # This needs to match Trainer(num_nodes=...)
#SBATCH --cpus-per-task=4 # total cpus = cpus-per-task*ntasks-per-node
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_predictor.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsanchoz@unav.es

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module load Python
source activate /data/jsanchoz/conda-env/DeepRBP
PYTHON_EXEC="/data/jsanchoz/conda-env/DeepRBP/bin/python"

# Check Python version
python --version

# Check active conda environment
conda info --envs

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
    -m deeprbp.training_module.main_predictor \
    --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_model_train.yaml" \
    --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor_SLURM_try9_06" \
    --epochs 100 \
    --num_workers 4 \
    --min_delta 0.001 \
    --patience 3