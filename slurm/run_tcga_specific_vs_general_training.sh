#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=test
#SBATCH --job-name=run_tcga_specific_vs_general_training
#SBATCH --gres=gpu:1 # prueba con 4
#SBATCH --ntasks-per-node=4 # prueba con 4 This needs to match Trainer(devices=...), must be number of gpus
#SBATCH --constraint=rtx3090 # --constraint=rtx3090; a100-sxm4
#SBATCH --mem=20gb
#SBATCH --nodes=1 # This needs to match Trainer(num_nodes=...)
#SBATCH --cpus-per-task=4 # total cpus = cpus-per-task*ntasks-per-node
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_tcga_specific_vs_general_training.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsanchoz@unav.es

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module load Miniforge3
conda activate /data/jsanchoz/conda-env/DeepRBP

# Check Python version
python --version

# Check active conda environment
conda info --envs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

# Run the training script with torchrun
torchrun \
    --nproc_per_node=$SLURM_NTASKS \
    --nnodes=$SLURM_JOB_NUM_NODES \
    --node_rank=$SLURM_NODEID \
    --master_addr=$(hostname) \
    --master_port=$(shuf -i 20000-30000 -n 1) \
    -m deeprbp.training_module.tumor_specific_training.main_specific_vs_general \
    --output_base_dir "/scratch/jsanchoz/DeepRBP/output/results/tumor_specific_training" \
    --epochs 100 \
    --num_workers 4 \
    --min_delta 0.001 \
    --patience 3 \
    --save_top_k 1 \
    --verbose 1 \
    --train_batch_size 8 \
    --val_batch_size 8
