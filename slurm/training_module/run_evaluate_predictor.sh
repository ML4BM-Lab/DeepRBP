#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=deeprbp_eval
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4   # must match Trainer(devices=...) in the evaluator
#SBATCH --constraint=a100-sxm4
#SBATCH --mem=60gb
#SBATCH --nodes=1             # must match Trainer(num_nodes=...)
#SBATCH --cpus-per-task=2     # total cpus = cpus-per-task * ntasks-per-node
#SBATCH -o /scratch/$USER/DeepRBP/output/logs/deeprbp_eval.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@domain.com

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module purge
module load Miniforge3
source activate /data/jsanchoz/conda-env/DeepRBP

python --version
conda info --envs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONPATH="/scratch/$USER/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

# Run evaluator with torchrun (DDP-style, consistent with training scripts)
torchrun \
  --nproc_per_node=$SLURM_NTASKS \
  --nnodes=$SLURM_JOB_NUM_NODES \
  --node_rank=$SLURM_NODEID \
  --master_addr=$(hostname) \
  --master_port=$(shuf -i 20000-30000 -n 1) \
  -m deeprbp.training_module.evaluate_predictor \
  --config_path "/scratch/$USER/DeepRBP/src/deeprbp/configs/config_model_eval.yaml" \
  --model_checkpoint "/scratch/$USER/DeepRBP/pretrained_model/model.ckpt" \
  --output_dir "/scratch/$USER/DeepRBP/output/results/eval_pretrained" \
  --num_workers $SLURM_CPUS_PER_TASK \
  --verbose 1

echo "✅ All processes are finished"



