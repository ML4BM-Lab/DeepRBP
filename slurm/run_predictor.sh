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

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.training_module.main_predictor \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_model_train.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor" \
  --epochs 100 \
  --num_workers 4 \
  --min_delta 0.001 \
  --patience 3


  