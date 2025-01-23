#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=run_predictor
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=90gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_predictor.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

# python /scratch/jsanchoz/DeepRBP/src/deeprbp/predictor_pipeline.py \
#   --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml" \
#   --external_config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml"

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.training_module.main_predictor \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml" \
  --external_config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_gtex.yaml"
