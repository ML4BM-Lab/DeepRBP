#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=run_explainer
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_explainer.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.explainability_module.main_explainer \
  --config_path_explain "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_model_explain_deeplift_knock_t_stat.yaml" \
  --config_path_train "/scratch/jsanchoz/DeepRBP/output/results/run_deeprbp_predictor/results/config.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainability_deeplift_knock_t_stat"

