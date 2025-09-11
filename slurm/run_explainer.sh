#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_explainer🚀
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_explainer.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module purge
module load Miniforge3

source activate /data/jsanchoz/conda-env/DeepRBP
#PYTHON_EXEC="/data/jsanchoz/conda-env/DeepRBP/bin/python"

# Check Python version
python --version

# Check active conda environment
conda info --envs

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

python3.9 -m deeprbp.explainability_module.main_explainer \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_explainer_dl_kout_t_stat.yaml" \
  --model_ckpt_path "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/checkpoint_model/deeprbp-predictor-epoch=124-validation_loss=0.08.ckpt" \
  --scaler_dir "/scratch/jsanchoz/DeepRBP/final_results/run_deeprbp_predictor/data" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/explainer_dl_kout_t_stat" \
  --select_category "Liver_Hepatocellular_Carcinoma,Acute_Myeloid_Leukemia,Kidney_Chromophobe"
