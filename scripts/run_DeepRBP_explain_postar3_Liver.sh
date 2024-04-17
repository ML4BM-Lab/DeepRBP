#!/bin/bash

#SBATCH --qos=test
#SBATCH --job-name=run_DeepRBP_explain_postar3_Liver
#SBATCH --cpus-per-task=2
#SBATCH --mem=35GB
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/ML4BM-Lab/DeepRBP/logs/run_DeepRBP_explain_postar3_Liver.log

# Set the variable
CELL_LINE="all"

echo ===================================
echo ===     Load the packages       ===
echo ===================================
echo $(date)

module load Python
conda activate /scratch/jsanchoz/envs/PyTorch
module unload Python
which python

# Explainability with Deeplift and differential expression
python /scratch/jsanchoz/ML4BM-Lab/DeepRBP/main_explain_postar3.py \
  --tumor_type 'Liver_Hepatocellular_Carcinoma' \
  --cell_line $CELL_LINE \
  --source_explain 'TCGA' \
  --model_selected '1024N_2HL_8f' \
  --path_result '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/results' \
  --path_data '/scratch/jsanchoz/data'
