#!/bin/bash

#SBATCH --qos=regular
#SBATCH --job-name=real_kds_lines_SRR296_RBM47_kd    
#SBATCH --cpus-per-task=2
#SBATCH --mem=35GB
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/ML4BM-Lab/DeepRBP/logs/execute_explain_real_kds_lines_SRR296_RBM47_kd.log

echo ===================================
echo ===     Load the packages       ===
echo ===================================
echo $(date)

module load Python
conda activate /scratch/jsanchoz/envs/PyTorch
module unload Python
which python

python /scratch/jsanchoz/ML4BM-Lab/DeepRBP/main_explain_real_kds.py \
    --readme_txt 'sample_info' \
    --path_data '/scratch/jsanchoz/ML4BM-Lab2/DeepRBP/data/validation_real_kds/GSE75491' \
    --experiment 'SRR296' \
    --rbp_interest 'RBM47' \
    --path_model '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/data/trained_models/model_1024N_2HL_8f' \
    --path_result '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/results'


