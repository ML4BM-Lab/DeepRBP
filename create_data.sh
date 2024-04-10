#!/bin/bash

#SBATCH --partition=regular
#SBATCH --job-name=cd_cancer_gns
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/ML4BM-Lab/DeepRBP/logs/create_data_cancer_gns.log

echo ===================================
echo ===     Load the packages       ===
echo ===================================
echo `date`

module load Python
conda activate /scratch/jsanchoz/envs/PyTorch
module unload Python
which python

python /scratch/jsanchoz/ML4BM-Lab/DeepRBP/create_data.py --chunksize 8000 --select_genes 'cancer_genes' --path_data '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/data'





