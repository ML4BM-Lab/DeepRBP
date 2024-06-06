#!/bin/bash
#SBATCH --job-name=create_data
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH -o /scratch/jsanchoz/DeepRBP/logs/out1.out

SCRIPT_DIR="$PWD"

path_deepRBP="/scratch/jsanchoz/DeepRBP" # Change this path to your DeepRBP folder
PATH_DATA="$path_deepRBP/data/input_create_model"
echo "Ruta de PATH_DATA: $PATH_DATA"

python "$SCRIPT_DIR/create_data.py" --chunksize 1000 --select_genes 'cancer_genes' --path_data "$PATH_DATA"
