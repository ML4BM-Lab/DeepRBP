#!/bin/bash
SCRIPT_DIR="$PWD"

path_deepRBP="/scratch/jsanchoz/DeepRBP" # Change this path to your DeepRBP folder
PATH_DATA="$path_deepRBP/data/input_create_model"
echo "Ruta de PATH_DATA: $PATH_DATA"

python "$SCRIPT_DIR/create_data.py" --chunksize 5000 --select_genes 'cancer_genes' --path_data "$PATH_DATA"
