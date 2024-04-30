#!/bin/bash
SCRIPT_DIR=$(dirname "$0")

path_deepRBP="" # Change this path to your DeepRBP folder

export PYTHONPATH="$path_deepRBP/model:$PYTHONPATH"
PATH_DATA="$path_deepRBP/data/input_create_model"

echo "Ruta de PATH_DATA: $PATH_DATA"
echo $PYTHONPATH
python "$SCRIPT_DIR/create_data.py" --chunksize 5000 --select_genes 'cancer_genes' --path_data "$PATH_DATA"
