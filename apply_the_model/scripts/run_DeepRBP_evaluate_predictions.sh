#!/bin/bash

# Script directory
SCRIPT_DIR=$(dirname "$0")

# DeepRBP paths and data
path_deepRBP="/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP"
export PYTHONPATH="$path_deepRBP/model:$PYTHONPATH"
PATH_DATA="$path_deepRBP/data/input_create_model/processed"

# Trained files paths and output folder
path_model="$path_deepRBP/model"
PATH_TRAIN_FILES="$path_model/output/example"
PATH_SAVE="$SCRIPT_DIR/../output/transcript_expression"

echo "PATH_DATA: $PATH_DATA"
echo "PATH_TRAIN_FILES: $PATH_TRAIN_FILES"
echo "PATH_SAVE: $PATH_SAVE"

# Execute python script
python "$SCRIPT_DIR/../main_evaluate_predictions.py" \
    --path_data "$PATH_DATA" \
    --path_train_files "$PATH_TRAIN_FILES" \
    --path_save "$PATH_SAVE" \
    --plot_results True
