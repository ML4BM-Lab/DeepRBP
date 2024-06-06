#!/bin/bash

# Script directory
SCRIPT_DIR="$PWD"

# DeepRBP paths and data
path_deepRBP="/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP"
export PYTHONPATH="$path_deepRBP/model:$PYTHONPATH"

# Paths and data
path_model="$path_deepRBP/model"
path_data="$path_deepRBP/data/input_create_model"
path_train_files="$path_model/output/example/"
path_exp="$path_deepRBP/data/data_real_kds/experiments/GSE136366"
experiment=$(basename "$path_exp")
rbp_interest="TARDBP"
path_save="../output/explain_real_kds"

# Arguments for Python script
echo "PATH_MODEL: $path_model"
echo "PATH_DATA: $path_data"
echo "PATH_TRAIN_FILES: $path_train_files"
echo "PATH_EXP: $path_exp"
echo "EXPERIMENT: $experiment"
echo "RBP_INTEREST: $rbp_interest"
echo "PATH_SAVE: $path_save"

echo "Executing Python script..."

# Execute python script
python "$SCRIPT_DIR/../main_explain_real_kds.py" \
    --path_exp "$path_exp" \
    --experiment "$experiment" \
    --rbp_interest "$rbp_interest" \
    --path_train_files "$path_train_files" \
    --path_data "$path_data" \
    --path_save "$path_save"