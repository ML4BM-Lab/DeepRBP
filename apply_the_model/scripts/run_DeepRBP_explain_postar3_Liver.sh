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
PATH_SAVE="$SCRIPT_DIR/../output/explain_postar"

# Arguments for Python script
TUMOR_TYPE="Liver_Hepatocellular_Carcinoma"  # Example value, replace with your desired tumor type
SOURCE_EXPLAIN="TCGA"  # Example value, replace with your source for explainability

echo "PATH_DATA: $PATH_DATA"
echo "PATH_TRAIN_FILES: $PATH_TRAIN_FILES"
echo "PATH_SAVE: $PATH_SAVE"
echo "TUMOR_TYPE: $TUMOR_TYPE"
echo "SOURCE_EXPLAIN: $SOURCE_EXPLAIN"

# Execute python script
python "$SCRIPT_DIR/../main_explain_postar3.py" \
    --path_save "$PATH_SAVE" \
    --path_data "$PATH_DATA" \
    --tumor_type "$TUMOR_TYPE" \
    --source_explain "$SOURCE_EXPLAIN" \
    --path_train_files "$PATH_TRAIN_FILES"