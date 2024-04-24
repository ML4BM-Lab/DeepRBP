#!/bin/bash

SCRIPT_DIR=$(dirname "$0")


path_deepRBP="/Users/joseba/Downloads/ML4BM-Lab2/DeepRBP"

export PYTHONPATH="$path_deepRBP/model:$PYTHONPATH"
PATH_DATA="$path_deepRBP/data/input_create_model/processed"

echo "PATH_DATA: $PATH_DATA"

python "$SCRIPT_DIR/main_predictor.py" \
    --sel_optimizer 'adamW' \
    --source_train 'TCGA' \
    --max_node 1024 \
    --num_hidden_layers 2 \
    --learning_rate 0.001 \
    --activation_layer 'relu' \
    --node_divition_factor 8 \
    --same_num_nodes True \
    --batch_size 128 \
    --epochs 1000 \
    --train_tumor_types 'all' \
    --test_tumor_types 'all' \
    --path_data "$PATH_DATA" \
    --plot_results True
