#!/bin/bash
SCRIPT_DIR=$(dirname "$0")

python "$SCRIPT_DIR/create_data.py" --chunksize 1500 --select_genes 'cancer_genes' --path_data "$SCRIPT_DIR/data/TCGA_GTeX"



