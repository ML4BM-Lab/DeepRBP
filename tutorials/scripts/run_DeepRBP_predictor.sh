#!/bin/bash

#SBATCH --qos=regular
#SBATCH --job-name=run_DeepRBP_predictor
#SBATCH --gres=gpu:rtx3090:1
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jsanchoz@unav.es
#SBATCH -o /scratch/jsanchoz/ML4BM-Lab/DeepRBP/logs/run_DeepRBP_predictor.log

echo ===================================
echo ===     Load the packages       ===
echo ===================================
echo `date`

module load Python
conda activate /scratch/jsanchoz/envs/PyTorch
module unload Python
which python

# Train the model and do predictions on Test data
python /scratch/jsanchoz/ML4BM-Lab/DeepRBP/main_predictor.py \
    --sel_optimizer 'adamW' \
    --source_train 'TCGA' \
    --activation_layer 'relu' \
    --learning_rate 0.001 \
    --max_node 1024 \
    --node_divition_factor 8 \
    --num_hidden_layers 2 \
    --same_num_nodes True \
    --batch_size 128 \
    --epochs 1000 \
    --model_name 'deepSF' \
    --train_tumor_types 'all' \
    --test_tumor_types 'all' \
    --path_result '/scratch/jsanchoz/ML4BM-Lab/DeepRBP/results' \
    --path_data '/scratch/jsanchoz/data' \
    --use_trained_model True \
    --save_model False