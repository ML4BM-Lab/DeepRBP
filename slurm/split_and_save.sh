#!/bin/bash
#SBATCH --job-name=split_and_save
#SBATCH --qos=test
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/split_and_save.out
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

python /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/split_data_and_save.py \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_data_split.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/splitted_datasets" \
                