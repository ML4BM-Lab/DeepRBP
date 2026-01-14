#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_realkds
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_realkds.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.explainability_module.real_knockdowns.main_real_knockdowns \
    --config_path '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_real_knockdowns.yaml' \
    --processed_data_dir '/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/GSE136366/processed' \
    --output_dir '/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns/GSE136366'