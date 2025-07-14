#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=xlong
#SBATCH --job-name=run_bench_ridge
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_bench_ridge.out
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Miniforge3
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

python -u -m deeprbp.training_module.benchmark_methods.train_baseline \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_benchmark_methods.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_model_benchmark" \
  --algorithm ridge

