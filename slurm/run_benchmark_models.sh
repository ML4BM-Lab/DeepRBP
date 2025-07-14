#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=run_bench_models
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_bench_models.out
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Miniforge3
conda activate /data/jsanchoz/conda-env/DeepRBP
#PYTHON_EXEC="/data/jsanchoz/conda-env/DeepRBP/bin/python"

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

python -u -m deeprbp.training_module.benchmark_methods.train_baselines \
  --config_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_benchmark_methods.yaml" \
  --output_dir "/scratch/jsanchoz/DeepRBP/output/results/run_model_benchmark"

