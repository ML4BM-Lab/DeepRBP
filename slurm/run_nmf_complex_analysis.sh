#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_nmf_complex_analysis
#SBATCH --ntasks-per-node=4
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_nmf_complex_analysis.out
#SBATCH --mail-type=ALL       
#SBATCH --mail-user=jsanchoz@unav.es

echo "########################################"
echo "Starting job at: $(date)"
echo "Current time in Hondarribia: $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module load Miniforge3
conda activate /data/jsanchoz/conda-env/DeepRBP

# Check Python version
python --version

# Check active conda environment
conda info --envs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
export PYTHONUNBUFFERED=1
python -c "import deeprbp; print('Package found')"

python -m deeprbp.explainability_module.complex_analysis.run_detect_complex_nmf_analysis \
    --corum_file_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/complex_analysis/corum_results.txt" \
    --scores_file_path "/scratch/jsanchoz/DeepRBP/output/results/stuff/explainability_deeplift_knock_t_stat/results/df_scores_GxRBP.csv" \
    --output_dir "/scratch/jsanchoz/DeepRBP/output/results/nmf_detect_complex_analysis" \
    

echo "✅ Process is finished"

