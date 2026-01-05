#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=run_corum_complex_analysis_1aml
#SBATCH --ntasks-per-node=4
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/jsanchoz/DeepRBP/final_results/explainability/explainer_dl_kout_t_stat/output_final_layer/Acute_Myeloid_Leukemia/corum_complex_analysis/logs/run_corum_complex_analysis_aml.out
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

python3.9 -m deeprbp.explainability_module.complex_analysis.run_corum_complex_analysis \
    --corum_file_path "/scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/complex_analysis/corum_results.txt" \
    --getBM_file_path "/scratch/jsanchoz/DeepRBP/data/annotation/getBM_gencode_v23.csv" \
    --scores_file_path "/scratch/jsanchoz/DeepRBP/final_results/explainability/explainer_dl_kout_t_stat/output_final_layer/Acute_Myeloid_Leukemia/df_scores_GxRBP.csv" \
    --output_dir "/scratch/jsanchoz/DeepRBP/final_results/explainability/explainer_dl_kout_t_stat/output_final_layer/Acute_Myeloid_Leukemia/corum_complex_analysis" \
    --include-non-family
    
echo "✅ Process is finished"

 
