#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=run_explainer_postar
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=90gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/run_explainer_postar.out
#SBATCH --mail-type=START,END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:$PYTHONPATH"
python -c "import deeprbp; print('Package found')"

python -m deeprbp.explainability_module.main_explainer \
  --config_path_explain "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml" \
  --config_path_train "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml"

# #### Results Visualization (EXECUTE THIS NOW!!)
# # Load R and run the visualization script
# module load R/4.3.2
# Rscript /scratch/jsanchoz/DeepRBP/src/deeprbp/explainability_module/results_visualization/run_postar_plot_generation.R \
#   --input_path /scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_all_2025-01-22_100_new_good_trained_model/explain_prediction_model/results/DeepLIFT_knockdown_reference_t-statistic_max_absolute_value \
#   --output_path /scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_all_2025-01-22_100_new_good_trained_model/explain_prediction_model/results/DeepLIFT_knockdown_reference_t-statistic_max_absolute_value/results_visualization \
#   --output_filename plot_score_results.pdf \
#   --results_filename df_results_summary.csv \
#   --list_rbps_postar_filename list_rbps_postar_ordered.csv \
#   --list_genes_postar_filename list_genes_postar_ordered.csv \
#   --getBM_filename getBM.csv \
#   --save_plot TRUE \
#   --index_start 1 \
#   --index_end 4 \
#   --max_iterations 7