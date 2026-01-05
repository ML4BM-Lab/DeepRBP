#!/bin/bash
#SBATCH --job-name=gen_model_inputs
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/generate_model_inputs.out
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

echo "########################################"
echo " TCGA comparative regulation job started "
echo " Started at: $(date)"
echo " Local time (Hondarribia): $(TZ='Europe/Madrid' date '+%Y-%m-%d %H:%M:%S')"
echo "########################################"

module purge
module load Miniforge3
conda activate /home/jsanchoz/.conda/envs/DeepRBP

python --version
conda info --envs

export PYTHONPATH="/scratch/jsanchoz/DeepRBP/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python -c "import deeprbp; print('[OK] deeprbp package found')"

python3.9 -m deeprbp.data_preprocessing.preprocess_data \
  --raw_data_dir "/scratch/jsanchoz/DeepRBP/data/training_module/raw" \
  --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/processed" \
  --feature_spec_file "/scratch/jsanchoz/DeepRBP/data/training_module/feature_specs/DeepRBP_feature_spec.xlsx" \
  --transcript_expression_file "TcgaTargetGtex_rsem_isoform_tpm.gz" \
  --gene_expression_file "TcgaTargetGtex_rsem_gene_tpm.gz" \
  --gene_counts_file "TcgaTargetGTEX_gene_expected_count.gz" \
  --transcript_counts_file "TcgaTargetGtex_isoform_expected_count.gz" \
  --phenotype_data_file "TcgaTargetGTEX_phenotype.txt" \
  --chunk_size 4500
