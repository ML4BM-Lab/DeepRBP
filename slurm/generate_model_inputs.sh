#!/bin/bash
#SBATCH --job-name=gen_model_inputs
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/generate_model_inputs.out
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

module load Python
conda activate /data/jsanchoz/conda-env/DeepRBP

python /scratch/jsanchoz/DeepRBP/src/deeprbp/data_preprocessing/prep_model_inputs.py \
  --raw_data_dir "/scratch/jsanchoz/DeepRBP/data/training_module/raw" \
  --selected_genes_dir "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps" \
  --output_dir "/scratch/jsanchoz/DeepRBP/data/training_module/processed" \
  --transcript_expression_file "TcgaTargetGtex_rsem_isoform_tpm.gz" \
  --gene_expression_file "TcgaTargetGtex_rsem_gene_tpm.gz" \
  --phenotype_data_file "TcgaTargetGTEX_phenotype.txt" \
  --chunk_size 4500 \
  --gene_selection True \
  --gene_transcript_mapping_file "getBM.csv" \
  --splicing_genes_file "Table_S5_Cancer_splicing_gene_eyras.xlsx" \
  --cancer_genes_file "Table_S6_Cancer_gene_eyras.xlsx" \
  --gene_census_file "Table_Cancer_Gene_Census.tsv" \
  --rbp_genes_file "Table_S2_list_RBPs_eyras.xlsx"
