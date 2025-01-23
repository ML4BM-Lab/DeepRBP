#!/bin/bash

OUTPUT_DIR=./data/training_module/raw

# Create the directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Download the files using the base path
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_isoform_tpm.gz -o "$OUTPUT_DIR/TcgaTargetGtex_rsem_isoform_tpm.gz"
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz -o "$OUTPUT_DIR/TcgaTargetGtex_rsem_gene_tpm.gz"
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz -o "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz -o "$OUTPUT_DIR/TcgaTargetGTEX_gene_expected_count.gz"

# Decompress the metadata file
gzip -d "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"
