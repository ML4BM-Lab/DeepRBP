#!/bin/bash

OUTPUT_DIR=./data/training_module/raw

# Crear el directorio si no existe
mkdir -p "$OUTPUT_DIR"

# Descargar los archivos usando la ruta base
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_isoform_tpm.gz -o "$OUTPUT_DIR/TcgaTargetGtex_rsem_isoform_tpm.gz"
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz -o "$OUTPUT_DIR/TcgaTargetGtex_rsem_gene_tpm.gz"
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz -o "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"

# Descomprimir el archivo
gzip -d "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"
