#!/bin/bash
set -euo pipefail

# Get the directory where this script lives
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Assume repo root is one level above slurm/
REPO_ROOT="$( dirname "$SCRIPT_DIR" )"

OUTPUT_DIR="$REPO_ROOT/data/training_module/raw"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"
echo "Downloading data to: $OUTPUT_DIR"

# Function to download a file only if it does not already exist
download_if_missing () {
    local url="$1"
    local outfile="$2"

    if [[ -f "$outfile" ]]; then
        echo "✔ File already exists, skipping: $(basename "$outfile")"
    else
        echo "⬇ Downloading: $(basename "$outfile")"
        curl -L "$url" -o "$outfile"
    fi
}

# =========================
# Gene expression (RSEM)
# =========================
download_if_missing \
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_isoform_tpm.gz" \
  "$OUTPUT_DIR/TcgaTargetGtex_rsem_isoform_tpm.gz"

download_if_missing \
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz" \
  "$OUTPUT_DIR/TcgaTargetGtex_rsem_gene_tpm.gz"

download_if_missing \
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz" \
  "$OUTPUT_DIR/TcgaTargetGTEX_gene_expected_count.gz"

# =========================
# Transcript expression RNAseq
# RSEM expected_count (n=19,109)
# =========================
download_if_missing \
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_expected_count.gz" \
  "$OUTPUT_DIR/TcgaTargetGtex_isoform_expected_count.gz"

# =========================
# Phenotype / metadata
# =========================
download_if_missing \
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz" \
  "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"

# Decompress phenotype metadata if needed
if [[ -f "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz" ]] && \
   [[ ! -f "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt" ]]; then
    echo "📦 Decompressing phenotype metadata"
    gzip -d "$OUTPUT_DIR/TcgaTargetGTEX_phenotype.txt.gz"
else
    echo "✔ Phenotype metadata already decompressed or missing"
fi

echo "✅ All requested files are present."