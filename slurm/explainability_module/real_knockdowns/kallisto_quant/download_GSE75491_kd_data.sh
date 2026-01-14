#!/bin/bash
#SBATCH --job-name=download_GSE75491_kd_data
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/download_GSE75491_kd_data.out   
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

set -e  # Exit immediately if a command exits with a non-zero status

source /scicomp/builds/Rocky/8.7/Common/software/Miniforge3/24.11.3-2/etc/profile.d/conda.sh
conda activate kallisto_env

echo "🔎 Kallisto version:"
kallisto version

# -------------------------
# Paths
# -------------------------
PROJECT_DATA="/scratch/jsanchoz/DeepRBP/data"

# Shared annotation
ANNOTATION_DIR="${PROJECT_DATA}/annotation"
TRANSCRIPTOME="${ANNOTATION_DIR}/gencode.v23.transcripts.fa.gz"
KALLISTO_INDEX="${ANNOTATION_DIR}/gencode.v23.transcripts.idx"

# Dataset-specific paths
RAW_DIR="${PROJECT_DATA}/explainability_module/real_kds"
OUTPUT_DIR="${RAW_DIR}/GSE75491"
FASTQ_DIR="${OUTPUT_DIR}/raw"

# Define quantification parameters
THREADS=2
BOOTSTRAPS=100

# Create directories
mkdir -p "${RAW_DIR}" "${OUTPUT_DIR}" "${FASTQ_DIR}"
echo "📁 Created directories: ${RAW_DIR}, ${OUTPUT_DIR}, ${FASTQ_DIR}"

# Download fastq files (only if not already downloaded)
echo "Checking and downloading FASTQ files if needed..."

declare -A FILE_URLS=(
  ["SRR2966446_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/006/SRR2966446/SRR2966446_1.fastq.gz"
  ["SRR2966446_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/006/SRR2966446/SRR2966446_2.fastq.gz"
  ["SRR2966447_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/007/SRR2966447/SRR2966447_1.fastq.gz"
  ["SRR2966447_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/007/SRR2966447/SRR2966447_2.fastq.gz"
  ["SRR2966448_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/008/SRR2966448/SRR2966448_1.fastq.gz"
  ["SRR2966448_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/008/SRR2966448/SRR2966448_2.fastq.gz"
  ["SRR2966449_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/009/SRR2966449/SRR2966449_1.fastq.gz"
  ["SRR2966449_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/009/SRR2966449/SRR2966449_2.fastq.gz"
  ["SRR2966450_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/000/SRR2966450/SRR2966450_1.fastq.gz"
  ["SRR2966450_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/000/SRR2966450/SRR2966450_2.fastq.gz"
  ["SRR2966451_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/001/SRR2966451/SRR2966451_1.fastq.gz"
  ["SRR2966451_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/001/SRR2966451/SRR2966451_2.fastq.gz"
)

for FILE in "${!FILE_URLS[@]}"; do
  FILE_PATH="${FASTQ_DIR}/${FILE}.fastq.gz"
  if [ -f "${FILE_PATH}" ]; then
    echo "✔️  ${FILE}.fastq.gz already exists, skipping."
  else
    echo "⬇️  Downloading ${FILE}.fastq.gz..."
    wget -c "${FILE_URLS[$FILE]}" -O "${FILE_PATH}"
    echo "✅  Downloaded ${FILE}.fastq.gz."
  fi
done

echo "All FASTQ files are present."

# Download GENCODE transcriptome if not already present
if [ ! -f "${TRANSCRIPTOME}" ]; then
  echo "⬇️  Downloading GENCODE v23 transcriptome..."
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.transcripts.fa.gz -O "${TRANSCRIPTOME}"
  echo "✅  Transcriptome downloaded at ${TRANSCRIPTOME}"
else
  echo "✔️  Transcriptome already exists at ${TRANSCRIPTOME}, skipping download."
fi

# Create Kallisto index if not already done
if [ ! -f "${KALLISTO_INDEX}" ]; then
    echo "🔧 Creating Kallisto index..."
    kallisto index -i "${KALLISTO_INDEX}" "${TRANSCRIPTOME}"
    echo "✅  Kallisto index created: ${KALLISTO_INDEX}"
else
    echo "✔️  Kallisto index already exists at ${KALLISTO_INDEX}, skipping creation."
fi

# Create base output dir for kallisto results
mkdir -p "${OUTPUT_DIR}/kallisto_output"

# Run Kallisto quantification
for SAMPLE in SRR2966446 SRR2966447 SRR2966448 SRR2966449 SRR2966450 SRR2966451; do
    SAMPLE_OUTPUT="${OUTPUT_DIR}/kallisto_output/${SAMPLE}"
    mkdir -p "${SAMPLE_OUTPUT}"

    echo "🚀 Running Kallisto for ${SAMPLE}..."
    kallisto quant -i "${KALLISTO_INDEX}" -o "${SAMPLE_OUTPUT}" -b ${BOOTSTRAPS} -t ${THREADS} \
      "${FASTQ_DIR}/${SAMPLE}_1.fastq.gz" "${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
    echo "✅  Kallisto completed for ${SAMPLE}."
    
done

echo "🎉 All processing completed successfully."
