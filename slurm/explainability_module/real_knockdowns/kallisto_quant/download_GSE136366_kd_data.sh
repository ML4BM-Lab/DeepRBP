#!/bin/bash
#SBATCH --job-name=download_GSE136366_kd_data
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/download_GSE136366_kd_data.out   
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
OUTPUT_DIR="${RAW_DIR}/GSE136366"
FASTQ_DIR="${OUTPUT_DIR}/raw"

# Define quantification parameters
THREADS=2
BOOTSTRAPS=100

# Create folders
mkdir -p "${RAW_DIR}" "${OUTPUT_DIR}" "${FASTQ_DIR}"
echo "📁 Created directories: ${RAW_DIR}, ${OUTPUT_DIR}, ${FASTQ_DIR}"

# Download fastq files (only if not already downloaded)
echo "Checking and downloading FASTQ files if needed..."

declare -A FILE_URLS=(
  ["SRR10045016_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/016/SRR10045016/SRR10045016_1.fastq.gz"
  ["SRR10045016_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/016/SRR10045016/SRR10045016_2.fastq.gz"
  ["SRR10045017_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/017/SRR10045017/SRR10045017_1.fastq.gz"
  ["SRR10045017_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/017/SRR10045017/SRR10045017_2.fastq.gz"
  ["SRR10045018_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/018/SRR10045018/SRR10045018_1.fastq.gz"
  ["SRR10045018_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/018/SRR10045018/SRR10045018_2.fastq.gz"
  ["SRR10045019_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/019/SRR10045019/SRR10045019_1.fastq.gz"
  ["SRR10045019_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/019/SRR10045019/SRR10045019_2.fastq.gz"
  ["SRR10045020_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/020/SRR10045020/SRR10045020_1.fastq.gz"
  ["SRR10045020_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/020/SRR10045020/SRR10045020_2.fastq.gz"
  ["SRR10045021_1"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/021/SRR10045021/SRR10045021_1.fastq.gz"
  ["SRR10045021_2"]="ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/SRR100/021/SRR10045021/SRR10045021_2.fastq.gz"
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

echo "📦 All FASTQ files are present."

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
for SAMPLE in SRR10045016 SRR10045017 SRR10045018 SRR10045019 SRR10045020 SRR10045021; do
    SAMPLE_OUTPUT="${OUTPUT_DIR}/kallisto_output/${SAMPLE}"
    mkdir -p "${SAMPLE_OUTPUT}"

    echo "🚀 Running Kallisto for ${SAMPLE}..."
    kallisto quant -i "${KALLISTO_INDEX}" -o "${SAMPLE_OUTPUT}" -b ${BOOTSTRAPS} -t ${THREADS} \
    "${FASTQ_DIR}/${SAMPLE}_1.fastq.gz" "${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
    echo "✅  Kallisto completed for ${SAMPLE}."
done

echo "🎉 All processing completed successfully."

