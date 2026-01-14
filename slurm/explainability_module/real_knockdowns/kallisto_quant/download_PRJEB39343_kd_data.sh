#!/bin/bash
#SBATCH --job-name=download_PRJEB39343_kd_data
#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/download_PRJEB39343_kd_data.out   
#SBATCH --mail-type=END,FAIL        
#SBATCH --mail-user=jsanchoz@unav.es

set -e  # Exit on error

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
OUTPUT_DIR="${RAW_DIR}/PRJEB39343"
FASTQ_DIR="${OUTPUT_DIR}/raw"

THREADS=2
BOOTSTRAPS=100

# Create directories
mkdir -p "${RAW_DIR}" "${OUTPUT_DIR}" "${FASTQ_DIR}"
echo "📁 Created directories: ${RAW_DIR}, ${OUTPUT_DIR}, ${FASTQ_DIR}"

# Download FASTQ files
echo "⬇️  Checking and downloading FASTQ files..."

declare -A FILE_URLS=(
  ["ERR4352445_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/005/ERR4352445/ERR4352445_1.fastq.gz"
  ["ERR4352445_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/005/ERR4352445/ERR4352445_2.fastq.gz"
  ["ERR4352446_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/006/ERR4352446/ERR4352446_1.fastq.gz"
  ["ERR4352446_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/006/ERR4352446/ERR4352446_2.fastq.gz"
  ["ERR4352447_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/007/ERR4352447/ERR4352447_1.fastq.gz"
  ["ERR4352447_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/007/ERR4352447/ERR4352447_2.fastq.gz"
  ["ERR4352448_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/008/ERR4352448/ERR4352448_1.fastq.gz"
  ["ERR4352448_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/008/ERR4352448/ERR4352448_2.fastq.gz"
  ["ERR4352449_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/009/ERR4352449/ERR4352449_1.fastq.gz"
  ["ERR4352449_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/009/ERR4352449/ERR4352449_2.fastq.gz"
  ["ERR4352450_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/000/ERR4352450/ERR4352450_1.fastq.gz"
  ["ERR4352450_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/000/ERR4352450/ERR4352450_2.fastq.gz"
  ["ERR4352451_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/001/ERR4352451/ERR4352451_1.fastq.gz"
  ["ERR4352451_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/001/ERR4352451/ERR4352451_2.fastq.gz"
  ["ERR4352452_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/002/ERR4352452/ERR4352452_1.fastq.gz"
  ["ERR4352452_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/002/ERR4352452/ERR4352452_2.fastq.gz"
  ["ERR4352453_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/003/ERR4352453/ERR4352453_1.fastq.gz"
  ["ERR4352453_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/003/ERR4352453/ERR4352453_2.fastq.gz"
  ["ERR4352454_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/004/ERR4352454/ERR4352454_1.fastq.gz"
  ["ERR4352454_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/004/ERR4352454/ERR4352454_2.fastq.gz"
  ["ERR4352455_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/005/ERR4352455/ERR4352455_1.fastq.gz"
  ["ERR4352455_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/005/ERR4352455/ERR4352455_2.fastq.gz"
  ["ERR4352456_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/006/ERR4352456/ERR4352456_1.fastq.gz"
  ["ERR4352456_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/006/ERR4352456/ERR4352456_2.fastq.gz"
  ["ERR4352457_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/007/ERR4352457/ERR4352457_1.fastq.gz"
  ["ERR4352457_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/007/ERR4352457/ERR4352457_2.fastq.gz"
  ["ERR4352458_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/008/ERR4352458/ERR4352458_1.fastq.gz"
  ["ERR4352458_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/008/ERR4352458/ERR4352458_2.fastq.gz"
  ["ERR4352459_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/009/ERR4352459/ERR4352459_1.fastq.gz"
  ["ERR4352459_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/009/ERR4352459/ERR4352459_2.fastq.gz"
  ["ERR4352460_1"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/000/ERR4352460/ERR4352460_1.fastq.gz"
  ["ERR4352460_2"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/000/ERR4352460/ERR4352460_2.fastq.gz"
)

# Download each file to FASTQ_DIR
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

# Download GENCODE transcriptome if needed
if [ ! -f "${TRANSCRIPTOME}" ]; then
  echo "⬇️  Downloading GENCODE v23 transcriptome..."
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.transcripts.fa.gz -O "${TRANSCRIPTOME}"
  echo "✅  Transcriptome downloaded."
else
  echo "✔️  Transcriptome already present."
fi

# Create Kallisto index if not exists
if [ ! -f "${KALLISTO_INDEX}" ]; then
  echo "🔧 Creating Kallisto index..."
  kallisto index -i "${KALLISTO_INDEX}" "${TRANSCRIPTOME}"
  echo "✅  Kallisto index created."
else
  echo "✔️  Kallisto index already exists."
fi

mkdir -p "${OUTPUT_DIR}/kallisto_output"

# Run kallisto quantification
for SAMPLE in \
  ERR4352445 ERR4352446 ERR4352447 ERR4352448 ERR4352449 ERR4352450 \
  ERR4352451 ERR4352452 ERR4352453 ERR4352454 ERR4352455 ERR4352456 \
  ERR4352457 ERR4352458 ERR4352459 ERR4352460; do

  SAMPLE_OUTPUT="${OUTPUT_DIR}/kallisto_output/${SAMPLE}"
  mkdir -p "${SAMPLE_OUTPUT}"

  echo "🚀 Running Kallisto for ${SAMPLE}..."
  kallisto quant -i "${KALLISTO_INDEX}" -o "${SAMPLE_OUTPUT}" -b ${BOOTSTRAPS} -t ${THREADS} \
    "${FASTQ_DIR}/${SAMPLE}_1.fastq.gz" "${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
  echo "✅  Kallisto completed for ${SAMPLE}."

done
echo "🎉 All processing completed successfully."
