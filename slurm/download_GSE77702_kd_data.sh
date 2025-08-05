#!/bin/bash
#SBATCH --job-name=download_GSE77702_kd_data
#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/download_GSE77702_kd_data.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -e  # Exit immediately on error

# ⚙️ Activar entorno con Kallisto
source /scicomp/builds/Rocky/8.7/Common/software/Miniforge3/24.11.3-2/etc/profile.d/conda.sh
conda activate kallisto_env

# 📁 Configuración de paths
RAW_DIR="/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds"
OUTPUT_DIR="${RAW_DIR}/GSE77702"
FASTQ_DIR="${OUTPUT_DIR}/raw"
TRANSCRIPTOME="${RAW_DIR}/gencode.v23.transcripts.fa.gz"
KALLISTO_INDEX="${RAW_DIR}/gencode.v23.transcripts.idx"
THREADS=2
BOOTSTRAPS=100

# 📂 Crear directorios si no existen
mkdir -p "${RAW_DIR}" "${OUTPUT_DIR}" "${FASTQ_DIR}"
echo "📁 Created directories: ${RAW_DIR}, ${OUTPUT_DIR}, ${FASTQ_DIR}"

# 📥 Lista de muestras a descargar (FASTQ single-end)
declare -A FILE_URLS=(
  ["SRR3153251"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/001/SRR3153251/SRR3153251.fastq.gz"
  ["SRR3153252"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/002/SRR3153252/SRR3153252.fastq.gz"
  ["SRR3153253"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/003/SRR3153253/SRR3153253.fastq.gz"
  ["SRR3153254"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/004/SRR3153254/SRR3153254.fastq.gz"
  ["SRR3153255"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/005/SRR3153255/SRR3153255.fastq.gz"
  ["SRR3153256"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/006/SRR3153256/SRR3153256.fastq.gz"
  ["SRR3153257"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/007/SRR3153257/SRR3153257.fastq.gz"
  ["SRR3153258"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/008/SRR3153258/SRR3153258.fastq.gz"
  ["SRR3153259"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/009/SRR3153259/SRR3153259.fastq.gz"
  ["SRR3153260"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/000/SRR3153260/SRR3153260.fastq.gz"
  ["SRR3153261"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/001/SRR3153261/SRR3153261.fastq.gz"
  ["SRR3153262"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/002/SRR3153262/SRR3153262.fastq.gz"
  ["SRR3153263"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/003/SRR3153263/SRR3153263.fastq.gz"
  ["SRR3153264"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/004/SRR3153264/SRR3153264.fastq.gz"
  ["SRR3153265"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/005/SRR3153265/SRR3153265.fastq.gz"
  ["SRR3153266"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/006/SRR3153266/SRR3153266.fastq.gz"
)

# ⬇️ Descarga de FASTQ (si no existen)
for SAMPLE in "${!FILE_URLS[@]}"; do
  FILE_PATH="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
  if [ -f "${FILE_PATH}" ]; then
    echo "✔️  ${SAMPLE}.fastq.gz already exists, skipping."
  else
    echo "⬇️  Downloading ${SAMPLE}.fastq.gz..."
    wget -c "${FILE_URLS[$SAMPLE]}" -O "${FILE_PATH}"
    echo "✅  Downloaded ${SAMPLE}.fastq.gz."
  fi
done

echo "📦 All FASTQ files are present."

# 📥 Descargar transcriptoma si no está
if [ ! -f "${TRANSCRIPTOME}" ]; then
  echo "⬇️  Downloading GENCODE v23 transcriptome..."
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.transcripts.fa.gz -O "${TRANSCRIPTOME}"
  echo "✅  Transcriptome downloaded at ${TRANSCRIPTOME}"
else
  echo "✔️  Transcriptome already exists at ${TRANSCRIPTOME}, skipping download."
fi

# 🔧 Crear índice de Kallisto si no existe
if [ ! -f "${KALLISTO_INDEX}" ]; then
  echo "🔧 Creating Kallisto index..."
  kallisto index -i "${KALLISTO_INDEX}" "${TRANSCRIPTOME}"
  echo "✅  Kallisto index created: ${KALLISTO_INDEX}"
else
  echo "✔️  Kallisto index already exists at ${KALLISTO_INDEX}, skipping creation."
fi

# 📁 Crear carpeta de salida para resultados de Kallisto
mkdir -p "${OUTPUT_DIR}/kallisto_output"

# 🚀 Correr Kallisto (single-end)
for SAMPLE in "${!FILE_URLS[@]}"; do
  INPUT_FASTQ="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
  SAMPLE_OUTPUT="${OUTPUT_DIR}/kallisto_output/${SAMPLE}"
  mkdir -p "${SAMPLE_OUTPUT}"

  echo "🚀 Running Kallisto for ${SAMPLE}..."
  kallisto quant -i "${KALLISTO_INDEX}" -o "${SAMPLE_OUTPUT}" -b ${BOOTSTRAPS} -t ${THREADS} \
    --single -l 200 -s 20 "${INPUT_FASTQ}"
  echo "✅  Kallisto completed for ${SAMPLE}."
done

echo "🎉 All processing completed successfully for GSE77702!"


# this is the only file with one-end data.
# https://www.ebi.ac.uk/ena/browser/view/PRJNA311234

# Documentación oficial de Kallisto
# 🔗 Fuente: https://pachterlab.github.io/kallisto/about

# For single-end reads, the average fragment length and the standard deviation must be provided using the -l and -s options.