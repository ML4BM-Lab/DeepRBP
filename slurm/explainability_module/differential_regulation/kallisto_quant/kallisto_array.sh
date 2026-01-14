#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-300%30
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/kallisto_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

# -------------------------
# Environment
# -------------------------
source /scicomp/builds/Rocky/8.7/Common/software/Miniforge3/24.11.3-2/etc/profile.d/conda.sh
conda activate kallisto_env

echo "🔎 Kallisto version:"
kallisto version

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch $0 <DATASET_ID>"
  exit 1
fi

DATASET="$1"

# -------------------------
# Paths
# -------------------------
PROJECT_DATA="/scratch/jsanchoz/DeepRBP/data"

ANNOTATION_DIR="${PROJECT_DATA}/annotation"
TRANSCRIPTOME="${ANNOTATION_DIR}/gencode.v23.transcripts.fa.gz"
KALLISTO_INDEX="${ANNOTATION_DIR}/gencode.v23.transcripts.idx"

RAW_DIR="${PROJECT_DATA}/explainability_module/differential_regulation"
DATASET_DIR="${RAW_DIR}/${DATASET}"
FASTQ_DIR="${DATASET_DIR}/raw_fastq"
KALLISTO_OUT="${DATASET_DIR}/kallisto_output"
SRR_LIST_FILE="${DATASET_DIR}/SRR_Acc_List.txt"

THREADS=2
BOOTSTRAPS=100

mkdir -p "${ANNOTATION_DIR}" "${KALLISTO_OUT}"

# -------------------------
# Transcriptome & index (SAFE IN ARRAY)
# -------------------------
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
  if [ ! -f "${TRANSCRIPTOME}" ]; then
    echo "⬇️ Downloading GENCODE v23 transcriptome"
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.transcripts.fa.gz \
         -O "${TRANSCRIPTOME}"
  fi

  if [ ! -f "${KALLISTO_INDEX}" ]; then
    echo "🔧 Building kallisto index"
    kallisto index -i "${KALLISTO_INDEX}" "${TRANSCRIPTOME}"
  fi
fi

# Wait until index exists (for other tasks)
while [ ! -f "${KALLISTO_INDEX}" ]; do
  echo "⏳ Waiting for kallisto index to be ready..."
  sleep 30
done

# -------------------------
# Resolve SRR
# -------------------------
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SRR_LIST_FILE}")

if [[ -z "${SRR}" ]]; then
  echo "❌ No SRR found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "🚀 Running kallisto for ${SRR}"

OUTDIR="${KALLISTO_OUT}/${SRR}"
mkdir -p "${OUTDIR}"

kallisto quant \
  -i "${KALLISTO_INDEX}" \
  -o "${OUTDIR}" \
  -b ${BOOTSTRAPS} \
  -t ${THREADS} \
  "${FASTQ_DIR}/${SRR}_1.fastq.gz" \
  "${FASTQ_DIR}/${SRR}_2.fastq.gz"

echo "✅ Kallisto completed for ${SRR}"