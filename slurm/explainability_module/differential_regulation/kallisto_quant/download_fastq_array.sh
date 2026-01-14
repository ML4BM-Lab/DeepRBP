#!/bin/bash
#SBATCH --job-name=dl_fastq_array
#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-300%20
#SBATCH -o /scratch/jsanchoz/DeepRBP/output/logs/dl_fastq_array_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jsanchoz@unav.es

set -euo pipefail

# -------------------------
# Environment
# -------------------------
source /scicomp/builds/Rocky/8.7/Common/software/Miniforge3/24.11.3-2/etc/profile.d/conda.sh
conda activate kallisto_env

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch $0 <DATASET_ID>"
  exit 1
fi

DATASET="$1"

# -------------------------
# Paths
# -------------------------
PROJECT_DATA="/scratch/jsanchoz/DeepRBP/data"
RAW_DIR="${PROJECT_DATA}/explainability_module/differential_regulation"
DATASET_DIR="${RAW_DIR}/${DATASET}"
FASTQ_DIR="${DATASET_DIR}/raw_fastq"
SRR_LIST_FILE="${DATASET_DIR}/SRR_Acc_List.txt"

mkdir -p "${FASTQ_DIR}"

# -------------------------
# Resolve SRR for this task
# -------------------------
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SRR_LIST_FILE}")

if [[ -z "${SRR}" ]]; then
  echo "❌ No SRR found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "📥 Processing ${SRR}"

# Skip if already done
if [[ -f "${FASTQ_DIR}/${SRR}_1.fastq.gz" && -f "${FASTQ_DIR}/${SRR}_2.fastq.gz" ]]; then
  echo "✔️ ${SRR} already downloaded"
  exit 0
fi

# -------------------------
# Local scratch
# -------------------------
LOCAL_SCRATCH="/lscratch/${SLURM_JOB_ID}_${SRR}"
mkdir -p "${LOCAL_SCRATCH}"

# -------------------------
# Download + convert
# -------------------------
echo "⬇️ Prefetching ${SRR} to local scratch"
prefetch "${SRR}" --output-directory "${LOCAL_SCRATCH}"

echo "🔄 Converting ${SRR} to FASTQ on local disk"
fasterq-dump "${LOCAL_SCRATCH}/${SRR}/${SRR}.sra" \
  --split-files \
  --threads 2 \
  --outdir "${LOCAL_SCRATCH}"

gzip -f "${LOCAL_SCRATCH}/${SRR}_1.fastq"
gzip -f "${LOCAL_SCRATCH}/${SRR}_2.fastq"

mv "${LOCAL_SCRATCH}/${SRR}_1.fastq.gz" "${FASTQ_DIR}/"
mv "${LOCAL_SCRATCH}/${SRR}_2.fastq.gz" "${FASTQ_DIR}/"

# -------------------------
# Cleanup
# -------------------------
rm -rf "${LOCAL_SCRATCH}"

echo "✅ ${SRR} download completed"