#!/bin/bash
set -euo pipefail

source /scicomp/builds/Rocky/8.7/Common/software/Miniforge3/24.11.3-2/etc/profile.d/conda.sh
conda activate DeepRBP

# ============================================
# Kallisto QC: summary + metadata + plots
# ============================================
# Usage:
#   kallisto_qc_full.sh \
#     <kallisto_output_dir> \
#     <SraRunTable.csv> \
#     <output_prefix>
# ============================================

if [[ $# -ne 3 ]]; then
  echo "❌ Usage: $0 <kallisto_output_dir> <SraRunTable.csv> <output_prefix>"
  exit 1
fi

KALLISTO_OUT="$1"
METADATA="$2"
PREFIX_NAME="$3"
DATASET_DIR="$(dirname "${KALLISTO_OUT}")"
PREFIX="${DATASET_DIR}/${PREFIX_NAME}"

if [[ ! -d "${KALLISTO_OUT}" ]]; then
  echo "❌ Directory not found: ${KALLISTO_OUT}"
  exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "❌ Metadata file not found: ${METADATA}"
  exit 1
fi

python3.9 <<EOF
import json
import csv
from pathlib import Path
import statistics
import sys
import matplotlib.pyplot as plt

base = Path("${KALLISTO_OUT}")
metadata_file = "${METADATA}"
prefix = "${PREFIX}"

# -------------------------
# Collect kallisto metrics
# -------------------------
qc = {}
pct = []

for run_info in base.glob("SRR*/run_info.json"):
    srr = run_info.parent.name
    try:
        with open(run_info) as f:
            d = json.load(f)
        p = float(d.get("p_pseudoaligned"))
        pct.append(p)

        if p >= 45:
            flag = "PASS"
        elif p >= 30:
            flag = "LOW"
        else:
            flag = "FAIL"

        qc[srr] = (p, flag)
    except Exception:
        qc[srr] = (None, "FAIL")

n = len(pct)
if n == 0:
    print("❌ No valid run_info.json files found")
    sys.exit(1)

mean_pct = statistics.mean(pct)
median_pct = statistics.median(pct)
low_map = sum(x < 60 for x in pct)

qc_status = "PASSED"
if mean_pct < 60 or low_map > 0.2 * n:
    qc_status = "WARNING"

# -------------------------
# Print + save summary
# -------------------------
summary = f"""
Kallisto QC summary
-------------------
Samples processed: {n}
Mean pseudoalignment rate: {mean_pct:.1f}%
Median pseudoalignment rate: {median_pct:.1f}%
Samples < 60% mapped: {low_map}
QC status: {qc_status}
"""

print(summary.strip())

with open(f"{prefix}_summary.txt", "w") as f:
    f.write(summary.strip() + "\\n")

# -------------------------
# Annotate metadata
# -------------------------
out_meta = f"{prefix}_with_kallisto_qc.csv"

with open(metadata_file, newline="") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames + ["p_pseudoaligned", "kallisto_qc_flag"]
    rows = []

    for row in reader:
        srr = row["Run"]
        p, flag = qc.get(srr, (None, "FAIL"))
        row["p_pseudoaligned"] = f"{p:.1f}" if p is not None else "NA"
        row["kallisto_qc_flag"] = flag
        rows.append(row)

with open(out_meta, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

# -------------------------
# QC plots
# -------------------------
plt.figure()
plt.hist(pct, bins=20)
plt.xlabel("Pseudoalignment rate (%)")
plt.ylabel("Number of samples")
plt.title("Kallisto pseudoalignment rate")
plt.tight_layout()
plt.savefig(f"{prefix}_hist.png")
plt.close()

plt.figure()
plt.violinplot(pct, showmeans=True)
plt.ylabel("Pseudoalignment rate (%)")
plt.title("Kallisto pseudoalignment rate distribution")
plt.tight_layout()
plt.savefig(f"{prefix}_violin.png")
plt.close()

print(f"✅ Metadata written to: {out_meta}")
print(f"📊 QC plots saved as: {prefix}_hist.png, {prefix}_violin.png")

if qc_status != "PASSED":
    sys.exit(2)
EOF
