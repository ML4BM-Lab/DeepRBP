## Complex Analysis: CORUM complex-level correlation of explainability scores
This submodule evaluates whether RBPs belonging to the same protein complex (annotated in **CORUM**) show more similar **DeepRBP explainability profiles across genes** than RBPs outside that complex.
In the paper, we apply this analysis to **spliceosome subcomplexes** (C, pre-B, B, A, E) and to multiple **TCGA tumor types**.

### Contents
- Overview  
- Inputs  
- Method  
- Outputs  
- How to run  
- Notes and common pitfalls  

---

## Overview
Given a **G×RBP explainability score matrix** (e.g. `df_scores_GxRBP.csv`) computed with the main explainability pipeline, this module:

1. Computes an **RBP–RBP Pearson correlation matrix** using gene-wise explainability vectors (one vector per RBP).
2. For each CORUM spliceosomal complex, compares:
  - **Within-complex correlations** (Complex pairs)
  - **Complex vs outside correlations** (Not-complex pairs)
3. Tests whether **within-complex correlations are higher** (one-sided test).
4. Combines per-complex significance within a tissue using **Stouffer’s method**  
  (and also reports **Fisher’s method** in the CLI output).

---

## Inputs
### 1) CORUM complexes file (`corum_results.txt`)
Tab-separated file (CORUM export / curated snapshot) including at least:

- `complex_id`
- `complex_name`
- `subunits_gene_name`
- `subunits_gene_name_synonyms`

> The script parses `subunits_gene_name` and `subunits_gene_name_synonyms` as **semicolon-separated lists**.

### 2) Gene mapping (`getBM.csv`)

CSV used across DeepRBP. Must include at least:

- `Gene_ID` (Ensembl gene ID, e.g. `ENSG...`)
- `Gene_name` (HGNC-like symbol)

This file is used to map **CORUM subunit names / synonyms → Ensembl `Gene_ID`**.

### 3) Explainability scores (`df_scores_GxRBP.csv`)
A **gene × RBP** matrix produced by the main explainability module.

- Rows: `Gene_ID` (Ensembl)
- Columns: RBP identifiers (typically Ensembl IDs)
- Values: explainability scores

**Important (code behavior):**  The script currently applies `abs(scores)` before computing correlations, i.e. it measures **similarity in magnitude of regulation patterns**, independent of direction.

---

## Method

### A) Correlation matrix
For a given tumor type / cohort:

1. Load `df_scores_GxRBP.csv`
2. Take absolute values: `abs(scores_GxRBP)`
3. Compute: `corr_scores = scores_GxRBP.corr()``
 → Pearson correlations between RBPs, using their gene-wise score vectors.

### B) Complex vs outside comparisons
We focus on CORUM spliceosome complexes (default order in code):

- **C complex**: `complex_id = 8369`
- **pre-B**: `8370`
- **B**: `8371`
- **A**: `8372`
- **E**: `8391`

For each complex *K*:

- **Complex pairs**  
  All unique RBP–RBP pairs within complex *K*  
  (upper triangle of the *K×K* block; diagonal excluded)

- **Not-complex pairs**  
  All RBP–RBP pairs between RBPs in complex *K* and RBPs outside *K*  
  (*K × notK* block; members of *K* excluded by definition)

#### Statistical test
- **Mann–Whitney U test**, one-sided (`alternative = "greater"`)
- Hypothesis:  
  *within-complex correlations > outside correlations*

#### Across complexes (per tissue)
- Combine per-complex p-values using **Stouffer’s method**  
  (reported and used in the manuscript)
- The script also reports the **Fisher combined p-value** for reference.

### C) Name → ID mapping
CORUM provides gene names and synonyms; these are mapped to `Gene_ID`
using `getBM.csv` as follows:

1. Direct match on **`Gene_name`**
2. If not found, try **comma-separated synonyms**
3. Unmapped genes are set to `None` and are ignored downstream

**Note:**  
Genes present in multiple CORUM complexes are **not removed** and may therefore
appear in multiple groups (see code note).

---

## Outputs
When `--output_dir` is provided, the script saves:

- **`correlation_overview_long.csv`** (or `correlation_overview_long_withF.csv`)  
  Long/tidy table used for overview plotting (*Complex vs Outside* distributions)

- **`per_pair_info.csv`** (or `per_pair_info_withF.csv`)  
  Per-complex p-values and plotting helper values (e.g., y-axis placement)

- **`correlation_boxplot_overview.png`** (or `correlation_boxplot_overview_withF.png`)  
  Overview plot comparing *Complex vs Outside* across complexes

Additionally, **per-complex plots** (one per complex) are generated if plotting
is enabled in your environment (see `plots_complex.py`).

### Console output
The console output includes:

- Per-complex **U statistic** and **p-value**
- **Combined p-values** (Stouffer and Fisher)

---

## How to run
### Run from Python module (recommended)

```bash
python -m deeprbp.explainability_module.complex_analysis.run_corum_complex_analysis \
  --corum_file_path "src/deeprbp/explainability_module/complex_analysis/corum_results.txt" \
  --getBM_file_path "/scratch/jsanchoz/DeepRBP/data/annotation/getBM_gencode_v23.csv" \
  --scores_file_path "<EXPLAINER_OUTPUT>/<TUMOR>/df_scores_GxRBP.csv" \
  --output_dir "<EXPLAINER_OUTPUT>/<TUMOR>/corum_complex_analysis" \
  --include-non-family
```

---

### What does --include-non-family do?
Adds a group F containing **all RBPs not assigned** to any of the A–E spliceosome complexes. This is useful because:

- without `--include-non-family`, the “outside” pool is limited to RBPs that appear in the A–E groups (more restrictive)
- with `--include-non-family`, the “outside” pool covers the broader RBP universe (recommended for the paper-style analysis)

## SLURM / HPC
This submodule is typically run on HPC using predefined SLURM job scripts, one per tumor type.
These scripts wrap the Python entry point with the correct paths, resources, and parameters.

### Available SLURM scripts
| Tumor type                     | TCGA code | SLURM script                         |
| ------------------------------ | --------- | ------------------------------------ |
| Acute Myeloid Leukemia         | LAML      | `run_corum_complex_analysis_LAML.sh` |
| Kidney Chromophobe             | KICH      | `run_corum_complex_analysis_KICH.sh` |
| Liver Hepatocellular Carcinoma | LIHC      | `run_corum_complex_analysis_LIHC.sh` |

All scripts are located under:

```bash
slurm/explainability_module/complex_analysis/
```

### How to submit a job
To run the analysis for a given tumor type:

```bash
sbatch slurm/explainability_module/complex_analysis/run_corum_complex_analysis_LAML.sh
```

## Notes and common pitfalls
- **Mapping failures (CORUM → getBM)**:
If many subunits are not found, confirm `getBM.csv` contains the expected `Gene_name` symbols and that CORUM synonyms match your naming convention.

- **Score matrix identifiers:**
Columns in `df_scores_GxRBP.csv` should match the identifiers used for RBP features (in your pipeline, typically Ensembl gene IDs).

- **Interpretation of correlation:**
This analysis tests whether RBPs in the same complex have more similar gene-wise explainability magnitude profiles (since abs() is applied). If you want signed similarity, remove abs() in the script.

- **Run per tissue/tumor type:**
This module is typically executed once per tumor type (e.g., LAML, KICH, LIHC), using the corresponding df_scores_GxRBP.csv.