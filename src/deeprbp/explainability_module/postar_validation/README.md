
# POSTAR validation: experimental support for DeepRBP explainability scores
This submodule validates **DeepRBP explainability scores** against **experimental RNA–protein binding evidence** from **POSTAR3** (Zhao et al., 2022), a comprehensive database of CLIP-seq–derived RBP binding sites.

The goal is to assess whether **gene–RBP** pairs predicted by **DeepRBP** to be strongly regulated are enriched for **experimentally observed binding events**, and to quantify this agreement using statistical tests and ROC-based metrics.

This module is used in the manuscript as the **core experimental validation** of DeepRBP-derived regulatory scores.

---

## Conceptual overview
### What is being validated?
DeepRBP produces **gene × RBP (G×RBP) explainability scores**, where larger absolute values indicate stronger regulatory influence of an RBP on a gene.

This module tests the hypothesis:

> Gene–RBP pairs labeled as **Binding** in POSTAR3 should receive **higher DeepRBP scores** than pairs labeled as **Not Binding**.

### What is POSTAR3?
POSTAR3 provides experimentally observed RBP–RNA binding sites derived from CLIP-seq experiments. Using POSTAR3 together with **EventPointer** annotations, we construct **binary gene–RBP matrices**:

- **1 (Binding)**: at least one splicing event within the gene is bound by the RBP
- **0 (Not Binding)**: the RBP was profiled but no binding was detected
- **NA (Unknown)**: the RBP was not profiled in the relevant tissue/cell line

Only **gene-level matrices** (G×RBP) are used, as event-to-transcript mapping is ambiguous.

---

## Tissue and cohort mapping (paper context)
POSTAR3 coverage is highly cell-line dependent. Based on availability of CLIP experiments, validation is performed using:

| TCGA tumor                            | Cell line(s) used in POSTAR3 |
| ------------------------------------- | ---------------------------- |
| Liver hepatocellular carcinoma (LIHC) | HepG2, Huh7                  |
| Acute myeloid leukemia (LAML)         | K562                         |
| Kidney chromophobe (KICH)             | HEK293                       |

For LAML and KICH, **Celligner** is used to map cell lines to the closest TCGA tumor samples.
Only **primary tumor samples** are used.

---

## Data access
Raw POSTAR3-derived data required for this module is distributed via Zenodo:

> **Zenodo link**: https://zenodo.org/uploads/15337302

Expected files:
- `human.txt`
POSTAR3 binding data (CLIP peaks), tissue/cell-line annotated

- `Events_Regions_gc23_400nt.RData`
Genomic regions defining splicing events (EventPointer)

- `EventsFound_gencode23.txt`
Event metadata (event type, genomic coordinates, IDs)

## Generating POSTAR gene × RBP matrices
POSTAR data must be preprocessed into tissue-specific **G×RBP binary matrices** using:
```bash
src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R
```

**Example: liver (HepG2 + Huh7)**
```bash
module load R/4.3.2

Rscript create_gene_rbp_postar_matrix.R \
  --input_path /data/DeepRBP/data/explainability_module/postar3 \
  --output_path /scratch/DeepRBP/data/explainability_module/postar3/processed \
  --output_file_name human_liver \
  --postar_file human.txt \
  --events_regions_file Events_Regions_gc23_400nt.RData \
  --events_gencode_file EventsFound_gencode23.txt \
  --selected_tissue_cell_line HepG2,Huh7 \
  --getBM_path /scratch/jsanchoz/DeepRBP/data/annotation/getBM_gencode_v23.csv
```
This produces: `human_liver_GxRBP.csv``

### Cell-line mapping used in the paper
- **LAML** → `K562`
- **KICH** → `HEK293`
- **LIHC** → `HepG2,Huh7`

**Argument summary:**

| Argument | Description |
|--------|-------------|
| `--input_path` | Directory containing the raw POSTAR3 and EventPointer files (`human.txt`, event region and metadata files). |
| `--output_path` | Directory where the processed POSTAR gene × RBP matrix will be written. |
| `--output_file_name` | Base name for the output file (suffix `_GxRBP.csv` is added automatically). |
| `--postar_file` | POSTAR3 file containing CLIP-seq RBP binding sites. |
| `--events_regions_file` | RData file defining genomic regions for splicing events (EventPointer output). |
| `--events_gencode_file` | Event metadata file linking events to genes and genomic coordinates. |
| `--selected_tissue_cell_line` | Comma-separated list of POSTAR cell lines to include (e.g. `HepG2,Huh7`). |
| `--getBM_path` | Path to `getBM.csv` mapping gene symbols to Ensembl Gene IDs. |

---

## Validation workflow
The validation consists of four main steps:

1. **Align DeepRBP scores with POSTAR matrix**
2. **Integrate POSTAR labels into the explainability result table**
3. **Statistical evaluation (Wilcoxon tests)**
4. **Threshold and ROC-based performance analysis**

All steps are handled by the `run-postar-validator` CLI.

## Running the POSTAR validator

### Basic command
```bash
run-postar-validator \
  --postar_matrix_dir "/path/to/postar/processed" \
  --postar_file "human_liver_GxRBP.csv" \
  --scores_result_dir "/path/to/explainer/output/Liver_Hepatocellular_Carcinoma" \
  --output_dir "/path/to/explainer/output/Liver_Hepatocellular_Carcinoma/run_postar_validator
```

### Required inputs

| Argument              | Description                                                                            |
| --------------------- | -------------------------------------------------------------------------------------- |
| `--postar_matrix_dir` | Directory containing POSTAR G×RBP matrices                                             |
| `--postar_file`       | POSTAR binary matrix (genes × RBPs)                                                    |
| `--scores_result_dir` | Explainer output directory (must contain `df_scores_GxRBP.csv` and `result_table.csv`) |
| `--output_dir`        | Output directory for validation results                                                |

### Optional inputs

| Argument        | Description                                                               |
| --------------- | ------------------------------------------------------------------------- |
| `--pvalues_csv` | CSV with per-RBP p-values (columns: `RBP_ID` + `p_adj` / `p_value` / `p`) |
| `--verbose`     | Logging verbosity (default: 1)                                            |


---

## What the validator computes
1. **Score–POSTAR alignment**
    - Matches genes and RBPs between DeepRBP scores and POSTAR matrix
    - Non-overlapping genes/RBPs are set to `NA`
    - Only overlapping entries are used downstream

2. **Statistical validation**
 For **each RBP** and **each gene**:
    - Compare score distributions between:
        - **Binding (1)**
        - **Not Binding (0)**
    - **Mann–Whitney–Wilcoxon test**
    - **Bonferroni correction**

3. **RBP-specific thresholds (Youden’s J)**
 For each RBP:
    - Compute ROC curve using absolute scores
    - Select threshold maximizing:
        $$
        J = \mathrm{TPR} - \mathrm{FPR}
        $$

    - Compute AUC
 This enables **binary calling of regulatory relationships** from continuous DeepRBP scores.

---

## Outputs
All outputs are written to:

```bash
<output_dir>/postar_validation/
```

### Main result files
| File                         | Description                                       |
| ---------------------------- | ------------------------------------------------- |
| `result_table_completed.csv` | Explainability table augmented with POSTAR labels |
| `optimal_thresholds.csv`     | Per-RBP score thresholds                          |
| `auc_results.csv`            | Per-RBP AUC values                                |
| `count_genes_per_rbp.csv`    | Number of bound genes per RBP                     |
| `count_rbps_per_gen.csv`     | Number of RBPs bound per gene                     |


### Diagnostic plots
- Per-RBP score distributions (Binding vs Not Binding)
- ROC curves with optimal threshold
- Optional annotation with external p-values

---

## Visualization (R)
To generate publication-ready boxplots and statistical summaries:

```bash
module load R/4.3.2

Rscript run_postar_plot_generation.R \
  --input_path "<output_dir>/postar_validation" \
  --results_filename "result_table_completed.csv" \
  --count_genes_per_rbp_file "count_genes_per_rbp.csv" \
  --count_rbps_per_gen_file "count_rbps_per_gen.csv" \
  --getBM_path "/path/to/getBM" \
  --getBM_filename "getBM.csv" \
  --output_path "<output_dir>/postar_validation" \
  --output_filename "plot_score_results.pdf" \
  --save_plot TRUE
```

### Interpretation of labels
- **Binding**: POSTAR3 evidence of binding (1)
- **Not Binding**: profiled but no binding detected (0)
- **Unknown**: RBP not profiled in the relevant context (NA)

---

## Important implementation notes
- **Absolute scores are used throughout**
Validation focuses on magnitude of regulatory influence, independent of sign.

- **Unknown POSTAR entries are excluded**
Only 0 vs 1 comparisons are used for statistics and ROC analysis.

- **Gene-level validation only**
Transcript-level POSTAR matrices are not computed due to ambiguous mapping.

- **One experiment per tumor type**
The validator is run independently for LIHC, LAML, and KICH.

---

## References
- **POSTAR3**: Zhao et al., 2022, Nucleic Acids Research
- **EventPointer**: Ferrer et al., 2022
- **POSTAR-based regulation analysis**: Carazo et al., 2019; Lobato-Fernandez et al., 2024
- **Celligner: Warren et al., 2021**