# DeepRBP Explainability Module
This module computes explainability scores on top of a trained **DeepRBP Predictor**, producing:

- **TxRBP** (transcript × RBP) explainability scores
- **GxRBP** (gene × RBP) explainability scores (collapsed from TxRBP)
- Optional:
    - **TxRBP per-sample tensor** (T×RBP×N) saved as a MultiIndex CSV (DeepLIFT-only)
    - **HL×RBP** (last hidden layer neurons × RBP) attributions (DeepLIFT-only)

It supports two complementary paradigms:
1. **DeepLIFT (Captum)**: attribution against a baseline (reference)
2. **Pseudoknockdown**: in-silico perturbations (Kout/Kup/half-Kout) + per-transcript statistics

This README includes the TCGA score recipe, category mode, CLI/YAML precedence, metadata filtering fallback, and 
reproduction checklist, while pointing advanced downstream analyses to sub-READMEs.

## Table of Contents
| Section                                                                                               | What you’ll find                                                 |
| ----------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| [1. Conceptual overview](#what-the-explainer-computes-conceptually)                                   | What is being explained, and how to interpret the sign           |
| [2. Outputs](#output-files)                                                                           | Produced CSVs (TxRBP, GxRBP, per-sample, hidden layer)           |
| [3. Quickstart](#quickstart-run-the-explainer)                                                        | CLI examples (ALL vs per-category, multi-category)               |
| [4. TCGA score recipe (DeepLIFT)](#tcga-score-computation-recipe-deeplift-default)                    | Reference choice, t-stat reduction, filtering, gene collapse     |
| [5. Alternative method: Pseudoknockdown](#alternative-method-pseudoknockdown-in-silico-perturbation)  | Knockout/knockup simulation + log2FC t-tests                     |
| [6. Category mode](#category-mode-tcga-like-metadata)                                                 | CLI/YAML precedence, required keys, folder layout                |
| [7. Condition filtering + fallback](#disease_condition--select_condition-filters--automatic-fallback) | Robust loading when labels don’t match Primary_Tumor conventions |
| [8. Hidden layer attributions](#hidden-layer-attributions-hlxrbp)                                     | DeepLIFT-only HL×RBP and constraints                             |
| [9. YAML cheat-sheet](#configuration-cheat-sheet-yaml)                                                | Key config parameters and meanings                               |
| [10. TCGA context + checklist](#tcga-experiments-tumor-types-and-disease-conditions-paper-context)    | Tumor types used and minimal reproducibility list                |
| [11. Next steps](#where-to-go-next-advanced-analyses-live-in-submodules)                              | Links to POSTAR / knockdowns / complex analysis                  |
| [References](#references-key)                                                                         | Core citations                                                   |

---

## What the explainer computes (conceptually)
**DeepRBP predicts transcript abundance** (in the model’s output space, aligned to training; typically `log2(TPM+1)`-like). The explainer asks:

> For a given cohort (e.g., TCGA tumor type), which RBPs contribute most to increasing/decreasing each transcript (and each **gene**)?

**Sign convention**
- **Positive score**: the RBP **positively** regulates the model prediction for that transcript/gene (increases predicted abundance).
- **Negative score**: the RBP **negatively** regulates the model prediction (decreases predicted abundance).

These are model-attribution signs (direction of predicted change), not necessarily causal regulation in vivo.

---

## Output files
Depending on flags/config:

- `df_scores_TxRBP.csv`
TxRBP matrix (rows: transcripts, columns: RBPs)

- `df_scores_GxRBP.csv`
GxRBP matrix (rows: genes, columns: RBPs)

- `result_table.csv`
For each (gene, RBP), records which transcript drove the gene score (max-|score|), plus metadata (gene/transcript names, biotype, #transcripts per gene).

Optional:
- `df_scores_TxRBP_per_sample.csv`
Only if `save_per_sample_scores: true`. MultiIndex columns: `(Sample, RBP)`.

- `df_scores_HLxRBP.csv`
Only with `--analyze_hidden_layer` (DeepLIFT-only). Rows are hidden neurons `(H0..H127)`, columns are RBPs.

---

## Quickstart: run the explainer
### Single-run mode (ALL samples)
If you do not provide categories, the pipeline runs once on all samples found at `test_path_files` (no TCGA metadata required).

```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/explainer_ALL
```

### Category mode (recommended for TCGA-like cohorts)
Compute per-category (e.g., per tumor type) into one subfolder per category:

```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/explainer_by_category \
  --select_category "Liver_Hepatocellular_Carcinoma"
```

### Multiple categories (comma-separated or repeat the flag):

```bash
run-deeprbp-explainer \
  --config_path src/deeprbp/configs/config_model_explain.yaml \
  --model_ckpt_path pretrained_model/model.ckpt \
  --scaler_dir pretrained_model \
  --output_dir output/explainer_by_category \
  --select_category "Liver_Hepatocellular_Carcinoma,Acute_Myeloid_Leukemia,Kidney_Chromophobe"
```

---

## TCGA score computation recipe (DeepLIFT default)
This is the exact logic we use to compute and aggregate scores for a TCGA cohort (tumor type), with alternatives noted.

1) **Per-sample DeepLIFT attributions (T×RBP×N)**
DeepLIFT attributes model outputs to RBP inputs by comparing each sample to a reference (baseline) and 
propagating reference-relative differences through the network.

Implementation notes:
- The model runs in **eval mode** to keep `BatchNorm` statistics stable during attribution.
- Attributions are computed only w.r.t. **RBP inputs**. Gene expression is passed as context (additional forward args) and not attributed by default.

The raw attribution tensor has shape:
- T transcripts × R RBPs × N samples

2) **Reference (baseline) choices (DeepLIFT)**
Configured by `reference_data`:

- **`knockout_reference` (recommended / used in TCGA experiments**)
baseline is an all-zero RBP vector (in scaled space as implemented). Interpretable as “contribution relative to RBP absence”.

Alternatives:
- **`median_reference`**
 baseline is the per-RBP median across the cohort.

- **`half_reference`**
 constant 0.5 baseline (mostly experimental / niche).

Why we prefer knockout reference for regulation maps:
- It tends to produce clearer RBP→target patterns when the goal is “regulatory influence relative to no RBP signal”.
- Median baselines can blur interpretability depending on cohort composition.

3) **Collapse across samples: batch_reduction_method (N → 1)**
After computing per-sample attributions, we collapse across samples to obtain a single **TxRBP** matrix.

Configured by `batch_reduction_method`:

- **`t-statistic` (default / recommended)**

  \[
  t = \frac{\mu}{\sigma / \sqrt{n}}
  \]

  where *n* is the number of samples.  
  This favors signals that are both **strong** and **consistent** across the cohort.

- **`sum_scores`**
  Sums scores across samples; can overemphasize cohort size and be less discriminative.

Practical interpretation:
- **t-statistic** behaves like a “cohort-level reliability-weighted effect”.
- **sum_scores** behaves like “total mass of attribution”.

4) **Expression-aware filtering (reduce artifacts)**
Before collapsing transcripts to genes, we apply two filters on **TxRBP** to avoid inflated noise from non-/low-expressed features:

1. **Never-expressed transcripts**
   If a transcript has total TPM = 0 across the cohort → set all its TxRBP scores to 0.

2. **Low-expressed genes**
   If a gene’s mean abundance is below a threshold (default: **1 TPM**) → set the corresponding transcript scores to 0 (attenuates background).

5) **Collapse transcripts to genes: `gene_collapse_method` (TxRBP → GxRBP)**
Configured by `gene_collapse_method`:

- **`max_absolute_value` (default / used in paper)**
  For each *(gene, RBP)* pair, select the transcript with **maximum |score|**, and keep the **sign** of that transcript.

### Outputs
- **`df_scores_GxRBP.csv`**  
  Final G×RBP matrix

- **`result_table.csv`**  
  Records the chosen transcript per *(gene, RBP)*, plus metadata and the number of transcripts per gene

### Interpretation note (important)
- The max-|score| strategy introduces a mild *“max-of-many”* effect: genes with more isoforms have a higher chance of showing a large absolute score.  
  This is often biologically plausible for splicing-rich genes, but should be kept in mind for ranking-based analyses.

---

## Alternative method: Pseudoknockdown (in-silico perturbation)
If `explanation_method: Pseudoknockdown`, we do not use DeepLIFT.  
Instead, for each RBP:

1. **Create perturbed cohorts** by modifying the corresponding RBP feature column:
   - **`Kout`** (set to 0)
   - **`half-Kout`** (set to 0.5)
   - **`Kup`** (set to 1)
   - **`control`** (no change)

2. **Run the model** for both conditions, convert predictions back to TPM scale if needed, and compute per-sample **log2 fold-change**.
3. For each transcript, run a **one-sample t-test** across samples to test whether mean log2FC ≠ 0.
4. Use the **transcript t-statistic** as the **TxRBP score** for that RBP.

This method is useful when a perturbation-flavored *“effect of changing one RBP”* view is desired, and is complementary to gradient/attribution-based **DeepLIFT**.

## Category mode (TCGA-like metadata)
Category mode is designed for datasets that include phenotype metadata (e.g., TCGA tumor types) and when you want **one run per category**.

## Precedence: CLI overrides YAML
Categories are decided in the following order:

1. **`--select_category` (CLI)**  
   - Can be repeated or comma-separated  
   - Parsed and deduplicated while preserving order

2. **`select_category` (YAML)**  
   - Can be a string or a list

3. **If none provided → ALL mode**  
   - Runs once on all samples (no metadata required)

---

## Required YAML key when using categories
If you provide categories (via CLI or YAML), you must define:

- **`sample_category`**: the metadata column used to select categories  
  (e.g., `detailed_category`)

If missing, the run fails with a clear error instructing you to either set  
`sample_category` or run in **ALL** mode.

---

## Output layout in category mode
For each category `<cat>`, results are written to:

```text
<output_dir>/<cat>/
  df_scores_TxRBP.csv
  df_scores_GxRBP.csv
  result_table.csv
  df_scores_TxRBP_per_sample.csv   # only if save_per_sample_scores: true
  df_scores_HLxRBP.csv              # only if --analyze_hidden_layer
```

In ALL mode, files are written directly under `<output_dir>/`.

---

## `disease_condition` / `select_condition` filters + automatic fallback
You can additionally filter **within a category** using:

- **`disease_condition`**: metadata column (e.g., `sample_type`)
- **`select_condition`**: allowed values  
  (e.g., `["Primary_Tumor"]`)

### Fallback behavior (robust to TCGA label quirks)
When running in **category mode**, the data loader proceeds as follows:

1. Load samples using **category + condition** filters
2. If **0 samples** are found → retry **without condition filters**
3. If still **0 samples** → raise an error suggesting to check metadata columns/labels

This is important for cohorts where standard  
*“Primary_Tumor / Solid_Tissue_Normal”* conventions do not apply cleanly  
(e.g., some blood-derived malignancies).

---

### Hidden-layer attributions (HL×RBP)
Enable with:

```bash
--analyze_hidden_layer
```
What it does:
- Computes attributions from **RBPs** to each neuron in the **last hidden layer**
(128-unit block), producing an **HL×RBP** matrix.

Hard constraint:
- **Only supported** when `explanation_method: "DeepLIFT"`.
If enabled with another method, the run raises an explicit error:

> Hidden-layer attributions are only supported with DeepLIFT. Set explanation_method: 'DeepLIFT' or disable --analyze_hidden_layer.

## Configuration cheat-sheet (YAML)
Below are the keys you’ll most often touch. Keep this file as your single source of truth for experiment reproducibility.

| Key                         |      Required | Typical values                                               | Notes                                                |
| --------------------------- | ------------: | ------------------------------------------------------------ | ---------------------------------------------------- |
| `test_path_files`           |             ✅ | path to Test split                                           | Cohort used to compute scores (e.g., held-out TCGA). |
| `getBM_path`                |             ✅ | path to `getBM.csv`                                          | Transcript↔gene mapping + annotations.               |
| `trans_col_name`            |             ✅ | `Transcript_ID`                                              | Column name in getBM for transcript IDs.             |
| `gene_col_name`             |             ✅ | `Gene_ID`                                                    | Column name in getBM for gene IDs.                   |
| `explanation_method`        |             ✅ | `DeepLIFT` / `Pseudoknockdown`                               | Main backend.                                        |
| `reference_data`            |      DeepLIFT | `knockout_reference` / `median_reference` / `half_reference` | Baseline for DeepLIFT.                               |
| `target_mode`               |      optional | `final` (default) / `logit`                                  | Attribute final output vs pre-sigmoid logit.         |
| `batch_reduction_method`    |             ✅ | `t-statistic` / `sum_scores`                                 | Collapse sample dimension.                           |
| `gene_collapse_method`      |             ✅ | `max_absolute_value`                                         | Collapse TxRBP → GxRBP.                              |
| `save_per_sample_scores`    |      optional | `false` (default)                                            | Writes Tx×RBP×Sample MultiIndex CSV (can be large).  |
| `sample_category`           | category mode | e.g. `detailed_category`                                     | Metadata column defining categories (tumor types).   |
| `select_category`           |      optional | string or list                                               | Used only if CLI does not provide categories.        |
| `disease_condition`         |      optional | e.g. `sample_type`                                           | Metadata column for condition filter.                |
| `select_condition`          |      optional | e.g. `["Primary_Tumor"]`                                     | Values kept under disease_condition.                 |
| `condition1` / `condition2` |   Pseudoknock | `Kout`, `control`, etc.                                      | Perturbation vs control settings.                    |

*Example (TCGA-style) YAML snippet:*
```bash
test_path_files: "/path/to/Test"
sample_category: "detailed_category"
disease_condition: "sample_type"
select_condition:
  - "Primary_Tumor"

explanation_method: "DeepLIFT"
reference_data: "knockout_reference"
target_mode: "final"

batch_reduction_method: "t-statistic"
gene_collapse_method: "max_absolute_value"
save_per_sample_scores: false

getBM_path: "/path/to/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
```

## TCGA experiments: tumor types and disease conditions (paper context)
In the paper workflow, we focused on TCGA samples that  
(i) were **not used in training** and  
(ii) which cohorts align with cell lines that have strong **CLIP/POSTAR** evidence:

- **Liver hepatocellular carcinoma**  
  (liver tumor context; HepG2 / Huh7)

- **Acute myeloid leukemia**  
  (mapped to K562 via Celligner)

- **Kidney chromophobe**  
  (mapped to HEK293 via Celligner)

#### Common condition choices
- **Primary tumor only** (when metadata is consistent), or
- Run **separate sweeps for tumor vs normal** to compare regulatory patterns.

---

## TCGA reproduction checklist (minimal + practical)
To reproduce a TCGA run end-to-end:

1. **Checkpoint**: `--model_ckpt_path /path/to/model.ckpt`
2. **Scaler directory**: `--scaler_dir /path/to/scaler_dir`
 (must include the scaler + sigma used during training)
3. **Held-out cohort**: `test_path_files: /path/to/Test`
4. **Transcript ↔ gene mapping**: `getBM_path: /path/to/getBM.csv`
5. **Category sweep (recommended)**: `--select_category "TumorA,TumorB,..."`
6. **Optional**
  - Per-sample tensor: `save_per_sample_scores: true`
  - hidden layer: add `--analyze_hidden_layer` (DeepLIFT-only)

## HPC / SLURM usage pattern (optional)
If your repo includes SLURM helpers, a typical pattern is to toggle hidden-layer analysis via an env var:

```bash
# default
cd slurm/explainability_module
sbatch run_explainer_tcga.sh

# with hidden layer (DeepLIFT-only)
ANALYZE_HL=true sbatch run_explainer_tcga.sh
```

## Where to go next (advanced analyses live in submodules)
This module produces the matrices; downstream validation/analysis is intentionally separated:

- **Complex / module analyses** (CORUM, NMF, higher-order structure)
→ `complex_analysis/README.md`

- **POSTAR3 validation** (binding vs not-binding comparisons, event→gene aggregation)
→ `postar_validation/README.md`

- **Real knockdowns** (full RNA-seq pipelines, quantification, DE/splicing comparisons)
→ `real_knockdowns/README.md`

- **"TCGA_normal vs tumor (work in progress!)**

## References (key)
- **DeepLIFT**: Shrikumar et al., 2017
- **POSTAR3**: Zhao et al., 2022
- **Celligner** (DepMap): Warren et al., 2021
- **EventPointer**: Ferrer et al., 2022
- POSTAR-based splicing regulation mapping: Carazo et al., 2019; Lobato-Fernandez et al., 2024
