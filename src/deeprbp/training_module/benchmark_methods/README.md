
# Benchmarking DeepRBP against classical machine learning methods
This module benchmarks the **DeepRBP predictor** against a set of traditional machine learning regression models.

The goal is to assess whether the deep learning–based approach provides improved transcript abundance prediction 
compared to classical baselines when using identical training and validation splits.

---

## Methods
We evaluate the following regression algorithms using a `MultiOutputRegressor` framework:

- Support Vector Regression (SVR)
- Decision Tree Regressor
- Elastic Net
- Ridge Regression

Each model predicts isoform-level abundances, which are subsequently scaled by the corresponding gene expression (TPM) values 
to obtain transcript-level TPMs. As in DeepRBP, predictions are evaluated on `log2(TPM + 1)` values.

---

## Evaluation metrics
Performance is assessed using:

- R²
- Mean Squared Error (MSE)
- Pearson correlation
- Spearman correlation

Metrics are computed:
- across all transcripts,
- and aggregated at the gene level.

---

## Experimental design
To ensure a fair comparison with DeepRBP, all benchmark methods reuse the
**same pre-defined training and validation splits** generated during the DeepRBP training pipeline.

Results for all benchmark methods are saved as CSV tables and used to generate manuscript figures.

---

## Configuration
Benchmark experiments are controlled via a YAML configuration file:

```yaml
# src/deeprbp/configs/config_benchmark_methods.yaml
pre_split_train_path: "..."
pre_split_val_path: "..."
getBM_path: "data/annotations/getBM.csv"
gene_col_name: "Gene_ID"
trans_col_name: "Transcript_ID"
plot_results: False
```

## Running the benchmark (HPC)
Benchmark experiments are typically executed on an HPC cluster using **a single, generic SLURM script**.
The specific machine learning algorithm is selected via a command-line argument.

### Available algorithms
The following algorithms are currently supported:

- `svr`
- `decision_tree`
- `random_forest`
- `gradient_boosting`
- `xgboost`
- `lightgbm`
- `knn`
- `elastic_net`
- `ridge`

### SLURM execution
A SLURM script example is provided:

```bash
sbatch slurm/training_module/benchmarking/run_benchmark_svr.sh
```

Inside the SLURM script, the algorithm is specified via the `--algorithm` argument:
```bash
python -u -m deeprbp.training_module.benchmark_methods.train_baselines \
  --config_path src/deeprbp/configs/config_benchmark_methods.yaml \
  --output_dir output/results/run_model_benchmark \
  --algorithm svr
```

To benchmark a different method, simply change the value of `--algorithm` (e.g., `ridge`, `elastic_net`, `random_forest`).

Each algorithm is executed independently and writes its results to a method-specific subdirectory under the specified output directory.