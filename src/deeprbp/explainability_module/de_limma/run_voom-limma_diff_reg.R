
# src/deeprbp/explainability_module/de_limma/run_voom-limma_diff_reg.R

# ================================================================
# run_voom-limma_diff_reg.R
# ================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# ------------------------------------------------
# Resolve paths
# ------------------------------------------------
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))

source(file.path(script_dir, "voom_limma_core.R"))
source(file.path(script_dir, "plot_volcano.R"))

# ------------------------------------------------
# CLI
# ------------------------------------------------
option_list <- list(
  make_option("--counts_path", type="character", help="TSV file: features x samples"),
  make_option("--metadata_path", type="character", help="TSV file: samples x covariates"),
  make_option("--condition_col", type="character", help="Column in metadata defining conditions"),
  make_option("--condition_a", type="character", help="Reference condition"),
  make_option("--condition_b", type="character", help="Comparison condition"),
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--feature_level", type="character", default="genes", help="genes | transcripts"),
  make_option("--getBM_path", type="character", default=NULL, help="Optional getBM annotation"),
  make_option("--logfc_thresh", type="numeric", default=1),
  make_option("--p_cut_type", type="character", default="fdr"),
  make_option("--p_cut_value", type="numeric", default=0.05),
  make_option("--show_legend", type="logical", default=FALSE),
  #make_option("--use_fpkm", type="logical", default=FALSE, help="Whether to transform counts to FPKM before voom")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("📥 Loading counts and metadata...\n")
counts <- fread(args$counts_path)
meta   <- fread(args$metadata_path)

rownames(counts) <- counts$Feature_ID
counts$Feature_ID <- NULL

rownames(meta) <- meta$Sample_ID
meta$Sample_ID <- NULL

# Same order
counts <- counts[, rownames(meta), drop=FALSE]

cat("🧪 Running voom-limma core...\n")

res <- voom_limma_core(
  counts_matrix = as.matrix(counts),
  metadata      = meta,
  condition_col = args$condition_col,
  condition_a   = args$condition_a,
  condition_b   = args$condition_b,
  feature_level = args$feature_level
)

top_table <- res$top_table

out_csv <- file.path(
  args$output_dir,
  paste0("DE_", args$feature_level, "_",
         args$condition_b, "_vs_", args$condition_a, ".csv")
)

write.csv(top_table, out_csv, row.names = FALSE)
cat("💾 DE table written to:", out_csv, "\n")

plot_volcano(
  df            = top_table,
  getBM         = if (!is.null(args$getBM_path)) fread(args$getBM_path) else NULL,
  path_dataset  = args$output_dir,
  p_cut_type    = args$p_cut_type,
  p_cut_value   = args$p_cut_value,
  logfc_thresh  = args$logfc_thresh,
  output_dir    = args$output_dir,
  file_stem     = paste0("volcano_", args$feature_level),
  show_legend   = args$show_legend
)

cat("🎉 Differential regulation DE completed successfully.\n")
