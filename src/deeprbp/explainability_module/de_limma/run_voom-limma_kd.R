# ================================================================
# run_voom-limma_kd.R
# ================================================================
# Differential expression analysis for REAL knockdown experiments
# using voom-limma.
#
# Responsibilities of THIS script:
#  - parse CLI arguments (KD-specific)
#  - load kallisto outputs
#  - build transcript-level expression matrix
#  - build metadata (control vs KD)
#  - call voom_limma_core()
#  - annotate + summarize (gene-level)
#  - save results
#  - generate volcano plots
#
# Statistical core is delegated to:
#   deeprbp/explainability_module/de_limma/voom_limma_core.R
# ================================================================
rm(list = ls())

# ------------------------------------------------
# CRAN + packages
# ------------------------------------------------
chooseCRANmirror(graphics = FALSE, ind = 1)
options(repos = c(CRAN = "https://cloud.r-project.org/"))

required_packages <- c(
  "optparse",
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2"
)

install_if_missing <- function(pkgs) {
  miss <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(miss)) install.packages(miss, dependencies = TRUE)
}
install_if_missing(required_packages)

library(optparse)
library(data.table)
library(dplyr)
library(tidyr)

# ------------------------------------------------
# Parse CLI arguments
# ------------------------------------------------
cat("🔧 Parsing command-line arguments...\n")

option_list <- list(
  make_option(c("--path"), type = "character", help = paste(
      "Path to the real knockdown dataset directory.",
      "This directory MUST contain:",
      "- kallisto_output/<sample_id>/abundance.tsv",
      "- info_samples.txt (with run_accession and sample_title columns)."
    ),
    metavar = "character"
  ),
  make_option(c("--output_dir"), type = "character", help = "Directory where all differential expression results will be written.",
    metavar = "character"
  ),
  make_option(c("--path_getBM"), type = "character", help = paste(
      "Path to getBM CSV file mapping Transcript_ID to Gene_ID.",
      "Used for gene-level aggregation and volcano plot annotation."
    ), metavar = "character"
  ),
  make_option(c("--condition_control"), type = "character", help = "Label of the CONTROL condition as it appears in info_samples.txt.",
    metavar = "character"
  ),
  make_option(c("--condition_kd"), type = "character", help = "Label of the KNOCKDOWN condition as it appears in info_samples.txt.",
    metavar = "character"
  ),
  make_option(c("--logfc_thresh"), type = "numeric", default = 1, help = "Absolute log2 fold-change threshold used for volcano plot coloring.",
    metavar = "numeric"
  ),
  make_option(c("--p_cut_type"), type = "character", default = "fdr", help = "Significance metric for volcano plots: 'fdr' or 'pvalue'.",
    metavar = "character"
  ),
  make_option(c("--p_cut_value"), type = "numeric", default = 0.05, help = "Threshold for the selected significance metric.",
    metavar = "numeric"
  ),
  make_option(c("--show_legend"), type = "logical", default = TRUE, help = "Whether to display the color legend in volcano plots."
  )
)

args <- parse_args(OptionParser(option_list = option_list))

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("🧪 Running KD DE:",
    args$condition_kd, "vs", args$condition_control, "\n")

# ------------------------------------------------
# Source shared modules
# ------------------------------------------------
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
source(file.path(script_dir, "voom_limma_core.R"))
source(file.path(script_dir, "plot_volcano.R"))

# ------------------------------------------------
# Load sample metadata (KD-specific)
# ------------------------------------------------
info_sample <- fread(
  file.path(args$path, "info_samples.txt"),
  sep = "\t",
  header = TRUE,           
  data.table = TRUE,
  check.names = FALSE
)

stopifnot(all(c("run_accession", "sample_title") %in% colnames(info_sample)))
info_sample[, sample_title := trimws(sample_title)]

# ------------------------------------------------
# Helper: clean kallisto abundance → transcript FPKM
# ------------------------------------------------
clean_abundance_df_fpkm <- function(df, sample_name) {

  total_est_counts <- sum(df$est_counts)

  df %>%
    tidyr::separate(
      target_id,
      into = c("Transcript_ID", paste0("V", 2:8)),
      sep = "\\|",
      extra = "drop"
    ) %>%
    mutate(
      Transcript_ID = sub("\\..*", "", Transcript_ID)
    ) %>%
    group_by(Transcript_ID) %>%
    summarise(
      est_counts = sum(est_counts),
      eff_length = mean(eff_length),
      .groups = "drop"
    ) %>%
    mutate(
      !!sample_name := (est_counts / total_est_counts / eff_length) * 1e9
    ) %>%
    select(Transcript_ID, !!sym(sample_name))
}

# ------------------------------------------------
# Load samples per condition
# ------------------------------------------------
load_condition_matrix <- function(condition_label) {
  samples <- info_sample$run_accession[
    info_sample$sample_title == condition_label
  ]

  if (length(samples) == 0) {
    stop("No samples found for condition: ", condition_label)
  }

  first <- fread(
    file.path(args$path, "kallisto_output", samples[1], "abundance.tsv")
  )
  mat <- clean_abundance_df_fpkm(first, samples[1])

  if (length(samples) > 1) {
    for (s in samples[-1]) {
      df <- fread(
        file.path(args$path, "kallisto_output", s, "abundance.tsv")
      )
      mat <- left_join(
        mat,
        clean_abundance_df_fpkm(df, s),
        by = "Transcript_ID"
      )
    }
  }

  return(list(
    matrix = mat,
    samples = samples
  ))
}

cat("📥 Loading control samples...\n")
ctrl <- load_condition_matrix(args$condition_control)

cat("📥 Loading KD samples...\n")
kd <- load_condition_matrix(args$condition_kd)


# ------------------------------------------------
# Build expression matrix (features x samples)
# ------------------------------------------------
Expression <- cbind(
  as.matrix(ctrl$matrix[, -1]),
  as.matrix(kd$matrix[, -1])
)

rownames(Expression) <- ctrl$matrix$Transcript_ID
colnames(Expression) <- c(ctrl$samples, kd$samples)

cat("✅ Expression matrix:",
    nrow(Expression), "features ×",
    ncol(Expression), "samples\n")

# ------------------------------------------------
# Build metadata for limma core
# ------------------------------------------------
metadata <- data.frame(
  sample_id = colnames(Expression),
  condition = c(
    rep("control", length(ctrl$samples)),
    rep("kd", length(kd$samples))
  ),
  row.names = colnames(Expression)
)

# ------------------------------------------------
# Run voom-limma core
# ------------------------------------------------
res <- voom_limma_core(
  counts_matrix = Expression,
  metadata = metadata,
  condition_col = "condition",
  condition_a = "control",
  condition_b = "kd",
  feature_level = "transcripts"
)

# ------------------------------------------------
# Annotation (getBM) + gene-level summary
# ------------------------------------------------
getBM <- fread(args$path_getBM)

annot <- merge(
  res$top_table,
  getBM,
  by.x = "Feature_ID",
  by.y = "Transcript_ID"
)

# For each gene, keep the transcript with the lowest adj.P.Val
gene_summary <- annot %>%
  group_by(Gene_ID) %>%
  arrange(adj.P.Val) %>%
  slice(1) %>% # select most significant transcript per gene
  ungroup()

# ------------------------------------------------
# Save results
# ------------------------------------------------
out_dir <- file.path(args$output_dir, "voom_limma_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  res$top_table,
  file.path(out_dir, "DE_transcripts_KD_vs_Control.csv"),
  row.names = FALSE
)

write.csv(
  gene_summary,
  file.path(out_dir, "DE_genes_KD_vs_Control.csv"),
  row.names = FALSE
)

# ------------------------------------------------
# Volcano plot
# ------------------------------------------------
plot_volcano(
  df = res$top_table,
  getBM = getBM,
  path_dataset = args$output_dir,
  p_cut_type = args$p_cut_type,
  p_cut_value = args$p_cut_value,
  logfc_thresh = args$logfc_thresh,
  output_dir = out_dir,
  file_stem = "volcano_transcripts_KD_vs_Control",
  show_legend = args$show_legend
)

cat("🎉 KD differential expression completed successfully.\n")

