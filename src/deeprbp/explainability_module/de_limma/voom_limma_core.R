
# src/deeprbp/explainability_module/de_limma/voom_limma_core.R

# ================================================================
# voom_limma_core.R
# ================================================================
# Core statistical engine for differential expression using voom-limma.
#
# This script MUST NOT:
#  - read kallisto outputs
#  - filter samples
#  - write files
#  - generate plots
#
# It only:
#  - builds the design matrix
#  - applies voom + limma
#  - returns DE tables
# ================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

# ------------------------------------------------
# Main core function
# ------------------------------------------------
voom_limma_core <- function(
  counts_matrix,        # matrix: features x samples
  metadata,             # data.frame: samples x covariates
  condition_col,        # string: column name in metadata
  condition_a,          # string: reference condition
  condition_b,          # string: comparison condition
  feature_level = c("genes", "transcripts")
) {
  feature_level <- match.arg(feature_level)

  # ------------------------------
  # Defensive checks
  # ------------------------------
  if (!condition_col %in% colnames(metadata)) {
    stop("condition_col not found in metadata: ", condition_col)
  }

  if (!all(c(condition_a, condition_b) %in% metadata[[condition_col]])) {
    stop("Conditions not found in metadata: ",
         paste(condition_a, condition_b, sep = ", "))
  }

  # Ensure same sample order
  samples <- colnames(counts_matrix)
  metadata <- metadata[samples, , drop = FALSE]

  # ------------------------------
  # Design matrix
  # ------------------------------
  group <- factor(metadata[[condition_col]], levels = c(condition_a, condition_b))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c(condition_a, condition_b)

  # Contrast: B - A
  contrast <- makeContrasts(
    contrasts = paste0(condition_b, "-", condition_a),
    levels = design
  )

  # ------------------------------
  # voom + limma
  # ------------------------------
  v <- voom(counts_matrix, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)

  # ------------------------------
  # Extract DE table
  # ------------------------------
  top_table <- topTable(
    fit,
    coef = 1,
    number = Inf,
    sort.by = "B",
    adjust.method = "fdr"
  )

  # Preserve feature IDs
  top_table$Feature_ID <- rownames(top_table)

  # ------------------------------
  # Return structured output
  # ------------------------------
  return(list(
    fit = fit,
    top_table = top_table,
    design = design,
    contrast = contrast,
    feature_level = feature_level,
    condition_a = condition_a,
    condition_b = condition_b
  ))
}
