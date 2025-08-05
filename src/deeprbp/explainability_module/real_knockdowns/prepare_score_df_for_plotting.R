# ================================================================
# Function: prepare_score_df_for_plotting
# Purpose: Prepare a melted and labeled dataframe for plotting
#          transcript or gene scores, distinguishing DE vs non-DE.
# Input:
#   - score_matrix: data.frame with rownames as gene/transcript IDs,
#                   columns as Gene_IDs (e.g., "ENSG000...")
#   - rbp: character string (e.g. "MBNL1")
#   - significant_ids: character vector with DE rownames (Gene_IDs or Tx_IDs)
#   - getBM: data.frame with Gene_ID ↔ Gene_name mapping
# Output:
#   - A melted data.frame with columns:
#       • RBP (Gene_name)
#       • RawScore
#       • DE_Status
#       • Scores (log10-transformed)
# ================================================================

prepare_score_df_for_plotting <- function(score_matrix, rbp, significant_ids, getBM) {
  # 1. Find matching Gene_ID for the given RBP (Gene_name)
  rbp_gene_ids <- unique(getBM$Gene_ID[getBM$Gene_name == rbp])
  
  if (length(rbp_gene_ids) == 0) {
    stop("❌ Gene_ID for RBP not found in getBM: ", rbp)
  }
  
  # 2. Filter those that are present in the score matrix
  rbp_gene_id <- rbp_gene_ids[rbp_gene_ids %in% colnames(score_matrix)][1]
  
  if (is.na(rbp_gene_id)) {
    stop("❌ No matching Gene_ID for RBP found in score matrix: ", rbp)
  }
  
  # 3. Extract column and rename to Gene_name (rbp)
  df_scores <- score_matrix[, rbp_gene_id, drop = FALSE]
  colnames(df_scores) <- rbp  # Rename column to Gene_name
  
  # 4. Add DE status
  df_scores$DE_Limma <- "No"
  df_scores[rownames(df_scores) %in% significant_ids, "DE_Limma"] <- "Yes"
  
  # 5. Melt and format
  df_melted <- reshape2::melt(
    df_scores,
    id.vars = "DE_Limma",
    variable.name = "RBP",
    value.name = "RawScore"
  ) %>%
    dplyr::mutate(
      DE_Status = factor(DE_Limma, levels = c("Yes", "No"), labels = c("DE", "non-DE")),
      Scores = log10(RawScore + 1)
    )
  
  return(df_melted)
}
