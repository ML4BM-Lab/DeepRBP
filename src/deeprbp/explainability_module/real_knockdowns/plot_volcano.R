#' plot_volcano
#'
#' Generate a volcano plot for transcript-level differential expression data.
#' This function highlights significantly regulated transcripts and optionally adds labels
#' using annotation data (e.g., from biomaRt).
#'
#' @param df Data frame with differential expression results. Must contain at least: logFC, adj.P.Val, and rownames or a column with transcript IDs.
#' @param id_col (Optional) String. Column name used as transcript identifier (e.g., "Transcript_ID"). If NULL, rownames are used.
#' @param getBM (Optional) Data frame with transcript annotation. Must contain "Transcript_ID" and optionally "Transcript_name".
#' @param label_column (Optional) String. Column name to label top-n significant transcripts. If NULL, uses "Transcript_name".
#' @param fdr_thresh Numeric. FDR threshold for significance (default = 0.05).
#' @param logfc_thresh Numeric. Absolute log2 fold change threshold (default = 1).
#' @param top_n Integer. Number of top significant transcripts to label (default = 10).
#' @param title String. Title of the plot (default = "Volcano Plot").
#' @param output_dir String. Folder path where the plot PDF will be saved.
#' @param save_plot Logical. Whether to save the plot as a PDF (default = TRUE).
#'
#' @return A ggplot2 volcano plot object.
#' @export
#'
#' @examples
#' plot_volcano(
#'   df = tT1,
#'   id_col = "Transcript_ID",
#'   getBM = getBM,
#'   title = "Volcano_Transcripts_KD_vs_Control",
#'   output_dir = file.path("results", "voom_limma_output"),
#'   save_plot = TRUE
#' )
plot_volcano <- function(df,
                         id_col = NULL,
                         getBM = NULL,
                         label_column = NULL,
                         fdr_thresh = 0.05,
                         logfc_thresh = 1,
                         top_n = 5,
                         title = "Volcano Plot",
                         output_dir = ".",
                         save_plot = TRUE) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Use rownames as ID if id_col is not provided
  if (is.null(id_col)) {
    df[["Transcript_ID"]] <- rownames(df)
    id_col <- "Transcript_ID"
  } else {
    if (!(id_col %in% colnames(df))) {
      df[[id_col]] <- rownames(df)
    }
  }
  
  # Ensure ID is character
  df[[id_col]] <- as.character(df[[id_col]])
  
  # Calculate -log10 FDR
  df$log10FDR <- -log10(df$adj.P.Val)
  
  # Determine significance
  df$significant <- with(df, adj.P.Val < fdr_thresh & abs(logFC) > logfc_thresh)
  
  # Merge with annotation if available
  if (!is.null(getBM)) {
    getBM$Transcript_ID <- as.character(getBM$Transcript_ID)
    
    if (!("Transcript_ID" %in% colnames(getBM))) {
      stop("getBM must contain a 'Transcript_ID' column.")
    }
    
    df <- merge(df, getBM, by.x = id_col, by.y = "Transcript_ID", all.x = TRUE)
    
    if (is.null(label_column)) {
      label_column <- "Transcript_name"
    }
    
    if (!(label_column %in% colnames(df))) {
      warning(paste("⚠️ Label column", label_column, "not found; using ID column instead."))
      label_column <- id_col
    }
    
    if (all(is.na(df[[label_column]]))) {
      warning("⚠️ No labels found after merging. Check your annotation table.")
    }
  } else {
    if (is.null(label_column)) {
      label_column <- id_col
    }
  }
  
  # Top labels
  top_labels <- df %>%
    filter(significant) %>%
    arrange(adj.P.Val) %>%
    slice(1:top_n)
  
  # Volcano plot
  p <- ggplot(df, aes(x = logFC, y = log10FDR)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_thresh), col = "blue", linetype = "dashed") +
    geom_text_repel(data = top_labels,
                    aes_string(label = label_column),
                    size = 3,
                    max.overlaps = 20,
                    box.padding = 0.4,
                    segment.color = "grey50",
                    na.rm = TRUE) +
    scale_color_manual(values = c("gray70", "firebrick")) +
    theme_minimal() +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10 FDR",
         color = "Significant")
  
  # Save plot
  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    title_clean <- gsub("[^a-zA-Z0-9]", "_", title)
    file_name <- paste0("volcano_", title_clean, ".pdf")
    full_path <- file.path(output_dir, file_name)
    
    ggsave(filename = full_path, plot = p, width = 8, height = 6)
    message("✅ Volcano plot saved to: ", full_path)
  }
  
  return(p)
}
