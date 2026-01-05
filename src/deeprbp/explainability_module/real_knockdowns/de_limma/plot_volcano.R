
# src/deeprbp/explainability_module/real_knockdowns/de_limma/plot_volcano.R

#' plot_volcano
#'
#' Generate a publication-ready volcano plot for transcript-level differential
#' expression (DE) results. Points that meet significance and effect-size
#' thresholds are colored by the *direction and magnitude* of change (log2FC);
#' non-significant points are shown in grey. The figure includes threshold
#' guide lines, an optional experiment subtitle inferred from a dataset path,
#' and an up/down count annotation.
#'
#' @param df Data frame with DE results. Must contain columns:
#'   \code{logFC}, \code{adj.P.Val}, \code{AveExpr}. Transcript IDs can be
#'   in rownames or in a dedicated column (see \code{id_col}).
#' @param id_col (Optional) Character. Name of the column with transcript IDs
#'   (e.g., \code{"Transcript_ID"}). If \code{NULL} or missing, rownames are
#'   moved to a new column \code{"Transcript_ID"} and used as IDs. Transcript
#'   versions (e.g., \code{.1}) are stripped automatically.
#' @param getBM (Optional) Data frame with transcript annotation. Must contain
#'   \code{Transcript_ID} and may include \code{Transcript_name}. If provided,
#'   it is joined to \code{df} to enable optional labeling.
#' @param path_dataset (Optional) Character. Path to the dataset root. The last
#'   path component (e.g., \code{"GSE136366"}) is shown as a centered subtitle.
#' @param fdr_thresh Numeric. FDR threshold for significance (default \code{0.05}).
#' @param logfc_thresh Numeric. Absolute log2 fold-change threshold (default \code{1}).
#'   Both \code{fdr_thresh} and \code{logfc_thresh} must be satisfied to color a point.
#' @param top_n Integer. Number of top significant transcripts to label
#'   (ranked by \code{adj.P.Val}). Default \code{0} (no labels).
#' @param output_dir Character. Directory where figures are saved. Created if absent.
#' @param file_stem Character. Basename (without extension) for the output files.
#' @param width,height Numeric. Figure dimensions in inches (defaults \code{7} × \code{5}).
#' @param dpi Integer. PNG resolution (default \code{300}).
#' @param save_pdf,save_png Logical. Whether to save PDF and/or PNG (defaults \code{TRUE}).
#'
#' @details
#' Coloring encodes \strong{Expression change}: a divergent gradient
#' (blue → grey → red) based on capped \code{logFC} with midpoint at 0.
#' Only significant points (\code{adj.P.Val < fdr_thresh} AND
#' \code{|logFC| > logfc_thresh}) are colored; others are grey.
#' Vertical dashed lines mark \code{±logfc_thresh} and the horizontal dashed
#' line marks \code{-log10(fdr_thresh)}. The annotation \code{"Count = N (↓n, ↑n)"}
#' summarizes significant down/up points.
#'
#' If \code{getBM} is supplied, the function attempts to use
#' \code{Transcript_name} for labels; otherwise, the ID column is used.
#'
#' @return A \code{ggplot} object (also saved to disk if requested).
#'
#' plot_volcano
#' (docstring igual que antes…)
#' plot_volcano  (docstring igual que antes…)

plot_volcano <- function(
    df,
    id_col = NULL,
    getBM = NULL,
    path_dataset = NULL,
    p_cut_type = c("fdr", "pvalue"),
    p_cut_value = 0.05,
    logfc_thresh = 1, 
    top_n = 0,
    output_dir = ".",
    file_stem = "volcano_transcripts_kd_vs_control",
    width = 4, height = 3, dpi = 450,           # << panel compacto por defecto
    save_pdf = TRUE, save_png = TRUE,
    show_legend = TRUE
){
  p_cut_type <- match.arg(p_cut_type)
  stopifnot(all(c("logFC","adj.P.Val","AveExpr") %in% colnames(df)))
  
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(ggrepel)
    library(stringr); library(rlang); library(tibble); library(grid)
  })
  
  if (is.null(id_col) || !(id_col %in% names(df))) {
    df <- tibble::rownames_to_column(df, var = "Transcript_ID")
    id_col <- "Transcript_ID"
  }
  df[[id_col]] <- sub("\\..*", "", as.character(df[[id_col]]))
  
  # 2) Elegir columna de p-valor según tipo
  p_cut_type <- match.arg(p_cut_type)
  p_col_sig  <- if (p_cut_type == "fdr") "adj.P.Val" else "P.Value"
  p_label_sig <- if (p_cut_type == "fdr") "FDR" else "P-value"
  
  # p_for_sig = columna usada para decidir significancia/color (FDR o P.Value)
  df$p_for_sig <- df[[p_col_sig]]
  df$log10P  <- -log10(df$P.Value)
  
  #df <- df %>%
  #  mutate(
  #    log10FDR   = -log10(adj.P.Val),
  #    significant = (adj.P.Val < fdr_thresh) & (abs(logFC) > logfc_thresh)
  #  )
  # 3) Significancia
  df <- df %>%
    mutate(
      significant = (p_for_sig < p_cut_value) & (abs(logFC) > logfc_thresh)
    )
  
  #df <- df %>%
  #  mutate(
  #    p_for_sig = ifelse(
  #      p_cut_type == "fdr",
  #      adj.P.Val,
  #      P.Value
  #    ),
  #    log10P    = -log10(p_for_sig),
  #    significant = (p_for_sig < p_cut_value) & (abs(logFC) > logfc_thresh)
  #  )
  
  # 4) Anotación opcional
  label_column <- id_col
  
  if (!is.null(getBM)) {
    
    if (id_col == "Gene_ID" && all(c("Gene_ID","Gene_name") %in% names(getBM))) {
      # Caso gene-level: usar Gene_name
      gene_annot <- getBM %>%
        dplyr::select(Gene_ID, Gene_name) %>%
        dplyr::distinct()
      
      df <- df %>%
        dplyr::left_join(gene_annot, by = "Gene_ID")
      
      if ("Gene_name" %in% names(df)) {
        label_column <- "Gene_name"
      }
      
    } else if ("Transcript_ID" %in% names(getBM)) {
      # Caso transcript-level (como antes)
      getBM$Transcript_ID <- sub("\\..*", "", as.character(getBM$Transcript_ID))
      df <- df %>%
        dplyr::left_join(
          getBM,
          by = dplyr::join_by( !!sym(id_col) == Transcript_ID )
        )
      if ("Transcript_name" %in% names(df)) {
        label_column <- "Transcript_name"
      }
    }
  }
  
  df_lab <- df %>%
    filter(significant) %>%
    arrange(p_for_sig) %>%
    slice_head(n = top_n)
  
  #df_lab <- df %>% filter(significant) %>% arrange(adj.P.Val) %>% slice_head(n = top_n)
  
  # 5) Subtítulo
  subtitle_exp <- if (!is.null(path_dataset)) {
    str_trim(basename(normalizePath(path_dataset, winslash = "/", mustWork = FALSE)))
  } else ""
  
  # 6) Conteos
  n_up   <- sum(df$significant & df$logFC >  0, na.rm = TRUE)
  n_down <- sum(df$significant & df$logFC <  0, na.rm = TRUE)
  n_sig  <- n_up + n_down
  
  cap <- max(2, logfc_thresh)
  df$logFC_capped <- pmax(pmin(df$logFC, cap), -cap)
  
  # 7) Límites y margen en Y
  y_hi <- max(df$log10P, na.rm = TRUE)
  
  if (any(df$significant)) {
    # Subtabla solo con los significativos
    sig_df <- df %>% filter(significant)
    
    # Fila con la "peor" significancia según p_for_sig (FDR o P.Value)
    worst_row <- sig_df %>%
      arrange(desc(p_for_sig)) %>%     # p_for_sig más alto = peor significativo
      slice(1)
    
    # Usamos su P.Value para posicionar la línea en el eje Y
    worst_pvalue_for_y <- worst_row$P.Value
    y_lo <- -log10(worst_pvalue_for_y)
  } else {
    # Si no hay significativos, usamos el umbral nominal
    y_lo <- -log10(p_cut_value) #y_lo <- -log10(p_cut_value) OLD
  }
  
  y_rng <- max(1e-6, y_hi - y_lo)
  
  headroom <- 0.12 * y_rng          # ~8% del rango # margen extra arriba
  y_annot  <- y_hi + 1.15 * headroom
  
  y_lab_expr <- expression(-log[10]("P-value"))
  
  # 8) Plot
  p <- ggplot(df, aes(x = logFC, y = log10P)) +
    geom_point(aes(color = ifelse(significant, logFC_capped, NA_real_)),
               size = 1.2, alpha = 0.75, show.legend = TRUE) +
    scale_color_gradient2(
      name = "log2FC",
      low = "steelblue", mid = "grey90", high = "firebrick",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      breaks = c(-cap, 0, cap), labels = c("Down", "0", "Up"),
      na.value = "grey90"
    ) +
    guides(color = guide_colorbar(direction = "vertical",
                                  title.position = "top", title.hjust = 0.5)) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = y_lo, linetype = "dashed") +
    geom_text_repel(
      data = df_lab, aes(label = !!sym(label_column)),
      size = 2.2, box.padding = 0.3, max.overlaps = 20,
      segment.color = "grey60", seed = 123, na.rm = TRUE
    ) +
    annotate("label",
             x = Inf, y = y_annot,
             label = paste0("Count = ", n_sig, " (↓", n_down, ", ↑", n_up, ")"),
             hjust = 1.05, vjust = 1, size = 3.6,
             fill = "white", color = "black",
             label.size = 0, label.padding = unit(0.12, "lines")) +
    labs(
      title = NULL,
      subtitle = subtitle_exp,
      x = "log2(Fold Change)",
      y = y_lab_expr
    ) +
    theme_classic(base_size = 10.5) +     # base menor para panel pequeño
    theme(
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x  = element_text(margin = margin(t = 4), size = 9.5),
      axis.title.y  = element_text(margin = margin(r = 4), size = 9.5),
      axis.text.x   = element_text(size = 8.5),
      axis.text.y   = element_text(size = 8.5),
      legend.title  = element_text(size = 9.5),
      legend.text   = element_text(size = 8),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(0.28, "cm"),
      plot.margin = margin(t = 6, r = 6, b = 4, l = 8)
    ) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(y_lo, y_hi + headroom), clip = "on") +
    scale_x_continuous(breaks = seq(-10, 10, 5))
  if (!isTRUE(show_legend)) p <- p + theme(legend.position = "none")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (save_png) ggsave(file.path(output_dir, paste0(file_stem, ".png")), p,
                       width = width, height = height, units = "in", dpi = dpi)
  if (save_pdf) ggsave(file.path(output_dir, paste0(file_stem, ".pdf")), p,
                       width = width, height = height, units = "in")
}
