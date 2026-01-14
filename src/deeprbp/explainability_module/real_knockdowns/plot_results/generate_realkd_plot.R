
# src/deeprbp/explainability_module/real_knockdowns/plot_results/generate_realkd_plot.R

# 🧹 Clear environment
rm(list = ls())

# 🌍 Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 📦 List of required packages
required_packages <- c("data.table", "dplyr", "tibble", "tidyr")

# 🔧 Function to install missing packages
# install_if_missing <- function(packages) {
#   missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(missing_packages)) {
#     install.packages(missing_packages, dependencies = TRUE)
#   }
# }

# ⚙️ Install required packages if not already installed
# install_if_missing(required_packages)

# 📚 Load libraries
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubr)
library(rstatix)
library(latex2exp)
library(optparse)

# --------------------- LOAD FUNCTION ----------------------------
cat("📄 Loading prepare score df for ploting \n")
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
source(file.path(script_dir, "prepare_score_df_for_plotting.R"))

generate_score_plot <- function(df_melted,
                                stat.test,
                                rbp,
                                experiment,
                                tipo = c("transcripts", "genes"),
                                transform = c("log10", "none")) {
  tipo      <- match.arg(tipo)
  transform <- match.arg(transform)
  
  # Leyenda según tipo
  legend_title <- if (tipo == "genes") {
    "Genes with\nDE transcript(s)"
  } else {
    "DE transcripts"
  }
  
  # Etiqueta de Y según transformación
  y_lab <- if (transform == "log10") {
    latex2exp::TeX(r'(Scores in $log_{10}(x+1)$ scale)')
  } else {
    "Scores"
  }
  
  # y.position: un poco por encima del máximo observado en 'Scores'
  max_score <- max(df_melted$Scores, na.rm = TRUE)
  bump      <- if (is.finite(max_score) && max_score != 0) 0.10 * abs(max_score) else 0.1
  stat.test$y.position <- max_score + bump
  
  # Plot con ggplot2
  p <- ggplot(df_melted, aes(x = DE_Status, y = Scores, fill = DE_Status)) +
    geom_boxplot(
      outlier.shape = NA,
      width         = 0.55,
      color         = "#333333",
      linewidth     = 0.4
    ) +
    geom_point(
      position = position_jitter(width = 0.18, height = 0),
      shape    = 21,
      size     = 1.6,
      alpha    = 0.6,
      stroke   = 0.2,
      color    = "#333333"
    ) +
    ylab(y_lab) +
    xlab(rbp) +
    ggtitle(experiment) +
    ggpubr::stat_pvalue_manual(
      stat.test,
      label      = "p.adj.signif",
      tip.length = 0.008,
      size       = 3.2     # más pequeño
    ) +
    # menos aire arriba/abajo
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    scale_fill_manual(
      values = c("DE" = "#d95f02", "non-DE" = "#7570b3"),
      breaks = c("DE", "non-DE"),
      name   = legend_title
    ) +
    theme_bw(base_size = 11) +   # base_size más pequeño
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line          = element_blank(),
      axis.line.x.bottom = element_line(color = "black"),
      axis.line.y.left   = element_line(color = "black"),
      
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      
      axis.title.x     = element_text(margin = margin(t = 6)),
      axis.title.y     = element_text(margin = margin(r = 6)),
      
      # título algo más pequeño para que no “estire” tanto
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      legend.key       = element_rect(fill = "white", colour = NA),
      
      # recortar márgenes exteriores
      plot.margin      = margin(t = 3, r = 3, b = 3, l = 3)
    )
  
  
  return(p)
  
}


# --------------------- PARSE ARGUMENTS --------------------------
cat("🔧 Parsing command-line arguments...\n")

library(optparse)

option_list <- list(
  make_option(c("--output_dir"), type = "character", default = "./output",
              help = "Directory where output results will be saved", metavar = "character"),
  
  make_option(c("--path_getBM"), type = "character", 
              default = "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv",
              help = "Path to getBM gene-to-transcript mapping file (CSV)", metavar = "character"),
  
  make_option(c("--path_de_limma"), type="character",
              default="/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns/de_limma"),
  
  make_option(c("--path_results"), type = "character", 
              default = "/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns",
              help = "Path to the base results folder", metavar = "character"),
  
  make_option(c("--experiment_list"), type = "character", 
              default = "PRJEB39343,GSE75491,GSE77702-FUS,GSE77702-TAF15,GSE77702-TARDBP,GSE136366",
              help = "Comma-separated list of experiments", metavar = "character"),
  
  make_option(c("--rbp_interest_list"), type = "character", 
              default = "MBNL1,RBM47,FUS,TAF15,TARDBP,TARDBP",
              help = "Comma-separated list of RBPs corresponding to the experiments", metavar = "character"),
  
  make_option(c("--logfc_thresh"), type = "numeric",
              default = 1,
              help = "Absolute log2 fold-change threshold used to define DE [default %default]",
              metavar = "numeric"),
  
  make_option(c("--p_cut_type"), type = "character",
              default = "fdr",
              help = "Type of p-value cut: 'fdr' (adj.P.Val) or 'pvalue' (P.Value) [default %default]",
              metavar = "character"),
  
  make_option(c("--p_cut_value"), type = "numeric",
              default = 0.05,
              help = "Numeric threshold for the chosen p_cut_type [default %default]",
              metavar = "numeric")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# --------------------- EXTRACT VALUES --------------------------
output_dir <- args$output_dir
path_getBM <- args$path_getBM
path_de_limma <- args$path_de_limma
path_results <- args$path_results
experiment_list <- unlist(strsplit(args$experiment_list, split = ","))
rbp_interest_list <- unlist(strsplit(args$rbp_interest_list, split = ","))


logfc_thresh <- args$logfc_thresh
p_cut_type   <- match.arg(args$p_cut_type, c("fdr","pvalue"))
p_cut_value  <- args$p_cut_value

cat("📊 logFC threshold:      |log2FC| >", logfc_thresh, "\n")
cat("📊 P-cut type:           ", p_cut_type, "\n")
cat("📊 P-cut value:          ", p_cut_value, "\n\n")

# --------------------- PRINT VALUES --------------------------
cat("📁 Output directory:     ", output_dir, "\n")
cat("📁 getBM file:           ", path_getBM, "\n")
cat("📁 Results path:         ", path_results, "\n")
cat("🧪 Experiment list:      ", paste(experiment_list, collapse = ", "), "\n")
cat("🧬 RBP interest list:    ", paste(rbp_interest_list, collapse = ", "), "\n\n")

# Create output directory if not exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---------- Helper: resolver la carpeta de limma por experimento ----------
# Soporta nombres con "/" (p.ej. GSE77702/FUS) y su versión con "_"
resolve_limma_dir <- function(base_limma, experiment){
  p1 <- file.path(base_limma, experiment, "voom_limma_output")
  if (dir.exists(p1)) return(p1)
  p2 <- file.path(base_limma, gsub("/", "_", experiment), "voom_limma_output")
  if (dir.exists(p2)) return(p2)
  p3 <- file.path(base_limma, gsub("_", "/", experiment), "voom_limma_output")
  if (dir.exists(p3)) return(p3)
  return(NA_character_)
}

trans_file <- "df_scores_TxRBP.csv"
genes_file <- "df_scores_GxRBP.csv" # para hacer el plot tendremos que coger una u otra.

getBM <- read.table(path_getBM, header = TRUE, sep = ",")

plotlist_trans <- list()
plotlist_genes <- list()

safe_wilcox <- function(df, formula, alternative = "greater") {
  y <- df$RawScore
  g <- df$DE_Status
  
  if (length(unique(g)) < 2) {
    warning("⚠️ Only one DE_Status group present, skipping test.")
    return(data.frame(
      group1 = "DE", group2 = "non-DE",
      p = NA, p.adj = NA, p.adj.signif = NA
    ))
  }
  
  if (length(unique(y[g == unique(g)[1]])) < 2 ||
      length(unique(y[g == unique(g)[2]])) < 2) {
    warning("⚠️ One group has constant scores, skipping Wilcoxon test.")
    return(data.frame(
      group1 = "DE", group2 = "non-DE",
      p = NA, p.adj = NA, p.adj.signif = NA
    ))
  }
  
  out <- tryCatch({
    rstatix::wilcox_test(df, formula, alternative = alternative)
  }, error = function(e) {
    warning("⚠️ Wilcoxon test failed: ", e$message)
    return(data.frame(
      group1 = "DE", group2 = "non-DE",
      p = NA, p.adj = NA, p.adj.signif = NA
    ))
  })
  
  return(out)
}

safe_add_xy <- function(stat.test, df, xvar) {
  out <- tryCatch({
    rstatix::add_xy_position(stat.test, x = xvar)
  }, error = function(e) {
    message("⚠️ Skipping add_xy_position(): ", e$message)
    stat.test$y.position <- max(df$Scores, na.rm = TRUE) * 1.1
    stat.test
  })
  return(out)
}


# ----------------------------------------------------------------
# Loop over experiments
# ----------------------------------------------------------------
for (i in seq_along(experiment_list)) {
  # i <- 1
  experiment <- experiment_list[i]
  rbp <- rbp_interest_list[i]
  
  cat("\n🧪 Processing experiment:", experiment, "| RBP:", rbp, "\n")
  
  # Paths
  base_path <- file.path(path_results, experiment)
  l_base_path <- file.path(path_de_limma, experiment)
  
  control_path <- file.path(base_path, "control")
  limma_path <- file.path(l_base_path, "voom_limma_output")
  
  # Files in control
  trans_path <- file.path(control_path, trans_file)
  genes_path <- file.path(control_path, genes_file)
  
  # Check inputs
  if (!file.exists(trans_path)) {
    warning("❌ Transcript scores file not found: ", trans_path)
    next
  }
  if (!file.exists(genes_path)) {
    warning("❌ Gene scores file not found: ", genes_path)
    next
  }
  
  cat("📄 Reading control scores in absolute values...\n")
  df_score_trans <- data.table::fread(trans_path) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("V1") %>%
    abs()
  
  # dim(df_score_trans)
  # df_score_trans[1:5,1:5]
  # max(df_score_trans)
  
  df_score_gns <- data.table::fread(genes_path) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Gene_ID") %>%
    abs()
  cat("✅ Loaded", nrow(df_score_trans), "transcript scores and", nrow(df_score_gns), "gene scores.\n")

  # min(df_score_gns)
    
  # List files in limma output
  limma_files <- list.files(limma_path, full.names = TRUE)
  
  # Find DE files
  de_tx_file <- limma_files[grepl("^DE_transcripts_.*\\.csv$", basename(limma_files))]
  de_gx_file <- limma_files[grepl("^DE_genes_.*\\.csv$", basename(limma_files))]
  
  # Check existence
  if (length(de_tx_file) == 0) {
    warning("❌ No DE_transcripts_*.csv file found in:", limma_path)
    next
  }
  if (length(de_gx_file) == 0) {
    warning("❌ No DE_genes_*.csv file found in:", limma_path)
    next
  }
  
  cat("📄 Reading limma DE results...\n")
  de_tx <- read.csv(de_tx_file[1], row.names = 1)  # assuming only one match
  de_gx <- read.csv(de_gx_file[1], row.names = 1)
  cat("✅ Loaded", nrow(de_tx), "transcripts from DE results and", nrow(de_gx), "genes.\n")
  
  # --------------------- FILTER SIGNIFICANT TRANSCRIPTS --------------------------
  p_col <- if (p_cut_type == "fdr") "adj.P.Val" else "P.Value"
  p_label <- if (p_cut_type == "fdr") "FDR" else "P-value"
  
  cat("🔎 Filtering significant transcripts (|log2FC| >", logfc_thresh,
      ", ", p_label, " <", p_cut_value, ")...\n")
  
  # Asegurarse de que rownames están limpios
  rownames(de_tx) <- sub("\\..*", "", rownames(de_tx))
  rownames(df_score_trans) <- sub("\\..*", "", rownames(df_score_trans))
  
  stopifnot(all(c("adj.P.Val","P.Value","logFC") %in% colnames(de_tx)))
  
  df_significant <- de_tx %>%
    dplyr::filter(rownames(de_tx) %in% rownames(df_score_trans)) %>%
    dplyr::filter(.data[[p_col]] < p_cut_value,
                  abs(logFC) > logfc_thresh) %>%
    dplyr::mutate(Transcript_ID = rownames(.))
  
  significant_samples <- df_significant$Transcript_ID
  cat("✅ Found", length(significant_samples),
      "significant transcripts (|log2FC| >", logfc_thresh, ", ",
      p_label, " <", p_cut_value, ") with score.\n")
  
  # --------------------- FILTER SIGNIFICANT GENES --------------------------
  cat("🔎 Filtering significant genes (|log2FC| >", logfc_thresh,
      ", ", p_label, " <", p_cut_value, ")...\n")
  
  rownames(de_gx)      <- sub("\\..*", "", rownames(de_gx))
  rownames(df_score_gns) <- sub("\\..*", "", rownames(df_score_gns))
  
  stopifnot(all(c("adj.P.Val","P.Value","logFC") %in% colnames(de_gx)))
  
  df_significant_gns <- de_gx %>%
    dplyr::filter(rownames(de_gx) %in% rownames(df_score_gns)) %>%
    dplyr::filter(.data[[p_col]] < p_cut_value,
                  abs(logFC) > logfc_thresh) %>%
    dplyr::mutate(Gene_ID = rownames(.))
  
  significant_samples_gns <- df_significant_gns$Gene_ID
  cat("✅ Found", length(significant_samples_gns),
      "significant genes (|log2FC| >", logfc_thresh, ", ",
      p_label, " <", p_cut_value, ") with score.\n")
  
  
  # --------------------- PREPARE MELTED SCORE DATAFRAMES -------------------
  cat("📦 Preparing melted score dataframes for plotting...\n")
  # old: con transformación log10(x+1)
      # df_melted_trans <- prepare_score_df_for_plotting(df_score_trans, "MBNL1",
                                                 #  significant_samples, getBM,
                                                 #  transform = "log10")
  
  # For transcripts
  df_melted_trans <- prepare_score_df_for_plotting(
    score_matrix = df_score_trans,
    rbp = rbp,
    significant_ids = significant_samples,
    getBM = getBM, transform = "log10"
  )
  
  # For genes
  df_melted_genes <- prepare_score_df_for_plotting(
    score_matrix = df_score_gns,
    rbp = rbp,
    significant_ids = significant_samples_gns,
    getBM = getBM, transform = "log10"
  )
  
  
  # --------------------- STATISTICAL TEST FOR TRANSCRIPTS ---------------------
  cat("🧪 Performing statistical test (transcripts)...\n")
  cat("ℹ️  Comparing raw RBP scores between DE and non-DE transcripts using Mann-Whitney-Wilcoxon test.\n")
  
  # Skip if only one DE_Status present
  if (length(unique(df_melted_trans$DE_Status)) < 2) {
    warning("❌ Only one DE_Status group present for transcripts in ", rbp, ". Skipping test.")
    next
  }
  
  # --- Transcripts ---
  stat.test.trans <- safe_wilcox(df_melted_trans, RawScore ~ DE_Status, alternative = "greater") %>%
    rstatix::adjust_pvalue(method = "bonferroni") %>%
    rstatix::add_significance("p.adj")
  
  stat.test.trans <- safe_add_xy(stat.test.trans, df = df_melted_trans, xvar = "DE_Status")
  
  
  cat("✅ Transcript-level Wilcoxon test complete for RBP:", rbp, "\n")
  cat("   → Adjusted p-value:", stat.test.trans$p.adj, "| Significance:", stat.test.trans$p.adj.signif, "\n\n")
  
  #stat.test.trans <- df_melted_trans %>%
   # rstatix::wilcox_test(RawScore ~ DE_Status, alternative = "greater") %>%  # <--- use RawScore here
  #  rstatix::adjust_pvalue(method = "bonferroni") %>%
  #  rstatix::add_significance("p.adj") %>%
  #  rstatix::add_xy_position(x = "DE_Status")
  
  # --------------------- STATISTICAL TEST FOR GENES --------------------------
  cat("🧪 Performing statistical test (genes)...\n")
  cat("ℹ️  Comparing raw RBP scores between DE and non-DE genes using Mann-Whitney-Wilcoxon test.\n")
  
  if (length(unique(df_melted_genes$DE_Status)) < 2) {
    warning("❌ Only one DE_Status group present for genes in ", rbp, ". Skipping test.")
    next
  }
  
  # --- Genes ---
  stat.test.genes <- safe_wilcox(df_melted_genes, RawScore ~ DE_Status, alternative = "greater") %>%
    rstatix::adjust_pvalue(method = "bonferroni") %>%
    rstatix::add_significance("p.adj")
  
  stat.test.genes <- safe_add_xy(stat.test.genes, df = df_melted_genes, xvar = "DE_Status")
  
  cat("✅ Gene-level Wilcoxon test complete for RBP:", rbp, "\n")
  cat("   → Adjusted p-value:", stat.test.genes$p.adj, "| Significance:", stat.test.genes$p.adj.signif, "\n\n")
  
  # --------------------- PLOT TRANSCRIPTS ---------------------
  cat("🎨 Plotting RBP scores for differentially expressed transcripts...\n")
  cat("   → RBP:", rbp, "\n")
  cat("   → Experiment:", experiment, "\n")
  plotlist_trans[[i]] <- generate_score_plot(df_melted_trans, 
                                             stat.test.trans, 
                                             rbp, experiment, 
                                             tipo = "transcripts", transform = "log10")
  cat("✅ Transcript plot added to list at index", i, "\n\n")
  
  # --------------------- PLOT GENES --------------------------
  cat("🎨 Plotting RBP scores for genes with DE transcripts...\n")
  cat("   → RBP:", rbp, "\n")
  cat("   → Experiment:", experiment, "\n")
  plotlist_genes[[i]] <- generate_score_plot(df_melted_genes, 
                                             stat.test.genes, 
                                             rbp, experiment,
                                             tipo = "genes", transform = "log10")
  cat("✅ Gene plot added to list at index", i, "\n\n")
  
}

# ----------------- PLOT ARRANGE: TRANSCRIPTS -----------------
cat("🧩 Arranging transcript-level plots...\n")
labels_trans <- c("A", "B", "C", "D", "E", "F")

arranged_trans <- ggpubr::ggarrange(
  plotlist = plotlist_trans,
  nrow = 1, ncol = 3,
  labels = labels_trans,
  common.legend = TRUE,
  legend = "right"
)

cat("✅ Transcript plots arranged.\n")

transcript_file <- file.path(output_dir, "transcripts_arranged.pdf")
ggsave(
  filename = transcript_file,
  plot = arranged_trans,
  height = 40, width = 100, units = "mm", scale = 2
)
cat("💾 Transcript figure saved to:", transcript_file, "\n\n")


# ----------------- PLOT ARRANGE: GENES -----------------
cat("🧩 Arranging gene-level plots...\n")
labels_genes <- c("D", "E", "F")

arranged_genes <- ggpubr::ggarrange(
  plotlist = plotlist_genes,
  nrow = 1, ncol = 3,
  labels = labels_genes,
  common.legend = TRUE,
  legend = "right"
)

cat("✅ Gene plots arranged.\n")

genes_file <- file.path(output_dir, "genes_arranged.pdf")
ggsave(
  filename = genes_file,
  plot = arranged_genes,
  height = 40, width = 100, units = "mm", scale = 2
)
cat("💾 Gene figure saved to:", genes_file, "\n\n")



# ----------------- PLOT ARRANGE: TRANSCRIPTS + GENES (STACKED) -----------------
cat("🧩 Building combined transcripts+genes stacked figure (no panel letters)...\n")

library(grid)

# Número de paneles (asumimos misma longitud para transcritos y genes)
n_panels <- length(plotlist_trans)

# 1) TRANSCRITOS: sin ejes X/Y (eje Y global), con títulos de panel,
#    y márgenes pequeños
plots_trans_clean <- lapply(plotlist_trans, function(p) {
  p +
    labs(x = NULL, y = NULL) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin  = margin(t = 2, r = 4, b = 4, l = 4)
    )
})

# 2) GENES: sin eje Y (global), sin títulos de panel,
#    con eje X activo y márgenes pequeños
plots_genes_clean <- lapply(plotlist_genes, function(p) {
  p +
    labs(y = NULL) +
    theme(
      axis.title.y = element_blank(),
      plot.title   = element_blank(),
      plot.margin  = margin(t = 4, r = 4, b = 2, l = 4)
    )
})

# 3) Fila superior (transcritos) con su leyenda, SIN labels de panel
row_trans <- ggpubr::ggarrange(
  plotlist      = plots_trans_clean,
  nrow          = 1,
  ncol          = n_panels,
  common.legend = TRUE,
  legend        = "right",
  align         = "v"
)

# 4) Fila inferior (genes) con su propia leyenda, SIN labels de panel
row_genes <- ggpubr::ggarrange(
  plotlist      = plots_genes_clean,
  nrow          = 1,
  ncol          = n_panels,
  common.legend = TRUE,
  legend        = "right",
  align         = "v"
)

# 5) Apilamos filas (transcritos arriba, genes abajo) con el mismo tamaño que te gustaba
combined_base <- ggpubr::ggarrange(
  row_trans,
  row_genes,
  ncol    = 1,
  heights = c(0.95, 0.95)
)

# 6) Añadimos eje Y global único
y_lab_global <- latex2exp::TeX(r'(Scores in $log_{10}(x+1)$ scale)')

combined_annotated <- ggpubr::annotate_figure(
  combined_base,
  left = grid::textGrob(
    label = y_lab_global,
    rot   = 90,
    vjust = 0.5,
    gp    = grid::gpar(cex = 1.0)
  )
)

# 7) Guardar la figura combinada
combined_file <- file.path(output_dir, "scores_transcripts_genes_stacked.pdf")
ggsave(
  filename = combined_file,
  plot     = combined_annotated,
  height   = 80,   # mantenemos el tamaño que te gustaba
  width    = 100,
  units    = "mm",
  scale    = 1.7
)

cat("💾 Combined transcripts+genes stacked figure saved to:", combined_file, "\n\n")
