
# 🧹 Clear environment
rm(list = ls())

# 🌍 Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 📦 List of required packages
required_packages <- c("data.table", "dplyr", "tibble", "tidyr")

# 🔧 Function to install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}

# ⚙️ Install required packages if not already installed
install_if_missing(required_packages)

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
setwd('/Users/joseba/Desktop/data:real_kds') # Update as needed
source(paste0(getwd(),"/prepare_score_df_for_plotting.R")) 

generate_score_plot <- function(df_melted, stat.test, rbp, experiment, tipo = c("transcripts", "genes")) {
  tipo <- match.arg(tipo)
  
  # Elegir nombre de la leyenda según tipo
  legend_title <- if (tipo == "genes") {
    "Genes with\nDE transcript(s)"
  } else {
    "DE transcripts"
  }
  
  # 📏 Calcular y.position como el máximo + margen
  max_score <- max(df_melted$Scores, na.rm = TRUE)
  stat.test$y.position <- max_score + 0.1 * max_score  # por ejemplo, 10% más arriba
  
  # 🎨 Plot
  p <- ggpubr::ggboxplot(
    df_melted,
    x = "DE_Status", y = "Scores", fill = "DE_Status",
    add = "jitter", add.params = list(shape = 21, alpha = 0.5, color = NULL)
  ) +
    ylab(latex2exp::TeX(r'(Scores in $log_{10}$ scale)')) +
    xlab(rbp) +
    ggtitle(experiment) +  
    ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  
    scale_fill_manual(
      values = c("#d95f02", "#7570b3"),
      breaks = c("DE", "non-DE"),
      name = legend_title          
    ) +
    theme(
      legend.position = "right",  
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
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
  
  make_option(c("--path_results"), type = "character", 
              default = "/scratch/jsanchoz/DeepRBP/output/results/real_knockdowns",
              help = "Path to the base results folder", metavar = "character"),
  
  make_option(c("--experiment_list"), type = "character", 
              default = "PRJEB39343,GSE75491,GSE77702-FUS,GSE77702-TAF15,GSE77702-TARDBP,GSE136366",
              help = "Comma-separated list of experiments", metavar = "character"),
  
  make_option(c("--rbp_interest_list"), type = "character", 
              default = "MBNL1,RBM47,FUS,TAF15,TARDBP,TARDBP",
              help = "Comma-separated list of RBPs corresponding to the experiments", metavar = "character")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# --------------------- EXTRACT VALUES --------------------------
output_dir <- args$output_dir
path_getBM <- args$path_getBM
path_results <- args$path_results
experiment_list <- unlist(strsplit(args$experiment_list, split = ","))
rbp_interest_list <- unlist(strsplit(args$rbp_interest_list, split = ","))

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

#path_results <- '/Users/joseba/Desktop/real_knockdowns' 
    #default: /scratch/jsanchoz/DeepRBP/output/results/real_knockdowns
  
#experiment_list <- c('PRJEB39343', 'GSE75491', 'GSE77702-FUS', 'GSE77702-TAF15', 'GSE77702-TARDBP', 'GSE136366') 
#rbp_interest_list <- c('MBNL1', 'RBM47', 'FUS', 'TAF15', 'TARDBP', 'TARDBP')  

trans_file <- "df_scores_TxRBP.csv"
genes_file <- "df_scores_GxRBP.csv" # para hacer el plot tendremos que coger una u otra.

getBM <- read.table(path_getBM, header = TRUE, sep = ",")

plotlist_trans <- list()
plotlist_genes <- list()

# ----------------------------------------------------------------
# Loop over experiments
# ----------------------------------------------------------------

for (i in seq_along(experiment_list)) {
  experiment <- experiment_list[i]
  rbp <- rbp_interest_list[i]
  
  cat("\n🧪 Processing experiment:", experiment, "| RBP:", rbp, "\n")
  
  # Paths
  base_path <- file.path(path_results, experiment)
  control_path <- file.path(base_path, "control")
  limma_path <- file.path(base_path, "voom_limma_output")
  
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
  
  df_score_gns <- data.table::fread(genes_path) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Gene_ID") %>%
    abs()
  cat("✅ Loaded", nrow(df_score_trans), "transcript scores and", nrow(df_score_gns), "gene scores.\n")
  
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
  cat("🔎 Filtering significant transcripts (P.Value < 0.05)...\n")
  
  # Asegurarse de que rownames están limpios
  rownames(de_tx) <- sub("\\..*", "", rownames(de_tx))
  rownames(df_score_trans) <- sub("\\..*", "", rownames(df_score_trans))
  
  # Filtrar por presencia y por adj p-value
  df_significant <- de_tx %>%
    dplyr::filter(rownames(de_tx) %in% rownames(df_score_trans)) %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::mutate(Transcript_ID = rownames(.))
  
  significant_samples <- df_significant$Transcript_ID
  cat("✅ Found", length(significant_samples), "significant transcripts with P.Value < 0.05 and score.\n")
  
  
  # --------------------- FILTER SIGNIFICANT GENES --------------------------
  cat("🔎 Filtering significant genes (P.Value < 0.05)...\n")
  
  rownames(de_gx) <- sub("\\..*", "", rownames(de_gx))
  rownames(df_score_gns) <- sub("\\..*", "", rownames(df_score_gns))
  
  df_significant_gns <- de_gx %>%
    dplyr::filter(rownames(de_gx) %in% rownames(df_score_gns)) %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::mutate(Gene_ID = rownames(.))
  
  significant_samples_gns <- df_significant_gns$Gene_ID
  cat("✅ Found", length(significant_samples_gns), "significant genes with P.Value < 0.05 and score.\n")
   
  # --------------------- PREPARE MELTED SCORE DATAFRAMES -------------------
  cat("📦 Preparing melted score dataframes for plotting...\n")
  # For transcripts
  df_melted_trans <- prepare_score_df_for_plotting(
    score_matrix = df_score_trans,
    rbp = rbp,
    significant_ids = significant_samples,
    getBM = getBM
  )
  
  # For genes
  df_melted_genes <- prepare_score_df_for_plotting(
    score_matrix = df_score_gns,
    rbp = rbp,
    significant_ids = significant_samples_gns,
    getBM = getBM
  )
  
  # --------------------- STATISTICAL TEST FOR TRANSCRIPTS ---------------------
  cat("🧪 Performing statistical test (transcripts)...\n")
  cat("ℹ️  Comparing raw RBP scores between DE and non-DE transcripts using Mann-Whitney-Wilcoxon test.\n")
  
  stat.test.trans <- df_melted_trans %>%
    rstatix::wilcox_test(RawScore ~ DE_Status) %>%  # <--- use RawScore here
    rstatix::adjust_pvalue(method = "bonferroni") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x = "DE_Status")
  
  cat("✅ Transcript-level Wilcoxon test complete for RBP:", rbp, "\n")
  cat("   → Adjusted p-value:", stat.test.trans$p.adj, "| Significance:", stat.test.trans$p.adj.signif, "\n\n")
  
  # --------------------- STATISTICAL TEST FOR GENES --------------------------
  cat("🧪 Performing statistical test (genes)...\n")
  cat("ℹ️  Comparing raw RBP scores between DE and non-DE genes using Mann-Whitney-Wilcoxon test.\n")
  
  stat.test.genes <- df_melted_genes %>%
    rstatix::wilcox_test(RawScore ~ DE_Status) %>%  # <--- use RawScore here too
    rstatix::adjust_pvalue(method = "bonferroni") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x = "DE_Status")
  
  cat("✅ Gene-level Wilcoxon test complete for RBP:", rbp, "\n")
  cat("   → Adjusted p-value:", stat.test.genes$p.adj, "| Significance:", stat.test.genes$p.adj.signif, "\n\n")
  
  # --------------------- PLOT TRANSCRIPTS ---------------------
  cat("🎨 Plotting RBP scores for differentially expressed transcripts...\n")
  cat("   → RBP:", rbp, "\n")
  cat("   → Experiment:", experiment, "\n")
  plotlist_trans[[i]] <- generate_score_plot(df_melted_trans, 
                                             stat.test.trans, 
                                             rbp, experiment, 
                                             tipo = "transcripts")
  cat("✅ Transcript plot added to list at index", i, "\n\n")
  
  # --------------------- PLOT GENES --------------------------
  cat("🎨 Plotting RBP scores for genes with DE transcripts...\n")
  cat("   → RBP:", rbp, "\n")
  cat("   → Experiment:", experiment, "\n")
  plotlist_genes[[i]] <- generate_score_plot(df_melted_genes, 
                                             stat.test.genes, 
                                             rbp, experiment,
                                             tipo = "genes")
  cat("✅ Gene plot added to list at index", i, "\n\n")
  
}

# ----------------- PLOT ARRANGE: TRANSCRIPTS -----------------
cat("🧩 Arranging transcript-level plots...\n")
labels_trans <- c("A", "B", "C", "D", "E", "F")

arranged_trans <- ggpubr::ggarrange(
  plotlist = plotlist_trans,
  nrow = 3, ncol = 2,
  labels = labels_trans,
  common.legend = TRUE,
  legend = "right"
)

cat("✅ Transcript plots arranged.\n")

transcript_file <- file.path(output_dir, "transcripts_arranged.pdf")
ggsave(
  filename = transcript_file,
  plot = arranged_trans,
  height = 140, width = 100, units = "mm", scale = 2
)
cat("💾 Transcript figure saved to:", transcript_file, "\n\n")


# ----------------- PLOT ARRANGE: GENES -----------------
cat("🧩 Arranging gene-level plots...\n")
labels_genes <- c("C", "D", "E", "F", "G", "H")

arranged_genes <- ggpubr::ggarrange(
  plotlist = plotlist_genes,
  nrow = 3, ncol = 2,
  labels = labels_genes,
  common.legend = TRUE,
  legend = "right"
)

cat("✅ Gene plots arranged.\n")

genes_file <- file.path(output_dir, "genes_arranged.pdf")
ggsave(
  filename = genes_file,
  plot = arranged_genes,
  height = 140, width = 100, units = "mm", scale = 2
)
cat("💾 Gene figure saved to:", genes_file, "\n\n")


