
# ================================================================
# Differential Expression Analysis for Explainability Validation
# ================================================================
# In the explainability module validation section, we use real 
# knockdown data from two conditions: control and knockdown.
# 
# Before proceeding, we identify differentially expressed transcripts 
# and genes between these two conditions.
# 
# The following code demonstrates how to perform differential 
# expression analysis using the voom-limma methodology.
#
# We then loaded the expression of the transcripts for each sample. On the one hand we load for the control samples and on the other hand for the knockout samples:

# Clean workspace
rm(list = ls())

# Set CRAN mirror
chooseCRANmirror(graphics = FALSE, ind = 1)
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Function to install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}

# Required packages
required_packages <- c(
  "optparse",     # for parsing command-line arguments
  "limma",        # for voom-limma pipeline
  "edgeR",        # required by voom
  "data.table",   # fast reading of tables (optional but useful)
  "dplyr",        # data manipulation
  "tidyr",        # data wrangling
  "ggplot2" 
)

# Install and load
install_if_missing(required_packages)

library(optparse)
library(limma)
library(edgeR)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define the command-line options
# --------------------- PARSE ARGUMENTS --------------------------
cat("🔧 Parsing command-line arguments...\n")
option_list <- list(
  make_option(c("--path"), type = "character", default = "/scratch/jsanchoz/DeepRBP/data/explainability_module/real_kds/raw/GSE136366",
              help = "Path to directory containing kallisto output and sample info", metavar = "character"),
  make_option(c("--output_dir"), type = "character",
              help = "Directory where output results will be saved", metavar = "character"),
  make_option(c("--path_getBM"), type = "character", default = "/scratch/jsanchoz/DeepRBP/data/training_module/selected_genes_rbps/getBM.csv",
              help = "Path to getBM gene-to-transcript mapping file (CSV)", metavar = "character"),
  make_option(c("--condition_control"), type = "character", default = "Rescued_tdp43",
              help = "Label for the control condition in the experiment", metavar = "character"),
  make_option(c("--condition_kd"), type = "character", default = "tdp43_ko",
              help = "Label for the knockdown condition in the experiment", metavar = "character")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

path <- args$path
output_dir <- args$output_dir
path_getBM <- args$path_getBM
condition_control <- args$condition_control
condition_kd <- args$condition_kd

cat("📂 Input path:           ", path, "\n")
cat("📁 Output directory:     ", output_dir, "\n")
cat("📁 getBM file:           ", path_getBM, "\n")
cat("🧪 Condition (control):  ", condition_control, "\n")
cat("🧬 Condition (KD):       ", condition_kd, "\n")

# Create output directory if not exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --------------------- LOAD FUNCTION ----------------------------
cat("📄 Loading volcano plot function...\n")
setwd('/Users/joseba/Desktop/data:real_kds') # Update as needed
source(paste0(getwd(),"/plot_volcano.R")) 

# --------------------- LOAD SAMPLE INFO -------------------------
cat("📑 Reading sample information...\n")
info_sample <- read.table(file = file.path(path, "info_samples.txt"), sep = "\t")

# --------------------- CLEANING FUNCTION ------------------------
clean_abundance_df <- function(df, sample_name) {
  library(dplyr)
  library(tidyr)
  
  df_clean <- df %>%
    tidyr::separate(target_id, into = c("Transcript_ID", paste0("V", 2:8)), sep = "\\|", extra = "drop", remove = FALSE) %>%
    dplyr::select(Transcript_ID, est_counts) %>%
    dplyr::mutate(Transcript_ID = sub("\\..*", "", Transcript_ID)) %>%  # remove version
    dplyr::group_by(Transcript_ID) %>%
    dplyr::summarise(!!sample_name := sum(est_counts), .groups = "drop")
  
  #colnames(df_clean)[2] <- sample_name
  return(df_clean)
}

# --------------------- LOAD CONTROL -----------------------------
cat("📥 Loading control samples...\n")
index_control <- which(info_sample$V2 == condition_control)
sample_names_control <- info_sample$V1[index_control]

# Print number and names
cat("🔎 Found", length(sample_names_control), "control sample(s):\n")
print(sample_names_control)

# start with first sample
first_df <- read.csv(file.path(path, "kallisto_output", sample_names_control[1], "abundance.tsv"), sep = "\t")
data_control <- clean_abundance_df(first_df, sample_names_control[1])

# concatenate rest of the samples
for (jjx in 2:length(sample_names_control)) {
  df <- read.csv(file.path(path, "kallisto_output", sample_names_control[jjx], "abundance.tsv"), sep = "\t")
  df_clean <- clean_abundance_df(df, sample_names_control[jjx])
  data_control <- dplyr::left_join(data_control, df_clean, by = "Transcript_ID")
}
cat("✅ Control expression matrix assembled.\n")

# --------------------- LOAD KNOCKDOWN ---------------------------
cat("📥 Loading knockdown samples...\n")
index_kd <- which(info_sample$V2 == condition_kd)
sample_names_kd <- info_sample$V1[index_kd]

# Print number and names
cat("🔎 Found", length(sample_names_kd), "knockdown sample(s):\n")
print(sample_names_kd)

# Leer la primera muestra
first_df_kd <- read.csv(file.path(path, "kallisto_output", sample_names_kd[1], "abundance.tsv"), sep = "\t")
data_kd <- clean_abundance_df(first_df_kd, sample_names_kd[1])

# Concatenar el resto
for (jjx in 2:length(sample_names_kd)) {
  df <- read.csv(file.path(path, "kallisto_output", sample_names_kd[jjx], "abundance.tsv"), sep = "\t")
  df_clean <- clean_abundance_df(df, sample_names_kd[jjx])
  data_kd <- dplyr::left_join(data_kd, df_clean, by = "Transcript_ID")
}
cat("✅ Knockdown expression matrix assembled.\n")

# --------------------- BUILD EXPRESSION MATRIX ------------------
cat("🧬 Merging expression matrices...\n")
# Final matrix
transcripts_names <- data_control$Transcript_ID
expression_matrix_control <- as.matrix(data_control[, -1])
rownames(expression_matrix_control) <- transcripts_names
colnames(expression_matrix_control) <- sample_names_control

transcripts_names_kd <- data_kd$Transcript_ID
expression_matrix_kd <- as.matrix(data_kd[, -1])
rownames(expression_matrix_kd) <- transcripts_names_kd
colnames(expression_matrix_kd) <- sample_names_kd

Expression <- cbind(expression_matrix_control, expression_matrix_kd)
all_sample_names <- c(sample_names_control, sample_names_kd)
cat("✅ Expression matrix ready with", nrow(Expression), "transcripts and", ncol(Expression), "samples.\n")

# --------------------- DESIGN & CONTRAST ------------------------
cat("📐 Creating design and contrast matrices...\n")
X <- model.matrix(~ 0 + factor(c(rep(c("control"), each=length(sample_names_control)), 
                                 rep(c("kd"), each=length(sample_names_kd)))))
colnames(X) <- c("control", "kd")
C <- makeContrasts(control-kd, levels = X)

# --------------------- VOOM + LIMMA -----------------------------
cat("🧪 Running voom transformation and linear modeling...\n")
# Apply voom transformation
y <- voom(Expression, X)

# voom() transforms your raw count matrix into log2-counts with precision weights.
# It estimates the mean-variance relationship, allowing the use of linear models on RNA-seq data.
# Internally, it performs a transformation similar to log2(counts + small_number), but weighted.
# y$E contains the transformed expression matrix

# Fit linear model
fit <- lmFit(y, X) # Fits a linear model for each gene (or transcript) using the design matrix X.

# Apply the contrast
fit2 <- contrasts.fit(fit, C) # Applies the contrast matrix C to the fitted model. In this case: compares control - kd.
 
# Empirical Bayesian moderation
fit2<- eBayes(fit2)
# Applies empirical Bayes shrinkage to the standard errors.
# This improves variance estimation, especially with small sample sizes.
cat("✅ Differential expression analysis complete.\n")

# --------------------- EXTRACT RESULTS --------------------------
cat("📊 Extracting top differentially expressed transcripts...\n")
tT1 <- topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=Inf)
tT1$Transcript_ID <- rownames(tT1)
# coef = 1 → selects the first contrast (control - kd)
# adjust = "fdr" → applies multiple testing correction (Benjamini-Hochberg)
# sort.by = "B" → sorts by B-statistic (log-odds of being differentially expressed)
# number = Inf → returns all genes (no filtering at this step)

# --------------------- ANNOTATION -------------------------------
cat("🔎 Reading transcript annotation (getBM)...\n")
# We grouped the results obtained by gene. As explained in the manuscript, we will consider genes that have at least one differentially expressed transcript. 
getBM <- read.table(path_getBM, header = TRUE, sep = ",")

cat("🔄 Merging DE results with annotations...\n")
annotated_tT1 <- merge(tT1, getBM, by = "Transcript_ID")

# Check if the merge yielded results
if (nrow(annotated_tT1) == 0) {
  stop("❌ No matches found between DE transcripts (tT1) and annotation (getBM). Check identifiers.")
}
cat("✅", nrow(annotated_tT1), "transcripts matched with annotation.\n")

# --------------------- GENE-LEVEL SUMMARY -----------------------
cat("📚 Summarizing DE results at the gene level...\n")

# For each gene, keep the transcript with the lowest adj.P.Val
gT1 <- annotated_tT1 %>%
  group_by(Gene_ID) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice(1) %>%  # select most significant transcript per gene
  ungroup() %>%
  select(Gene_ID, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Convert to data frame and set rownames
gT1 <- as.data.frame(gT1)
rownames(gT1) <- gT1$Gene_ID
gT1$Gene_ID <- NULL
cat("✅ Gene-level summary complete with", nrow(gT1), "genes.\n")

tT1 <- tT1[,1:6]

# --------------------- EXPORT RESULTS ---------------------------
cat("💾 Saving results to disk...\n")
# The results are stored in the folder: "voom_limma_output".
# Save DE transcripts
tx_csv <- file.path(
    output_dir, "voom_limma_output",
    paste0("DE_transcripts_", condition_kd, "_vs_", condition_control, ".csv"))
dir.create(dirname(tx_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(tT1, file = tx_csv, row.names = TRUE)

# Save DE genes
gx_csv <- file.path(
    output_dir, "voom_limma_output", 
    paste0("DE_genes_", condition_kd, "_vs_", condition_control, ".csv"))
write.csv(gT1, file = gx_csv, row.names = TRUE)

cat("✅ DE tables saved:\n")
cat("   • Transcripts: ", tx_csv, "\n")
cat("   • Genes:       ", gx_csv, "\n")

# --------------------- VOLCANO PLOT -----------------------------
cat("📈 Generating volcano plot...\n")
plot_volcano(
  df = tT1,
  id_col = "Transcript_ID",
  getBM = getBM,
  title = "Transcripts_DE_KD_vs_Control",
  output_dir = file.path(output_dir, "voom_limma_output"),
  save_plot = TRUE
)

cat("🎉 All done! Volcano plot and DE results are ready.\n")


 

