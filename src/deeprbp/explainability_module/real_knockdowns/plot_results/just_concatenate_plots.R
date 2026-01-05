# src/deeprbp/explainability_module/real_knockdowns/plot_results/just_concatenate_plots.R

# install.packages(c("cowplot","magick")) # si no los tienes
library(cowplot)
library(magick)   # opcional si quieres controlar DPI/lectura


# leer el csv
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
})

## ===================== PARÁMETROS =====================
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(readxl)
  library(tibble)
  library(stringr)
})

## ===================== PARÁMETROS (CLI) =====================

option_list <- list(
  make_option(
    c("--base_dir"), type = "character",
    default = "/Users/joseba/Desktop/output_real_knockdowns/output/results/real_knockdowns/de_limma",
    help = "Root folder where de_limma results are stored (usually .../de_limma[/logFCx_FDRy])",
    metavar = "character"
  ),
  make_option(
    c("--excel_path"), type = "character",
    default = "/Users/joseba/Downloads/Supplementary File 0.xlsx",
    help = "Path to the Excel file with transcript/gene annotations",
    metavar = "character"
  ),
  make_option(
    c("--p_cut_type"), type = "character",
    default = "fdr",
    help = "Type of p-value cut: 'fdr' (adj.P.Val) or 'pvalue' (P.Value) [default %default]",
    metavar = "character"
  ),
  make_option(
    c("--p_cut_value"), type = "numeric",
    default = 0.05,
    help = "Numeric threshold for the chosen p_cut_type [default %default]",
    metavar = "numeric"
  ),
  make_option(
    c("--logfc_thresh"), type = "numeric",
    default = 1,
    help = "Absolute log2 fold-change threshold used for significance [default %default]",
    metavar = "numeric"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

base_dir     <- opt$base_dir
excel_path   <- opt$excel_path
p_cut_type   <- match.arg(opt$p_cut_type, c("fdr","pvalue"))
p_cut_value  <- opt$p_cut_value
logfc_thresh <- opt$logfc_thresh

cat("📁 Base dir:      ", base_dir, "\n")
cat("📄 Excel path:    ", excel_path, "\n")
cat("📊 Cut type:      ", p_cut_type, "\n")
cat("📊 Cut value:     ", p_cut_value, "\n")
cat("📊 logFC thr:     |log2FC| >", logfc_thresh, "\n\n")

# Carpeta para dejar salidas del batch
out_dir <- file.path(base_dir, "batch_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ===================== HELPERS =====================
strip_version <- function(x) sub("\\..*$", "", as.character(x))

read_de_csv <- function(path, id_name) {
  stopifnot(file.exists(path))
  df <- read.csv(path, check.names = FALSE, row.names = 1)
  df %>%
    rownames_to_column(var = id_name)
}

# Extrae un "label" legible del nombre del CSV (p. ej., "tdp43_ko_vs_Rescued_tdp43" o "TAF15 KD_vs_Control")
label_from_file <- function(f) {
  b <- basename(f)
  b <- sub("^DE_(transcripts|genes)_", "", b, ignore.case = TRUE)
  b <- sub("\\.csv$", "", b, ignore.case = TRUE)
  b
}

# Devuelve tibble con paths a transcripts y genes por cada carpeta voom_limma_output
discover_experiments <- function(root) {
  out_dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  out_dirs <- out_dirs[grepl("voom_limma_output$", out_dirs)]
  tibble(
    voom_dir = out_dirs,
    transcripts = map_chr(out_dirs, ~{
      cands <- list.files(.x, pattern = "^DE[_ ]?transcripts.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
      if (length(cands) == 0) NA_character_ else cands[1]
    }),
    genes = map_chr(out_dirs, ~{
      cands <- list.files(.x, pattern = "^DE[_ ]?genes.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
      if (length(cands) == 0) NA_character_ else cands[1]
    })
  ) %>%
    mutate(
      dataset = basename(dirname(voom_dir)),                      # p. ej., GSE136366
      label_t = ifelse(is.na(transcripts), NA, label_from_file(transcripts)),
      label_g = ifelse(is.na(genes),       NA, label_from_file(genes)),
      label   = coalesce(label_t, label_g)                        # intentamos unificar etiqueta
    ) %>%
    filter(!is.na(transcripts) & !is.na(genes))
}

## ===================== LECTURA DE EXCEL =====================
xls_trans <- read_excel(excel_path, sheet = 3)   # columnas: Transcript_name, Transcript_ID
xls_genes <- read_excel(excel_path, sheet = 2)   # columnas: Gene_name, Gene_ID

stopifnot(all(c("Transcript_ID") %in% names(xls_trans)))
stopifnot(all(c("Gene_ID", "Gene_name") %in% names(xls_genes)))

ids_excel_base <- xls_trans$Transcript_ID %>%
  as.character() %>%
  strip_version() %>%
  unique()

genes_list <- xls_genes %>%
  transmute(
    Gene_ID_base = strip_version(Gene_ID),
    Gene_name    = as.character(Gene_name)
  ) %>%
  distinct()

cat("📄 Excel: ", length(ids_excel_base), " Transcript_ID únicos (sin versión)\n", sep = "")

## ===================== PROCESADO POR EXPERIMENTO =====================
process_one <- function(tr_path, ge_path, dataset, label) {
  
  # columna de p segun tipo
  p_col <- if (p_cut_type == "fdr") "adj.P.Val" else "P.Value"
  
  # Lee DE - transcripts
  tT1 <- read_de_csv(tr_path, "Transcript_ID") %>%
    mutate(
      Transcript_ID_base = strip_version(Transcript_ID),
      significant = (.data[[p_col]] < p_cut_value) & (abs(logFC) > logfc_thresh)
    )
  
  # Lee DE - genes
  tG1 <- read_de_csv(ge_path, "Gene_ID") %>%
    mutate(
      Gene_ID_base = strip_version(Gene_ID),
      significant  = (.data[[p_col]] < p_cut_value) & (abs(logFC) > logfc_thresh)
    )
  
  # Lee DE
  #tT1 <- read_de_csv(tr_path, "Transcript_ID") %>%
  #  mutate(
  #    Transcript_ID_base = strip_version(Transcript_ID),
  #    significant = (adj.P.Val < fdr_thresh) & (abs(logFC) > logfc_thresh)
  #  )
  
  #tG1 <- read_de_csv(ge_path, "Gene_ID") %>%
  #  mutate(
  #    Gene_ID_base = strip_version(Gene_ID),
  #    significant  = (adj.P.Val < fdr_thresh) & (abs(logFC) > logfc_thresh)
  #  )
  
  # Transcript-level
  tT1_in_list <- tT1 %>% filter(Transcript_ID_base %in% ids_excel_base)
  n_total_sig   <- sum(tT1$significant, na.rm = TRUE)
  n_in_list     <- nrow(tT1_in_list)
  n_in_list_sig <- sum(tT1_in_list$significant, na.rm = TRUE)
  
  # Gene-level
  overlap <- tG1 %>% inner_join(genes_list, by = "Gene_ID_base")
  n_total_in_de    <- nrow(tG1)
  n_total_sig_de   <- sum(tG1$significant, na.rm = TRUE)
  n_list_in_de     <- nrow(overlap)
  n_list_sig_in_de <- sum(overlap$significant, na.rm = TRUE)
  
  # Guardar overlaps (opcional pero útil)
  safe_label <- gsub("[^A-Za-z0-9._-]+", "_", paste(dataset, label, sep = "_"))
  out_tr_csv <- file.path(out_dir, paste0("overlap_transcripts__", safe_label, ".csv"))
  out_ge_csv <- file.path(out_dir, paste0("overlap_genes__",       safe_label, ".csv"))
  suppressWarnings(write.csv(tT1_in_list, out_tr_csv, row.names = FALSE))
  suppressWarnings(write.csv(overlap,      out_ge_csv, row.names = FALSE))
  
  tibble(
    dataset,
    contrast_label = label,
    transcripts_total_sig = n_total_sig,
    transcripts_in_list   = n_in_list,
    transcripts_in_list_sig = n_in_list_sig,
    genes_total           = n_total_in_de,
    genes_total_sig       = n_total_sig_de,
    genes_from_excel_in_DE = n_list_in_de,
    genes_from_excel_sig   = n_list_sig_in_de,
    overlap_transcripts_csv = out_tr_csv,
    overlap_genes_csv       = out_ge_csv
  )
}

## ===================== DISCOVER + RUN =====================
exp_tbl <- discover_experiments(base_dir)

if (nrow(exp_tbl) == 0) {
  stop("No se encontraron carpetas 'voom_limma_output' con pares DE_transcripts/DE_genes.")
}

batch_summary <- pmap_dfr(
  list(exp_tbl$transcripts, exp_tbl$genes, exp_tbl$dataset, exp_tbl$label),
  process_one
)

# Guarda resumen combinado
summary_csv <- file.path(out_dir, "batch_summary.csv")
write.csv(batch_summary, summary_csv, row.names = FALSE)

cat("\n✅ Hecho. Resumen batch en:\n - ", summary_csv,
    "\n📁 Overlaps por experimento en:\n - ", out_dir, "\n", sep = "")

print(batch_summary, n = nrow(batch_summary))



################################################################################################
# Rutas a tus PNG (en el orden A, B, C)
# --- Input: rutas en orden A..F ---
# Rutas a tus PNG (en el orden A..)
library(magick)

# --- INPUT (A.. en orden) ---
datasets <- c("PRJEB39343", "GSE75491", "GSE136366")
# datasets <- c("PRJEB39343", "GSE75491", "GSE77702-FUS", "GSE77702-TAF15", "GSE77702-TARDBP", "GSE136366")

pngs <- file.path(base_dir, datasets, "voom_limma_output", "volcano_Transcripts_DE_KD_vs_Control.png")

missing <- pngs[!file.exists(pngs)]
if (length(missing)) stop("No se encontraron los PNG:\n", paste(missing, collapse = "\n"))

# --- Carga y normalización ---
imgs <- image_read(pngs)

# Altura fija + canvas base
target_h <- 1600   # un poco más alto
target_w <- 2000   # más ancho por panel
row_gap  <- 60     # separador entre filas más visible
legend_pad <- 300  

# Escalado base
imgs <- image_scale(imgs, paste0("x", target_h))                                   # ajusta altura
imgs <- image_extent(
  imgs,
  paste0(target_w, "x", target_h),                         # canvas base
  gravity = "center",
  color = "white"
)

# --- Padding extra para paneles con leyenda (3=C, 6=F si existen) ---
n <- length(imgs)
legend_idx_all <- c(3, 6)                 # posibles paneles con barra de color
legend_idx <- legend_idx_all[legend_idx_all <= n]  # sólo los que existen

legend_pad <- 260                          # píxeles extra a la derecha
for (i in legend_idx) {
  imgs[i] <- image_extent(
    imgs[i],
    paste0(target_w + legend_pad, "x", target_h),
    gravity = "west",                      # ancla a la izquierda → añade espacio a la derecha
    color = "white"
  )
}

# --- Etiquetas A.. (arriba-izquierda) ---
labs <- LETTERS[1:n]
for (i in seq_along(imgs)) {
  imgs[i] <- image_annotate(
    imgs[i], labs[i],
    gravity = "northwest",
    size = 80, weight = 700,
    location = "+25+18", color = "black"
  )
}

# --- Construye layout (1 o 2 filas) ---
if (n <= 3) {
  # Solo una fila
  combo_core <- image_append(imgs, stack = FALSE)
} else {
  # Primera fila: hasta 3 paneles
  row1 <- image_append(imgs[1:3], stack = FALSE)
  # Segunda fila: del 4 al n
  row2 <- image_append(imgs[4:n], stack = FALSE)
  
  # Ensancha a mismo ancho (por si leyenda añadió asimetría)
  wmax <- max(image_info(row1)$width, image_info(row2)$width)
  row1 <- image_extent(
    row1,
    paste0(wmax, "x", image_info(row1)$height),
    gravity = "center", color = "white"
  )
  row2 <- image_extent(
    row2,
    paste0(wmax, "x", image_info(row2)$height),
    gravity = "center", color = "white"
  )
  
  # Inserta gap entre filas
  row_gap_img <- image_blank(width = wmax, height = row_gap, color = "white")
  combo_core <- image_append(c(row1, row_gap_img, row2), stack = TRUE)
}

# Marco blanco y guardar
combo <- image_border(combo_core, color = "white", geometry = "20x20")

out_dir <- base_dir
image_write(combo, path = file.path(out_dir, "Figure_SX_volcano.png"),
            format = "png", density = "400x400")
# image_write(combo, file.path(out_dir, "Figure_SX_volcano.pdf"), format = "pdf")
