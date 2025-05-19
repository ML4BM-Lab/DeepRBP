# Load necessary libraries
rm(list = ls())

# Set a default CRAN mirror
chooseCRANmirror(graphics = FALSE, ind = 1)  # Selecciona el primer espejo disponible
options(repos = c(CRAN = "https://cloud.r-project.org/"))


packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

required_packages <- c("data.table", "dplyr",  "rstatix", "ggplot2", "ggpubr","optparse")

install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)
install.packages("dplyr")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("ggplot2")

library(data.table)
library(dplyr)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(optparse)
library(tidyr)

# Load the create_postar_plots function
source(paste0(getwd(),"/generate_postar_plots.R"))

# Define the command-line options
option_list <- list(
  make_option(c("--input_path"), type = "character", default = NULL, help = "Path to results summary file and additional postar files", metavar = "character"),
  make_option(c("--results_filename"), type = "character", default = "df_results_summary.csv", help = "Results file name (CSV)", metavar = "character"),
  make_option(c("--count_genes_per_rbp_file"), type = "character", default = "count_genes_per_rbp.csv", help = "Number of genes per RBP in Postar ordered by positive rbp-gene interactions(CSV) in decreasing order", metavar = "character"),
  make_option(c("--count_rbps_per_gen_file"), type = "character", default = "count_rbps_per_gen.csv", help = "Number of rbps per Gene in Postar ordered by positive rbp-gene interactions(CSV) in decreasion order", metavar = "character"),
  make_option(c("--getBM_path"), type = "character", default = NULL, help = "Path to the Gene ID mapping file", metavar = "character"),
  make_option(c("--getBM_filename"), type = "character", default = "getBM.csv", help = "Gene ID mapping file name (CSV)", metavar = "character"),
  make_option(c("--output_path"), type = "character", default = NULL, help = "Path to save output files", metavar = "character"),
  make_option(c("--output_filename"), type = "character", default = "plot_score_results.pdf", help = "Name of the output file (with extension)", metavar = "character"),
  make_option(c("--save_plot"), type = "logical", default = FALSE, help = "Whether to save the plot as a PDF", metavar = "logical"),
  make_option(c("--index_start"), type = "integer", default = 1, help = "Starting index for slicing the RBP and gene lists for plotting", metavar = "integer"),
  make_option(c("--index_end"), type = "integer", default = 4, help = "Ending index for slicing the RBP and gene lists for plotting", metavar = "integer"),
  make_option(c("--max_iterations"), type = "integer", default = 1, help = "Maximum number of iterations for processing RBPs and genes for plotting", metavar = "integer")
)

# Parse the command-line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assign arguments to variables
input_path <- args$input_path
results_filename <- args$results_filename
count_genes_per_rbp_file <- args$count_genes_per_rbp_file
count_rbps_per_gen_file <- args$count_rbps_per_gen_file
getBM_path <- args$getBM_path
getBM_filename <- args$getBM_filename
output_path <- args$output_path
output_filename <- args$output_filename
save_plot <- args$save_plot
index_start <- args$index_start
index_end <- args$index_end
max_iterations <- args$max_iterations

# Call the function with the main arguments
plotlist <- create_postar_plots(
      input_path = input_path, 
      results_filename = results_filename, 
      count_genes_per_rbp_file = count_genes_per_rbp_file,
      count_rbps_per_gen_file = count_rbps_per_gen_file,
      getBM_path = getBM_path,
      getBM_filename = getBM_filename, 
      output_path = output_path, 
      output_filename = output_filename,
      save_plot = save_plot,
      index_start = index_start,  
      index_end = index_end,    
      max_iterations = max_iterations
)
