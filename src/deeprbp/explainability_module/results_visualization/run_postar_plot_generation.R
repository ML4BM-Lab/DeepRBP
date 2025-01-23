# Load necessary libraries
rm(list = ls())

chooseCRANmirror(graphics = FALSE, ind = 1)  # Selecciona el primer espejo disponible
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("readr")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("rstatix")
#install.packages("latex2exp")
#install.packages("dplyr")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(latex2exp)
library(dplyr)

# Load necessary libraries
chooseCRANmirror(graphics = FALSE, ind = 1)  # Selecciona el primer espejo disponible
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("readr")
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(readr)
library(GenomicRanges)


# Load the create_postar_plots function
source("postar_plot_generator.R")

# Set paths and filenames
input_path <- '/Users/joseba/Desktop'
output_path <- '/Users/joseba/Desktop'
output_filename <- 'plot_score_results.pdf'
results_filename <- 'df_results_summary.csv'
list_rbps_postar_filename <- 'list_rbps_postar_ordered.csv'  
list_genes_postar_filename <- 'list_genes_postar_ordered.csv'  
getBM_filename <- 'getBM.csv'
save_plot <- TRUE
index_start <- 1
index_end <- 4
max_iterations <- 5

# Get command-line arguments if needed
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables if running from the command line
if (length(args) > 0) {
  input_path <- args[1]
  results_filename <- args[2]
  list_rbps_postar_filename <- args[3]
  list_genes_postar_filename <- args[4]
  getBM_filename <- args[5]
  output_path <- args[6]
  output_filename <- args[7]
  save_plot <- as.logical(args[8])  
  index_start <- args[9]
  index_end <- args[10]
  max_iterations <- args[11]
} 

# Call the function with the main arguments
plotlist <- create_postar_plots(
      input_path = input_path, 
      results_filename = results_filename, 
      list_rbps_postar_filename = list_rbps_postar_filename, 
      list_genes_postar_filename = list_genes_postar_filename, 
      getBM_filename = getBM_filename, 
      output_path = output_path, 
      output_filename = output_filename,
      save_plot = save_plot,
      index_start = index_start,  # Default value for index_start
      index_end = index_end,    # Default value for index_end
      max_iterations = max_iterations
)


                