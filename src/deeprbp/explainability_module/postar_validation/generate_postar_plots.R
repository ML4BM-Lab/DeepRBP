#' Create Visualizations for POSTAR Data Using Explainability Scores
#'
#' This function generates visualizations to compare scores derived from 
#' explainability techniques (such as DeepLIFT) against the POSTAR experimental data. 
#' It reads data from specified input files, processes the information, and produces 
#' box plots that illustrate the distribution of scores across different RNA Binding 
#' Proteins (RBPs) and genes, utilizing POSTAR labels for classification.
#' Optionally, the generated plots can be saved to a specified output file.
#'
#' @param input_path A character string specifying the path to the input files.
#' @param results_filename A character string indicating the name of the results file (CSV) 
#'                        that contains the calculated scores and postar classification.
#' @param count_genes_per_rbp_file A character string indicating the name of the file (CSV) 
#'                        that contains the number of genes per RBP in Postar ordered by positive rbp.
#'                        gene interactions in decreasing order.
#' @param count_rbps_per_gen_file A character string indicating the name of the file (CSV) 
#'                        that contains the number of rbps per Gene in Postar ordered by positive rbp.
#'                        gene interactions in decreasing order.
#' @param getBM_path A character string indicating the path of the getBM file (CSV) 
#'                       that contains gene ID mappings.
#' @param getBM_filename A character string indicating the name of the getBM file (CSV) 
#'                       that contains gene ID mappings.
#' @param output_path A character string specifying the path where the output files will be saved.
#' @param output_filename A character string indicating the name of the output file (with extension).
#' @param save_plot A logical value indicating whether to save the plot as a PDF. Default is FALSE.
#' @param index_start An integer specifying the starting index for slicing the RBP and gene lists for plotting. Default is 1.
#' @param index_end An integer specifying the ending index for slicing the RBP and gene lists  for plotting. Default is 4.
#' @param max_iterations An integer specifying the maximum number of iterations for processing RBPs and genes for plotting. Default is 1.
#'
#' @return A ggplot object containing the generated plots that visualize the comparison 
#'         of DeepLIFT scores to POSTAR labels across RBPs and genes. RBPs and genes are plotted in groups of five until 
#'         the entire list of RBPs has been displayed.
#'
#' @import data.table
#' @import ggpubr
#' @import dplyr
#' @import rstatix
#' @import tidyr
#'
#' @examples
#' create_postar_plots(input_path = "path/to/data", 
#'                     results_filename = "results.csv",
#'                     count_genes_per_rbp_file = "count_genes_per_rbp.csv",
#'                     count_rbps_per_gen_file = "count_rbps_per_gen.csv",
#'                     getBM_path = "path/to/getBM", 
#'                     getBM_filename = "getBM.csv", 
#'                     output_path = "path/to/output", 
#'                     output_filename = "plots", 
#'                     save_plot = TRUE,
#'                     index_start = 1,
#'                     index_end = 4,
#'                     max_iterations = 1)

create_postar_plots <- function(
    input_path, 
    results_filename, 
    count_genes_per_rbp_file,
    count_rbps_per_gen_file,
    getBM_path,
    getBM_filename,
    output_path,
    output_filename,
    save_plot = FALSE,
    index_start = 1,
    index_end = 4,
    max_iterations = 4
  ){
  
  # Load input data
  path_summary_results <- file.path(input_path, results_filename)
  path_count_genes_per_rbp <- file.path(input_path, count_genes_per_rbp_file)
  path_count_rbps_per_gen <- file.path(input_path, count_rbps_per_gen_file)
  path_getbm <- file.path(getBM_path, getBM_filename)
  
  df_results_summary <- data.table::fread(path_summary_results) %>% as.data.frame() 
  df_count_genes_per_rbp <- data.table::fread(path_count_genes_per_rbp) %>% as.data.frame()
  df_count_rbps_per_gen <- data.table::fread(path_count_rbps_per_gen) %>% as.data.frame()
  getBM <- data.table::fread(path_getbm) %>% as.data.frame()
  
  list_rbps_postar <- df_count_genes_per_rbp$RBPs
  list_genes_postar <- df_count_rbps_per_gen$Genes
  
  # Calculate the absolute scores
  df_results_summary$Score <- abs(df_results_summary$Score)
  
  # Change the list name of rbps from gene_id to gene_name
  indices <- match(list_rbps_postar, getBM$Gene_ID)
  new_names <- getBM$Gene_name[indices]
  list_rbps_postar <- new_names
  
  # Change dtype of Postar column to character and factorize
  #- **Binding**: Indicates that the gene is regulated by the RBP (Postar_Score = 1).
  #- **Not Binding**: Indicates that the gene is not regulated by the RBP (Postar_Score = 0).
  #- **Unknown**: Indicates that there is insufficient data to determine the binding status (Postar_Score = N/A).
  
  df_results_summary$Postar_Score <- as.character(df_results_summary$Postar_Score)
  df_results_summary$Postar_Score[is.na(df_results_summary$Postar_Score)] <- "Unknown"
  df_results_summary <- df_results_summary %>% 
    mutate(Postar_Score = factor(Postar_Score, 
                                 levels = c("1","0", "Unknown"), 
                                 labels = c("Binding","Not Binding","Unknown")))
  head(df_results_summary)
  
  # Initialize an empty list to hold all plot figures
  all_plot_figures <- list()
  label_counter <- 0  # Counter for labels (starting from 0)
  
  # Filter dataframes for RBPs and Genes
  df_results_summary1 <- df_results_summary %>% 
    filter(RBP_name %in% unique(list_rbps_postar)) %>% 
    group_by(RBP_ID) %>%
    filter(all(c("Not Binding", "Binding") %in% Postar_Score)) %>%  # Keep only RBP_ID with both values
    ungroup() %>%  # Ungroup to avoid issues in subsequent steps
    mutate(RBP_name = factor(RBP_name, levels = unique(list_rbps_postar)))
  
  # For genes just take first 100 genes with more RBP 1s. 
  df_results_summary2 <- df_results_summary %>%
    filter(Gene_ID %in% list_genes_postar) %>%  # Filter for the genes of interest
    group_by(Gene_ID) %>%
    filter(sum(Postar_Score == "Binding") >= 5,  # Ensure at least 5 in group 1
           sum(Postar_Score == "Not Binding") >= 5) %>%  # Ensure at least 5 in group 0
    ungroup() %>%  # Ungroup to avoid issues in subsequent steps
    filter(Gene_ID %in% list_genes_postar[index_start:100]) %>%  # Filter by genes in the specified range
    mutate(Gene_ID = factor(Gene_ID, levels = list_genes_postar[index_start:100]))  # Ensure Gene_ID is a factor

  # Calculate p-values onto the bar plots across RBPs
  # Add p-values onto the bar plots
  stat.test_rbps <- df_results_summary1 %>%
    group_by(RBP_name) %>%
    rstatix::wilcox_test(Score ~ Postar_Score, comparisons = list(c("Binding", "Not Binding"))) %>%
    add_xy_position(fun = "max", x = "RBP_name") %>% 
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  stat.test_genes <- df_results_summary2 %>%
    group_by(Gene_ID) %>%
    rstatix::wilcox_test(Score ~ Postar_Score, comparisons = list(c("Binding", "Not Binding"))) %>%
    add_xy_position(fun = "max", x = "Gene_ID") %>% 
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")  
  #%>% filter(Gene_ID %in% actual_genes) 
  
  # Save the statistical test results as CSV files
  write.csv(stat.test_rbps %>% unnest(cols = c(groups)), file = file.path(output_path, paste0("stat_test_rbps.csv")), row.names = FALSE)
  write.csv(stat.test_genes %>% unnest(cols = c(groups)), file = file.path(output_path, paste0("stat_test_genes.csv")), row.names = FALSE)
  
  stat.test_rbps <- stat.test_rbps %>%
    mutate(
      adjustment = (row_number() - 1) %/% 4 * 4,  
      x = x - adjustment,                          
      xmin = xmin - adjustment,                    
      xmax = xmax - adjustment                    
    ) %>%
    select(-adjustment) 
  
  stat.test_genes <- stat.test_genes %>%
    mutate(
      adjustment = (row_number() - 1) %/% 4 * 4,  
      x = x - adjustment,                          
      xmin = xmin - adjustment,                    
      xmax = xmax - adjustment                    
    ) %>%
    select(-adjustment)
  
  stat.test_genes <- merge(stat.test_genes, 
                           df_results_summary2 %>% filter(Postar_Score %in% c("Binding", "Not Binding")) %>% 
                             group_by(Gene_ID) %>% summarise(y.position.2 = max(Score)+2))
  
  while (label_counter < max_iterations) {
    # Increment the label counter for each iteration
    label_counter <- label_counter + 1
    
    # Use index_start and index_end to slice the lists
    actual_rbps <- list_rbps_postar[index_start:index_end]
    actual_rbps2 <- paste0("                  ", actual_rbps)
    actual_genes <- list_genes_postar[index_start:index_end]
    
    # Filter dataframes for RBPs and Genes
    actual_df_results_summary1 <- df_results_summary %>% 
      filter(RBP_name %in% actual_rbps) %>% 
      mutate(RBP_name = factor(RBP_name, levels = actual_rbps, labels = actual_rbps2))
  
    actual_df_results_summary2 <- df_results_summary %>% 
      filter(Gene_ID %in% actual_genes) %>% 
      mutate(Gene_ID = factor(Gene_ID, levels = actual_genes))
    
    # Filter stat.test_rbps to include only RBPs and Genes present in actual_df_results_summary1
    filtered_stat_rbps <- stat.test_rbps %>%
      filter(RBP_name %in% actual_rbps)
    
    filtered_stat_rbps <- filtered_stat_rbps %>%
      mutate(y.position = y.position - 0.5)
    
    filtered_stat_genes <- stat.test_genes %>%
      filter(Gene_ID %in% actual_genes)
    
    # Plotting
    plotlist <- list()
    
    plotlist[[1]] <- ggboxplot(
      actual_df_results_summary1, 
      notch = FALSE,
      x = "RBP_name", 
      y = "Score", 
      fill = "Postar_Score") + xlab("RBP name") + 
      scale_fill_manual(
        values = c("#7fc97f", "#beaed4"),  
        breaks = c("Binding", "Not Binding"),  
        labels = c("Binding", "Not Binding")     
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      scale_y_continuous(expand = c(0,0.25)) + 
      stat_pvalue_manual(filtered_stat_rbps,  label = "p.adj.signif") +
      guides(fill = "none") 
    
    plotlist[[2]] <- ggboxplot(
      actual_df_results_summary2,
      x = "Gene_ID", 
      y = "Score", 
      fill = "Postar_Score") + xlab("Gene ID") + 
      scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086"),  
                        breaks = c("Binding", "Not Binding", "Unknown"),
                        labels = c("Binding", "Not Binding", "Unknown")) +  # Labels for the legend) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      scale_y_continuous(expand = c(0,0.25)) +
      stat_pvalue_manual(
        filtered_stat_genes,  label = "p.adj.signif", 
        y.position = "y.position.2")
    
    # Create a vector to store labels dynamically
    labels <- letters[(label_counter - 1) * 2 + 1:2]  # Create labels A, B, C, D, etc.
    
    # Concatenate the figures
    plot_figure <- ggarrange(plotlist = plotlist, nrow = 1, ncol = 2, 
                             labels = labels, #labels = c("A", "B"), 
                             common.legend = T, legend = "right")
    
    # Add the plot figure to the all_plot_figures list
    all_plot_figures[[label_counter]] <- plot_figure
    
    # Save each figure as a separate PDF
    pdf(file = file.path(output_path, paste0(output_filename, "_", label_counter, ".pdf")), width = 8, height = 6)
    print(plot_figure)
    dev.off()
    
    # Update indices for the next iteration
    if (index_end < length(list_rbps_postar)) {
      # Increment the indices for the next batch
      index_start <- index_start + 4
      index_end <- index_end + 4
    } else {
      break  # Exit the loop if there are no more RBPs to process
    }
  }
  return(all_plot_figures)
}

