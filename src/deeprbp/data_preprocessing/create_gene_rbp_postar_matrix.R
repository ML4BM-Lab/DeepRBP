# src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R

# Load necessary libraries
chooseCRANmirror(graphics = FALSE, ind = 1)  # Selecciona el primer espejo disponible
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("readr")
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(readr)
library(GenomicRanges)

create_gxrbp <- function(input_path, output_path, output_file_name, postar_file, events_regions_file, events_gencode_file, selected_tissue_cell_line, getBM_file) {
  
  # Construct file paths based on input arguments
  path_events <- file.path(input_path, events_gencode_file)
  path_regions <- file.path(input_path, events_regions_file)
  path_postar <- file.path(input_path, postar_file)
  path_getbm <- file.path(input_path, getBM_file)
  
  # Load the EventsFound file and add EventID column
  EventsFound <- read.delim(file = path_events, stringsAsFactors = FALSE)
  EventsFound$EventID <- paste0(EventsFound$GeneName, "_", EventsFound$EventNumber)
  
  # Load the regions file (GRseq3 object)
  load(path_regions) # GRseq3
  
  # Read the POSTAR file
  postar_txt <- read_delim(
    path_postar, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE
  )
    #n_max = 10000 # Remove this argument to read the entire file
  
  # Assign column names to the POSTAR data
  colnames(postar_txt) <- c(
    "seqname", "start", "end", "experiment_metadata", "strand",
    "RBP_name", "technique", "raw_tissue", "experiment_accession", "value1", "value2"
  )
  
  # Apply tissue filtering
  postar_txt <- postar_txt[postar_txt$raw_tissue %in% selected_tissue_cell_line, ]
  
  # Convert POSTAR data to a data frame
  POSTAR <- as.data.frame(postar_txt)
  
  # Create a list where each RBP is a separate data frame
  POSTAR_L <- split.data.frame(x = POSTAR, f = POSTAR$RBP_name)
  print(paste0("Number of unique RBPs: ", length(POSTAR_L)))
  
  # List of unique RBP names, excluding specific RBPs
  mySF <- unique(c(names(POSTAR_L)))
  mySF <- mySF[!mySF %in% c("AGO2MNASE", "HURMNASE", "PAN-ELAVL")]
  nSF <- length(mySF)
  
  # Initialize the ExRBP matrix
  ExRBP <- matrix(0, nrow = nrow(EventsFound), ncol = nSF)
  rownames(ExRBP) <- EventsFound$EventID
  colnames(ExRBP) <- mySF
  
  # Iterate through each RBP to populate the ExRBP matrix
  for (i in seq_len(nSF)) {
    SF <- mySF[i]
    
    # a) Get the index position of the current RBP
    jjx <- match(SF, names(POSTAR_L))
    
    # b) Retrieve the data for the current RBP
    peaks <- POSTAR_L[[jjx]]
    
    # c) Filter peaks to include only those on valid chromosomes
    iD <- which(!(peaks$seqname %in% c(paste0("chr", c(1:22)), "chrX", "chrY")))
    if (length(iD) > 0) {
      peaks <- peaks[-iD, ]
    }
    
    # d) Convert peaks to a GRanges object
    peaks_GR <- GRanges(peaks)
    
    # e) Remove "chr" prefix from chromosome identifiers
    seqlevels(peaks_GR) <- sub("chr", "", seqlevels(peaks_GR))
    
    # f) Reduce the GRanges object to a non-redundant set of genomic intervals
    peaks_GR <- reduce(peaks_GR)
    
    # g) Identify overlaps between POSTAR intervals and Event Regions (GRseq3)
    Overlaps <- findOverlaps(peaks_GR, GRseq3)
    
    # h) Get the EventIDs where the RBP overlaps
    EvMatch <- as.character(elementMetadata(GRseq3)$EventID[subjectHits(Overlaps)])
    if (length(EvMatch) > 0) {
      ExRBP[EvMatch, i] <- 1
    }
  }
  
  # Create the GxRBP matrix from ExRBP
  GxRBP <- ExRBP
  
  # a) Remove the event number from row names to get Gene_ID
  rownames(GxRBP) <- gsub("\\..*", "", rownames(GxRBP))
  
  # b) Aggregate values by Gene_ID
  GxRBP <- aggregate(GxRBP, by = list(rownames(GxRBP)), sum)
  rownames(GxRBP) <- GxRBP$Group.1
  GxRBP$Group.1 <- NULL
  
  # c) Set values greater than 1 to 1
  GxRBP[GxRBP > 1] <- 1
  
  # d) Merge the information for sex chromosomes
  GxRBP$gene_unique_id <- rownames(GxRBP)
  GxRBP$gene_unique_id <- gsub("R", "0", GxRBP$gene_unique_id)
  GxRBP <- aggregate(. ~ gene_unique_id, data = GxRBP, sum)
  rownames(GxRBP) <- GxRBP$gene_unique_id
  GxRBP$gene_unique_id <- NULL
  
  # Set values greater than 1 to 1 again
  GxRBP[GxRBP > 1] <- 1
  
  # e) Change gene name for identify RBPs for gene_id
  # Read the getBM file
  getBM <- read.csv(path_getbm, header = TRUE, stringsAsFactors = FALSE)
  
  # Create a mapping vector between Gene_name and Gene_ID
  gene_map <- setNames(getBM$Gene_ID, getBM$Gene_name)
  
  # Replace column names in the original dataframe
  colnames(GxRBP) <- ifelse(colnames(GxRBP) %in% names(gene_map),
                                   gene_map[colnames(GxRBP)],
                                   colnames(GxRBP))
  
  # Add axis names (rownames and colnames)
  attr(GxRBP, "row_axis_name") <- "Gene_ID"
  attr(GxRBP, "col_axis_name") <- "RBP_ID"
  
  # Verify the result
  attributes(GxRBP)

  # f) Write the GxRBP matrix to a CSV file
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  write.csv(GxRBP, file = file.path(output_path, paste0(output_file_name, "_GxRBP.csv")))
  #saveRDS(GxRBP, file = file.path(output_path, paste0(output_file_name, "_GxRBP.rds")))
  
  # Clean up the workspace
  rm(postar_txt, POSTAR, POSTAR_L, mySF, ExRBP, peaks, peaks_GR, Overlaps, EvMatch)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_path <- args[1]
output_path <- args[2]
output_file_name <- args[3]
postar_file <- args[4]
events_regions_file <- args[5]
events_gencode_file <- args[6]
selected_tissue_cell_line <- unlist(strsplit(args[7], ","))
getBM_file <- args[8]

# Call the function with the parsed arguments
create_gxrbp(input_path, output_path, output_file_name, postar_file, events_regions_file, events_gencode_file, selected_tissue_cell_line, getBM_file)
