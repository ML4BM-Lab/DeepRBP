#src/deeprbp/data_preprocessing/create_gene_rbp_postar_matrix.R

# Load necessary libraries
chooseCRANmirror(graphics = FALSE, ind = 1)  # Selecciona el primer espejo disponible
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("readr")
install.packages("BiocManager")
install.packages("optparse")
BiocManager::install("GenomicRanges")
library(readr)
library(GenomicRanges)
library(optparse)

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

# Define command-line options
option_list <- list(
  make_option(c("--input_path"), type = "character", help = "Directory containing the input files."),
  make_option(c("--output_path"), type = "character", help = "Directory where the processed output will be saved."),
  make_option(c("--output_file_name"), type = "character", help = "The name of the output file."),
  make_option(c("--postar_file"), type = "character", help = "The POSTAR file containing RBP binding information."),
  make_option(c("--events_regions_file"), type = "character", help = "File specifying the genomic regions of the events."),
  make_option(c("--events_gencode_file"), type = "character", help = "File with metadata on events, including IDs and positions."),
  make_option(c("--selected_tissue_cell_line"), type = "character", help = "Specifies the cell lines from POSTAR experiments to include in the matrix (comma-separated)."),
  make_option(c("--getBM_file"), type = "character", help = "Name of the getBM file with info to relate gene_name with gene_ids.")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign arguments to variables
input_path <- opt$input_path
output_path <- opt$output_path
output_file_name <- opt$output_file_name
postar_file <- opt$postar_file
events_regions_file <- opt$events_regions_file
events_gencode_file <- opt$events_gencode_file
selected_tissue_cell_line <- unlist(strsplit(opt$selected_tissue_cell_line, ","))
getBM_file <- opt$getBM_file

# Call the function with the parsed arguments
create_gxrbp(input_path, output_path, output_file_name, postar_file, events_regions_file, events_gencode_file, selected_tissue_cell_line, getBM_file)
