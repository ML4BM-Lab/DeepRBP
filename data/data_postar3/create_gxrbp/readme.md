# Events x RBPs Matrix Creation in R

This R script is designed to generate an Events x RBP (Ribonucleoprotein) matrix for splicing factor enrichment analysis. The matrix helps in understanding the relationship between RNA binding proteins (RBPs) and splicing events across different tissues.

## Prerequisites:

### A:
- Data from POSTAR3:
  * Mouse
  * Human

### B:
- Information on splicing events.

## Steps:

### Load Packages and Set Arguments
- Libraries `readr` and `GenomicRanges` are loaded.
- Path to input data and tissue names are defined.

### Load Data
- **EventsFound_gencode23.txt**: This file contains information about splicing events such as position, event type, name, ID, etc.
- **Events_Regions_gc23_400nt.RData**: This file contains genomic regions of the splicing events.

### Process Data
- Iterate through the defined tissues (Liver, Myeloid, Kidney_embryo).
- Load POSTAR data for each tissue and preprocess it:
  * Remove index from Python.
  * Set column names.
  * Remove the first row with column names.

### Create POSTAR Dataframe for Each RBP
- Split the POSTAR dataframe based on the RBP column.
- Identify unique RBPs and exclude specific RBPs like "AGO2MNASE", "HURMNASE", and "PAN-ELAVL".

### Generate ExS Matrix
- Initialize an Events x RBPs matrix (`ExS`).
- For each RBP:
  * Filter genomic peaks based on chromosome (autosomes and sex chromosomes).
  * Convert peaks into a `GRanges` object.
  * Find overlaps between `POSTAR` and `Event Regions`.
  * Populate `ExS` matrix with 1s where overlaps occur.

### Generate GxS Matrix
- Duplicate the `ExS` matrix to create `GxS`.
- Group by Gene_ID to aggregate values.
- Set values greater than 1 to 1.
- Write the `GxS` matrix to a CSV file for each tissue.

### Cleanup
- Remove unnecessary variables from memory.

## Output:

The script will produce a CSV file for each tissue in the `./results` directory containing the GxS matrix. The GxS matrix represents the relationship between genes (rows) and RBPs (columns), with 1s indicating a presence of overlap between the gene and RBP in splicing events.

Ensure that the input data files (`EventsFound_gencode23.txt` and `Events_Regions_gc23_400nt.RData`) are correctly formatted and located in the specified path before running the script.
