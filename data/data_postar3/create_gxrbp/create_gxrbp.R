# Para crear la ExS necesitamos:
# 
# A:
#   - datos de POSTAR3.
#     * raton
#     * humano

# B
#   - Información de los eventos de splicing.

#' Events X RBPS matrix creation
#'
#' Generates the Events x RBP matrix for the splicing factor enrichment analysis.
#'

#### 1) Load packages and write args
library(readr)
library(GenomicRanges)
path = './../input_data'
tissues <- c("Liver", "Myeloid", "Kidney_embryo")

#### 2) Load Data
# 2.1) <- EventsFound_gencode23.txt: info de los eventos: posici�n, tipo de evento, nombre, id, etc...
# Path to eventsFound.txt with the information of all the events
pathtoeventstable <- paste0(path,"/EventsFound_gencode23.txt")

EventsFound <- read.delim(file=pathtoeventstable, stringsAsFactors=FALSE)
EventsFound$EventID <- paste0(EventsFound$GeneName, "_", EventsFound$EventNumber)

# 2.2) <- EventsRegions <- : las regiones de los eventos
load(paste0(path,"/Events_Regions_gc23_400nt.RData")) #GRseq3

for(tissue in tissues){
  print(tissue)
  
  # 2.4) <- Postar3 of the selected tissue: information of the RBPs that have been attatched in some positions of the genome 
  # Table with peaks of POSTAR
  human_txt <- read_delim(paste0(path,"/results/human_",tissue, ".txt"), 
                          delim = ",", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)
  
    # 2.4.1) Remove index from Python
  human_txt <- human_txt[,-1]
    # 2.4.2) Set colnames
  colnames(human_txt) <- human_txt[1,]
    # 2.4.3) Remove first row with the name of the columns
  human_txt <- human_txt[-1,]
  # dim(human_txt)
  # head(human_txt)
  POSTAR <- as.data.frame(human_txt) ##esto requiero un tipo largo y ocupa un total de 5.5 Gb
  
  # -) Type of events
  #typeA <- c("Retained Intron","Cassette Exon","Alternative 3' Splice Site","Alternative 5' Splice Site")
     
  # 3) Create a POSTAR dataframe for each of the RBPs
  POSTAR_L <- split.data.frame(x = POSTAR, f = POSTAR$RBP)
  print(paste0('Number of unique RBPs in ',tissue,': ',length(POSTAR_L)))
     
  # 4) List with the unique RBP names
  mySF <- unique(c(names(POSTAR_L)))
  mySF <- mySF[!mySF %in% c("AGO2MNASE", "HURMNASE", "PAN-ELAVL")]
  nSF <- length(mySF)
     
  # 5) Create the ExS matrix
  ExS <- matrix(0, nrow = nrow(EventsFound), ncol = nSF) 
  rownames(ExS) <- EventsFound$EventID
  colnames(ExS) <- mySF
     
  for(i in seq_len(nSF)){
    SF<- mySF[i]
    # a) index position of the actual SF
    jjx <- match(SF, names(POSTAR_L))
    # b) Take the Postar Info of the actual RBP
    peaks <- POSTAR_L[[jjx]]
    # c) filter and keep only genomic peaks (represented by rows in the data frame 
      # "peaks") found on autosomal chromosomes and human sex chromosomes. Remove chr different to "chr", 
      # c(1:22)), "chrX" and "chrY"
    iD <- which(!(peaks$seqname %in% c(paste0("chr", c(1:22)), "chrX", "chrY")))
    if(length(iD) > 0){
      peaks <- peaks[-iD, ]
    }
    # d) Convert the peaks object into a GRanges object so to perform later the olverlapping of regions
    peaks_GR <- GRanges(peaks)
    # e) Remove the "chr" prefix from each seqlevel which represent the chromosome identifiers
    seqlevels(peaks_GR) <- sub("chr", "", seqlevels(peaks_GR))
    # f) The peaks_GR object is reduced to a non-redundant set of genomic intervals
    peaks_GR<-reduce(peaks_GR)
    # g) identify overlaps between the intervals in the peaks_GR object (POSTAR) and the genomic intervals in 
    # the GRseq3 object (Event Regions)
    Overlaps <- findOverlaps(peaks_GR, GRseq3)
    # h) Event names where the RBP has overlapped
    EvMatch <- as.character(elementMetadata(GRseq3)$EventID[subjectHits(Overlaps)])
    if(length(EvMatch)>0){
      ExS[EvMatch,i] <- 1
    }
  }
 
  # 6) Create the GxS matrix
  GxS <- ExS
  # a) Drop the event number from the rownames so to get the Gene_ID 
  rownames(GxS) <- gsub("\\..*", "", rownames(GxS))
  # b) Aggrupate values for each Gene_ID
  GxS <- aggregate(GxS, by = list(rownames(GxS)), sum)
  rownames(GxS) <- GxS$Group.1
  GxS$Group.1 <- NULL
  # c) Values greater than 1 set to 1
  GxS[GxS>1] <- 1
  
  # d) Merge the information of the sex chr 
  GxS$gene_unique_id <- rownames(GxS)
  GxS$gene_unique_id <- gsub("R", "0", GxS$gene_unique_id)
  GxS <- aggregate(. ~ gene_unique_id, data = GxS, sum)
  rownames(GxS) <- GxS$gene_unique_id
  GxS$gene_unique_id <- NULL
  # Values greater than 1 set to 1
  GxS[GxS>1] <- 1
                     
  # e) Write the GxS matrix in a csv
  path_save <- ("./results")
  if (!dir.exists(path_save)) {
    dir.create(path_save, recursive = TRUE)
  }
  write.csv(GxS, file = paste0(path_save,'/',tissue,'_GxRBP.csv'))
  
  rm(human_txt, POSTAR, POSTAR_L, mySF, ExS, peaks, peaks_GR, Overlaps, EvMatch)
}


