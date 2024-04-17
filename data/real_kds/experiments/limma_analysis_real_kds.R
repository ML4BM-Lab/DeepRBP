
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("Glimma")
install.packages("Mus.musculus")

library(BiocManager)
library(limma)
library(Glimma)
library(edgeR)

# Data packaging
# Reading in count-data

path <- '/Users/joseba/Downloads'
experiment <- 'GSE77702-kd3'   #'SRR296'  #'HFE' # 'AGS'
rbp_interest <- 'TARDBP'   #'RBM47'  #'MBNL1' # 'ESRP2'
setwd(path)

# haz con sleuth: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html


########################
getBM <- read.table(paste0(path, '/getBM_total.csv'), header = TRUE, sep = ",")

# 1) Read data counts
data_control <- read.table(paste0(path,'/',experiment,'_control/',experiment,'_control_trans_est_counts.csv'), header = TRUE, sep = ",")
data_kd <- read.table(paste0(path,'/',experiment,'_', rbp_interest, '_kd/',experiment,'_', rbp_interest, '_kd_trans_est_counts.csv'), header = TRUE, sep = ",")

# Extrae los nombres de las columnas de muestras para ambos grupos
sample_names_control <- colnames(data_control)[-1]
sample_names_kd <- colnames(data_kd)[-1]

# Extrae las matrices de expresión para ambos grupos
expression_matrix_control <- as.matrix(data_control[, -1])
expression_matrix_kd <- as.matrix(data_kd[, -1])

# Combina los nombres de las muestras de control y kd
all_sample_names <- c(sample_names_control, sample_names_kd)
# Design Matrix: Crea la matriz de diseño con una sola columna indicando el grupo (control o kd)
# X <- model.matrix(~ 0+factor(rep(c("control", "kd"), each = length(all_sample_names)/2)))
X <- model.matrix(~ 0 + factor(c(rep(c("control"), each=length(sample_names_control)), rep(c("kd"), each=length(sample_names_kd)))))
colnames(X) <- c("control", "kd")

C <- makeContrasts(control-kd, levels = X)


# Usa solo la primera columna de la matriz de diseño
# design <- design[, "control", drop = FALSE]

# Crea un objeto DGEList para ambos grupos
#dge <- DGEList(counts = cbind(expression_matrix_control, expression_matrix_kd), genes = data_control$sample, group = design)

# Normaliza los datos utilizando TMM para ambos grupos
#dge <- calcNormFactors(dge)

# Realiza el análisis diferencial para ambos grupos usando limma
Expression <- cbind(expression_matrix_control, expression_matrix_kd)
rownames(Expression) <- data_control$sample

y <- voom(Expression, X)
fit <- lmFit(y, X)

fit2 <- contrasts.fit(fit, C)
fit2<- eBayes(fit2)

tT1 <- topTable(fit2,coef=1, adjust="fdr", sort.by="B", number=Inf)
tT1$Transcript_ID <- rownames(tT1)

merged_df <- merge(tT1, getBM, by="Transcript_ID")
merged_df <- merged_df[,1:8]

# Ruta del archivo CSV
csv_file_path <- paste0(getwd(), "/limma_", experiment, "_", rbp_interest, "_results.csv")

# Guardar los datos en el archivo CSV
write.csv(merged_df, file = csv_file_path, row.names = TRUE)
cat("Los datos se han guardado correctamente en", csv_file_path, "\n")










