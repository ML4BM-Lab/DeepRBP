 
import os
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
from sanbomics.tools import id_map
from sanbomics.plots import volcano
from config_loader import ConfigParser
from processing import DataImporter, DatasetLoader
from utils import ensure_directory_exists

class DifferentialExpressionAnalysis:
    def __init__(self, config_path, verbose=1):
        self.logger = Logger(verbose)
        self.config_path = config_path
        self.base_config = self.load_config()
        self.path_save_results = os.path.join(self.base_config['output_dir'])
        ensure_directory_exists(self.path_save_results)
        self.data_loader = self.initialize_data_loader()
        self.results_df = None

    def load_config(self):
        self.logger.log("Loading configuration.", level=1)
        config_parser = ConfigParser(self.config_path)
        return config_parser.get_base_config()

    def initialize_data_loader(self):
        self.logger.log("Initializing data loader.", level=1)
        data_importer = DataImporter(self.base_config['data_paths'])
        return DatasetLoader(data_importer, self.base_config)

    def process_data(self):
        self.logger.log("Processing data.", level=1)
        # Load data
        data = self.data_loader.load_data()
        counts_df = data['rbp_counts_df']

        # Transform counts
        #counts_df = np.power(2, counts_df) - 1
        #counts_df = counts_df.round().astype(int)
        # aqui poner raw counts mejor # aqui poner raw counts mejor # aqui poner raw counts mejor 
    
        # Prepare metadata
        metadata = data['metadata_df'].copy()
        metadata.rename(columns={'sample_type': 'condition'}, inplace=True)
        metadata['condition'] = metadata['condition'].str.replace('_', ' ', regex=False)

        # Order and filter metadata and counts
        condition_order = {'Primary Tumor': 1, 'Solid Tissue Normal': 2}
        metadata['order'] = metadata['condition'].map(condition_order)
        metadata = metadata.sort_values(by='order')
        counts_df = counts_df.loc[metadata.index]
        metadata = metadata.drop(columns=['order'])

        # Filter samples
        samples_to_keep = metadata['condition'].isin(['Primary Tumor', 'Solid Tissue Normal'])
        counts_df = counts_df.loc[samples_to_keep]
        metadata = metadata.loc[samples_to_keep]

        # Filter genes
        genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
        counts_df = counts_df[genes_to_keep]
        self.logger.log("Data processing completed.", level=1)
        return counts_df, metadata
    
    def run_analysis(self):
        self.logger.log("Running differential expression analysis.", level=1)
        counts_df, metadata = self.process_data()

        # DESeq2 setup
        inference = DefaultInference(n_cpus=8)
        ref_level = ["condition", "Solid Tissue Normal"]

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            inference=inference,
            ref_level=ref_level
        )

        # Run DESeq2
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["condition", "Primary Tumor", "Solid Tissue Normal"], inference=inference)
        self.results_df = stat_res.results_df.copy()

        # Sort results
        self.results_df = self.results_df.sort_values(by='stat', ascending=False)
        self.logger.log("Differential expression analysis completed.", level=1)

    def save_results(self):
        if self.results_df is not None:
            results_file_path = os.path.join(self.path_save_results, "results.csv")
            self.results_df.to_csv(results_file_path)
            self.logger.log(f"Results saved to {results_file_path}", level=1)
            print(f"Results saved to {results_file_path}")
        else:
            self.logger.warn("No results to save.", level=1)
        
    def generate_plots(self, logfc_threshold=2, pval_threshold=0.05):
        self.logger.log("Generating plots.", level=1)
        gene_id_to_name = pd.read_csv(self.base_config['data_paths']['getBM_path'])[['Gene_ID', 'Gene_name']].drop_duplicates()
        mapper = id_map(species='human')

        # Filter and order gene names
        filtered_gene_names = gene_id_to_name[gene_id_to_name.Gene_ID.isin(self.results_df.index.tolist())]
        ordered_gene_names = filtered_gene_names.set_index('Gene_ID').reindex(self.results_df.index).reset_index()
        
        self.results_df['Gene_name'] = ordered_gene_names['Gene_name'].tolist()
        self.results_df['Symbol'] = self.results_df.index.map(mapper.mapper)

        # Volcano plot
        volcano(self.results_df, symbol='Gene_name')
        plt.savefig(os.path.join(self.path_save_results, 'volcano_plot.png'), dpi=300)
        plt.close()

        # Filter significant results
        sign_df = self.results_df[(self.results_df.padj < pval_threshold) & 
                                  (self.results_df.log2FoldChange.abs() > logfc_threshold)]

        # Clustermap for significant results
        dds_sigs = dds[:, sign_df.index]
        grapher = pd.DataFrame(dds_sigs.layers['log1p'].T, index=dds_sigs.var_names, columns=dds_sigs.obs_names)
        sns.clustermap(grapher, z_score=0, cmap='RdYlBu_r', figsize=(8, 6))
        plt.savefig(os.path.join(self.path_save_results, 'clustermap_significant.png'), dpi=300)
        self.logger.log("Plots generated successfully.", level=1)
        plt.close()

# Uso de la clase
# rbps_of_interest = ['SLU7', 'SRFBP1', 'SRSF1', 'SRSF2', 'SRSF3', 'UPF1', 'IGF2BP1', 'IGF2']
config_path = "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_deg.yaml"
dea = DifferentialExpressionAnalysis(config_path)
dea.run_analysis()
dea.save_results()
dea.generate_plots()

















############################################################################################
# based on: https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_05_dge.html

# Differential expression is performed with the function rank_genes_group. The default method to compute differential 
# expression is the t-test_overestim_var. Other implemented methods are: logreg, t-test and wilcoxon.

# By default, the .raw attribute of AnnData is used in case it has been initialized, it can be changed by setting use_raw=False.

# As you can see, the X matrix only contains the variable genes, while the raw matrix contains all genes.

# Printing a few of the values in adata.raw.X shows that the raw matrix is normalized.

# For DGE analysis we would like to run with all genes, on normalized values, so we will have to revert back to the raw matrix. 
# In case you have raw counts in the matrix you also have to renormalize and logtransform.

# steps:
# 1) sacar las cuentas de TCGA (download RSEM expected_count) -> log2(expected_count+1)
        # https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz; Full metadata

# -	De la sección ‘Identifying potential RBP-Gene candidates’, hacer:
# o	1) DEG Control vs Cancer (all liver TCGA – see JIND).
# o	2) Test data     Expl. Control
#           Expl. Cancer Live
# o	Ver si para los DEG para los que pasan el theshold (se puede considerar como un 1) cambian en control vs liver cáncer. Y luego comparar con Postar. Y Luego mirar los RBPs de interés del CIMA.


# necesito acceder a all liver tcga (y coger DEG Control vs Cancer (necesito las cuentas?))

# analizando esto: 
# https://github.com/mousepixels/sanbomics_scripts/blob/main/PyDeseq2_DE_tutorial.ipynb

import matplotlib.pyplot as plt
from sanbomics.tools import id_map
from sanbomics.plots import volcano

from config_loader import ConfigParser
from processing import DataImporter, DatasetLoader
from utils import ensure_directory_exists

config_path = "/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_deg.yaml"

#class DifferentialExpressionAnalysis:
#def __init__(self, config_path, external_config_path=None):
        
config_parser = ConfigParser(config_path)
base_config = config_parser.get_base_config()
getBM = pd.read_csv(base_config['data_paths']['getBM_path'])

# Define paths for saving data and results
path_save_results = os.path.join(base_config['output_dir'])
ensure_directory_exists(path_save_results)

# Initialize DataImporter and DatasetLoader for raw data
data_importer = DataImporter(base_config['data_paths'])
data_loader = DatasetLoader(data_importer, base_config) 
data = data_loader.load_data()

# Data Loading 
# To perform differential expression analysis of genes (DEG), PyDESeq2 requires two types of inputs:
# A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
counts_df = data['rbp_counts_log2p_df']

############### steps que meter en el generate_input_data:
counts_df = np.power(2, counts_df) - 1
counts_df  = counts_df.round().astype(int)
##############

# Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.
metadata = data['metadata_df'].copy() 
metadata.rename(columns={'sample_type': 'condition'}, inplace=True)
metadata['condition'] = metadata['condition'].str.replace('_', ' ', regex=False)

# Crear un mapeo para las condiciones
condition_order = {
    'Primary Tumor': 1,
    'Solid Tissue Normal': 2
}
metadata['order'] = metadata['condition'].map(condition_order)


metadata = metadata.sort_values(by='order')
counts_df = counts_df.loc[metadata.index]

metadata = metadata.drop(columns=['order']) # Opcional: eliminar la columna de orden si no la necesitas

# Data filtering
#Before proceeding it is good practice to preprocess your data, e.g. to remove samples for which annotations are missing and exclude genes with very low levels of expression. 
#We start by removing samples that we dont want 
samples_to_keep = metadata['condition'].isin(['Primary Tumor', 'Solid Tissue Normal'])
counts_df = counts_df.loc[samples_to_keep]
metadata = metadata.loc[samples_to_keep]

# Filter genes with less that 10 read counts in total
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

# Single factor analysis: That is, we compare gene expressions of samples that have condition B to those that have condition A.

# Read counts modeling with the DeseqDataSet class
# We start by creating a DeseqDataSet object from the count and metadata data. A DeseqDataSet fits dispersion and log-fold change (LFC) 
# parameters from the data, and stores them.
inference = DefaultInference(n_cpus=8)
ref_level = ["condition", "Solid Tissue Normal"] # 

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
    inference=inference,
    ref_level=ref_level
)

# Once a DeseqDataSet was initialized, we may run the deseq2() method to fit dispersions and LFCs.
dds.deseq2()

# The DeseqDataSet class extends the AnnData class.
print(dds)

# Parameters are stored according to the AnnData data structure, with key-based data fields. In particular, 
# - X stores the count data,
# - obs stores design factors,
# - obsm stores sample-level data, such as "design_matrix" and "size_factors",
# - varm stores gene-level data, such as "dispersions" and "LFC".

# As an example, here is how we would access dispersions and LFCs (in natural log scale):

# Statistical analysis with the DeseqStats class
# Now that dispersions and LFCs were fitted, we may proceed with statistical tests to compute p-values and adjusted 
# p-values for differential expresion. This is the role of the DeseqStats class.
stat_res = DeseqStats(dds, contrast=["condition", "Primary Tumor", "Solid Tissue Normal"], inference=inference)
stat_res.summary()

# Obtenemos los resultados necesarios
# El volcano plot necesita:

# - Log2FoldChange (log2FoldChange): Diferencia de expresión en escala log2.
# - p-valor ajustado (padj): Para determinar la significancia estadística.
results_df = stat_res.results_df.copy()
# Ordenar results_df por la columna 'stat' en orden descendente
results_df = results_df.sort_values(by='stat', ascending=False)
results_df.to_csv(path_save_results, "results.csv")

gene_id_to_name = getBM[['Gene_ID', 'Gene_name']].drop_duplicates() # we are adding two columns with gene names to have more synonims. For most of them will be equal.
mapper = id_map(species = 'human')

# Filtrar gene_id_to_name para que contenga solo los genes que están en results_df
filtered_gene_names = gene_id_to_name[gene_id_to_name.Gene_ID.isin(results_df.index.tolist())]
# Reordenar filtered_gene_names para que siga el orden de results_df
ordered_gene_names = filtered_gene_names.set_index('Gene_ID').reindex(results_df.index).reset_index()

results_df['Gene_name'] = ordered_gene_names['Gene_name'].tolist()
results_df['Symbol'] = results_df.index.map(mapper.mapper) 

# Filtramos los significativos y los que tengan más menos del threshold
logfc_threshold = 2 # porque hay muchos
pval_threshold = 0.05
sign_df = results_df[(results_df.padj < pval_threshold) & (results_df.log2FoldChange.abs() > logfc_threshold)] # quedarnos con el score +10/-10 por encima, p

# The results are then stored in the results_df attribute (stat_res.results_df). As with as DeseqDataSet, the whole DeseqStats object may be 
# saved using pickle. However, it is often more convenient to have the results as a CSV. Hence, we may export stat_res.results_df as CSV, using 
# pandas.DataFrame.to_csv().

volcano(results_df, symbol='Gene_name')
plt.savefig('/scratch/jsanchoz/DeepRBP/stuff/example.png')

import numpy as np
import seaborn as sns

dds.layers['normed_counts']
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

sign_df

dds_sigs = dds[:, sign_df.index]
dds_sigs

grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                       index=dds_sigs.var_names, columns=dds_sigs.obs_names)


sns.clustermap(grapher, z_score=0, cmap = 'RdYlBu_r')
plt.savefig('/scratch/jsanchoz/DeepRBP/stuff/example.png', dpi=300) 
plt.tight_layout() 
plt.close()

## RBPs of interest
rbps_of_interest = ['SLU7', 'SRFBP1', 'SRSF1', 'SRSF2', 'SRSF3', 'UPF1', 'IGF2BP1', 'IGF2']
results_df[results_df.Gene_name.isin(rbps_of_interest)].index

dds_sub = dds[:, results_df[results_df.Gene_name.isin(rbps_of_interest)].index]

grapher = pd.DataFrame(dds_sub.layers['log1p'].T,
                       index=dds_sub.var_names, columns=dds_sub.obs_names)

grapher.index = grapher.index.map(mapper.mapper)
grapher.index.name = 'RBP name'  # Cambia 'Gene' por el nombre que desees
grapher.columns.name = 'Sample'

sns.clustermap(grapher, z_score=0, cmap = 'RdYlBu_r', figsize=(4,4))
plt.savefig('/scratch/jsanchoz/DeepRBP/stuff/example.png', dpi=300) 
plt.tight_layout() 
plt.close()

# usa los raw o los log2 pero simplemente redondenda porque es expected_counts.
#####################################################################
#####################################################################

# Result Columns Explanation

# baseMean:
# This column represents the mean of the normalized read counts for each gene
# across all samples in both conditions (Primary Tumor and Solid Tissue Normal).
# In other words, it is the average expression of the gene in both conditions.
# A higher value indicates that the gene is generally more expressed.

# log2FoldChange:
# This value indicates the change in gene expression between the two conditions 
# on a logarithmic scale (base 2).
# A positive value (e.g., 0.125828 for ENSG00000188976) indicates that the gene 
# is more expressed in "Primary Tumor" compared to "Solid Tissue Normal".
# A negative value (e.g., -0.253558 for ENSG00000242485) indicates that the gene 
# is less expressed in "Primary Tumor" compared to "Solid Tissue Normal".
# The log2 of the fold change is useful because it allows easy interpretation 
# of changes of two times (2^1 = 2) or more (2^2 = 4).

# lfcSE (Standard Error of the Log2 Fold Change):
# This is the standard deviation of the log2 fold change. 
# A lower value suggests that the estimation of the change in expression is 
# more reliable. For ENSG00000188976, the standard error is 0.085534, 
# indicating that the change estimation is quite precise.

# stat:
# This is the Wald test statistic used to calculate the p-value. 
# It measures how far the log2 fold change is from zero in terms of standard error. 
# A higher absolute value indicates a more significant change in expression.

# pvalue:
# This is the p-value associated with the Wald test. 
# It indicates the probability of observing a change as extreme as the observed 
# (or more extreme) if there is no true effect (i.e., if the true log2 fold change is zero).
# A low p-value (commonly < 0.05) suggests that the change in expression is 
# statistically significant. For ENSG00000242485, the p-value is 0.0081, 
# indicating that the change in expression is significant.

# padj (Adjusted p-value):
# This is the p-value adjusted for multiple tests using a method such as Benjamini-Hochberg.
# This adjustment is important in differential expression analysis where multiple gene 
# comparisons are made. A low adjusted value (also commonly < 0.05) indicates that 
# the change in expression is significant after considering the total number of tests.
# For example, ENSG00000130764 has a padj of 1.306835e-10, indicating that it is highly significant.

# Result Interpretation
# Significant Gene: A gene is considered significantly differentially expressed 
# if its p-value or padj is less than a threshold (commonly 0.05).
# Example: ENSG00000242485 has a p-value of 0.0081 and a padj of 0.0136, 
# indicating that this gene is differentially expressed between the conditions.

# Expression Change: The log2FoldChange indicates the direction and magnitude 
# of the change in expression.
# Example: ENSG00000130764 with a log2FoldChange of -0.442297 indicates that 
# the gene is significantly less expressed in "Primary Tumor" compared to 
# "Solid Tissue Normal".