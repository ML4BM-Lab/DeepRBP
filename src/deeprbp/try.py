### TRYING THE NEW CODE!!
# aqui hay que poner los puntos relativos.

from models import ExplainerModel
from deeprbp.explainer_validator_postar import ExplainerValidatorPostar

############################################## llamada ##############################################
# Uso de la clase
config_path_explain = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_explain.yaml'
config_path_train = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_tcga_train.yaml'

explainer_model = ExplainerModel(config_path_explain=config_path_explain, config_path_train=config_path_train)
data = explainer_model.load_and_process_data()
model = explainer_model.load_trained_predictor_model()
outputs = explainer_model.perform_explainer()

# Acceso a los resultados
df_transcript_scores = outputs['df_scores_TxRBP']
df_gene_scores = outputs['df_scores_GxRBP']
df_results_summary = outputs['result_table']

# Imprimir resultados
print("Transcript Scores DataFrame (TxRBP):")
print(df_transcript_scores)
print("Gene Scores DataFrame (GxRBP):")
print(df_gene_scores)
print("Result Table:")
print(df_results_summary)



# VALIDATOR POSTAR
validator = ExplainerValidatorPostar(explainer_model, df_results_summary)
# Cargar datos de POSTAR
df_postar_scores = validator.load_postar_data()
validator.process_postar_data(df_postar_scores, df_gene_scores)

# Ahora llama al método para contar y ordenar la matriz de POSTAR
df_rbps_per_gene_count, df_genes_per_rbp_count = validator.count_and_sort_postar_matrix()

# este orden lo voy a necesitar para los plots:  # Store the ordered lists of genes and RBPs
        #self.list_ordered_genes_postar = df_count_rbps_per_gen['Genes'].values.tolist()
        #self.list_ordered_rbps_postar = df_count_genes_per_rbp['RBPs'].values.tolist()

# Devolver los resultados summary.
df_results_summary = validator.return_summary_results()
# Calculate thresholds
df_optimal_thresholds = validator.calculate_rbp_thresholds(explainer_model.path_save_results)
# Save results
validator.save_results(explainer_model.path_save_results)

## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 
## ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  V ##  ##  ##  V V 



# code to yet develop 

# 3) Plot Explainer computed scores vs Postar

# 3) Plot DeepLIFT scores vs Postar (esto ayudarme de un x.py - piensa un nombre guay)
# ME HE QUEDADO AQUI BROTHER !!!

#def analyze_results_per_rbp_or_gene(df_score_GxRBP, df_val_nan_included_GxRBP, path_save, path_data, getBM): 
# Analyze the results per RNA binding protein (RBP) or gene, including plotting distributions, 
# ROC curves, and identifying NaN candidates.

# 


# (PARA ESTOS SINO BUSCA LA VERSION DE R Y LISTO. AHORA LO UNICO A SABER DONDE ESTÁ!!!!)
