from data_loader import CustomDataset
from config_loader import get_config

### main de otro archivo ####################################
config = get_config()
output_dir = config['train_prediction_model']['output_dir']
batch_size = config['train_prediction_model']['batch_size']
train_config = {k: config[k] for k in ('train_prediction_model', 'data')}
val_config = {k: config[k] for k in ('val_prediction_model', 'data')}
select_tumor_types = train_config['train_prediction_model']['train_tumor_types']

# opcion 1.1 - (primera ejecución con el config) [PROBADO Y FUNCIONA]
dataset = CustomDataset(
                    config=train_config, 
                    select_tumor_types = select_tumor_types,
                    output_dir=output_dir, 
                    save_files=True
                    )

# opcion 1.2 - (primera ejecución con paths) [PROBADO Y FUNCIONA]
paths = { # esto sobre todo para volver a ejecutar un pre-scaling data ya utilizado o unos paths q no estén en la config.
     "rbp_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/RBPs_log2p_tpm.csv",
     "gene_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/gn_expr_each_iso_tpm.csv",
     "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/trans_log2p_tpm.csv",
     "metadata_path": "/scratch/jsanchoz/DeepRBP/data/training_module/processed/TCGA/phenotype_metadata.csv"
    }

dataset = CustomDataset( 
                    config=train_config, 
                    select_tumor_types = select_tumor_types,
                    paths=paths, 
                    output_dir=output_dir,
                    save_files=True
                    )
###
# opcion 2 - (usar una ejecución ya creada anteriormente - reusar pre-scaling data con scaler e idx ya calculados anteriormente)  [PROBADO Y FUNCIONA]
paths = { # esto sobre todo para volver a ejecutar un pre-scaling data ya utilizado o unos paths q no estén en la config.
     "rbp_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-10-30/train_prediction_model/data/pre-scaling/rbp_expr_log2p_tpm.csv",
     "gene_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-10-30/train_prediction_model/data/pre-scaling/gene_expr_tpm.csv",
     "isoform_expr_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-10-30/train_prediction_model/data/pre-scaling/trans_expr_log2p_tpm.csv",
     "metadata_path": "/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-10-30/train_prediction_model/data/pre-scaling/metadata_df.csv"
    }

scaler, sigma, train_idx, valid_idx, test_idx = CustomDataset.load_scaler_and_idx(
                    path_save_files='/scratch/jsanchoz/DeepRBP/output/results/analysis/TCGA_Lung-Breast_2024-10-30/train_prediction_model/data'
                    )

dataset = CustomDataset(
                    config=train_config, 
                    paths=paths, 
                    output_dir=output_dir,
                    scaler=scaler,
                    sigma=sigma,
                    train_idx=train_idx,
                    valid_idx=valid_idx,
                    test_idx=test_idx
                    )

# opcion 3 - usar un scaler hecho en TCGA para GTEX (Me faltaría este)
# config

dataset = CustomDataset(
                     config=config, 
                     source_name=val_config, 
                     scaler=scaler,
                     sigma=sigma
                     )

dataloaders = create_dataloaders(dataset, batch_size=batch_size)

# Acceder a los loaders que se hayan creado
train_loader = dataloaders.get('train_loader', None)
valid_loader = dataloaders.get('validation_loader', None)
test_loader = dataloaders.get('test_loader', None)

# esto ya es parte del train.py
num_epochs = 10
for epoch in range(num_epochs):
    for batch_index, (rbp_exp, gene_exp, trans_exp) in enumerate(train_loader):
        print(type(rbp_exp), type(gene_exp), type(trans_exp))  # Verificar tipos
        print(rbp_exp, gene_exp, trans_exp)  # Inspeccionar contenido
        print(rbp_exp.shape, gene_exp.shape, trans_exp.shape) 
        break
    break

    
