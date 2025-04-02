# /src/deeprbp/training_module/hyperparameter_optimization/grid_search_optuna.py

import os
import logging
import sys
import argparse
import pandas as pd
import optuna
from optuna.samplers import TPESampler
import optuna.visualization.matplotlib as optuna_plt
import optuna.logging

from ..util.logger import Logger
from torch.utils.data import DataLoader
from ..data_loading.config_loader import ConfigParser
from ..data_loading.data_loader import DataImporter, DataSplitter
from .train_model import TrainPredictor
from .model import PredictorModel
from .plots import plot_loss_curve
from .evaluation import calculate_metrics, save_metrics_summary, calculate_metrics_per_category

from ..util.utils import CustomTensorDataset, filter_data_by_sample_ids, save_data, adjust_batch_size

# Define the default configuration path
config_path = '/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_hyper_optimization.yaml'
output_dir = 'results/'

#args = parse_args()

# Load configuration settings from the specified YAML file
print("Loading configuration from:", config_path)
config_parser = ConfigParser(config_path)
config = config_parser.load_config()
print("Configuration loaded successfully.")

# Load the data using the data importer
print("Loading data...")
data_importer = DataImporter(config['data_paths'])
data = data_importer.load()
print("Data loaded successfully")

# Filter a portion of the samples to optimize time and computational resources
subset_idx, _ = DataSplitter.split_data_class(
                                            data=data, 
                                            config=config, 
                                            sample_category=config['sample_category'], 
                                            test_size=config['sample_fraction']
                                            )
data_subset = filter_data_by_sample_ids(data, subset_idx)

# Perform a train/validation split with the data subset
splitter = DataSplitter(data_subset, config)
train_data, valid_data = splitter.split_data_sets(test_name='validation')

# Save the training and val samples used in the optimization to the specified output directory
print(f"Saving training data to: {os.path.join(output_dir, 'hyper_opt_data/Train')}")
save_data(train_data, 
            os.path.join(output_dir, 'hyper_opt_data/Train'), 
            custom_names = {
                        'rbp_df': 'train_RBPs_log2p_tpm.csv',
                        'isoform_df': 'train_trans_log2p_tpm.csv',
                        'gene_df': 'train_gn_tpm.csv',
                        'metadata_df': 'train_phenotype_metadata.csv'}
                        )
print("Training data saved successfully.")

print(f"Saving test data to: {os.path.join(output_dir, 'hyper_opt_data/Val')}")
save_data(valid_data, 
            os.path.join(output_dir, 'hyper_opt_data/Val'), 
            custom_names = {
                        'rbp_df': 'val_RBPs_log2p_tpm.csv',
                        'isoform_df': 'val_trans_log2p_tpm.csv',
                        'gene_df': 'val_gn_tpm.csv',
                        'metadata_df': 'val_phenotype_metadata.csv'}
                        )
print("Test data saved successfully.")

# Scale data
scaler = Scaler()
train_data['scaled_rbp_df'] = scaler.fit_transform(train_data['rbp_df'])
valid_data['scaled_rbp_df'] = scaler.transform(valid_data['rbp_df'])

# Create data loaders

# el gene_df de alguna forma tienes que convertirse en gn_expr_each_iso_tpm dentro de la función de forward que ahora está reducida!
# tendré que poner de alguna manera los nombres de los genes y los nombres de los transcritos en el tensorDataset y hacer a continuación 
# la operación
getBM = pd.read_csv(config['getBM_path'])
train_dataset, valid_dataset = [
        CustomTensorDataset(
            data,
            getBM,
            rbp_data_key='scaled_rbp_df', 
            gene_data_key='gene_df', 
            transcript_data_key='isoform_df',
            trans_col_name=config['trans_col_name'],
            gene_col_name=config['gene_col_name']
        ) for data in [train_data, valid_data]]

## objective
def objective(trial, config, train_dataset, valid_dataset, val_batch_size=512):
    ### Suggest Optuna: Sample hyperparameters for this Trial 
    # Select hyperparameters
    num_hidden_layers = trial.suggest_int('num_hidden_layers', 0, 4)
    hidden1_nodes = trial.suggest_categorical('hidden1_nodes', [64, 128, 256, 512, 1024, 2048, 4096])  
    uniform_nodes = trial.suggest_categorical('uniform_nodes', [True, False])
    node_shrink_factor = trial.suggest_categorical('node_shrink_factor', [2, 4, 8])
    activation_func = trial.suggest_categorical('activation_func', ["relu", "tanh", "sigmoid"])
    learning_rate = trial.suggest_categorical('learning_rate', [0.0001, 0.001, 0.01])
    optimizer_name = trial.suggest_categorical('optimizer_name', ['sgd90', 'asgd', 'adam', 'adagrad', 'adadelta', 'adamW'])
    train_batch_size = trial.suggest_categorical('batch_size', [32, 64, 128, 256, 512, 1024, 2048, 4096]) 
    num_epochs = trial.suggest_categorical('num_epochs', [50, 100, 500, 1000, 2000, 3000])  # Sugerir valores categóricos para epochs

    config["num_hidden_layers"] = num_hidden_layers
    config["hidden1_nodes"] = hidden1_nodes
    config["uniform_nodes"] = uniform_nodes
    config["node_shrink_factor"] = node_shrink_factor
    config["activation_func"] = activation_func
    config["learning_rate"] = learning_rate
    config["optimizer_name"] = optimizer_name
    config["batch_size"] = train_batch_size
    config["num_epochs"] = num_epochs

    # Create loaders
    train_loader, val_loader = [
        DataLoader(
            dataset,
            batch_size=adjust_batch_size(dataset, (train_batch_size if idx == 0 else val_batch_size)),
            shuffle=(idx == 0),  # Solo hacer shuffle en el conjunto de entrenamiento
            drop_last=(idx == 0) # Solo drop_last en el conjunto de entrenamiento
        )
        for idx, dataset in enumerate([train_dataset, valid_dataset])   
    ]

    try:
        model = PredictorModel(
                input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1],
                output_size=next(iter(train_loader))['isoform_df'].shape[1],
                config=config
            )
            
        # Proceed with training the model
        trainer = TrainPredictor(
                model=model,
                config=config,
                input_features=('scaled_rbp_df', 'gene_df'), 
                output_features=('isoform_df',)
        )
        train_history, val_history, _ = trainer.fit(train_loader, 
                                                    val_loader, 
                                                    epochs=config["num_epochs"],
                                                    path_save_results='/scratch/jsanchoz/DeepRBP/stuff')


        ## pseudo try code:
        trainer.generate_predictions(val_loader)
            
        #preds_labels = [trainer.generate_predictions(loader) for loader in [train_loader, val_loader, test_loader]]
        #metrics = [calculate_metrics(pred, label) for pred, label, _ in preds_labels]


        plot_loss_curve(train_history, val_history, output_dir=self.path_save_results)
        self.logger.log("✅ Model and history saved.", level=1)

        # Here calculate other metrics with the trained model
        

        return val_history[-1]
    
    except ValueError as e:
        # Catch the ValueError raised in the PredictorModel
        print(f"Training failed due to configuration error: {e}")
        return float('inf')  # Return a high value to indicate this trial was unsuccessful
    

# The authors should provide the full table of results for the hyperparameter optimization runs
# to be able to validate the claim that more complex models (more hidden layers) are necessary.
 


    # Verify that the suggested trial is elegible
    # 1) si num_hidden_layers == 0 -> hidden1_nodes, uniform_nodes, node_shrink_factor y activation_func no aplican "not_used" (bien)
    # 2) si num_hidden_layers == 1 -> uniform_nodes y node_shrink_factor no aplican "not_used" (bien)
    # 3) si num_hidden_layers == 2 -> uniform_nodes no aplica "not_used" (bien)
    # 4) si num_hidden_layers == 3 -> no hay error (con uniform_nodes = False)
    # 5) si num_hidden_layers == 4 -> si hidden1_nodes es 64 y node_shrink_factor es 8, 128, 256, 512, 1024, 2048 error! (con uniform_nodes = False)
    # return model parameters (esto mejorar luego)
   
config["save_best_model" ] = True
config["num_hidden_layers"] = 3
config["hidden1_nodes"] = 128
config["uniform_nodes"] = False
config["node_shrink_factor"] = 2
config["activation_func"] = "relu"
config["learning_rate"] = 0.0001
config["optimizer_name"] = 'adamW'
config["batch_size"] = train_batch_size = 256
config["num_epochs"] = 300
val_batch_size = 256  









     
## hacer aqui las mias y mirar en el paper y en la revision cuantas tengo que pedir!
    # Put the actual suggestion in Config (aqui en el config tengo que comprobar que lo que está sugeriendo "
    # es plausible y si no lo es pasar a un nuevo trial o modificar las variables que no se vayan a UserWarning
    # como por ejemplo num_hidden_layers = 0 con todo lo demás. )"


## grid-Search with Optuna TPESampler
path_save_results = f'../results/optuna/{sample_id}'
study_name = "DeepRBP_Predictor-optimization"
optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
study = optuna.create_study(study_name=study_name, sampler=TPESampler(seed=config["seed"]), direction='minimize') # Recommended budgets with this sampler (#trials: 100-1000)
study.optimize(lambda trial: objective(trial, ), n_trials=n_trials) 








    print(f"Best trial: {study.best_trial}, with parameters: {study.best_params} and objective value:{study.best_value}")
    trials_df = study.trials_dataframe()
    trials_df.to_csv(path_save_results+f'/trials_df_{config["loss_func"]}_loss.csv', index=False)

    ### Visualize study history to analayze the hyperparams-performance relationship
    plt.rcParams['figure.figsize'] = (16*3.3, 9*3.3)
    plt.rcParams['figure.dpi'] = 300
    visualize_study_history(study, path_save_results)
    torch.cuda.empty_cache()


# Definir batch sizes
train_batch_size = 256   
val_batch_size = 256  

# Crear loaders
train_loader, val_loader = [
    DataLoader(
        dataset,
        batch_size=adjust_batch_size(dataset, (train_batch_size if idx == 0 else val_batch_size)),
        shuffle=(idx == 0),  # Solo hacer shuffle en el conjunto de entrenamiento
        drop_last=(idx == 0) # Solo drop_last en el conjunto de entrenamiento
    )
    for idx, dataset in enumerate([train_dataset, valid_dataset])   
]








model = PredictorModel(
            input_size=next(iter(train_loader))['scaled_rbp_df'].shape[1],
            output_size=next(iter(train_loader))['isoform_df'].shape[1],
            config=config
        )

 # Train the model
       # train_history, val_history, trainer = self.train_model(train_loader, val_loader)
       # self.save_model_and_history(trainer, train_history, val_history)

trainer = TrainPredictor(model=model, config=config)
train_history, val_history = trainer.fit(train_loader, val_loader, epochs=self.training_config['epochs'])
       

## ESTARIA GUAY QUE LA CLASE DeepRBPredictorPipeline fuera lo que usara aquí directamente, porque sino vamos 
# a repetir el mismo código dos veces

