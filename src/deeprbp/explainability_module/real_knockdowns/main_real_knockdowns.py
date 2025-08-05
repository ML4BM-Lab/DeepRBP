
# src/deeprbp/explainability_module/real_knockdowns/main_real_knockdowns.py

import argparse
from ...data_loading.config_loader import ConfigParser
from ...training_module.model import PredictorModel
from ...training_module.evaluation import evaluate_and_visualize

from ...explainability_module.main_explainer import run_explainability_pipeline, save_results
from .utils_kds import build_datamodule, build_trainer

def parse_args():
    parser = argparse.ArgumentParser(
        description='Execute the evaluation of the DeepRBP trained predictor on both control and knockdown datasets. ' \
        'Additionally, compute explainability scores specifically for the control dataset.'
    )
    parser.add_argument('--config_path', type=str, required=False, 
        default='/scratch/jsanchoz/DeepRBP/src/deeprbp/configs/config_real_knockdowns.yaml',
        help='Specifies the path to the YAML configuration file containing all necessary parameters for processing.'
    )
    parser.add_argument('--processed_data_dir', type=str, required=True, 
        help='Indicates the directory where processed input data and corresponding labels for control and knockdown conditions are stored.'
    )
    parser.add_argument('--output_dir', type=str, required=True, help='Designates the directory where the output results, ' \
          'including evaluation metrics and explainability scores, will be saved.'
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # Load the base configuration from YAML
    print("\n[main_real_knockdowns] 📄 Loading base configuration from YAML...")
    config_base = ConfigParser(config_path=args.config_path)
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Crear DataModules
    print("\n[main_real_knockdowns] 🔧 Creating DataModules for control and knockdown conditions...")
    dm_control = build_datamodule(config_base, "control", args.processed_data_dir, args.output_dir)
    dm_knockdown = build_datamodule(config_base, "knockdown", args.processed_data_dir, args.output_dir)
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Crear Trainers usando el mismo output_dir que el datamodule
    print("\n[main_real_knockdowns] Creating Trainers for each condition...")
    trainer_control = build_trainer(config_base, dm_control.output_dir)
    trainer_knockdown = build_trainer(config_base, dm_knockdown.output_dir)
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Setup test data for both DataModules
    print("\n[main_real_knockdowns] 🔄 Setting up test data for both DataModules...")
    dm_control.setup(stage='test')     # Loads, scales, and prepares control test data
    dm_knockdown.setup(stage='test')   # Loads, scales, and prepares knockdown test data
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Access the final processed tensor datasets
    print("\n[main_real_knockdowns] 📊 Accessing final processed tensor datasets...")
    dataset_control = dm_control.test_dataset
    #dataset_knockdown = dm_knockdown.test_dataset
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Load trained model
    print('\n[main_real_knockdowns] 🚀 Loading trained model...')
    model = PredictorModel.load_from_checkpoint(dm_control.config.get('model_checkpoint_path'))
    print('\n[main_real_knockdowns] ──────────────────────────────────────')
        
    # Evaluate on CONTROL
    print("\n[main_real_knockdowns] 🧪 Testing model on CONTROL data...")
    trainer_control.test(model, dm_control)
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Evaluate on KNOCKDOWN
    print("\n[main_real_knockdowns] 🧪 Testing model on KNOCKDOWN data...")
    trainer_knockdown.test(model, dm_knockdown)
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Evaluate model performance in detail (here all samples is just one category)
    print("\n[main_real_knockdowns] 📊 Evaluate model performance in detail ...")
    datasets = [
        {
            "name": "control",
            "dm": dm_control,
            "trainer": trainer_control
        },
        {
            "name": "knockdown",
            "dm": dm_knockdown,
            "trainer": trainer_knockdown
        }
    ]

    for entry in datasets:
        print(f"\n[main_real_knockdowns] 📈 Evaluating on the {entry['name'].upper()} set...")
        evaluate_and_visualize(
            trainer=entry["trainer"],
            model=model,
            dm=entry["dm"],
            output_dir=entry["dm"].output_dir,   
            set_name=entry["name"]
        )
    print('\n[main_real_knockdowns] ──────────────────────────────────────')

    # Calculate explainability scores using control data
    results = run_explainability_pipeline(config_base, dataset_control, model, dm_control.getBM)

    # Save results as CSV files
    save_results(dm_control.output_dir, results['df_scores_TxRBP'], results['df_scores_GxRBP'], results['result_table']) 
    
if __name__ == "__main__":
    main()
