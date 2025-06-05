# src/deeprbp/explainability_module/main_explainer.py

import os
import argparse
import pandas as pd

from ..util.utils import print_if_main
from ..data_loading.config_loader import ConfigParser
from ..training_module.model import PredictorModel
from ..data_preparation.prepare_data import PrepareData
 
from .utils_explainer import *

def data_preparation(config, args):
    getBM = pd.read_csv(config.get('getBM_path'))
    
    prep_data = PrepareData(
                getBM, 
                config.get('trans_col_name'), 
                config.get('gene_col_name'), 
                config.get('sample_category'), 
                args.output_dir
    )
    
    data = prep_data.load_data(
                path=config.get('test_path_files'), 
                select_category=config.get('select_category'),
                disease_condition=config.get('disease_condition'),
                select_condition=config.get('select_condition')
    )
    prep_data.load_scaler(args.scaler_dir)
    data_scaled = prep_data.scale_data(data)
    dataset = prep_data.create_tensor_dataset(data_scaled)
    return dataset, getBM

def main():
    args = parse_args()

    # Load configuration and auxiliary file
    print_if_main('\n[main_explainer] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path) 
    print_if_main('\n[main_explainer] ──────────────────────────────────────')
    
    # Data preparation
    print_if_main('\n[main_explainer] 🚀 Preparing data...')
    dataset, getBM = data_preparation(config, args)
    print_if_main('\n[main_explainer] ──────────────────────────────────────')

    # Load trained model
    print_if_main(f"\n[main_explainer] 🚀 Loading trained model. Using checkpoint: {args.model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(args.model_ckpt_path)
    print_if_main('\n[main_explainer] ──────────────────────────────────────')

    # Running explainability pipeline
    print_if_main(f"\n[main_explainer] 🚀 Running explainability pipeline...")
    results = run_explainability_pipeline(config, dataset, model, getBM)
    print_if_main('\n[main_explainer] ──────────────────────────────────────')

    # Save results as CSV files
    save_results(args.output_dir, results['df_scores_TxRBP'], results['df_scores_GxRBP'], results['result_table']) 
    
def run_explainability_pipeline(config, dataset, model, getBM):
    print_if_main("\n [run_explainability_pipeline] 🚀 Starting the explanation process...")
    # Initialize explainer
    print_if_main("\n[run_explainability_pipeline] Initializing explainer handler...")
    explainer_handler = initialize_explainer_handler(config, dataset, model)
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")

    # Calculate scores at transcript level
    print_if_main("\n[run_explainability_pipeline] Calculating scores at the transcript level...")
    df_scores_TxRBP = explainer_handler.calculate_scores_transcript_level()
    print_if_main(f"\n[run_explainability_pipeline] ── Calculated transcript-level scores: {df_scores_TxRBP.shape[0]} transcripts.") # si no es el 0 es el 1
    
    # Filter scores for low-expressed transcripts
    print_if_main("\n[run_explainability_pipeline] Filtering scores for low-expressed transcripts...")
    df_scores_TxRBP = filter_scores_for_low_expressed_transcripts(df_scores_TxRBP, dataset)   
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")

    # Filter scores for low-expressed genes
    print_if_main("\n[run_explainability_pipeline] Filtering scores for low-expressed genes...")
    df_scores_TxRBP = filter_scores_for_low_expressed_genes(df_scores_TxRBP, dataset) 
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")

    # Collapse scores to genes (RBP x G)
    print_if_main("\n[run_explainability_pipeline] Collapsing scores from transcript level to gene level...")
    results = collapse_transcript_scores_to_genes(df_scores_TxRBP, getBM, config.get('gene_collapse_method'), dataset)  
    print_if_main("\n[run_explainability_pipeline] ── Collapsed scores to gene level successfully.")
    
    print_if_main("\n[run_explainability_pipeline] ✅ Explanation process completed successfully.")
    return {
        'df_scores_TxRBP': df_scores_TxRBP,
        'df_scores_GxRBP': results.df_scores_GxRBP,
        'result_table': results.result_table
    }

def save_results(path_save_results, df_scores_TxRBP, df_scores_GxRBP, result_table):
    """Save the results as CSV files."""
    try:
        print_if_main("💾 Saving results...")
        os.makedirs(path_save_results, exist_ok=True)
        df_scores_TxRBP.to_csv(os.path.join(path_save_results, 'df_scores_TxRBP.csv'), index=True)
        df_scores_GxRBP.to_csv(os.path.join(path_save_results, 'df_scores_GxRBP.csv'), index=True)
        result_table.to_csv(os.path.join(path_save_results, 'result_table.csv'), index=True)
        print_if_main("✅ Results saved successfully.")
    except Exception as e:
        raise RuntimeError(f"❌ Error saving results: {e}")

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP explainer to calculate transcript x RBP and genes x RBP explainability scores.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file for explainability.')
    parser.add_argument('--model_ckpt_path', type=str, required=True, help='Path to a trained DeepRBP predictor model checkpoint file')
    parser.add_argument('--scaler_dir', type=str, required=True, help='Directory where the trained scaler and sigma is located')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory where the results will be saved.')  
    return parser.parse_args()

if __name__ == "__main__":
    main()
