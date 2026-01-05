# src/deeprbp/explainability_module/main_explainer.py

import os
import re
import argparse
import pandas as pd

from ..util.utils import print_if_main
from ..data_loading.config_loader import ConfigParser
from ..training_module.model import PredictorModel
from ..data_preparation.prepare_data import PrepareData
 
from .utils_explainer import (initialize_explainer_handler, filter_scores_for_low_expressed_transcripts, 
                              filter_scores_for_low_expressed_genes, collapse_transcript_scores_to_genes, 
                              _sanitize, _parse_categories)

def parse_args():   
    parser = argparse.ArgumentParser(description='Run the DeepRBP explainer to calculate transcript x RBP and genes x RBP explainability scores.')
    parser.add_argument('--config_path', type=str, required=True, help='Path to the configuration file for explainability.')
    parser.add_argument('--model_ckpt_path', type=str, required=True, help='Path to a trained DeepRBP predictor model checkpoint file')
    parser.add_argument('--scaler_dir', type=str, required=True, help='Directory where the trained scaler and sigma is located')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory where the results will be saved.')  
    parser.add_argument('--select_category', action='append',
        help='Category or categories to select (repeat flag or comma-separated). Overrides YAML select_category if provided.'
    )
    parser.add_argument('--analyze_hidden_layer', action='store_true', help='Also compute attribution scores for the last hidden layer (HL x RBP).')
    return parser.parse_args()

def main():
    args = parse_args()

    # Load configuration  
    print_if_main('\n[main_explainer] 🚀 Loading configuration...')
    config = ConfigParser(args.config_path) 
    print_if_main('\n[main_explainer] ──────────────────────────────────────')
    
    # Decide categories (CLI overrides YAML).
    # If none is provided -> single-run mode on ALL samples (no metadata/categories required).
    categories = _parse_categories(args.select_category)
    if not categories:
        cfg_cat = config.get('select_category')
        if isinstance(cfg_cat, (list, tuple)):
            categories = [str(x) for x in cfg_cat if str(x).strip()]
        elif isinstance(cfg_cat, str) and cfg_cat.strip():
            categories = [cfg_cat.strip()]

    if not categories:
        print_if_main("[main_explainer] No categories provided. Running in single-run mode (ALL samples).")
        categories = ["ALL"]

    # If user requested categories, we MUST have sample_category to filter by.
    if categories != ["ALL"] and not config.get("sample_category"):
        raise ValueError(
            "You provided --select_category / select_category but 'sample_category' is missing in the YAML. "
            "Set sample_category to the metadata column used for categories (e.g., 'detailed_category'), "
            "or run without categories (ALL samples)."
        )

    # Load trained model once
    print_if_main(f"\n[main_explainer] 🚀 Loading trained model. Using checkpoint: {args.model_ckpt_path}")
    model = PredictorModel.load_from_checkpoint(args.model_ckpt_path)
    print_if_main('\n[main_explainer] ──────────────────────────────────────')

    # Loop over categories
    for cat in categories:
        print_if_main(f"\n[main_explainer] 🚀 Preparing data... (category: {cat})")
        dataset, getBM = data_preparation(config, args, select_category=cat)
        print_if_main('\n[main_explainer] ──────────────────────────────────────')

        print_if_main(f"\n[main_explainer] 🚀 Running explainability pipeline... (category: {cat})")
        results = run_explainability_pipeline(config, dataset, model, getBM, analyze_hidden_layer=args.analyze_hidden_layer)
        print_if_main('\n[main_explainer] ──────────────────────────────────────')

        # Save results into a subfolder per category
        # In ALL mode, save directly under output_dir to keep paths clean.
        out_dir_cat = args.output_dir if cat == "ALL" else os.path.join(args.output_dir, _sanitize(cat))
        save_results(
            out_dir_cat, 
            results['df_scores_TxRBP'], 
            results['df_scores_GxRBP'], 
            results['result_table'],
            results['df_scores_HLxRBP'],
            results['df_per_sample']
        )
        
def data_preparation(config, args, select_category):
    getBM = pd.read_csv(config.get('getBM_path'))
    
    prep_data = PrepareData(
                getBM, 
                config.get('trans_col_name'), 
                config.get('gene_col_name'), 
                config.get('sample_category'), 
                args.output_dir
    )
    
    # If running without categories, disable filtering by category/condition
    # This makes the pipeline work for non-TCGA datasets without phenotype_metadata.csv.
    if select_category in (None, "ALL"):
        select_category = None
        dis_cond = None
        sel_cond = None
    else:
        dis_cond = config.get("disease_condition")
        sel_cond = config.get("select_condition")


    def _load(path, cat, dis_cond, sel_cond):
        return prep_data.load_data(
            path=path,
            select_category=cat,
            disease_condition=dis_cond,
            select_condition=sel_cond
        )

    # 1) First try with filter coming from yaml
    data = _load(
        path=config.get('test_path_files'),
        cat=select_category,
        dis_cond=dis_cond,
        sel_cond=sel_cond
    )

    #  (caso AML u otros “blood-derived” que no encajan con Primary_Tumor/Solid_Tissue_Normal)
    def _n_samples(d):
        return 0 if d is None or 'rbp_df' not in d or d['rbp_df'] is None else d['rbp_df'].shape[0]

    if _n_samples(data) == 0:
        print_if_main(
            f"[data_preparation] ⚠️ No samples found for category '{select_category}' "
            f"after applying disease_condition/select_condition "
            f"({config.get('disease_condition')}, {config.get('select_condition')}). "
            "Retrying without condition filters..."
        )
        data = _load(
            path=config.get('test_path_files'),
            cat=select_category,
            dis_cond=None,          # <<< desactiva el filtro por columna
            sel_cond=None           # <<< desactiva el filtro por valores
        )
    
    if _n_samples(data) == 0:
        raise ValueError(
            f"No samples available for category '{select_category}' "
            "even without disease/condition filters. Check phenotype_metadata.csv, sample_category and labels."
        )

    prep_data.load_scaler(args.scaler_dir)
    data_scaled = prep_data.scale_data(data)
    dataset = prep_data.create_tensor_dataset(data_scaled)
    return dataset, getBM

def run_explainability_pipeline(config, dataset, model, getBM, analyze_hidden_layer=False):
    print_if_main("\n [run_explainability_pipeline] 🚀 Starting the explanation process...")
    # 1) Initialize explainer
    print_if_main("\n[run_explainability_pipeline] Initializing explainer handler...")
    explainer_handler = initialize_explainer_handler(config, dataset, model)
    explain_method = config.get('explanation_method')
    print_if_main(f"[run_explainability_pipeline] ▶️ Explainer method: {explain_method}")
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")
    target_mode = config.get('target_mode', 'final')
    print_if_main(f"[run_explainability_pipeline] 🎯 target_mode: {target_mode}\n")
    
    # 2) Calculate scores at transcript level (output)
    print_if_main("\n[run_explainability_pipeline] Calculating scores at the transcript level...")
    df_scores_TxRBP, df_per_sample = explainer_handler.calculate_scores_transcript_level()
    print_if_main(f"\n[run_explainability_pipeline] ── Calculated transcript-level scores: {df_scores_TxRBP.shape[0]} transcripts.")  
    
    # Optional: save per-sample attributions (Tx x RBP x Samples).
    # If disabled, we keep only the collapsed TxRBP matrix.
    save_per_sample = config.get("save_per_sample_scores", default=False)

    if not save_per_sample:
        df_per_sample = None

    # 3) Filter output scores for low-expressed transcripts & low-expressed genes
    print_if_main("\n[run_explainability_pipeline] Filtering scores for low-expressed transcripts...")
    df_scores_TxRBP = filter_scores_for_low_expressed_transcripts(df_scores_TxRBP, dataset)   
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")
    print_if_main("\n[run_explainability_pipeline] Filtering scores for low-expressed genes...")
    df_scores_TxRBP = filter_scores_for_low_expressed_genes(df_scores_TxRBP, dataset) 
    print_if_main("\n[run_explainability_pipeline] ──────────────────────────────────────")
    
    # 4) Collapse output scores to genes (RBP x G)
    print_if_main("\n[run_explainability_pipeline] Collapsing scores from transcript level to gene level...")
    results = collapse_transcript_scores_to_genes(df_scores_TxRBP, getBM, config.get('gene_collapse_method'), dataset)  
    print_if_main("\n[run_explainability_pipeline] ── Collapsed scores to gene level successfully.")
    
    # 5) (Opcional) Last hidden layer explanation
    df_scores_HLxRBP = None
    if analyze_hidden_layer:
        if explain_method != "DeepLIFT":
            raise ValueError(
                "[run_explainability_pipeline] Hidden-layer attributions are only supported with DeepLIFT. "
                f"Current explainer is '{explain_method}'. Set explanation_method: 'DeepLIFT' or disable --analyze_hidden_layer."
            )
        print_if_main("\n[run_explainability_pipeline] 🔬 Calculating scores at the last hidden layer...")
        df_scores_HLxRBP = explainer_handler.calculate_scores_hidden_layer()
        print_if_main("\n[run_explainability_pipeline] ✅ Hidden layer scores computed.")

    print_if_main("\n[run_explainability_pipeline] ✅ Explanation process completed successfully.")
    return {
        'df_scores_TxRBP': df_scores_TxRBP,
        'df_scores_GxRBP': results.df_scores_GxRBP,
        'result_table': results.result_table,
        'df_scores_HLxRBP': df_scores_HLxRBP, 
        'df_per_sample': df_per_sample   
    } 

def save_results(path_save_results, df_scores_TxRBP, df_scores_GxRBP, result_table, df_scores_HLxRBP=None, df_per_sample=None):
    """Save the results as CSV files."""
    print_if_main("💾 Saving results...")
    os.makedirs(path_save_results, exist_ok=True)
    df_scores_TxRBP.to_csv(os.path.join(path_save_results, 'df_scores_TxRBP.csv'), index=True)
    df_scores_GxRBP.to_csv(os.path.join(path_save_results, 'df_scores_GxRBP.csv'), index=True)
    result_table.to_csv(os.path.join(path_save_results, 'result_table.csv'), index=True)
    
    if df_scores_HLxRBP is not None:
        df_scores_HLxRBP.to_csv(os.path.join(path_save_results, 'df_scores_HLxRBP.csv'), index=True)

    if df_per_sample is not None:
        df_per_sample.to_csv(os.path.join(path_save_results, 'df_scores_TxRBP_per_sample.csv'), index=True)
        print_if_main("✅ Saved per-sample attribution matrix (df_scores_TxRBP_per_sample.csv).")

    print_if_main("✅ Results saved successfully.")

if __name__ == "__main__":
    main()