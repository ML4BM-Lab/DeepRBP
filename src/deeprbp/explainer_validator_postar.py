# # /scratch/jsanchoz/DeepRBP/src/deeprbp/postar_utils.py
                              
# import pandas as pd
# import numpy as np

# def match_scores_and_postar_data(deeplift_scores_genes, postar_score_genes):
#     """
#     Parameters:
#         deeplift_scores_genes (pd.DataFrame): DataFrame containing DeepLIFT scores at the gene level,
#                                         where rows correspond to genes and columns 
#                                         correspond to RBPs. 
#         postar_score_genes (pd.DataFrame): DataFrame containing Postar experiment results at the gene level,
#                                         where rows correspond to genes and columns 
#                                         correspond to RBPs. 
        
#     Returns:
#         pd.DataFrame: Modified DeepLIFT scores DataFrame with matching genes and RBPs.
#         pd.DataFrame: Modified Postar scores DataFrame with matching genes and RBPs, with NaNs where applicable.
#     """
    
#     ### Find the matching and non-matching genes and RBPs between the Scores dataframe and Postar
#     genes_match = [x for x in deeplift_scores_genes.index if x in postar_score_genes.index]
#     genes_not_match = [x for x in deeplift_scores_genes.index if x not in postar_score_genes.index]  # Genes present in deepLIFT but not in Postar
#     rbps_match = [x for x in deeplift_scores_genes.columns if x in postar_score_genes.columns]
#     rbps_not_match = [x for x in deeplift_scores_genes.columns if x not in postar_score_genes.columns]  # RBPs present in deepLIFT but not in Postar
#     print('[match_scores_and_postar_data] Find the matching and non-matching genes and RBPs between the Scores dataframe and Postar')
    
#     ### Process postar_score_genes matrix to match deeplift_scores_genes' shape, with NaN values for the RBPs-genes not presented in Postar
#     # Create the postar dataframe with NaN values for the RBPs and genes where we have scores and reorder rows and columns
#     # Concatenate with NaN values
#     nan_df = pd.DataFrame(index=genes_not_match, columns=rbps_not_match, dtype=np.float32).fillna(np.nan)
#     postar_score_genes_with_nan = pd.concat([postar_score_genes, nan_df], axis=0)  # Concatenate along rows (axis=0)
    
#     # Reindex the dataframes to align their shapes
#     deeplift_scores_genes = deeplift_scores_genes.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
#     postar_score_genes_with_nan = postar_score_genes_with_nan.loc[genes_match + genes_not_match, rbps_match + rbps_not_match]
    
#     # Renaming indices for the postar_score_genes_with_nan DataFrame
#     postar_score_genes_with_nan.index.name = 'Gene_ID'
#     postar_score_genes_with_nan.columns.name = 'RBP_ID'
#     return deeplift_scores_genes, postar_score_genes_with_nan

# def count_and_sort_postar_matrix(matched_postar_scores):
#     """
#     Analyze the POSTAR matrix to count the number of RBPs per gene and the number of genes per RBP,
#     and sort the results based on the number of Class 1 occurrences.

#     Parameters:
#         matched_postar_scores (pd.DataFrame): DataFrame containing POSTAR scores with Gene_ID as index and RBP_ID as columns.

#     Returns:
#         Tuple[pd.DataFrame, pd.DataFrame]: DataFrames containing counts of RBPs per gene and genes per RBP,
#                                             both sorted by the number of Class 1 occurrences.
#     """
#     # Create a filtered copy of the matched POSTAR scores
#     matched_postar_scores_filtered = matched_postar_scores.copy()
#     # Number of RBPs per Gene
#     df_count_rbps_per_gen = pd.DataFrame()
#     df_count_rbps_per_gen['Class 0'] = matched_postar_scores_filtered.apply(lambda x: (x == 0).sum(), axis=1)
#     df_count_rbps_per_gen['Class 1'] = matched_postar_scores_filtered.apply(lambda x: (x == 1).sum(), axis=1)
#     df_count_rbps_per_gen['Class NaN'] = matched_postar_scores_filtered.apply(lambda x: x.isna().sum(), axis=1)
#     df_count_rbps_per_gen['Genes'] = df_count_rbps_per_gen.index
#     df_count_rbps_per_gen = df_count_rbps_per_gen.reset_index(drop=True)
#     # Sort the RBPs per Gene DataFrame by Class 1 in descending order
#     df_count_rbps_per_gen = df_count_rbps_per_gen.sort_values(by='Class 1', ascending=False).reset_index(drop=True)
#     # Number of genes per RBP
#     df_count_genes_per_rbp = pd.DataFrame()
#     df_count_genes_per_rbp['Class 0'] = matched_postar_scores_filtered.apply(lambda x: (x == 0).sum(), axis=0)
#     df_count_genes_per_rbp['Class 1'] = matched_postar_scores_filtered.apply(lambda x: (x == 1).sum(), axis=0)
#     df_count_genes_per_rbp['Class NaN'] = matched_postar_scores_filtered.apply(lambda x: x.isna().sum(), axis=0)
#     df_count_genes_per_rbp['RBPs'] = df_count_genes_per_rbp.index
#     df_count_genes_per_rbp = df_count_genes_per_rbp.reset_index(drop=True)
#     # Sort the Genes per RBP DataFrame by Class 1 in descending order
#     df_count_genes_per_rbp = df_count_genes_per_rbp.sort_values(by='Class 1', ascending=False).reset_index(drop=True)
#     return df_count_rbps_per_gen, df_count_genes_per_rbp
