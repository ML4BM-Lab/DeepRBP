import pandas as pd
import os

filename = './human.txt'

# Define a dictionary to map unique values to their corresponding morphology, tissue, disease values and extra infor
mapping_dict = {
    'HeLa': {'morphology': 'epithelial', 'tissue':'Uterus_Cervix', 'disease': 'Adenocarcinoma', 
             'extra_info':'https://www.atcc.org/products/ccl-2'},
    
    'HEK293T': {'morphology': 'epithelial', 'tissue': 'Kidney_embryo', 'disease': None, 
                'extra_info':'Embryo, https://www.atcc.org/products/crl-3216'},
    
    'HEK293': {'morphology': 'epithelial', 'tissue': 'Kidney_embryo', 'disease': None, 
               'extra_info':'Embryo, https://www.atcc.org/products/crl-1573'},
    
    'Brain': {'morphology': None, 'tissue': 'Brain', 'disease': None, 'extra_info':None},
    
    'K562': {'morphology': 'lymphoblast', 'tissue': 'Myeloid', 'disease': 'Chronic_Myelogenous_Leukemia', 
             'extra_info':'https://www.atcc.org/products/ccl-243#:~:text=K%2D562%20are%20lymphoblast%20cells,system%20disorder%20and%20immunology%20research.&text=Discounts%20may%20be%20available%20for%20our%20fellow%20nonprofit%20organizations'},
    
    'HEK_293_FRT': {'morphology': 'epithelial', 'tissue': 'Kidney_embryo', 'disease': None, 
                    'extra_info':'Embryo, T antigen, which refers to the large T-antigen of the simian virus 40 (SV40) that is expressed in these cells'},
    
    'HepG2': {'morphology': 'epithelial', 'tissue': 'Liver', 'disease': '', 'extra_info':None},
    
    'TREX_FLP-In_293T_cells': {'morphology': None, 'tissue': 'Kidney_embryo', 'disease': None, 'extra_info':'Embryo'},
    
    'T_cell': {'morphology': None, 'tissue': 'T_cell', 'disease': None, 'extra_info':None},
    
    'HEK293T_no_SRRM4': {'morphology': None, 'tissue': 'Kidney_embryo', 'disease': None, 
                        'extra_info':'Embryo, a variant or subclone of the HEK293T cell line that has been genetically modified to knock down or remove the expression of the SRRM4 gene (protein involved in splicing of RNA)'},
    
    'MGG8': {'morphology': None, 'tissue': 'Brain_gliobastoma', 'disease': 'Glioblastoma', 
             'extra_info':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1324083'},
    
    'Huh7': {'morphology': 'epithelial', 'tissue': 'Liver', 'disease': 'Hepatocellular_Carcinoma', 
             'extra_info':'https://huh7.com/#:~:text=Huh%2D7%20is%20an%20immortal,typically%20grow%20as%202D%20monolayers'},
    
    'Cardiac': {'morphology': None, 'tissue': 'myocardium', 'disease': None, 'extra_info':None},
    
    'HEK293T_plus_SRRM4': {'morphology': 'epithelial', 'tissue': 'Kidney_embryo', 'disease': None, 
                           'extra_info':'Embryo, Overexpression or increased expression of SRRM4 can have various effects on cellular processes and gene expression, depending on the specific context'},
    
    'SH-SY5Y': {'morphology': 'epithelial', 'tissue': 'Bone_Marrow', 'disease': 'Neuroblastoma', 
                'extra_info':'https://www.atcc.org/products/crl-2266'},
    
    '22RV1': {'morphology': 'epithelial', 'tissue': 'Prostate_carcinoma', 'disease': 'Carcinoma', 
              'extra_info':'https://www.atcc.org/products/crl-2505'},
    
    'HME': {'morphology': 'epithelial', 'tissue': 'Breast', 'disease': None, 'extra_info':None},
    
    'adrenal gland': {'morphology': None, 'tissue': 'adrenal_gland', 'disease': None, 'extra_info':None},
    
    'HeLa_treated_with_puromycoin': {'morphology': 'epithelial', 'tissue': 'Uterus_Cervix', 'disease': 'Adenocarcinoma', 
                                      'extra_info':None},
    
    'T98G': {'morphology': 'fibroblast', 'tissue': 'Brain_gliobastoma', 'disease': 'Glioblastoma_Multiforme', 
              'extra_info':'https://www.atcc.org/products/crl-1690'},
    
    'PC3': {'morphology': 'epithelial', 'tissue': 'Prostate_adenocarcinoma', 'disease': 'Adenocarcinoma', 
            'extra_info':'Grade IV, https://www.atcc.org/products/crl-1435'},
    
    'SK-MEL-103': {'morphology': None, 'tissue': 'Skin', 'disease': 'Malignant_Melanoma', 'extra_info':None},
    
    'MDA-LM2': {'morphology': 'epithelial', 'tissue': 'Breast_adenocarcinoma', 'disease': 'Adenocarcinoma', 'extra_info':'Mammary gland'},
    
    'LAPC4': {'morphology': None, 'tissue': 'Prostate_carcinoma', 'disease': 'Carcinoma', 'extra_info':None},
    
    'HK-2': {'morphology': None, 'tissue': 'Kidney_papilloma', 'disease': 'Papilloma', 'extra_info':'Cortex; Proximal tubule'},
    
    'DU145': {'morphology': 'epithelial', 'tissue': 'Prostate_carcinoma', 'disease': 'Carcinoma', 
              'extra_info':'https://www.atcc.org/products/htb-81'},
    
    'LNCaP': {'morphology': None, 'tissue': 'Prostate_carcinoma', 'disease': 'Carcinoma', 'extra_info':None},
    
    'Cortical': {'morphology': None, 'tissue': 'Cortical', 'disease': None, 'extra_info':None},
    
    'H9': {'morphology': None, 'tissue': 'Blastocyst_Embryo', 'disease': None, 'extra_info':None},
    
    'frontalcortex': {'morphology': None, 'tissue': 'Brain', 'disease': None, 'extra_info':None},
    
    'Panc1': {'morphology': None, 'tissue': 'Pancreas_epith_carcinoma', 'disease': 'Epithelioid_Carcinoma', 'extra_info':'Duct'},
    
    'A2780': {'morphology': None, 'tissue': 'Ovarian', 'disease': 'Cancer', 'extra_info':None},
    
    'MDA-MB-231': {'morphology': 'epithelial', 'tissue': 'Breast_adenocarcinoma', 'disease': 'Adenocarcinoma', 'extra_info':'Mammary gland'},
    
    'hippocampus': {'morphology': None, 'tissue': 'Brain', 'disease': None, 'extra_info':None},
    
    'PL45': {'morphology': 'epithelial', 'tissue': 'Pancreas_adenocarcinoma', 'disease': 'Adenocarcinoma', 'extra_info':'Ductal'},
    
    'HCT116': {'morphology': None, 'tissue': 'Colon', 'disease': 'Carcinoma', 'extra_info':'Colorectal, Large intestine'}
}

def replace_comma(filename):
    with open(filename, 'r', encoding='ISO-8859-1') as file:
        lines = file.readlines()
    with open(filename, 'w') as file:
        for line in lines:
            line = line.replace(',', '_')
            file.write(line)  

def process_postar3(filename, mapping_dict):
    # 0) Replace the commas in file with "_"
    replace_comma(filename)
    # 1) Read the Postar3
    df_postar3 = pd.read_csv(filename, sep='\t', header=None)
    # 2) Add the colnames
    df_postar3.columns = ["chromosome", "start", "end", "RBP_specif", "strand", "RBP_name", "technique", "raw_tissue", "experiment", "PhastCons_Score", "PhyloP_Score"]
    # 3) Reorder columns
    df_postar3 = df_postar3[['RBP_name', 'chromosome', "start", "end", "raw_tissue", "technique", "experiment", "strand", "PhastCons_Score", "PhyloP_Score"]]
    # 4) Remove RBP_name == 'RBP_occupancy'
    df_postar3 = df_postar3[df_postar3['RBP_name'] != 'RBP_occupancy']
    # 5) Remove the samples with raw_tissue equal to NaN
    df_postar3 = df_postar3.dropna(subset=['raw_tissue'])
    # 6) Create a Table to descrive each of the raw_tissue and add 'Tissue' column to df_postar3
    # Create 'result' directory if it doesn't exist
    if not os.path.exists('./results'):
        os.makedirs('./results')
    df_tissue_description = pd.DataFrame.from_dict(mapping_dict, orient='index')
    df_tissue_description.to_csv('./results/df_tissue_description.txt')
    df_postar3['tissue'] = df_postar3['raw_tissue'].apply(lambda x: mapping_dict[x]['tissue'])
    # 7) Create the tissue-specific Posta3 files
    [df_postar3[df_postar3.tissue == tissue].reset_index(drop='Index').to_csv(f'./results/human_{tissue}.txt') for tissue in df_postar3.tissue.unique()]

if __name__ == '__main__':
    process_postar3(filename, mapping_dict)