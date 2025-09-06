# tcga_codes.py
from __future__ import annotations

# Mapa canónico categoría → código TCGA
TCGA_CODE = {
    "Breast_Invasive_Carcinoma": "BRCA",
    "Kidney_Clear_Cell_Carcinoma": "KIRC",
    "Lung_Adenocarcinoma": "LUAD",
    "Thyroid_Carcinoma": "THCA",
    "Head_&_Neck_Squamous_Cell_Carcinoma": "HNSC",
    "Lung_Squamous_Cell_Carcinoma": "LUSC",
    "Prostate_Adenocarcinoma": "PRAD",
    "Brain_Lower_Grade_Glioma": "LGG",
    "Skin_Cutaneous_Melanoma": "SKCM",
    "Stomach_Adenocarcinoma": "STAD",
    "Bladder_Urothelial_Carcinoma": "BLCA",
    "Ovarian_Serous_Cystadenocarcinoma": "OV",
    "Liver_Hepatocellular_Carcinoma": "LIHC",
    "Colon_Adenocarcinoma": "COAD",
    "Kidney_Papillary_Cell_Carcinoma": "KIRP",
    "Cervical_&_Endocervical_Cancer": "CESC",
    "Sarcoma": "SARC",
    "Uterine_Corpus_Endometrioid_Carcinoma": "UCEC",
    "Esophageal_Carcinoma": "ESCA",
    "Pheochromocytoma_&_Paraganglioma": "PCPG",
    "Pancreatic_Adenocarcinoma": "PAAD",
    "Acute_Myeloid_Leukemia": "LAML",
    "Glioblastoma_Multiforme": "GBM",
    "Testicular_Germ_Cell_Tumor": "TGCT",
    "Thymoma": "THYM",
    "Rectum_Adenocarcinoma": "READ",
    "Kidney_Chromophobe": "KICH",
    "Mesothelioma": "MESO",
    "Adrenocortical_Cancer": "ACC",
    "Uveal_Melanoma": "UVM",
    "Uterine_Carcinosarcoma": "UCS",
    "Diffuse_Large_B_Cell_Lymphoma": "DLBC",
    "Cholangiocarcinoma": "CHOL",
}

ORDER_POS = {code: i for i, code in enumerate(TCGA_CODE.values())}

def get_tcga_code(label: str) -> str:
    """Return the short TCGA code if label is a long name or already a code."""
    if label in ORDER_POS:          
        return label
    return TCGA_CODE.get(label, label)
