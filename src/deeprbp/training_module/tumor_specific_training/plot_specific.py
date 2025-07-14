
import matplotlib.pyplot as plt
import seaborn as sns


ABBREVIATION_DICT = {
    'Thyroid_Carcinoma': 'THCA',
    'Testicular_Germ_Cell_Tumor': 'TGCT',
    'Prostate_Adenocarcinoma': 'PRAD',
    'Skin_Cutaneous_Melanoma': 'SKCM',
    'Sarcoma': 'SARC',
    'Mesothelioma': 'MESO',
    'Uterine_Corpus_Endometrioid_Carcinoma': 'UCEC',
    'Pheochromocytoma_&_Paraganglioma': 'PCPG',
    'Uterine_Carcinosarcoma': 'UCS',
    'Lung_Adenocarcinoma': 'LUAD',
    'Stomach_Adenocarcinoma': 'STAD',
    'Uveal_Melanoma': 'UVM',
    'Thymoma': 'THYM',
    'Lung_Squamous_Cell_Carcinoma': 'LUSC',
    'Rectum_Adenocarcinoma': 'READ',
    'Ovarian_Serous_Cystadenocarcinoma': 'OV',
    'Pancreatic_Adenocarcinoma': 'PAAD',
    'Kidney_Clear_Cell_Carcinoma': 'KIRC',
    'Glioblastoma_Multiforme': 'GBM',
    'Head_&_Neck_Squamous_Cell_Carcinoma': 'HNSC',
    'Liver_Hepatocellular_Carcinoma': 'LIHC',
    'Colon_Adenocarcinoma': 'COAD',
    'Cervical_&_Endocervical_Cancer': 'CESC',
    'Diffuse_Large_B_Cell_Lymphoma': 'DLBC',
    'Breast_Invasive_Carcinoma': 'BRCA',
    'Esophageal_Carcinoma': 'ESCA',
    'Kidney_Chromophobe': 'KICH',
    'Kidney_Papillary_Cell_Carcinoma': 'KIRP',
    'Cholangiocarcinoma': 'CHOL',
    'Acute_Myeloid_Leukemia': 'LAML',
    'Bladder_Urothelial_Carcinoma': 'BLCA',
    'Brain_Lower_Grade_Glioma': 'LGG',
    'Adrenocortical_Cancer': 'ACC',
    'all': 'all'  # para mantener coherencia si "all" es una categoría especial
}


def plot_metric_confusion_matrix(
    df,
    metric_name,
    output_path="confusion_matrix.png",
    title=None
):
    metric_name = metric_name.lower()
    error_metrics = ['mse']
    if title is None:
        if metric_name in error_metrics:
            title = f"{metric_name.replace('_', ' ').upper()} across tumor types"
            cbar_label = metric_name.replace('_', ' ').upper()
        else:
            title = f"{metric_name.replace('_', ' ').capitalize()} across tumor types"
            cbar_label = metric_name.replace('_', ' ').capitalize()
    df_plot = df.copy()
    df_plot.index = [ABBREVIATION_DICT.get(idx, idx) for idx in df_plot.index]
    df_plot.columns = [ABBREVIATION_DICT.get(col, col) for col in df_plot.columns]
    df_plot.index.name = 'Train Tumor Type'
    df_plot.columns.name = 'Eval Tumor Type'
    fig_width = min(max(4, len(df_plot.columns) * 0.6), 20)   # desde 6 hasta 20
    fig_height = min(max(3, len(df_plot) * 0.45), 15)  
    plt.figure(figsize=(fig_width, fig_height))
    base_fontsize = 10
    scaling_factor = min(fig_width / 10, fig_height / 6)
    font_size = max(6, base_fontsize * scaling_factor)
    cmap = "YlGnBu"
    if metric_name in error_metrics:
        cmap += "_r"
    vmin, vmax = df_plot.min().min(), df_plot.max().max()
    ax = sns.heatmap(
        df_plot,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        cbar_kws={"label": cbar_label},
        annot_kws={"fontsize": font_size * 0.6},
        linewidths=0.5,
        linecolor='lightgrey',
        vmin=vmin,
        vmax=vmax
    )
    ax.set_title(title, fontsize=font_size+2, pad=20)
    plt.xticks(rotation=90, fontsize=font_size)
    plt.yticks(rotation=0, fontsize=font_size)
    plt.xlabel(df_plot.columns.name, fontsize=font_size, labelpad=15)
    plt.ylabel(df_plot.index.name, fontsize=font_size, labelpad=15)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=font_size)
    cbar.set_label(cbar_label, fontsize=font_size, labelpad=6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


