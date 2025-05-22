import setuptools
from setuptools import setup

setup(
    name='deeprbp',
    version='1.0.0',
    description='Deep learning model for inferring RBP-Gene regulations',
    url='https://github.com/ML4BM-Lab/DeepRBP',
    author='Joseba Sancho Zamora',
    packages=setuptools.find_packages(),
    install_requires=[
        "h5py",
        "matplotlib",
        "pandas",
        "scanpy",
        "scikit-learn",
        "numpy",
        "scipy",
        "seaborn",
        "statsmodels",
        "torch",
        "torchvision",
        "tornado",
        "tqdm",
        "plotly",
        "notebook",
        "ipykernel",
        "openpyxl",
        "captum",        
        "optuna",        
        "joblib",
        "pydeseq2",
        "sanbomics",
        "colorama",
        "torchinfo",
        "psutil",
        "lightning",
        "GPUtil",
        "optuna-integration[pytorch_lightning]",
        "tabulate"
    ],
    entry_points={
        'console_scripts': [
            'preprocess-data=src.deeprbp.data_preprocessing.preprocess_data:main',
            'split-and-save=src.deeprbp.data_preprocessing.split_data_and_save:main',
            'run-hyper-optimization-optuna=src.deeprbp.training_module.hyperparameter_optimization.grid_search_optuna:main',
            'run-deeprbp-predictor=src.deeprbp.training_module.main_predictor:main',
            'run-deeprbp-explainer=src.deeprbp.explainability_module.main_explainer:main',
            'run-postar-validator=src.deeprbp.explainability_module.postar_validation.postar_validator:main'
        ],
    },
    zip_safe=False,
)

 
