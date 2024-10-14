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
        "joblib"  
    ],

    entry_points={
        'console_scripts': [
            'prepare-model-inputs=src.deeprbp.data_preprocessing.prep_model_inputs:main',
        ],
    },
    zip_safe=False
)
 