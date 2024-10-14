import pandas as pd
import os

def compare_patients_with_excel(output_dir: str, patient_ids: set):
    """Compare processed patients with patient IDs and log results."""
    # Load the processed patients from the Excel file
    processed_patients_file = os.path.abspath(f'{output_dir}/processed_patients.xlsx')
    df_processed_patients = pd.read_excel(processed_patients_file)
    # Create a set of processed patient IDs for comparison
    processed_patient_ids = set(df_processed_patients['Patient_ID'])
    # Filter out patient IDs that start with 'TARGET-*' from the original patient_ids
    patient_ids = {pid for pid in patient_ids if not pid.startswith('TARGET-')}
    # Find the differences
    missing_from_processed = patient_ids - processed_patient_ids
    print('\n')
    # Log results
    if missing_from_processed:
        print(f"Missing from processed: {missing_from_processed}")
    else:
        print("All original patients are present in the processed data.")

def save_processed_patients_to_excel(processed_patients: dict, output_dir: str):
    """Save processed patient information to an Excel file."""
    df_patients = pd.DataFrame(processed_patients)
    excel_file = os.path.abspath(f'{output_dir}/processed_patients.xlsx')
    df_patients.to_excel(excel_file, index=False)

# def sort_indices(output_dir: str) -> None:
#     """
#     Sort the indices (patient IDs) of all CSV files in the 'TCGA' and 'GTEX' directories.

#     Parameters:
#     - output_dir (str): The directory containing 'TCGA' and 'GTEX' directories with CSV files to align.

#     This function will modify the CSV files in place to ensure their indices are aligned.
#     """
#     for dataset in ["TCGA", "GTEX"]:
#         dataset_dir = os.path.join(output_dir, dataset)
#         csv_files = [file for file in os.listdir(dataset_dir) if file.endswith('.csv')]
#         dataframes = {file: pd.read_csv(os.path.join(dataset_dir, file), index_col=0) for file in csv_files}
#         common_ids = None
#         for df in dataframes.values():
#             if common_ids is None:
#                 common_ids = df.index
#             else:
#                 common_ids = common_ids.intersection(df.index)
#         for file, df in dataframes.items():
#             aligned_df = df.loc[common_ids].sort_index()
#             aligned_df.to_csv(os.path.join(dataset_dir, file))
#         print(f"Sorted DataFrames in {dataset_dir}")