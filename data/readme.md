## Model Directory (`./data`)

The `./data` directory contains code and utilities necessary for downloading, generating, and processing data for the DeepRBP model.

### Directory Structure

- **input_create_model**: 
    - **extra**: Contains the `get bm` script with information linking each transcript to each gene. After processing the data, a reduced version is generated with the genes and transcripts used for modeling.
    - **processed**: Includes the script for generating datasets (`create_data.sh`) and processing them for DeepRBP input. The files are found in `splitted_datasets`.
    - **raw**: Holds the download script for fetching TCGA and GTEx data (`download_script.sh`). The downloaded files are stored in this folder.
    - **selected_genes_rbps**: Tables listing the genes and RBPs used in the modeling.

- **data_real_kds**: 
    - Contains folders for experiments and tutorials explaining how to download fastq files, run kallisto, and perform differential expression analysis with voom-limma. In the `experiments` folder, for each of the 4 experiments, the following information is saved:
        - A folder for each sample with its abundances.
        - A `limma` folder where a CSV table with the differential expression analysis results should be saved.
        - A `datasets` folder with subfolders for control and knockdown samples processed for DeepRBP input.
        - A TXT file describing each sample as either knockdown or control.

- **data_postar3**: 
    - Contains an `input_data` folder with:
        - The general postar file detailing RBP binding sites in the genome.
        - `EventsFound_gencode23.txt` with event information.
        - `Events_Regions_gc23_400nt.txt` with event region details.
        - A shell script `process_script.sh` that processes the postar matrix. This script generates:
            - A postar file for each tissue in the `./results` folder.
            - A TXT file describing each tissue or cell line.
    - Contains a `create_gxrbp` folder with:
        - An R script with the same name. This script, using:
            - Specific Postar3 matrices for each tissue.
            - Event information.
        - Creates the binary RBPxgene matrix of Postar3 with regulations. These results can be found in `./results`.
