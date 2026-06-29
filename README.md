# TCR Similarity Project
This repository is code for a project investigating how various representations of TCR similarity encode biological signals. It handles data pre-processing (ambient RNA removal, doublet detection, merger of scRNA with scVDJ data), multi-sample integration, as well as defining three different similarity representations to investigate.

## Installation and Requirements
This project requires R for SoupX ambient RNA removal (`runsoupx.r`) and Python for everything else. To install, clone the repository, and create two environments:

For `process_singlesample_TCR.ipynb` the required Python environment can be installed via:
```
conda env create -f singlesample_env.yml
```

For all the other `.ipynb` files, the required Python environment can be installed via:
```
conda env create -f multisample_env.yml
```

You will also need CellRanger 10x 5' scRNA and scVDJ objects to analyze. Each object should be of a single patient.

## Usage

### Pre-processing
The preprocessing comes in three steps. 

Firstly, run `runsoupx.r` in your R environment of choice. This handles ambient RNA removal, part of the quality control steps in data preprocessing.

After ambient RNA removal is complete, each patient sample can be processed (such as doublet removal and mitochondrial cutoffs) via running `process_singlesample_TCR.ipynb`; this will produce a `.h5ad` file for that patient. Each new sample requires a full execution of the notebook from start to finish.

Once all the samples are individually processed, `process_multisample_TCR.ipynb` can be ran once from start to finish to integrate all the patient samples' data contained in the generated `.h5ad` files. The notebook can also annotate and subset to T cells only, and save a new `.h5ad` object containing only T cells.