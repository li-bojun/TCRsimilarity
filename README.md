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