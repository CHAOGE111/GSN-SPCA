# GSN-SPCA

## Introduction

This project is a novel sparse PCA method specifically designed for high-throughput sequencing.  
It aims to perform sparse and gene extraction in clusters to obtain more complete pathways and a more comprehensive set of genes.

- **Python Preprocessing Module**  
  Cleans expression matrix → Builds interaction networks → Identifies core gene cliques  

- **R GSN-SPCA Model**  
  Feature extraction → GSN-SPCA  



# Data Set
• The following GEO datasets used in this study are accessible through the following links:
GSE174330: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174330
GSE224449: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224449
GSE34053: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34053

# Pathway Data
The gene pathway network data used in this study was obtained from the Pathway Commons database:
Pathway Commons:​ http://www.pathwaycommons.org/




## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Experiment](#simulation-experiment)
- [File Structure](#file-structure)
- [Critical Notes](#critical-notes)


## Installation

### Python Requirements
```bash
pip install pandas numpy tqdm scikit-learn ast
```

### R Requirements

```R
install.packages("reticulate")
```

## Usage

### Input Files for R GSN-SPCA

- `result-1_p2.txt` : Expression matrix
- `gene_new_o.txt` : Gene-gene relationships
- `data.csv` : Sample labels

### Steps

1. For the **simulation experiment**, first run the files in the `ready` folder to obtain the initial data matrix and edge relationships.
2. Run the **Python for GSN-SPCA** scripts to construct the clique relationships.
3. Input the matrix and clique/edge relationships into the various methods in the `PCAs` folder for operation and testing.
4. For the `GSN-SPCA` folder, first use **Python for GSN-SPCA** to construct the clique relationships, then replace the input files in `R for GSN-SPCA` and run.

## Simulation Experiment

The simulation experiment validates the GSN-SPCA model:

- Generates initial data matrix and edge relationships (`ready` folder)
- Constructs clique relationships (Python preprocessing)
- Inputs data into multiple PCA methods (`PCAs` folder) for comparison
- Runs GSN-SPCA workflow with updated inputs (`R for GSN-SPCA` folder)

## File Structure

```
project_root/
│
├── GSN-SPCA_model/         # Main workflow
├── simulation_experiment/  # Validation scripts
│   └── ready/              # Initial data and edge generation
```

## Critical Notes

- Always configure the Python path in R before running scripts:

  ```R
  reticulate::use_python("/actual/python/path")
  ```


#Special Reminder

When using the full version of R for GSN-SPCA, please note that if you need the model to classify correctly, you need to modify lines 43 and 44 of the 1.py file.

