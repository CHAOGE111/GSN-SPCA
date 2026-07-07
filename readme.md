# GSN-SPCA

## Introduction

This project is a novel sparse PCA method specifically designed for high-throughput sequencing.  
It aims to perform sparse and gene extraction in clusters to obtain more complete pathways and a more comprehensive set of genes.

- **Python Preprocessing Module**  
  Cleans expression matrix → Builds interaction networks → Identifies core gene cliques  

- **R GSN-SPCA Model**  
  Feature extraction → GSN-SPCA  



# Data Set
The following GEO datasets used in this study are accessible through the following links:
- [GSE174330](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174330)
- [GSE224449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224449)  
- [GSE34053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34053)

# Pathway Data
The gene pathway network data used in this study was obtained from the Pathway Commons database:
Pathway Commons:​ http://www.pathwaycommons.org/




## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Experiment](#simulation-experiment)
- [Batch Simulation Engine](#batch-simulation-engine)
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

## Batch Simulation Engine

The `gsn-spca-simulation-engine/` directory provides a Python-first batch benchmarking framework that compares all four sparse PCA methods across 10 systematically varied parameter configurations.

### Parameter Families

| Family | Experiments | Purpose |
|--------|------------|---------|
| baseline | `baseline_current` | Reference setting |
| density | `density_60_40`, `70_30`, `80_20`, `90_10` | Varying edge network sparsity |
| noise_strength | `noise_c10`, `noise_c20` | Varying noise magnitude |
| gene_count | `genes_100`, `1000`, `10000` | Scalability stress tests |

### Methods Compared
ESPCA · DM-ESPCA · AWGE-ESPCA · GSN-SPCA

### Estimated Runtime
~1 hour on a single CPU core (varies by machine). Pre-computed example outputs are provided in `example output.zip`.

### Quick Start

```bash
pip install -r requirements.txt
python run.py --check-env
python run.py --mode run --parameter-sets all --methods all
```

See [gsn-spca-simulation-engine/README.md](gsn-spca-simulation-engine/README.md) for full setup instructions, environment requirements, and output format details.

## File Structure

```
project_root/
│
├── GSN-SPCA/                     # Main workflow
├── Simulation experiment/        # Validation scripts
│   ├── ready/                    # Initial data and edge generation
│   └── PCAs/                     # Multiple PCA method comparisons
├── gsn-spca-simulation-engine/   # Batch benchmarking framework
│   ├── methods/                  # 4 method implementations
│   ├── configs/                  # 10 parameter set definitions
│   ├── scripts/                  # Engine & runner scripts
│   └── example output.zip        # Pre-computed results
```

## Critical Notes

- Always configure the Python path in R before running scripts:

  ```R
  reticulate::use_python("/actual/python/path")
  ```


#Special Reminder

When using the full version of R for GSN-SPCA, please note that if you need the model to classify correctly, you need to modify lines 43 and 44 of the 1.py file.

