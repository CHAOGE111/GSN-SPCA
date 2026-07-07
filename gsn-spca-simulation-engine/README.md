# GSN-SPCA Simulation Engine

A Python-first runnable wrapper around four sparse PCA methods (ESPCA, DM-ESPCA, AWGE-ESPCA, GSN-SPCA) for comparing their performance across 10 simulated parameter sets.

## Quick Start

```powershell
pip install -r requirements.txt
```

Edit `configs/run_settings.txt` to set your local paths, then:

```powershell
# Check environment
python run.py --check-env

# Collect existing results and build report (Python only, no R needed)
python run.py --mode collect --parameter-sets all --methods all

# Run methods from scratch (requires R + reticulate)
python run.py --mode run --parameter-sets all --methods all
```

## Requirements

### Collect Mode
- Python 3.9+

### Run Mode (also requires collect-mode items)
- R 4.x with packages: `reticulate`, `readr`, `dplyr`
- Python packages: `pandas`, `scikit-learn`

```powershell
pip install -r requirements.txt
```

```r
install.packages(c("reticulate", "readr", "dplyr"))
```

See [ENVIRONMENT.md](ENVIRONMENT.md) for full setup instructions.

## Configuration

All user-facing settings are in a single file:

```text
configs/run_settings.txt
```

Edit this file to control:
- **Mode**: `collect` (gather existing outputs) or `run` (execute R methods)
- **Parameter sets**: which of the 10 simulation configs to include
- **Methods**: which of the 4 methods to run
- **Paths**: source data, Rscript, Python for reticulate
- **Seeds**: random seeds for reproducibility
- **Method parameters**: k, k_group, niter, etc.

## Structure

```
gsn-spca-simulation-engine/
├── run.py                    # Python entry point
├── requirements.txt          # pip dependencies
├── ENVIRONMENT.md            # Setup instructions
├── EXPERIMENTS.md            # 10 simulation experiments explained
├── configs/
│   ├── run_settings.txt      # User-editable control file
│   ├── default_config.json   # Internal base config
│   └── parameter_sets.csv    # 10 parameter set definitions
├── scripts/
│   ├── engine.py             # Run engine (prepares, executes, collects, reports)
│   ├── run_method.R          # R method runner
│   └── check_environment.py  # Environment validation
└── methods/
    ├── ESPCA/
    ├── DM-ESPCA/
    ├── AWGE-ESPCA/
    └── GSN-SPCA/
```

## Output

Every run creates:

```
outputs/
  run_YYYYMMDD_HHMMSS/
    run_result.html              # Full report with tables + SVG figures
    manifest.json                # Run metadata
    tables/
      method_run_summary.csv     # 40-row results table
      classified_loadings.csv    # Per-gene class labels
      default_heatmap_values.csv # Panel values for heatmaps
      stacked_bar_counts.csv     # Selected-gene counts by class
    figures/
      overview_10x4_progress_table.svg
      <parameter_set>__<method>__classes.svg
      <parameter_set>__method_class_counts.svg
      method_<method>_by_parameter_set.svg
```

## Methods

| Method | Description |
|---|---|
| **ESPCA** | Eigen-Sparse PCA with edge-based grouping |
| **DM-ESPCA** | Difference-Modulated ESPCA with t-test based power weighting |
| **AWGE-ESPCA** | Adaptive Weighted Group-Energy ESPCA |
| **GSN-SPCA** | GSN-SPCA (calls Python via reticulate for classification) |

## Parameter Sets

10 experiments organized into 4 families:

| Family | Experiments | Purpose |
|---|---|---|
| baseline | baseline_current | Reference setting |
| density | density_60_40, 70_30, 80_20, 90_10 | Varying edge network sparsity |
| noise_strength | noise_c10, noise_c20 | Varying noise magnitude |
| gene_count | genes_100, 1000, 10000 | Scalability stress tests |

## License

[MIT](LICENSE) (if applicable)
