# Environment Requirements

This package is Python-first, but the full method runner also calls R and, for GSN-SPCA, R calls Python through `reticulate`.

There are two practical levels of environment support.

## 1. Collect Mode

Use this when you want to collect existing final `U/D/V` outputs and build the HTML report, CSV tables, class visualizations, and stacked bar charts.

Required:

- Windows or another system that can access the paths configured in `configs/run_settings.txt`.
- Python 3.9 or newer.
- No third-party Python packages are required for `collect` mode.
- Existing result folders and summary CSV configured by:
  - `parameter_set_collect_summary`
  - `source_root`

Default command:

```powershell
python .\run.py --mode collect --parameter-sets all --methods all
```

## 2. Run Mode

Use this when you want to start the R methods from copied parameter-set inputs.

Required:

- Python 3.9 or newer.
- R 4.x with `Rscript.exe` available.
- R packages:
  - `reticulate` (required — bridges R to Python)
  - `readr` (optional, used by legacy helper)
  - `dplyr` (optional, used by legacy helper)
- Python packages for the GSN-SPCA classification helper:
  - `pandas`
  - `scikit-learn`

Install Python packages:

```powershell
pip install -r requirements.txt
```

Install R packages:

```r
install.packages(c("reticulate", "readr", "dplyr"))
```

## Configuration

Set these paths in `configs/run_settings.txt`:

```ini
[paths]
# Root directory containing all simulation data (e.g., F:\hggz)
source_root = <SOURCE_ROOT>

# Full path to Rscript.exe (e.g., C:\Program Files\R\R-4.4.1\bin\Rscript.exe)
rscript = <RSCRIPT_PATH>

# Python interpreter used by R's reticulate package.
# Must have pandas and scikit-learn installed.
reticulate_python = <RETICULATE_PYTHON_PATH>

# Parent directory containing one subfolder per parameter set.
# Each subfolder: matrix_x.txt, edges_1based.txt, metadata.json, etc.
parameter_set_input_root = <PARAMETER_SET_INPUT_ROOT>

# CSV summary from a previous run, used by collect mode.
parameter_set_collect_summary = <PARAMETER_SET_COLLECT_SUMMARY>
```

### reticulate_python: Important

When GSN-SPCA runs, R's `reticulate` invokes the Python interpreter specified by `reticulate_python`. There are three ways to handle this:

1. **(Recommended) Share the same Python** — Set `reticulate_python` to the same Python interpreter used to launch `run.py` (i.e. the one from which you ran `pip install -r requirements.txt`). This avoids version/package mismatches between R's Python and the engine's Python.
2. **Dedicated Python environment** — Create a separate conda/virtualenv for `reticulate`, install `pandas` and `scikit-learn` there, and point `reticulate_python` to it.
3. **Environment variable fallback** — If `reticulate_python` is left empty, `fun_GSN-SPCA.R` falls back to `Sys.getenv("RETICULATE_PYTHON")`. You can set this environment variable before running the engine.

In all cases, the Python pointed to must have `pandas` and `scikit-learn` installed.

## Environment Check

Run:

```powershell
python .\run.py --check-env
```

or:

```powershell
python .\scripts\check_environment.py --settings .\configs\run_settings.txt
```

The check reports:

- Current Python version.
- Whether configured paths exist.
- Whether `Rscript` is available.
- Whether required Python modules are installed.
- Whether R packages are installed when `Rscript` is available.

## Reproducibility Notes

The random seeds used by the experiment are declared in `configs\run_settings.txt`. The package keeps these seeds in generated manifests, but exact numerical reproduction also depends on compatible R/Python package versions and on the same input files.
