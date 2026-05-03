# spec4ml R Package Documentation

`spec4ml` is an R package for importing, preprocessing, normalizing, and modeling spectral data for machine-learning and chemometric workflows.

The package is intended for spectroscopy datasets where rows represent spectra or samples, non-spectral columns contain metadata or targets, and numeric spectral columns contain wavelengths, wavenumbers, or derived spectral features.

## Repository

- Package: `spec4ml`
- Language: R
- Repository: `https://github.com/elkadi/spec4ml`
- License: GPL >= 2
- Minimum R version: 3.6.0

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("elkadi/spec4ml")
```

Load the package:

```r
library(spec4ml)
```

## Main dependency groups

The package imports tools for spectral data handling, preprocessing, and modeling:

- `hyperSpec` for spectral data objects.
- `prospectr`, `pracma`, and package functions for preprocessing.
- `EMSC` for Extended Multiplicative Signal Correction workflows.
- `caret` for model training and cross-validation helpers.
- `dplyr`, `fs`, `rlang`, and `utils` for data manipulation and file operations.

## Typical workflow

A common workflow is:

1. Import spectral files or merged spectral datasets.
2. Identify metadata, target, and spectral feature columns.
3. Convert spectra into a suitable R object, commonly a `hyperSpec` object where relevant.
4. Apply one or more preprocessing strategies.
5. Save or select preprocessed spectra for modeling.
6. Run regression or PLS cross-validation.
7. Compare preprocessing/model combinations using validation metrics.

Example skeleton:

```r
library(spec4ml)

# 1. Import a semicolon-separated spectral file
spectra_df <- SpecImp("raw_spectra.csv", y = 8, delimiter = ";")

# 2. Convert to or prepare a hyperSpec object according to your project layout
# raw_hyper <- hyperSpec::read.*(...) or project-specific conversion

# 3. Generate preprocessing variants
# pp <- SpectraPreProcess(raw_hyper)

# 4. Evaluate regression models
# result <- RegCV(M = raw_hyper, V = "target_column", comp = 20, model = "pls")
```

The exact object-construction step depends on how your instrument export is represented.

## Data expectations

Most workflows assume:

- metadata columns before spectral columns,
- a target variable for modeling,
- numeric spectral columns representing wavelengths/wavenumbers or spectral channels,
- sample identifiers when repeated spectra or technical replicates are present.

Functions that operate on `hyperSpec` objects expect valid `hyperSpec` inputs. File-import helpers expect CSV-style instrument exports.

## Exported functions

The package exports the following functions.

### Import and file-management helpers

| Function | Purpose |
|---|---|
| `SpecImp(x, y = 8, delimiter = ";")` | Imports a delimited CSV spectral file after skipping header lines. |
| `avasoft_import(...)` | Imports Avantes/AvaSoft-style spectral exports. |
| `Import_and_Merge_Spectra(...)` | Imports and merges spectra across files. |
| `Merge_subfiles(...)` | Merges subfiles into a combined dataset. |
| `copy_from_subfolders(...)` | Copies files from nested folders into a target structure. |
| `Add_FolderName_to_Files(...)` | Adds folder-name metadata to files. |
| `Add_NestedFolderNames_to_Files(...)` | Adds nested folder-name metadata. |
| `Add_NestedFolderNames_to_Files2(...)` | Alternate nested-folder metadata helper. |
| `tnName(...)` | Utility for filename or target-name handling. |

### Spectral preprocessing and normalization

| Function | Purpose |
|---|---|
| `SpectraPreProcess(...)` | Generates a broad list of preprocessing variants for a `hyperSpec` object. |
| `SpectraPreProcessE(...)` | Generates preprocessing variants including EMSC-based preprocessing; use cautiously in CV because EMSC may introduce leakage if references are computed across folds. |
| `SpectraPreProcessOptimized(...)` | Optimized preprocessing workflow. |
| `SpectraSmoothing(...)` | Applies smoothing to spectra. |
| `emsc_normalize(...)` | Applies EMSC normalization to NIR spectra using grouped control spectra. |
| `emsch(...)` | EMSC helper. |
| `snvh(...)` | Standard normal variate helper. |
| `uvnormalize(...)` | Unit-vector or UV-style normalization helper. |
| `find_first_numeric_column(...)` | Detects the first spectral/numeric column. |

### Modeling and validation

| Function | Purpose |
|---|---|
| `PLSCV(M, V, ncompMax = 20, validation = "repeatedcv")` | Performs PLS regression cross-validation for a `hyperSpec` object and target. |
| `RegCV(M, V, comp = 20, model = "pls", validation = "repeatedcv")` | Performs regression cross-validation using a selected model type. |

## Important function notes

### `SpecImp()`

```r
spectra <- SpecImp(x = "spectra.csv", y = 8, delimiter = ";")
```

Use this for CSV files that contain non-data header rows before the table starts.

Arguments:

- `x`: CSV file path.
- `y`: number of lines to skip before reading the data table.
- `delimiter`: column delimiter; default is `;`.

Returns a data frame.

### `SpectraPreProcess()`

```r
pp <- SpectraPreProcess(
  rawspectra,
  SmLn = 3,
  SmMn = 5,
  SmHn = 11,
  ExtraDerivativeLSm = 15,
  ExtraDerivativeHSm = 17
)
```

Creates a list of many preprocessed versions of the same spectral dataset. The package documentation describes the return value as a list of 66 preprocessed spectra plus their names.

Use this when comparing multiple preprocessing strategies before machine-learning evaluation.

### `SpectraPreProcessE()`

```r
pp_e <- SpectraPreProcessE(rawspectra, SmLn = 3, SmMn = 5, SmHn = 11)
```

Creates preprocessing variants including EMSC-related variants. The package documentation describes the return value as a list of 74 preprocessed spectra plus their names.

Caution: EMSC preprocessing can cause data leakage if reference spectra are computed using information across training and validation folds. For strict cross-validation, prefer `SpectraPreProcess()` unless you explicitly control the fold-wise reference calculation.

### `emsc_normalize()`

```r
normalized <- emsc_normalize(
  spectra_file = "NIR_spectra.csv",
  controls_file = "Control_spectra.csv",
  output_file = "NIR_spectra.csv",
  group_col = "Day",
  spectra_dir = "../InputSpectra",
  output_dir = "../Normalized_Spectra",
  save_output = FALSE
)
```

This workflow:

1. reads spectral and control files,
2. identifies the first numeric spectral column,
3. computes grouped mean control/reference spectra,
4. applies EMSC normalization per group,
5. preserves metadata columns,
6. optionally writes a normalized CSV file.

Use `group_col` for batch/day/run-specific controls.

### `PLSCV()`

```r
pls_result <- PLSCV(M = spectra_hyper, V = "target", ncompMax = 20)
```

Performs cross-validation for PLS regression. The package documentation currently notes that this function needs further revision, so validate outputs carefully before publication or production use.

### `RegCV()`

```r
reg_result <- RegCV(M = spectra_hyper, V = "target", comp = 20, model = "pls")
```

Performs regression cross-validation using a selected model. The default model is `pls`. The package documentation currently notes that this function needs further revision.

## Technical replicates

The R package contains several utilities useful for spectral preprocessing and modeling, but it does not currently expose the same higher-level replicate-handling API used in the Python/Studio workflow, such as named modes for:

- no replicate handling,
- averaging spectra before modeling,
- averaging predictions after modeling.

For R workflows, handle technical replicates explicitly by:

1. keeping a sample/group identifier column,
2. averaging spectra by group before model training when appropriate,
3. using group-aware validation designs where possible,
4. avoiding train/test leakage between technical replicates of the same biological or physical sample.

Example group averaging:

```r
library(dplyr)

metadata_cols <- c("sample_id", "target")
spectral_cols <- setdiff(names(df), metadata_cols)

averaged <- df %>%
  group_by(sample_id) %>%
  summarise(
    target = mean(target, na.rm = TRUE),
    across(all_of(spectral_cols), ~ mean(.x, na.rm = TRUE)),
    n_replicates = dplyr::n(),
    .groups = "drop"
  )
```

For classification targets, do not average labels blindly. First check that all replicate rows in each group have the same label.

## Validation and leakage guidance

Spectral preprocessing can introduce leakage when parameters are learned from the full dataset before cross-validation. This is especially important for normalization, EMSC, feature selection, scaling, and any preprocessing that estimates dataset-level statistics.

Recommended practice:

- Split data by sample/group before computing fold-specific preprocessing parameters.
- Keep technical replicates from the same sample in the same fold.
- Use preprocessing functions that do not learn from held-out test data, or recompute preprocessing inside each training fold.
- Report whether metrics are row-level, spectrum-level, or sample-level.

## Output conventions

For reproducible analyses, save:

- raw input files,
- preprocessing settings,
- selected preprocessing names,
- model settings,
- validation splits,
- prediction tables,
- metrics tables,
- package versions and `sessionInfo()`.

Example:

```r
sessionInfo()
packageVersion("spec4ml")
```

## Troubleshooting

### `hyperSpec` object errors

Confirm that your object is a valid `hyperSpec` object and that spectral data are stored in the expected slot/format.

### Missing or non-numeric spectral columns

Use `find_first_numeric_column()` or inspect column names manually. Spectral columns should be numeric data, even if wavelength column names are character strings.

### EMSC gives unexpectedly high validation scores

Check for leakage. If EMSC references were computed using the complete dataset before validation, results may be optimistic.

### Replicates appear in both training and testing folds

Create folds by sample/group identifier rather than by individual spectrum rows.

## Relationship to the Python package and Studio

This repository is the R package. The Python package lives separately at `elkadi/spec4ml_py`. The Streamlit application lives separately at `elkadi/SpecML-Studio` and depends on the Python package.

Do not assume every feature exists in all three repositories. Keep package-specific documentation and examples synchronized only where the implementation actually matches.

## Documentation maintenance checklist

When adding or changing functions:

1. Update roxygen comments in the corresponding `R/*.R` file.
2. Run `devtools::document()` to regenerate `man/*.Rd` and `NAMESPACE`.
3. Add or update examples in this document.
4. Add a small reproducible example where possible.
5. Note any leakage-sensitive behavior clearly.
