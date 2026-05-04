# Spec4ML for R

Spec4ML for R is a package for handling and analyzing spectral data with special considerations for machine-learning and chemometric workflows.

## Documentation

Complete package documentation is available here:

- [Full documentation](docs/README.md)

The full documentation covers:

- installation,
- data expectations,
- import helpers,
- spectral preprocessing,
- EMSC normalization,
- PLS/regression cross-validation,
- technical replicate guidance,
- leakage prevention,
- troubleshooting,
- relationship to the Python package and Studio app.

## Installation

```r
install.packages("remotes")
remotes::install_github("elkadi/spec4ml")
```

```r
library(spec4ml)
```

## Related repositories

- Python package: `https://github.com/elkadi/spec4ml_py`
- Streamlit app: `https://github.com/elkadi/SpecML-Studio`
