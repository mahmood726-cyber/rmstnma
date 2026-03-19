# rmstnma 0.1.0

## Initial Release

This is the first release of the rmstnma package for CRAN.

### Main Features

* Network meta-analysis for restricted mean survival time (RMST) outcomes
* Bayesian inference using Stan through cmdstanr
* Support for both Kaplan-Meier and individual patient data
* Flexible baseline hazard models (piecewise and spline)
* Random effects with normal or Student-t distributions

### Key Functions

* `rmst_network()`: Create network meta-analysis dataset from survival data
* `rmst_nma()`: Fit Bayesian RMST network meta-analysis model
* `diagnose_tau()`: Diagnostic tools for selecting restriction times
* `select_tau()`: Automated tau selection using cross-validation
* `plot.rmst_nma_fit()`: Comprehensive plotting including forest plots, rankograms, and RMST curves

### Data

* `example_network`: Example dataset with three studies comparing three treatments

### Export Functions

* `export_cinema()`: Export results in CINeMA-compatible format for evidence synthesis

### Notes

* Requires optional dependency cmdstanr for Stan model fitting
* Stan models use new array syntax (requires CmdStan >= 2.33.0)
