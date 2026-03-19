# rmstnma: Operations & Workflow

This document outlines the workflow for conducting Bayesian Restricted Mean Survival Time (RMST) Network Meta-Analysis using the `rmstnma` package.

## Environment Setup

The package requires a working installation of **CmdStan**.

1.  **Install CmdStanR**:
    ```r
    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    ```
2.  **Install CmdStan**:
    ```r
    library(cmdstanr)
    check_cmdstan_toolchain()
    install_cmdstan(cores = 4)
    ```

## Standard Workflow

### 1. Data Ingestion and Networking
Use `rmst_network()` to create a standardized network object from survival data (Kaplan-Meier or IPD).

```r
data(example_network)
net <- rmst_network(
  data = example_network,
  study = "study_id",
  trt = "treatment",
  time = "time",
  surv = "survival"
)
```

### 2. Bayesian Model Fitting
Fit the NMA model using `rmst_nma()`. Choose an appropriate restriction time ($\tau$).

```r
fit <- rmst_nma(
  network = net,
  tau = 24,           # Restriction time
  baseline = "rp_spline", # Royston-Parmar splines
  chains = 4,
  iter = 2000
)
```

### 3. Convergence Diagnostics
Check for MCMC convergence using trace plots and R-hat values.

```r
# Trace plots
plot(fit, type = "trace", pars = c("d[2]", "sigma"))

# Detailed summary
fit$fit$summary()
```

### 4. Results and Ranking
Extract treatment contrasts and calculate SUCRA scores.

```r
# Forest plot
plot(fit, type = "forest")

# Rankogram and SUCRA
plot(fit, type = "rank")
```

## Sensitivity Analysis

It is critical to evaluate the sensitivity of results to the choice of $\tau$.

```r
# Fit across multiple tau values
fit_multi <- rmst_nma(network = net, tau = c(12, 24, 36))

# Plot RMST curves
plot(fit_multi, type = "rmst")
```

## Troubleshooting

- **Divergent Transitions**: If you encounter divergent transitions, try increasing `adapt_delta` (e.g., `adapt_delta = 0.99`) and `iter`.
- **Baseline Model Selection**: If spline models fail to converge, try the simpler `baseline = "piecewise"` model.
- **CmdStan Path**: Ensure `cmdstanr` can find CmdStan using `cmdstanr::set_cmdstan_path()`.
