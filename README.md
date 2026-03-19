# rmstnma: Restricted Mean Survival Time Network Meta-Analysis

The `rmstnma` package implements a Bayesian framework for conducting network meta-analysis (NMA) using Restricted Mean Survival Time (RMST) as the primary outcome.

## Key Features

- **Bayesian Inference**: Utilizes Stan (via `cmdstanr`) for robust parameter estimation and full posterior distributions.
- **Flexible Baseline Modeling**: Supports multiple baseline hazard specifications, including Royston-Parmar splines and piecewise constant models.
- **Uncertainty Propagation**: Handles uncertainty from Kaplan-Meier curve reconstruction (Guyot's method).
- **Ranking and SUCRA**: Automatically calculates treatment rankings and SUCRA scores across multiple restriction times (tau).
- **Visualization**: Built-in support for forest plots, rankograms, network diagrams, and RMST curves.

## Installation

```r
# Install from GitHub
# remotes::install_github("mahmood726-cyber/rmstnma")
```

Note: This package requires [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) to be installed.

## Quick Start

```r
library(rmstnma)

# Load example data
data(example_network)

# Create network object
net <- rmst_network(
  data = example_network,
  study = "study_id",
  trt = "treatment",
  time = "time",
  surv = "survival"
)

# Fit Bayesian NMA model
# fit <- rmst_nma(network = net, tau = 24)

# Visualize results
# plot(fit, type = "forest")
# plot(fit, type = "rank")
```

## Manuscript

A manuscript draft describing the methodology and implementation can be found in the `manuscript/` directory.

## License

MIT © Mahmood Ahmad
