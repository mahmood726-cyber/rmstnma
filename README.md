# rmstnma: Restricted Mean Survival Time Network Meta-Analysis

The `rmstnma` package implements a Bayesian framework for conducting network meta-analysis (NMA) using Restricted Mean Survival Time (RMST) as the primary outcome.

## Key Features

- **Bayesian Inference**: Stan (via `cmdstanr`) for posterior estimation of RMST contrasts.
- **Flexible Baseline Modeling**: Royston–Parmar splines and piecewise-constant hazards (with exponential as a sensitivity case).
- **Reconstruction-uncertainty hook**: API for propagating Kaplan–Meier reconstruction error from Guyot et al. (2012). The current helper detects reconstructed arms but requires arm-level draws; see Limitations.
- **Ranking and SUCRA**: Treatment rankings and SUCRA scores reported across multiple restriction times (τ).
- **Visualization**: Forest plots, rankograms, network diagrams, and RMST curves.

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

## Methods

`rmstnma` fits arm-based Bayesian network meta-analysis on the RMST scale. The Stan model (`inst/stan/rmst_arm_based.stan`) parameterises the survival function over time and integrates to the user-specified restriction time τ. Three baseline-hazard specifications are supported:

- **Royston–Parmar splines** on the log-cumulative-hazard scale (configurable knot count).
- **Piecewise-constant** hazards on a user-defined or data-driven time grid.
- **Exponential** as a sensitivity special case.

Random study effects on baseline are on by default; random treatment effects are optional. Setting `inconsistency = TRUE` adds design-by-treatment interaction terms so closed-loop networks can be checked.

When individual-patient times are not available, Kaplan–Meier curves are reconstructed using the Guyot et al. (2012) algorithm and entered as a reconstructed arm. The current `add_reconstruction_uncertainty()` helper is a stub: it detects reconstructed arms but requires arm-level reconstructed survival draws that the helper does not yet generate, so fits with `propagate_recon_uncertainty = TRUE` will stop with an explanatory message until that path is wired through.

`diagnose_tau()` returns a maturity grid (`mean_survival`, `min_survival`, `prop_mature`) over 20 τ values between 20%–90% of the minimum study follow-up. Pick τ from data, not from convenience.

## Limitations

- **CmdStan dependency.** Inference goes through `cmdstanr`; CmdStan must be installed separately. The package errors out cleanly if it cannot find a compiled model.
- **Consistency assumption.** Arm-based NMA assumes consistency. Closed-loop diagnostics are exposed but not run automatically — comparing `inconsistency = FALSE` and `inconsistency = TRUE` fits is the user's responsibility.
- **Reconstruction uncertainty is not yet propagated.** See the note above; current behaviour is to error out rather than silently absorb the error.
- **Single outcome per network.** Surrogate or composite-outcome mixing must be handled before calling `rmst_network()`.
- **Computational cost.** Default 2 chains × 2000 iterations on a 10-treatment network takes minutes on a laptop and scales with network size × τ-grid points.

## Conclusions

Use `rmstnma` when (a) treatment effects on survival are clearly time-varying so a single HR is misleading, (b) RMST at a clinically meaningful τ is the intended decision quantity, and (c) at least one arm has either IPD or reconstructable KM curves. For networks where proportional hazards is defensible and τ is not a decision variable, HR-based NMA remains simpler and runs faster.

## Manuscript

A manuscript draft describing the methodology and implementation is in the `manuscript/` directory.

## License

MIT © Mahmood Ahmad
