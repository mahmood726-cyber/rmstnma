# powerNMA: Comprehensive Network Meta-Analysis Suite

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)

> A powerful, unified R package combining time-varying NMA methods with advanced meta-analysis techniques

##  Features

### Core Capabilities

- **Standard Network Meta-Analysis** - Frequentist and Bayesian approaches using `netmeta` and `gemtc`
- **Time-Varying Methods** - RMST and milestone survival analyses for time-to-event data
- **Transportability** - Advanced weighting schemes for target population inference
- **Publication Bias** - PET-PEESE, selection models, trim-and-fill
- **Meta-Regression** - With spline smoothing and cluster-robust inference
- **Sensitivity Analyses** - Leave-one-out, leave-one-treatment-out, multiverse analysis
- **ML Heterogeneity** - Machine learning exploration of effect modifiers
- **Interactive Dashboard** - Comprehensive Shiny interface

### Methodological Innovations

1. **Transportability Weighting**
   - Mahalanobis and Euclidean distance metrics
   - Gaussian and tricube kernels
   - GRADE-informed and design-adjusted weighting

2. **Time-Varying Effects**
   - Restricted Mean Survival Time (RMST) NMA
   - Milestone survival analysis with continuity correction
   - Trajectory visualization

3. **Comprehensive Diagnostics**
   - Global and local inconsistency testing
   - Effective sample size calculations
   - Network geometry metrics
   - Convergence diagnostics for Bayesian models

### Advanced Methods (Phases 12-14)

#### Phase 12: Advanced Bayesian & External Integration
- **Stan-based Bayesian NMA** - Full HMC/NUTS sampling with comprehensive diagnostics
- **IPD Network Meta-Analysis** - One-stage and two-stage approaches with personalized predictions
- **Diagnostic Test Accuracy NMA** - Bivariate models, HSROC curves, sensitivity/specificity
- **Interactive HTML Reports** - Self-contained R Markdown reports with Plotly/DT
- **External API Integration** - Automated PubMed and ClinicalTrials.gov data retrieval
- **Advanced Meta-Regression** - Restricted cubic splines and fractional polynomials

#### Phase 13: Cutting-Edge Statistical Methods (2024-2025)
- **Multivariate NMA** - Joint efficacy + safety modeling with within-study correlations
- **Benefit-Risk MCDA** - Multi-criteria decision analysis with SMAA
- **ROPE Analysis** - Region of Practical Equivalence (Bayesian and frequentist)
- **TOST & Equivalence** - Two One-Sided Tests for non-inferiority
- **Transitivity Testing** - Network coherence and assumption validation
- **RWE Integration** - RCT + real-world evidence synthesis with power priors
- **Sequential NMA** - Trial Sequential Analysis and living network meta-analysis
- **Treatment Selection** - Precision medicine algorithms for optimal treatment choice

#### Phase 14: Deep Learning & Causal Inference (REVOLUTIONARY!)
- **Graph Neural Networks** - GNN/GCN/GAT/GraphSAGE for network structure learning
  - Treatment embeddings and representation learning
  - Link prediction for missing comparisons
  - Transfer learning across disease areas
  - GNN-based treatment clustering
- **Causal Inference Framework** - Rigorous causal effect estimation
  - G-formula for causal effects
  - Targeted Maximum Likelihood Estimation (TMLE)
  - Doubly robust estimation
  - E-values for unmeasured confounding sensitivity
- **Distributional NMA** - Beyond mean effects to full distributions
  - Quantile regression network meta-analysis
  - Distribution reconstruction from summary statistics
- **Competing Risks NMA** - Multi-event time-to-event analysis
  - Cause-specific hazards models
  - Fine-Gray subdistribution hazards
  - Multi-state transition models
- **Special Trial Designs** - Crossover and N-of-1 trials
  - Within-patient correlation modeling
  - Carryover effect adjustment
  - N-of-1 trial aggregation
- **Finite Mixture Models** - Subgroup discovery without pre-specification
  - Latent class analysis for treatment response
  - Growth mixture models for longitudinal outcomes
  - Responder vs non-responder identification

#### Phase 15: Production-Grade Quality & Performance (ENTERPRISE!)
- **Comprehensive Input Validation** - Defensive programming best practices
  - assert_not_null(), assert_data_frame(), assert_numeric() validators
  - assert_character(), assert_logical(), assert_nma_object() type checkers
  - assert_probability(), assert_positive_integer() specialized validators
  - assert_network_structure() for network data validation
  - Standardized error messages with function context
  - Safe try-catch wrappers with error recovery
- **Performance Optimization Infrastructure** - Production-grade performance
  - Intelligent memoization with automatic cache management
  - Parallel processing with automatic core detection
  - Progress tracking for long-running operations
  - Batch processing for large datasets
  - Memory-efficient data structures
  - Execution time profiling and benchmarking
  - System resource monitoring
- **Enterprise Quality Assurance** - Comprehensive testing
  - 25+ edge case tests (empty data, single study, extreme values)
  - Disconnected network detection
  - Large network stress testing (1000+ studies)
  - Negative value validation
  - Memory usage estimation
  - Integration testing across modules

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("your-org/powerNMA")

# Or install with all optional dependencies
devtools::install_github("your-org/powerNMA", dependencies = TRUE)
```

###  Required Packages

```r
install.packages(c(
  "netmeta", "survRM2", "survival",
  "ggplot2", "dplyr", "tidyr", "purrr", "tibble",
  "shiny", "DT", "plotly", "readr"
))
```

###  Optional (for full features)

```r
install.packages(c(
  "gemtc", "coda",        # Bayesian NMA
  "multinma",             # Advanced Bayesian methods
  "metafor",              # Meta-regression
  "clubSandwich",         # Robust inference
  "ranger",               # ML heterogeneity
  "weightr", "metasens",  # Publication bias
  "rmarkdown", "knitr"    # Reports
))
```

###  For Enterprise bs4Dash Dashboard

```r
install.packages(c(
  "bs4Dash",              # Bootstrap 4 dashboard framework
  "shinyWidgets",         # Enhanced UI widgets
  "waiter",               # Loading screens and progress bars
  "fresh",                # Custom themes
  "htmlwidgets",          # Interactive HTML widgets
  "webshot2"              # High-resolution static exports from HTML
))
```

## Quick Start

### Example 1: Standard Network Meta-Analysis

```r
library(powerNMA)

# Generate example data
data <- simulate_nma_data(n_studies = 40)

# Run comprehensive analysis
results <- run_powernma(
  data = data,
  data_type = "pairwise",
  config = setup_powernma(
    sm = "HR",
    use_bayesian = TRUE,
    export_results = TRUE
  )
)

# View results
print(results)
summary(results)

# Access components
results$results$main_nma        # Main NMA
results$results$ranking         # Treatment ranking
results$results$loo             # Leave-one-out
```

### Example 2: Time-Varying Analysis with IPD

```r
# Generate IPD
ipd <- generate_example_ipd(n_trials = 10, n_per_arm = 100)

# Configure for time-varying analysis
config <- setup_powernma(
  use_timevarying = TRUE,
  tau_list = c(90, 180, 365),
  milestone_times = c(90, 180, 365),
  export_plots = TRUE
)

# Run analysis
results <- run_powernma(
  data = ipd,
  data_type = "ipd",
  config = config
)

# RMST results
print(results$results$rmst_nma)

# Milestone results
print(results$results$milestone_nma)
```

### Example 3: Advanced Analysis with Transportability

```r
# Your trial data
data <- simulate_nma_data(n_studies = 50)

# Define target population
target_pop <- list(
  age_mean = 70,
  female_pct = 0.50,
  bmi_mean = 29,
  charlson = 2.0
)

# Comprehensive configuration
config <- setup_powernma(
  sm = "HR",
  use_transport = TRUE,
  use_bayesian = TRUE,
  run_sensitivity = TRUE,
  run_metareg = TRUE,
  export_results = TRUE,
  export_plots = TRUE
)

# Run full analysis
results <- run_powernma(
  data = data,
  data_type = "pairwise",
  target_population = target_pop,
  config = config
)
```

##  Interactive Dashboard

### Standard Dashboard (shinydashboard)

Launch the Shiny dashboard for interactive exploration:

```r
# Standard NMA
data <- simulate_nma_data()
results <- run_powernma(data, data_type = "pairwise")
launch_powernma_app(results)

# IPD with time-varying methods
ipd <- generate_example_ipd()
results <- run_powernma(ipd, data_type = "ipd",
                        config = setup_powernma(use_timevarying = TRUE))
launch_powernma_app(results)
```

### Enterprise bs4Dash Dashboard (NEW!)

**Modern, professional Bootstrap 4 interface with high-resolution export capabilities:**

```r
# Launch enterprise dashboard with high-resolution exports
launch_bs4dash_app(
  theme = "default",        # "default", "dark", or "light"
  high_res_dpi = 600,       # Up to 1200 DPI for publication-quality
  max_upload_mb = 100
)

# Dark theme with maximum resolution
launch_bs4dash_app(
  theme = "dark",
  high_res_dpi = 1200       # Maximum resolution for publications
)
```

**Key Features:**
- Modern Bootstrap 4 UI with responsive design
- High-resolution plot downloads (300-1200 DPI)
- Multiple export formats: PNG, PDF, SVG, TIFF, HTML
- Publication-quality outputs suitable for journals
- Integrated Phase 12-14 advanced methods
- Real-time progress tracking with waiter animations
- Customizable export dimensions and resolution
- Batch export functionality

**High-Resolution Downloads:**
Every visualization includes download buttons for:
- PNG at 300, 600, and 1200 DPI
- PDF (vector graphics)
- SVG (vector graphics)
- TIFF (lossless compression)
- HTML (interactive)

Configure export settings globally via the control bar (gear icon in top-right) or per-plot using the download buttons.

##  Data Format

### Pairwise Data (Standard NMA)

```r
# Required columns:
data.frame(
  studlab = c("Study_01", "Study_01", "Study_02"),
  treat1 = c("Control", "Control", "DrugA"),
  treat2 = c("DrugA", "DrugB", "DrugB"),
  TE = c(-0.15, -0.25, -0.10),      # Treatment effect
  seTE = c(0.10, 0.12, 0.09)        # Standard error
)

# Optional covariates for meta-regression:
# age_mean, female_pct, bmi_mean, charlson, is_rct, grade, etc.
```

### Individual Patient Data (IPD)

```r
# Required columns:
data.frame(
  trial = c("Trial_01", "Trial_01", "Trial_01"),
  treatment = c("Control", "Control", "DrugA"),
  time = c(150, 200, 180),          # Follow-up time
  status = c(1, 0, 1)               # Event indicator (1=event, 0=censored)
)
```

##  Configuration Options

```r
config <- setup_powernma(
  # Effect measure
  sm = "HR",                       # "HR", "OR", "RR", "MD", "SMD"

  # Methods
  use_bayesian = TRUE,             # Run Bayesian NMA
  use_timevarying = TRUE,          # Run RMST/milestone (requires IPD)
  run_sensitivity = TRUE,          # LOO, LOTO analyses
  run_metareg = TRUE,              # Meta-regression

  # Time-varying settings
  tau_list = c(90, 180, 365),      # RMST time horizons
  milestone_times = c(90, 180, 365),  # Milestone time points

  # Transportability
  use_transport = TRUE,            # Enable transportability weighting

  # Export
  export_results = TRUE,
  export_plots = TRUE,
  output_dir = "results",
  plot_dir = "plots",

  # Reproducibility
  seed = 42
)
```

##  Output Structure

```r
results <- run_powernma(...)

# Main components:
results$results$main_nma          # netmeta object
results$results$rmst_nma          # RMST analyses (if IPD)
results$results$milestone_nma     # Milestone analyses (if IPD)
results$results$bayesian          # Bayesian results
results$results$ranking           # Treatment ranking
results$results$loo               # Leave-one-out sensitivity
results$results$geometry          # Network metrics
results$results$global_inconsistency  # Design-by-treatment test
results$results$local_inconsistency   # Node-splitting

results$summary                   # Summary statistics
results$config                    # Configuration used
results$ref_treatment             # Reference treatment
```

##  Development Tools

```r
# Validate data
validate_ipd(ipd_data)
validate_nma_input(pairwise_data)

# Clean data
clean_data <- clean_nma_data(raw_data)

# Network geometry
geometry <- network_geometry(data)

# Clear cache
clear_powernma_cache()
```

##  Visualization

```r
# Coming soon: comprehensive plotting functions
# plot(results, type = "forest")
# plot(results, type = "network")
# plot(results, type = "sucra")
# plot(results, type = "trajectory")  # For time-varying
```

##  Examples

See `vignettes/` for detailed examples:
- `introduction.Rmd` - Getting started
- `time_varying.Rmd` - RMST and milestone analysis
- `advanced.Rmd` - Transportability and meta-regression
- `bayesian.Rmd` - Bayesian methods

##  Citation

If you use powerNMA in your research, please cite:

```
@misc{powernma2025,
  title = {powerNMA: Comprehensive Network Meta-Analysis Suite},
  author = {NMA Research Group},
  year = {2025},
  url = {https://github.com/your-org/powerNMA}
}
```

##  Contributing

Contributions are welcome! Please see `CONTRIBUTING.md` for guidelines.

##  License

MIT License - see `LICENSE` file

##  References

### Time-Varying Methods
- Guyot P, et al. (2012). Enhanced secondary analysis of survival data. *BMC Med Res Methodol* 12:9
- Wei Y, Royston P (2017). Reconstructing time-to-event data from published KM curves. *Stata J* 17(4):786-802

### Network Meta-Analysis
- Rücker G, Schwarzer G (2015). Ranking treatments in frequentist network meta-analysis works without resampling methods. *BMC Med Res Methodol* 15:58
- Dias S, et al. (2013). Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials. *Med Decis Making* 33:641-656

### Transportability
- Dahabreh IJ, et al. (2020). Extending inferences from a randomized trial to a target population. *Eur J Epidemiol* 35:719-722

##  Support

- **Issues**: https://github.com/your-org/powerNMA/issues
- **Discussions**: https://github.com/your-org/powerNMA/discussions
- **Email**: research@example.com

---

**Built with** | `netmeta` • `survRM2` • `gemtc` • `shiny` • `ggplot2` • `tidyverse`
