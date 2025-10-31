# Validation Datasets for powerNMA Benchmarking

**Purpose:** Real-world NMA datasets for validating powerNMA against published results

**Date:** October 31, 2025

---

## Strategy

We will validate powerNMA by:
1. Using datasets from R packages (netmeta, gemtc, meta)
2. Replicating published NMA papers
3. Comparing powerNMA results to published findings

---

## Available Datasets in R Packages

### 1. netmeta Package Datasets

**Installation:**
```r
install.packages("netmeta")
library(netmeta)
```

**Key Datasets:**

#### a) Senn2013 - Antidiabetic Drugs
```r
data(Senn2013)
```
- **Description:** Network of 26 treatments for diabetes
- **Outcome:** Change in HbA1c
- **Studies:** 37 trials
- **Publication:** Senn et al. (2013) Statistical issues in drug development
- **Use for:** Complex network validation, many treatments

#### b) Woods2010 - Thrombolytic Drugs
```r
data(Woods2010)
```
- **Description:** Network comparing thrombolytic drugs
- **Outcome:** Mortality
- **Studies:** 50 trials
- **Publication:** Woods et al. (2010) Network meta-analysis
- **Use for:** Multi-arm trial validation (has 3-arm trials)

#### c) Dong2013 - Stroke Prevention
```r
data(Dong2013)
```
- **Description:** Stroke prevention interventions
- **Outcome:** Stroke events
- **Studies:** 25 trials, 8 treatments
- **Publication:** Dong et al. (2013)
- **Use for:** Binary outcome validation

#### d) Franchini2012 - Depression Treatments
```r
data(Franchini2012)
```
- **Description:** Depression treatments
- **Outcome:** Response rate (binary)
- **Studies:** 24 trials
- **Publication:** Franchini et al. (2012)
- **Use for:** Odds ratio validation

#### e) Linde2016 - St John's Wort
```r
data(Linde2016)
```
- **Description:** St John's Wort vs antidepressants
- **Outcome:** Response rate
- **Studies:** 18 trials
- **Publication:** Linde et al. (2016)
- **Use for:** Simple network validation

### 2. gemtc Package Datasets

**Installation:**
```r
install.packages("gemtc")
library(gemtc)
```

**Key Datasets:**

#### a) smoking - Smoking Cessation
```r
data(smoking)
```
- **Description:** Smoking cessation interventions
- **Outcome:** Abstinence at 6-12 months
- **Studies:** 24 trials
- **Treatments:** 4 interventions
- **Use for:** Bayesian NMA validation

#### b) parkinson - Parkinson's Disease
```r
data(parkinson)
```
- **Description:** Parkinson's disease treatments
- **Outcome:** Mean off-time reduction
- **Studies:** 7 trials
- **Treatments:** 5 treatments
- **Use for:** Continuous outcome validation

#### c) blocker - Beta Blockers
```r
data(blocker)
```
- **Description:** Beta blockers for heart failure
- **Outcome:** Mortality
- **Studies:** 22 trials
- **Use for:** Classic validation (widely cited)

### 3. Published Benchmark Datasets

#### a) Lu & Ades (2004) - Thrombolytics
**Paper:** "Combination of direct and indirect evidence in mixed treatment comparisons"
**Journal:** Statistics in Medicine 23:3105-3124

```r
# Data structure:
studies <- c("GUSTO-I", "ISIS-3", "EMERAS", ...)
treatments <- c("SK", "tPA", "rtPA", "SKtPA", "PTCA")
outcome <- "30-day mortality"
n_studies <- 50
has_multi_arm <- TRUE  # Several 3-arm trials
```

**Expected Results (from paper):**
- SK (reference): OR = 1.00
- tPA vs SK: OR = 0.83 (95% CI: 0.73-0.94)
- PTCA vs SK: OR = 0.74 (95% CI: 0.55-1.00)

**Validation Test:**
```r
# powerNMA should match these results within rounding error
```

#### b) Caldwell et al. (2005) - Diabetes Drugs
**Paper:** "Simultaneous comparison of multiple treatments"
**Journal:** BMJ 331:897

```r
treatments <- c("Placebo", "Metformin", "Sulfonylurea", "Thiazolidinedione")
outcome <- "Mean change in HbA1c"
n_studies <- 20+
```

**Expected Results:**
- Metformin vs Placebo: MD = -0.97% (95% CI: -1.25 to -0.70)
- All treatments reduce HbA1c vs placebo

#### c) Cipriani et al. (2009) - Antidepressants
**Paper:** "Comparative efficacy and acceptability of 12 new-generation antidepressants"
**Journal:** Lancet 373:746-758

```r
treatments <- 12 antidepressants + placebo
outcome <- "Response rate" (binary)
n_studies <- 117 trials
network_size <- "large"
```

**Famous For:** One of the largest published NMAs
**Expected:** Detailed results in supplementary materials

---

## Validation Protocol

### Step 1: Exact Replication (Primary Validation)

For each dataset:

```r
# 1. Load dataset
data(Senn2013, package = "netmeta")

# 2. Run with netmeta directly
nm_published <- netmeta::netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = Senn2013,
  sm = "MD",
  reference.group = "plac"
)

# 3. Run with powerNMA
pnma_result <- powerNMA::run_powernma(
  data = Senn2013,
  data_type = "pairwise",
  mode = "standard"
)

# 4. Compare ALL results
testthat::test_that("powerNMA matches netmeta on Senn2013", {
  # Treatment effects
  testthat::expect_equal(
    pnma_result$network$TE.random,
    nm_published$TE.random,
    tolerance = 1e-10
  )

  # Standard errors
  testthat::expect_equal(
    pnma_result$network$seTE.random,
    nm_published$seTE.random,
    tolerance = 1e-10
  )

  # Heterogeneity
  testthat::expect_equal(
    pnma_result$network$tau,
    nm_published$tau,
    tolerance = 1e-10
  )

  # I²
  testthat::expect_equal(
    pnma_result$network$I2,
    nm_published$I2,
    tolerance = 1e-10
  )
})
```

**Expected:** EXACT match (tolerance < 1e-8)

### Step 2: Published Paper Replication

For Lu & Ades (2004) thrombolytics:

```r
# 1. Manually enter or load data
lu_ades_data <- read_csv("lu_ades_2004_data.csv")

# 2. Run analysis
result <- run_powernma(lu_ades_data, mode = "standard")

# 3. Compare to published Table 2
published_results <- data.frame(
  treatment = c("SK", "tPA", "PTCA"),
  OR = c(1.00, 0.83, 0.74),
  CI_lower = c(NA, 0.73, 0.55),
  CI_upper = c(NA, 0.94, 1.00)
)

# 4. Extract powerNMA results
pnma_OR <- exp(result$network$TE.random["tPA", "SK"])
pnma_CI <- exp(result$network$lower.random["tPA", "SK"]),
                exp(result$network$upper.random["tPA", "SK"]))

# 5. Compare (should match within rounding)
testthat::expect_equal(pnma_OR, 0.83, tolerance = 0.02)
```

---

## Validation Test Suite to Create

**File:** `powerNMA/tests/testthat/test-published-datasets.R`

```r
# Test Suite: Published Dataset Validation

test_that("Senn2013: Matches netmeta exactly", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nm <- netmeta::netmeta(...)
  pnma <- run_powernma(Senn2013, mode = "standard")

  expect_equal(pnma$network$TE.random, nm$TE.random, tolerance = 1e-10)
})

test_that("Woods2010: Multi-arm trials handled correctly", {
  skip_if_not_installed("netmeta")

  data(Woods2010, package = "netmeta")

  # This dataset has multi-arm trials
  nm <- netmeta::netmeta(...)
  pnma <- run_powernma(Woods2010, mode = "standard")

  # Should match exactly
  expect_equal(pnma$network$TE.random, nm$TE.random, tolerance = 1e-10)
})

test_that("Dong2013: Binary outcomes (OR) match", {
  skip_if_not_installed("netmeta")

  data(Dong2013, package = "netmeta")

  nm <- netmeta::netmeta(..., sm = "OR")
  pnma <- run_powernma(Dong2013, mode = "standard")

  expect_equal(pnma$network$TE.random, nm$TE.random, tolerance = 1e-10)
})

test_that("smoking (gemtc): Bayesian mode validates", {
  skip_if_not_installed("gemtc")
  skip_if_not_installed("coda")

  data(smoking, package = "gemtc")

  # Run gemtc directly
  network <- gemtc::mtc.network(smoking)
  model <- gemtc::mtc.model(network)
  mcmc <- gemtc::mtc.run(model, n.iter = 10000)

  # Run powerNMA Bayesian
  pnma <- run_powernma(smoking, mode = "standard", use_bayesian = TRUE)

  # Compare posterior medians (should be close)
  # Note: Exact match not expected due to MCMC randomness
})

test_that("Lu & Ades (2004): Replicates published results", {
  # Load manually entered data from the paper
  lu_ades <- read_published_data("lu_ades_2004")

  result <- run_powernma(lu_ades, mode = "standard")

  # Compare to Table 2 from paper
  # tPA vs SK: OR = 0.83 (95% CI: 0.73-0.94)
  tpa_vs_sk_OR <- exp(result$network$TE.random["tPA", "SK"])
  expect_equal(tpa_vs_sk_OR, 0.83, tolerance = 0.02)  # Allow rounding
})
```

---

## How to Access These Datasets

### Option 1: Directly from R Packages
```r
# Install packages
install.packages(c("netmeta", "gemtc", "meta"))

# List all datasets
data(package = "netmeta")$results[, "Item"]
data(package = "gemtc")$results[, "Item"]

# Load specific dataset
data(Senn2013, package = "netmeta")
head(Senn2013)
```

### Option 2: Download from Package Documentation
```bash
# netmeta datasets: https://github.com/guido-s/netmeta
# gemtc datasets: https://github.com/gertvv/gemtc
```

### Option 3: From Published Papers
- Lu & Ades (2004): Supplementary data available
- Caldwell et al. (2005): Data in paper appendix
- Cipriani et al. (2009): Supplementary materials

---

## Recommended Validation Workflow

### Phase 1: R Package Datasets (Quick Validation)
**Time:** 1-2 days

1. **Senn2013** - Large network, complex
2. **Woods2010** - Multi-arm trials
3. **Dong2013** - Binary outcomes
4. **smoking** - Bayesian validation

**Goal:** Confirm exact match with netmeta/gemtc

### Phase 2: Published Papers (Credibility)
**Time:** 1 week

1. **Lu & Ades (2004)** - Classic thrombolytics
2. **Caldwell et al. (2005)** - Diabetes drugs
3. **One time-to-event dataset** - For RMST validation

**Goal:** Show powerNMA replicates published findings

### Phase 3: Create Validation Vignette
**Time:** 3-5 days

```r
vignette("validation-benchmarks", package = "powerNMA")
```

Contents:
- All datasets tested
- Comparison tables (powerNMA vs published)
- Interpretation
- Conclusion: "powerNMA validated"

---

## Expected Outcomes

After completing validation:

1. **Standard Mode:**
   - ✓ Matches netmeta to 1e-10 on 5+ real datasets
   - ✓ Matches gemtc on 2+ real datasets
   - ✓ Replicates 3+ published papers

2. **Publications:**
   - Can cite: "Validated against netmeta using 5 published datasets"
   - Can cite: "Replicated results from Lu & Ades (2004), Caldwell et al. (2005)"

3. **User Confidence:**
   - Reviewers will trust results
   - Suitable for systematic reviews
   - Cochrane-acceptable

---

## Data Sources Summary

| Dataset | Source | N Studies | N Treatments | Outcome | Use For |
|---------|--------|-----------|--------------|---------|---------|
| **Senn2013** | netmeta | 37 | 26 | HbA1c (cont.) | Complex networks |
| **Woods2010** | netmeta | 50 | 8 | Mortality (bin.) | Multi-arm validation |
| **Dong2013** | netmeta | 25 | 8 | Stroke (bin.) | Binary outcomes |
| **Franchini2012** | netmeta | 24 | 9 | Response (bin.) | Odds ratios |
| **smoking** | gemtc | 24 | 4 | Abstinence (bin.) | Bayesian NMA |
| **parkinson** | gemtc | 7 | 5 | Off-time (cont.) | Continuous outcomes |
| **Lu & Ades** | Published | 50 | 5 | Mortality (bin.) | Classic benchmark |

---

## Implementation Plan

### Create Validation Test File

**File:** `powerNMA/tests/testthat/test-real-datasets.R`

```r
context("Validation with Real Published Datasets")

# Helper function to load datasets safely
load_if_available <- function(dataset, package) {
  if (requireNamespace(package, quietly = TRUE)) {
    data(list = dataset, package = package, envir = environment())
    get(dataset)
  } else {
    testthat::skip(paste(package, "not available"))
  }
}

test_that("Senn2013 (netmeta): Exact match", {
  ...
})

test_that("Woods2010 (netmeta): Multi-arm correct", {
  ...
})

# ... add all datasets
```

### Create Validation Vignette

**File:** `powerNMA/vignettes/validation-benchmarks.Rmd`

```rmd
---
title: "Validation Benchmarks: powerNMA vs Published Results"
author: "powerNMA Development Team"
output: rmarkdown::html_vignette
---

## Introduction

This vignette demonstrates that powerNMA produces identical results to
validated packages (netmeta, gemtc) and replicates published findings.

## Dataset 1: Senn2013 (Diabetes Drugs)

### Published Results
- Source: Senn et al. (2013)
- 37 trials, 26 treatments
- Outcome: HbA1c change

### Validation
[Code and results comparing powerNMA to netmeta]

### Conclusion
✓ Exact match (difference < 1e-10)

## Dataset 2: Lu & Ades (2004) (Thrombolytics)
...

## Summary Table
[Table showing all validations passed]
```

---

## Next Steps

1. **Install Required Packages**
   ```r
   install.packages(c("netmeta", "gemtc", "coda"))
   ```

2. **Create Test File**
   - `powerNMA/tests/testthat/test-real-datasets.R`
   - Add tests for each dataset

3. **Run Validation**
   ```r
   devtools::test()
   ```

4. **Document Results**
   - Create validation vignette
   - Update README with "Validated against X datasets"

5. **Commit**
   ```bash
   git add tests/testthat/test-real-datasets.R
   git commit -m "Add validation tests with real published datasets"
   ```

---

**Status:** Ready to implement
**Time Required:** 2-3 days for complete validation suite
**Impact:** Greatly increases credibility and trust in powerNMA

---

**Document Version:** 1.0
**Last Updated:** October 31, 2025
