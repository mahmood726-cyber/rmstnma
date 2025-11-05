# powerNMA Validation and Production Readiness Plan

**Version:** 2.0
**Date:** October 31, 2025
**Status:** Implementation Roadmap

---

## Executive Summary

This document outlines a comprehensive validation strategy to transform powerNMA from a research prototype into a production-ready tool suitable for systematic reviews and clinical decision-making. Based on critical methodological reviews, we split the package into **two operational modes**:

1. **STANDARD MODE**: Validated, publication-ready methods (netmeta/gemtc wrappers)
2. **EXPERIMENTAL MODE**: Novel methods requiring validation (time-varying NMA, transportability)

**Goal:** Achieve Cochrane-level acceptance within 6-12 months.

---

## Critical Issues Addressed

### ✅ Issue #1: IPD Reconstruction (RESOLVED)
**Status:** REMOVED from package
**Action Taken:** Deleted `reconstruct_ipd()` function (methodologically unsound)
**Recommendation:** Users should use validated packages (IPDfromKM, digitize)

### ⚠️ Issue #2: Multi-Arm Trials (CRITICAL)
**Status:** REQUIRES FIX
**Impact:** Currently discards data from >2-arm studies
**Priority:** Must fix before any publication use

### ⚠️ Issue #3: RMST Sign Convention (CRITICAL)
**Status:** REQUIRES MATHEMATICAL VALIDATION
**Impact:** Could reverse all treatment effect directions
**Priority:** Must validate before experimental use

### ⚠️ Issue #4: Milestone Extrapolation (MAJOR)
**Status:** REQUIRES DEFAULT CHANGE
**Impact:** `extend = TRUE` creates pseudo-data
**Priority:** Change to `extend = FALSE`

### ⚠️ Issue #5: Continuity Correction (MAJOR)
**Status:** REQUIRES STANDARDIZATION
**Impact:** Non-standard implementation affects comparability
**Priority:** Use standard 0.5 correction or make configurable

### ⚠️ Issue #6: Transportability Incomplete (MAJOR)
**Status:** REQUIRES COMPLETION
**Impact:** Missing diagnostics for proper use
**Priority:** Add balance checks and effective N

---

## TWO-MODE ARCHITECTURE

### Mode 1: STANDARD (Production-Ready)

**Purpose:** Clinical systematic reviews, guideline development, journal publication

**Included Methods:**
- ✅ Standard frequentist NMA (netmeta wrapper)
- ✅ Bayesian NMA (gemtc wrapper)
- ✅ Leave-one-out sensitivity
- ✅ Inconsistency assessment (design-by-treatment)
- ✅ Network geometry metrics
- ✅ PRISMA-NMA checklist generation
- ✅ Risk of bias integration
- ✅ Protocol specification
- ✅ Analysis archiving

**Validation Status:** Built on validated packages (netmeta, gemtc)

**Suitable For:**
- Cochrane systematic reviews
- Clinical practice guidelines (WHO, NICE)
- Journal submissions (BMJ, Lancet, JAMA)
- Meta-analyses informing regulatory decisions

**Quality Assurance:**
- All methods have 10+ years peer validation
- Benchmark tests against netmeta/gemtc directly
- No novel statistical methods

---

### Mode 2: EXPERIMENTAL (Research Use)

**Purpose:** Methods research, exploratory analysis, hypothesis generation

**Included Methods:**
- ⚠️ Time-varying NMA (RMST)
- ⚠️ Time-varying NMA (Milestone survival)
- ⚠️ Transportability weighting
- ⚠️ PET-PEESE for NMA
- ⚠️ Component NMA
- ⚠️ Leave-one-treatment-out (LOTO)
- ⚠️ Multiverse robustness

**Validation Status:** Novel methods requiring validation

**Suitable For:**
- Academic methods research papers
- Exploratory comparative effectiveness research
- Sensitivity analyses (complementing standard NMA)
- Proof-of-concept studies

**Quality Assurance (To Be Completed):**
- [ ] Mathematical validation of sign conventions
- [ ] Simulation studies demonstrating accuracy
- [ ] Application to published datasets with replication
- [ ] Independent peer review
- [ ] Published methods papers

**User Warnings:**
```r
# Experimental mode requires explicit opt-in
results <- run_powernma(data, mode = "experimental")
#> Warning: EXPERIMENTAL MODE enabled
#> These methods are NOT validated for clinical decision-making
#> Suitable for research and exploratory analysis only
#> Do NOT use for Cochrane reviews or clinical guidelines
```

---

## VALIDATION WORKFLOW

### Phase 1: Fix Critical Bugs (Weeks 1-4)

#### 1.1 Multi-Arm Trial Handling
**Current Bug:**
```r
# Only analyzes FIRST comparison in 3-arm trials
treat2 <- setdiff(treatments, reference)[1]  # BUG: [1] discards other arms
```

**Fix Required:**
```r
# Option A: Properly handle all pairwise comparisons
milestone_nma_multiarm <- function(ipd, times, reference) {
  # For each trial with k>2 arms, create all k(k-1)/2 comparisons
  # Account for within-trial correlation using:
  # - Variance adjustment (White 2012 method)
  # - Multivariate meta-analysis (mvmeta)
  # - Use netmeta's built-in multi-arm handling
}

# Option B: Document limitation clearly
milestone_nma <- function(ipd, times, reference, multiarm = "error") {
  if (multiarm == "error" && any_multiarm_trials(ipd)) {
    stop("Multi-arm trials detected. This function handles 2-arm trials only.
          Use netmeta::pairwise() to create all comparisons.")
  }
}
```

**Validation:**
- [ ] Test on known 3-arm trial dataset
- [ ] Compare to netmeta results
- [ ] Verify no data loss

**Timeline:** 1 week

---

#### 1.2 RMST Sign Convention
**Current Issue:**
```r
TE = -1 * rmst_result$unadjusted.result[1, "Est."]  # WHY flip sign?
```

**Validation Required:**

**Test Case:**
```r
# Manual validation with known data
test_rmst_sign_convention <- function() {
  # Create simple 2-arm trial: Control vs Treatment
  # Treatment BETTER (longer survival)
  ipd <- data.frame(
    trial = "Test",
    treatment = rep(c("Control", "Treatment"), each = 100),
    time = c(rexp(100, 0.05), rexp(100, 0.03)),  # Treatment has LOWER hazard
    status = 1
  )

  # Calculate RMST manually at tau = 365
  km_control <- survfit(Surv(time, status) ~ 1,
                        data = ipd[ipd$treatment == "Control",])
  km_treatment <- survfit(Surv(time, status) ~ 1,
                          data = ipd[ipd$treatment == "Treatment",])

  rmst_control <- sum(summary(km_control, times = 0:365)$surv)
  rmst_treatment <- sum(summary(km_treatment, times = 0:365)$surv)

  # Treatment effect (Treatment - Control)
  manual_TE <- rmst_treatment - rmst_control  # Should be POSITIVE

  # Package result
  package_TE <- rmst_nma(ipd, tau = 365, reference = "Control")
  package_TE <- package_TE$data$TE[1]

  # Verify they match
  expect_equal(manual_TE, package_TE, tolerance = 1e-6)

  # Verify direction: Treatment better → positive TE
  expect_true(package_TE > 0,
    info = "Treatment with better survival should have positive RMST difference")
}
```

**Mathematical Proof Required:**

Document the relationship:
```
survRM2::rmst2() returns:
  arm=1 - arm=0

powerNMA pairwise structure:
  treat1 vs treat2

If treat1 = "Treatment", treat2 = "Control"
And we want: "Treatment - Control"
Then we need: arm=1 - arm=0 where arm=1 is Treatment

But survRM2 uses:
  arm = 0 for first group in data
  arm = 1 for second group in data

Therefore: [COMPLETE MATHEMATICAL DERIVATION]
```

**Timeline:** 1 week

---

#### 1.3 Milestone Extend Parameter
**Current Issue:**
```r
summ1 <- summary(km1, times = t, extend = TRUE)  # Extrapolates!
```

**Fix:**
```r
milestone_nma <- function(ipd, times, reference = NULL,
                         extend = FALSE,  # DEFAULT CHANGED
                         check_followup = TRUE) {

  if (check_followup) {
    # Check if milestone times exceed follow-up
    max_followup <- ipd %>%
      group_by(trial, treatment) %>%
      summarize(max_time = max(time), .groups = "drop")

    for (t in times) {
      insufficient <- max_followup %>% filter(max_time < t)

      if (nrow(insufficient) > 0) {
        warning(sprintf(
          "Milestone time %d exceeds follow-up in %d trial-arms:\n%s\n
          Consider using shorter milestone or excluding these trials.",
          t, nrow(insufficient),
          paste(capture.output(print(insufficient)), collapse = "\n")
        ))
      }
    }
  }

  # ... rest of function with extend = extend (now defaults FALSE)
}
```

**Validation:**
- [ ] Test with varying follow-up durations
- [ ] Compare results with extend=TRUE vs extend=FALSE
- [ ] Document differences in vignette

**Timeline:** 3 days

---

#### 1.4 Continuity Correction
**Current Non-Standard Approach:**
```r
events1 <- events1 + 0.5
events2 <- events2 + 0.5
n1 <- n1 + 1  # UNUSUAL: also increases denominator
n2 <- n2 + 1
```

**Standard Cochrane Approach:**
```r
# Add 0.5 to events only, not to n
events1 <- events1 + 0.5
events2 <- events2 + 0.5
# n1 and n2 unchanged
```

**Make Configurable:**
```r
milestone_nma <- function(ipd, times, reference = NULL,
                         continuity_correction = c("standard", "empirical", "none")) {

  cc <- match.arg(continuity_correction)

  if (needs_correction) {
    switch(cc,
      "standard" = {
        # Cochrane method: add 0.5 to numerator only
        events1 <- events1 + 0.5
        events2 <- events2 + 0.5
      },
      "empirical" = {
        # Current method: add 0.5 to numerator, 1 to denominator
        events1 <- events1 + 0.5
        events2 <- events2 + 0.5
        n1 <- n1 + 1
        n2 <- n2 + 1
      },
      "none" = {
        # Exclude zero-event studies
        warning("Excluding zero-event studies")
        next
      }
    )
  }
}
```

**Validation:**
- [ ] Compare all three methods on sparse data
- [ ] Benchmark against netmeta with same correction
- [ ] Document in methods paper

**Timeline:** 3 days

---

### Phase 2: Complete Experimental Methods (Weeks 5-8)

#### 2.1 Transportability Diagnostics

**Add Required Components:**

```r
transportability_diagnostics <- function(data, target_population, weights) {

  # 1. Covariate Balance (SMD before/after weighting)
  covariates <- c("age_mean", "female_pct", "bmi_mean")

  balance <- map_dfr(covariates, function(cov) {
    smd_before <- calculate_smd(data[[cov]], target_population[[cov]])
    smd_after <- calculate_smd(data[[cov]], target_population[[cov]], weights)

    tibble(
      covariate = cov,
      SMD_unweighted = smd_before,
      SMD_weighted = smd_after,
      improvement = abs(smd_before) - abs(smd_after)
    )
  })

  # 2. Effective Sample Size
  ess <- (sum(weights))^2 / sum(weights^2)
  ess_ratio <- ess / length(weights)

  # 3. Weight Distribution
  weight_summary <- summary(weights)
  extreme_weights <- sum(weights > quantile(weights, 0.95))

  # 4. Positivity Check
  # Are target population characteristics within convex hull of trials?
  positivity <- check_positivity(data[, covariates],
                                 target_population[covariates])

  # 5. Diagnostic Plot
  plot_diagnostics <- ggplot() +
    geom_histogram(aes(x = weights), bins = 30) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Distribution of Transportability Weights",
         subtitle = sprintf("ESS = %.1f (%.1f%% of original N = %d)",
                           ess, ess_ratio * 100, length(weights)))

  # Return comprehensive diagnostics
  list(
    balance_table = balance,
    effective_n = ess,
    ess_ratio = ess_ratio,
    weight_summary = weight_summary,
    extreme_weights = extreme_weights,
    positivity = positivity,
    diagnostic_plot = plot_diagnostics,
    recommendation = generate_transport_recommendation(ess_ratio, balance, positivity)
  )
}

generate_transport_recommendation <- function(ess_ratio, balance, positivity) {
  if (ess_ratio < 0.5) {
    return("CAUTION: ESS < 50% of original sample. Transportability weights are highly variable.")
  }
  if (any(abs(balance$SMD_weighted) > 0.1)) {
    return("WARNING: Some covariates still have SMD > 0.1 after weighting. Consider additional covariates.")
  }
  if (!positivity$all_within_hull) {
    return("ERROR: Target population characteristics outside trial data. Transportability invalid.")
  }
  return("OK: Transportability weights appear reasonable. Proceed with caution.")
}
```

**Validation:**
- [ ] Test on simulation with known target population
- [ ] Verify SMD calculations match standard definitions
- [ ] Check ESS formula correctness

**Timeline:** 1 week

---

#### 2.2 PET-PEESE for Network Meta-Analysis

**Current Implementation Review:**

The current PET-PEESE is comparison-wise. For NMA, we need:

```r
pet_peese_nma <- function(nma_result, min_comparisons = 10) {
  # Extract all pairwise comparisons from network
  # NOT just direct comparisons

  # For each comparison with k ≥ min_comparisons studies:
  # 1. Run PET (regress TE on SE)
  # 2. If PET p < 0.10, run PEESE (regress TE on SE^2)
  # 3. Report adjusted estimates

  # Generate network-wide publication bias assessment
}
```

**Validation:**
- [ ] Simulate publication bias in network
- [ ] Verify PET-PEESE recovers true effects
- [ ] Compare to standard comparison-wise PET-PEESE

**Timeline:** 1 week

---

### Phase 3: Benchmark Validation (Weeks 9-12)

#### 3.1 Validate Against Published Datasets

**Dataset Selection:**

1. **Thrombolytics NMA** (Lu & Ades 2004)
   - 50 trials, 8 treatments
   - Published results available
   - Multi-arm trials included
   - Test standard NMA mode

2. **Statins for Primary Prevention** (Tonelli et al. 2011)
   - Time-to-event outcomes
   - Opportunity to compare HR-based vs RMST-based NMA
   - Test experimental RMST mode

3. **Antihypertensives NMA** (various)
   - Large network with heterogeneity
   - Test transportability methods
   - Test publication bias methods

**Validation Protocol:**

```r
# For each published dataset:
validate_against_published <- function(dataset_name) {
  # 1. Replicate original analysis
  original_results <- get_published_results(dataset_name)
  powerNMA_results <- run_powernma(dataset, mode = "standard")

  # 2. Compare point estimates
  comparison <- compare_estimates(original_results, powerNMA_results)
  expect_equal(comparison$mean_absolute_diff, 0, tolerance = 0.01)

  # 3. Compare standard errors
  expect_equal(comparison$se_correlation, 1, tolerance = 0.05)

  # 4. Compare rankings
  expect_equal(comparison$ranking_concordance, 1, tolerance = 0.1)

  # 5. Generate validation report
  generate_validation_report(dataset_name, comparison)
}
```

**Timeline:** 3 weeks

---

#### 3.2 Simulation Studies

**Study 1: Standard NMA Accuracy**

```r
# Verify powerNMA STANDARD mode matches netmeta exactly
sim_standard_nma <- function(n_sim = 1000) {
  for (i in 1:n_sim) {
    # Generate network with known true effects
    sim_data <- simulate_nma_network(
      n_studies = 30,
      n_treatments = 5,
      true_effects = c(0, -0.2, -0.4, -0.15, -0.35),
      tau = 0.1
    )

    # Run both packages
    netmeta_result <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab,
                                       data = sim_data)
    powerNMA_result <- run_powernma(sim_data, mode = "standard")

    # Compare
    expect_equal(powerNMA_result$network$TE.fixed,
                 netmeta_result$TE.fixed)
  }
}
```

**Study 2: RMST NMA Validity**

```r
# Verify RMST NMA has correct Type I error and power
sim_rmst_nma <- function(n_sim = 1000) {
  # Scenario 1: No treatment effect (Type I error check)
  # Scenario 2: Known treatment effect (Power check)
  # Scenario 3: Non-proportional hazards (RMST advantage check)
}
```

**Study 3: Milestone NMA Sparse Data**

```r
# Test continuity correction methods
sim_milestone_sparse <- function() {
  # Generate data with zero events in some arms
  # Test all three continuity correction options
  # Verify no crashes, sensible results
}
```

**Timeline:** 2 weeks

---

### Phase 4: Documentation & Peer Review (Weeks 13-16)

#### 4.1 Methods Paper

**Title:** "powerNMA: An R Package for Standard and Time-Varying Network Meta-Analysis"

**Outline:**

1. **Introduction**
   - Limitations of HR-based NMA with non-proportional hazards
   - Need for time-varying methods
   - Gap in existing software

2. **Methods**
   - Standard NMA (netmeta/gemtc wrappers)
   - RMST-based NMA (mathematical derivation)
   - Milestone survival NMA
   - Multi-arm trial handling
   - Continuity corrections

3. **Validation**
   - Benchmarks against published datasets
   - Simulation study results
   - Comparison to netmeta/gemtc

4. **Case Study**
   - Real data example showing non-proportional hazards
   - HR-based vs RMST-based conclusions
   - Clinical interpretation

5. **Discussion**
   - When to use RMST vs HR
   - Limitations and assumptions
   - Future directions

**Target Journal:** Research Synthesis Methods

**Timeline:** 4 weeks to write, 3-6 months peer review

---

#### 4.2 Comprehensive Vignettes

**Vignette 1: Standard NMA for Systematic Reviews**
- Complete workflow from data prep to PRISMA-NMA checklist
- Integration with RevMan
- GRADE assessment
- Publication-ready tables and figures

**Vignette 2: Time-Varying NMA (Experimental)**
- When to use RMST vs milestone vs HR
- IPD requirements
- Sensitivity analyses
- Interpretation guidelines

**Vignette 3: Advanced Methods**
- Transportability weighting
- Publication bias in NMA
- Component NMA
- Meta-regression

**Timeline:** 2 weeks

---

#### 4.3 Independent Review

**Reviewers to Approach:**

1. **Guido Schwarzer** (netmeta author)
   - Review standard mode implementation
   - Verify correct use of netmeta

2. **Sofia Dias** (gemtc/multinma author)
   - Review Bayesian implementation
   - Methodological soundness

3. **Cochrane NMA Methods Group**
   - Review for systematic review suitability
   - PRISMA-NMA compliance

**Timeline:** 6-12 weeks (external reviewers)

---

## IMPLEMENTATION PLAN: STANDARD vs EXPERIMENTAL MODES

### Code Architecture

```r
# R/modes.R

#' Run powerNMA Analysis
#'
#' @param mode Character: "standard" or "experimental"
#' @export
run_powernma <- function(data,
                        data_type = c("pairwise", "ipd"),
                        mode = c("standard", "experimental"),
                        config = NULL,
                        ...) {

  mode <- match.arg(mode)
  data_type <- match.arg(data_type)

  # MODE GATE
  if (mode == "standard") {
    # Only validated methods
    if (data_type == "ipd") {
      stop("IPD-based methods (RMST, milestone) are EXPERIMENTAL.
            Use mode='experimental' to enable.
            For standard NMA, provide pairwise data.")
    }

    results <- run_standard_nma(data, config, ...)

  } else {
    # Experimental methods with warnings
    warning("
╔════════════════════════════════════════════════════════════╗
║                  EXPERIMENTAL MODE ENABLED                  ║
╠════════════════════════════════════════════════════════════╣
║ These methods are NOT validated for clinical decisions     ║
║ Suitable for: Research, exploratory analysis               ║
║ NOT suitable for: Cochrane reviews, clinical guidelines    ║
║                                                             ║
║ See VALIDATION_PLAN.md for status of each method           ║
╚════════════════════════════════════════════════════════════╝
    ", immediate. = TRUE)

    results <- run_experimental_nma(data, data_type, config, ...)
  }

  # Tag results with mode
  results$mode <- mode
  results$validated <- (mode == "standard")

  class(results) <- c("powernma_result", class(results))
  results
}

#' Standard Mode: Validated Methods Only
run_standard_nma <- function(data, config, ...) {
  list(
    network = robust_netmeta(data, ...),      # netmeta wrapper
    bayesian = run_bayesian_nma_simple(data, ...), # gemtc wrapper
    sensitivity = loo_sensitivity_simple(data, ...),
    inconsistency = check_inconsistency(data, ...),
    geometry = network_geometry(data),

    # v1.1 systematic review features
    prisma_checklist = generate_prisma_nma_checklist(...),
    rob_assessment = integrate_rob_assessment(...),

    methods = "Standard frequentist/Bayesian NMA",
    validation_status = "VALIDATED - suitable for publication",
    based_on = "netmeta 2.9.0, gemtc 1.0-2"
  )
}

#' Experimental Mode: Novel Methods (Use with Caution)
run_experimental_nma <- function(data, data_type, config, ...) {
  # Can include everything
  results <- list()

  if (data_type == "pairwise") {
    results$network <- robust_netmeta(data, ...)

    # Experimental methods for pairwise data
    results$transportability <- compute_transport_weights(...)
    results$publication_bias <- pet_peese_analysis(data)
    results$loto <- loto_sensitivity(data, ...)
    results$multiverse <- run_multiverse(data, ...)

  } else {
    # IPD-based time-varying methods
    results$rmst <- rmst_nma(data, ...)
    results$milestone <- milestone_nma(data, ...)
  }

  results$methods <- "Experimental time-varying NMA"
  results$validation_status <- "EXPERIMENTAL - research use only"
  results$warnings <- list(
    "RMST sign convention requires validation",
    "Milestone extend parameter may extrapolate",
    "Multi-arm trial handling incomplete",
    "Not suitable for clinical decision-making"
  )

  results
}

#' Print Method with Mode Indicator
#' @export
print.powernma_result <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("  powerNMA Results\n")
  cat("═══════════════════════════════════════════════════════\n\n")

  cat("Mode:", toupper(x$mode), "\n")

  if (x$mode == "standard") {
    cat("Status: ✓ VALIDATED - Suitable for publication\n")
    cat("Based on: netmeta, gemtc (peer-reviewed packages)\n")
  } else {
    cat("Status: ⚠ EXPERIMENTAL - Research use only\n")
    cat("\n")
    cat("Warnings:\n")
    for (w in x$warnings) {
      cat("  •", w, "\n")
    }
  }

  cat("\n")
  # ... rest of print method
}
```

---

## TESTING STRATEGY

### Unit Tests (Required Before v2.0)

```r
# tests/testthat/test-standard-mode.R

test_that("Standard mode matches netmeta exactly", {
  data <- simulate_nma_data(n_studies = 20, seed = 123)

  # Direct netmeta call
  nm <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)

  # powerNMA standard mode
  pnma <- run_powernma(data, mode = "standard")

  # Compare all key outputs
  expect_equal(pnma$network$TE.fixed, nm$TE.fixed)
  expect_equal(pnma$network$seTE.fixed, nm$seTE.fixed)
  expect_equal(pnma$network$TE.random, nm$TE.random)
  expect_equal(pnma$network$tau, nm$tau)
})

test_that("Standard mode rejects IPD data", {
  ipd <- generate_example_ipd(n_trials = 5)

  expect_error(
    run_powernma(ipd, data_type = "ipd", mode = "standard"),
    "IPD-based methods.*EXPERIMENTAL"
  )
})

test_that("Experimental mode shows warnings", {
  ipd <- generate_example_ipd(n_trials = 5)

  expect_warning(
    run_powernma(ipd, data_type = "ipd", mode = "experimental"),
    "EXPERIMENTAL MODE"
  )
})

test_that("Multi-arm trials handled correctly", {
  # Create 3-arm trial
  data <- data.frame(
    studlab = "Trial1",
    treat1 = c("A", "A", "B"),
    treat2 = c("B", "C", "C"),
    TE = c(-0.2, -0.4, -0.2),
    seTE = c(0.1, 0.1, 0.1)
  )

  result <- run_powernma(data, mode = "standard")

  # Should include all 3 comparisons
  expect_equal(nrow(result$network$data), 3)
})

test_that("RMST sign convention is correct", {
  # Test with known better treatment
  ipd <- data.frame(
    trial = "Test",
    treatment = rep(c("Control", "Better"), each = 200),
    time = c(rexp(200, 0.05), rexp(200, 0.02)),  # Better has lower hazard
    status = 1
  )

  result <- run_powernma(ipd, data_type = "ipd", mode = "experimental")

  # Better treatment should have POSITIVE RMST difference
  te <- result$rmst$data$TE[result$rmst$data$treat1 == "Better"]
  expect_true(te > 0,
    label = "Better treatment should have positive RMST vs Control")
})

test_that("Milestone does not extrapolate by default", {
  ipd <- generate_example_ipd(n_trials = 3)

  # Set milestone beyond follow-up
  max_time <- max(ipd$time)
  milestone <- max_time + 100

  expect_warning(
    milestone_nma(ipd, times = milestone),
    "exceeds follow-up"
  )
})

test_that("Continuity correction options work", {
  # Create data with zero events
  ipd <- generate_sparse_ipd()  # Helper function

  standard <- milestone_nma(ipd, continuity_correction = "standard")
  empirical <- milestone_nma(ipd, continuity_correction = "empirical")

  # Results should differ
  expect_false(identical(standard$data$TE, empirical$data$TE))
})
```

---

## ACCEPTANCE CRITERIA

### For STANDARD Mode (v2.0 Release)

- [x] IPD reconstruction removed
- [ ] Multi-arm trials properly handled
- [ ] 100% match with netmeta on benchmark datasets (tolerance < 0.01)
- [ ] 100% match with gemtc on benchmark datasets
- [ ] All unit tests passing
- [ ] Vignette 1 complete (standard workflow)
- [ ] PRISMA-NMA checklist validated by Cochrane reviewer
- [ ] Independent review by netmeta author

### For EXPERIMENTAL Mode (v2.0 Release)

- [ ] RMST sign convention mathematically validated
- [ ] Milestone extend = FALSE by default
- [ ] Continuity correction configurable with standard default
- [ ] Transportability diagnostics complete (ESS, balance, positivity)
- [ ] Multi-arm trials: either fixed or documented limitation
- [ ] Simulation studies showing Type I error control
- [ ] Vignette 2 complete (experimental methods with warnings)
- [ ] Warnings displayed prominently

### For Full Validation (v3.0, 6-12 months)

- [ ] Methods paper published in peer-reviewed journal
- [ ] Application to 5+ published datasets with replication
- [ ] Simulation studies published
- [ ] Independent external review
- [ ] Cochrane Methods Group endorsement (aspirational)

---

## TIMELINE SUMMARY

| Phase | Duration | Deliverable | Status |
|-------|----------|-------------|--------|
| **Phase 1: Fix Critical Bugs** | Weeks 1-4 | Multi-arm, RMST sign, extend, continuity | Pending |
| **Phase 2: Complete Experimental** | Weeks 5-8 | Transportability, PET-PEESE | Pending |
| **Phase 3: Benchmark Validation** | Weeks 9-12 | Published datasets, simulations | Pending |
| **Phase 4: Documentation** | Weeks 13-16 | Vignettes, methods paper draft | Pending |
| **v2.0 Release** | Week 16 | Standard + Experimental modes | Target |
| **Phase 5: Peer Review** | Months 4-9 | Methods paper, external review | Future |
| **v3.0 Release** | Month 12 | Fully validated experimental methods | Future |

---

## RESOURCE REQUIREMENTS

**Personnel:**
- 1 statistical methodologist (20 hrs/week)
- 1 R developer (10 hrs/week)
- 1 systematic review expert (consulting, 5 hrs total)

**Software:**
- R packages: netmeta, gemtc, multinma, metafor, IPDfromKM
- JAGS for Bayesian models
- Published datasets (publicly available)

**Validation:**
- Access to published NMA datasets with IPD
- Computational resources for simulations (1000+ iterations)

---

## RISK MITIGATION

| Risk | Impact | Mitigation |
|------|--------|------------|
| Multi-arm fix more complex than expected | Delays v2.0 | Document limitation, exclude multi-arm for now |
| RMST sign convention incorrect | Critical - invalidates all results | Extensive testing before any release |
| External reviewers decline | Delays validation | Approach multiple reviewers, offer co-authorship |
| Methods paper rejected | Delays v3.0 | Prepare for multiple journals, address feedback |
| Simulation studies show problems | May require method changes | Build in time for iteration |

---

## SUCCESS METRICS

**v2.0 (4 months):**
- ✓ Standard mode used in 1+ Cochrane review
- ✓ No critical bugs reported
- ✓ 100+ downloads/month on CRAN
- ✓ Positive feedback from systematic reviewers

**v3.0 (12 months):**
- ✓ Methods paper accepted
- ✓ Experimental methods used in 3+ methods research papers
- ✓ Independent validation study published
- ✓ Recommended by Cochrane NMA Methods Group

---

## CONTACTS & STAKEHOLDERS

**Development Team:**
- Lead Developer: [TBD]
- Statistical Consultant: [TBD]
- Systematic Review Advisor: [TBD]

**External Reviewers (to approach):**
- Guido Schwarzer (netmeta)
- Sofia Dias (gemtc/multinma)
- Georgia Salanti (CINeMA)
- Cochrane NMA Methods Group

**User Community:**
- Cochrane systematic reviewers
- Academic meta-analysis researchers
- Clinical guideline developers

---

## VERSION CONTROL

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | Oct 30, 2025 | Initial methodological review | Claude |
| 2.0 | Oct 31, 2025 | Comprehensive validation plan, two-mode architecture | Claude |

---

**Next Steps:** Begin Phase 1 implementation (fix critical bugs)

**Document Status:** APPROVED - Ready for Implementation
