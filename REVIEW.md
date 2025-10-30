# cardioTVNMA R Package - Comprehensive Code Review

**Review Date:** 2025-10-30
**Reviewer:** Claude
**Package Version:** 1.0.0

---

## Executive Summary

This R package provides a framework for time-varying network meta-analysis (NMA) of clinical trial data with a focus on cardiovascular outcomes. The package implements RMST and milestone survival analyses, Bayesian NMA, sensitivity analyses, and an interactive Shiny dashboard.

**Overall Assessment:** The package shows good structure and addresses important statistical methodologies, but several critical issues need attention before production use.

**Rating:** ⚠️ **Major Revision Required**

---

## 1. PACKAGE STRUCTURE ASSESSMENT

### ✅ Strengths

- **Well-organized file structure** following R package conventions
- **Comprehensive DESCRIPTION file** with appropriate dependencies
- **Proper use of roxygen2** for documentation
- **Test suite included** using testthat
- **Interactive Shiny dashboard** for results exploration
- **MIT License** clearly specified

### ⚠️ Issues

1. **Missing `.Rbuildignore`** - The structure mentions it but no content provided
2. **Missing example data file** - `inst/extdata/example_data.csv` referenced but not provided
3. **No vignettes** - Would greatly improve usability
4. **Incomplete `.gitignore`** - Not mentioned in the structure

---

## 2. STATISTICAL METHODS REVIEW

### 2.1 RMST Network Meta-Analysis (`rmst_nma.R`)

#### ✅ Strengths
- Uses established `survRM2` package for RMST calculation
- Properly handles pairwise comparisons
- Appropriate use of `netmeta` for network synthesis

#### 🔴 Critical Issues

**Issue #1: Sign Convention Inconsistency**
```r
# Line with the sign flip
TE = -1 * rmst_result$unadjusted.result[1, "Est."]
```
**Problem:** The comment says "we need to flip the sign because netmeta expects treat2 vs treat1", but this logic is potentially flawed. The `rmst2()` function returns arm=1 vs arm=0, and the pairwise data structure should naturally align without manual sign flipping if the comparison order is correct.

**Recommendation:**
- Verify the comparison direction more carefully
- Add explicit validation tests that check if results match expected direction
- Document the rationale more clearly with references

**Issue #2: Limited Error Handling**
```r
tryCatch({
  rmst_result <- survRM2::rmst2(...)
  # ...
}, error = function(e) {
  warning(paste("RMST calculation failed..."))
})
```
**Problem:** Warnings are issued but execution continues. If many trials fail, the final NMA could be based on incomplete data without the user realizing.

**Recommendation:**
- Track and report the number of failed calculations
- Consider throwing an error if too many trials fail (e.g., >20%)
- Return diagnostic information about which trials failed

**Issue #3: No Minimum Trial Validation**
```r
if(nrow(pairwise) > 0) {
  nma_result <- netmeta::netmeta(...)
}
```
**Problem:** Network meta-analysis requires at least 2-3 trials to be meaningful, but this only checks for >0.

**Recommendation:** Add validation:
```r
if(nrow(pairwise) < 2) {
  warning(paste("Insufficient data for NMA at tau =", tau,
                ": only", nrow(pairwise), "comparisons available"))
  next
}
```

### 2.2 Milestone Survival Analysis (`milestone_nma.R`)

#### ✅ Strengths
- **Correct statistical approach** using event counts and odds ratios
- Much improved over the described "flawed log-odds calculation"
- Uses proper `survival::survfit()` for KM estimation

#### ⚠️ Issues

**Issue #4: Potential Edge Cases with Zero Events**
```r
events1 <- n1 - summ1$n.risk
events2 <- n2 - summ2$n.risk
```
**Problem:** If there are zero events in either arm, the odds ratio calculation will fail or produce extreme estimates.

**Recommendation:** Add continuity correction:
```r
# Check for zero events and add continuity correction if needed
if (events1 == 0 || events2 == 0 || events1 == n1 || events2 == n2) {
  warning(paste("Zero events detected in trial", trial_id, "at time", t))
  events1 <- events1 + 0.5
  events2 <- events2 + 0.5
  n1 <- n1 + 1
  n2 <- n2 + 1
}
```

**Issue #5: Time Extension Behavior**
```r
summ1 <- summary(km1, times = t, extend = TRUE)
```
**Problem:** The `extend = TRUE` parameter will extrapolate beyond the observed follow-up time, which can lead to misleading results if the milestone time exceeds the study duration.

**Recommendation:**
- Check maximum follow-up time per trial
- Warn users if milestone time exceeds observed data
- Consider making `extend` a parameter with default `FALSE`

### 2.3 IPD Reconstruction (`data_functions.R`)

#### 🔴 Critical Issues

**Issue #6: Oversimplified Algorithm**
```r
# Calculate number of events and censored
n_events <- round(n_start * (s_start - s_end) / s_start)
n_censored <- n_start - n_end - n_events
```

**Problems:**
1. The Guyot et al. (2012) algorithm is considerably more sophisticated
2. This simplified version doesn't properly account for:
   - The relationship between survival curve coordinates and at-risk numbers
   - Iterative refinement needed for accuracy
   - Constraints from reported summary statistics
3. Uniform distribution of event times within intervals is a crude approximation

**Recommendation:**
- Clearly state in documentation that this is a **simplified approximation**
- Add prominent warnings about accuracy limitations
- Consider implementing or wrapping a proper IPD reconstruction package (e.g., `IPDfromKM`)
- Add validation against known test cases

**Issue #7: No Validation of Reconstructed Data**
```r
# The function returns IPD without checking if it's reasonable
return(ipd)
```

**Recommendation:** Add validation:
```r
# Before returning, validate
total_reconstructed_events <- sum(ipd$status == 1)
if (abs(total_reconstructed_events - total_events) > 5) {
  warning(paste("Reconstructed events (", total_reconstructed_events,
                ") differ substantially from reported (", total_events, ")"))
}
```

### 2.4 Bayesian Analysis (`bayes_nma.R`)

#### ✅ Strengths
- Proper integration with `multinma` package
- Reasonable default priors
- Good MCMC settings (adapt_delta, iterations, chains)

#### ⚠️ Issues

**Issue #8: Limited Prior Options**
```r
prior_trt = multinma::normal(scale = prior_sd)
```
**Problem:** Only normal priors are supported, and the default SD of 10 might be too vague or too informative depending on the outcome scale.

**Recommendation:**
- Allow users to specify custom priors
- Provide guidance on prior selection based on outcome type (RMST vs OR)
- Consider scale-dependent defaults

**Issue #9: No Convergence Diagnostics**
```r
model <- multinma::nma(...)
return(model)
```
**Problem:** The function doesn't check or report convergence diagnostics (Rhat, effective sample size).

**Recommendation:**
```r
# Check convergence
if (requireNamespace("posterior", quietly = TRUE)) {
  diagnostics <- posterior::summarise_draws(model)
  if (any(diagnostics$rhat > 1.1, na.rm = TRUE)) {
    warning("Some parameters have Rhat > 1.1, indicating convergence issues")
  }
}
```

---

## 3. CODE QUALITY ANALYSIS

### 3.1 Error Handling

**Grade: C+**

Issues:
- Inconsistent error handling across functions
- Many functions fail silently or with warnings when they should error
- Limited input validation

Example of good practice to follow:
```r
validate_ipd <- function(ipd) {
  required_cols <- c("trial", "treatment", "time", "status")

  if(!all(required_cols %in% names(ipd))) {
    missing <- setdiff(required_cols, names(ipd))
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }
  # ... more validation
}
```

Recommendation: Apply this validation pattern to all main functions.

### 3.2 Code Duplication

**Issue #10: Repeated Treatment Assignment Logic**

The logic for determining which treatments to compare appears in multiple places:
```r
# This pattern repeats in multiple functions
if (reference %in% treatments) {
  treat1 <- reference
  treat2 <- setdiff(treatments, reference)
} else {
  treat1 <- treatments[1]
  treat2 <- treatments[2]
}
```

**Recommendation:** Extract to a helper function:
```r
assign_comparison_treatments <- function(treatments, reference) {
  if (length(treatments) < 2) {
    stop("Need at least 2 treatments for comparison")
  }

  if (reference %in% treatments) {
    list(
      treat1 = reference,
      treat2 = setdiff(treatments, reference)[1]
    )
  } else {
    list(
      treat1 = treatments[1],
      treat2 = treatments[2]
    )
  }
}
```

### 3.3 Performance Considerations

**Issue #11: Inefficient Data Frame Building**
```r
pairwise <- data.frame()
for(trial_id in unique(ipd$trial)) {
  # ...
  pairwise <- rbind(pairwise, data.frame(...))
}
```

**Problem:** Growing data frames with `rbind()` in loops is inefficient for large datasets.

**Recommendation:**
```r
# Use list accumulation
pairwise_list <- list()
for(trial_id in unique(ipd$trial)) {
  # ...
  pairwise_list[[trial_id]] <- data.frame(...)
}
pairwise <- do.call(rbind, pairwise_list)
```

### 3.4 Documentation Quality

**Grade: B+**

Strengths:
- All exported functions have roxygen2 documentation
- Examples provided for most functions
- Parameter descriptions are clear

Issues:
- Missing details on statistical methodology
- No discussion of assumptions
- Limited guidance on result interpretation
- References to source papers would strengthen credibility

---

## 4. SHINY DASHBOARD REVIEW

### 4.1 UI/UX Assessment

#### ✅ Strengths
- Clean, organized layout
- Multiple visualization types (forest plots, league tables, trajectories)
- Download functionality included
- Model selection (fixed/random effects) included

#### ⚠️ Issues

**Issue #12: Limited Error Messages to User**
```r
current_nma <- shiny::reactive({
  req(input$analysis_type, input$timepoint)
  key <- paste0(...)
  analysis_results[[input$analysis_type]][[key]]
})
```
**Problem:** If the key doesn't exist, the app might crash without a helpful error message.

**Recommendation:**
```r
current_nma <- shiny::reactive({
  req(input$analysis_type, input$timepoint)
  key <- paste0(...)

  result <- analysis_results[[input$analysis_type]][[key]]

  validate(
    need(!is.null(result), "No results available for this selection")
  )

  result
})
```

**Issue #13: Trajectory Plot Limited to RMST**
```r
output$trajectory_plot <- plotly::renderPlotly({
  req(input$analysis_type == "rmst") # Only for RMST
  # ...
})
```
**Problem:** The tab is visible even when milestone analysis is selected, leading to confusion.

**Recommendation:**
```r
# In UI, conditionally show the tab
shiny::conditionalPanel(
  condition = "input.analysis_type == 'rmst'",
  shiny::tabPanel("Time Trajectory", ...)
)
```

### 4.2 Security Considerations

**Issue #14: No Input Sanitization**
If the app is deployed publicly, consider adding input validation to prevent injection attacks or crashes from malformed inputs.

---

## 5. TESTING ASSESSMENT

**Grade: C**

Current tests are minimal:
- Basic IPD generation test
- Smoke tests for main functions
- No tests for edge cases
- No tests for error conditions
- No tests comparing results to known values

### Recommended Additional Tests

```r
test_that("RMST NMA handles insufficient data appropriately", {
  ipd <- generate_example_ipd(n_trials = 1, n_per_arm = 5)
  expect_warning(rmst_nma(ipd))
})

test_that("Milestone NMA handles zero events correctly", {
  # Create IPD with no events
  ipd <- data.frame(
    trial = rep("Trial1", 100),
    treatment = rep(c("A", "B"), each = 50),
    time = runif(100, 100, 200),
    status = rep(0, 100)  # All censored
  )

  result <- milestone_nma(ipd, times = 180)
  expect_warning(result, "Zero events")
})

test_that("Results match known benchmark", {
  # Use published data with known NMA results
  # ...
})
```

---

## 6. DEPENDENCY MANAGEMENT

**Issue #15: Heavy Dependency Load**

The package imports 12 packages, some with large dependency trees:
- `multinma` requires Stan (very large)
- Multiple network meta-analysis packages (`netmeta`, `multinma`)

**Recommendations:**
1. Consider making `multinma` a "Suggests" rather than "Imports" since Bayesian analysis is optional
2. Use `requireNamespace()` with helpful error messages:
```r
if (!requireNamespace("multinma", quietly = TRUE)) {
  stop("Package 'multinma' is required for Bayesian analysis.\n",
       "Install with: install.packages('multinma')")
}
```

---

## 7. MISSING ELEMENTS

### Critical Missing Files

1. **`.Rbuildignore`** - Should include:
```
^.*\.Rproj$
^\.Rproj\.user$
^\.github$
^README\.Rmd$
^LICENSE\.md$
```

2. **`inst/extdata/example_data.csv`** - Referenced but not provided

3. **`.gitignore`** - Should include:
```
.Rproj.user
.Rhistory
.RData
.Ruserdata
inst/doc
```

### Strongly Recommended Additions

1. **Vignette** - Create `vignettes/introduction.Rmd`:
   - Step-by-step tutorial
   - Real-world example
   - Interpretation guidance
   - Comparison of RMST vs Milestone approaches

2. **NEWS.md** - Track changes across versions

3. **CONTRIBUTING.md** - Guidelines for contributors

4. **Code of Conduct** - Standard `CODE_OF_CONDUCT.md`

5. **pkgdown configuration** - For website generation

---

## 8. RECOMMENDATIONS SUMMARY

### Priority 1: Critical (Must Fix Before Release)

1. ✅ **Fix IPD reconstruction** - Either implement properly or clearly state limitations
2. ✅ **Add zero-event handling** in milestone analysis with continuity correction
3. ✅ **Validate RMST sign convention** - Verify treatment effect direction is correct
4. ✅ **Add missing package files** - `.Rbuildignore`, example data
5. ✅ **Improve error handling** - Fail explicitly rather than silently

### Priority 2: Important (Should Fix)

6. ✅ **Expand test coverage** - Add edge case and validation tests
7. ✅ **Add convergence diagnostics** for Bayesian models
8. ✅ **Refactor duplicated code** - Extract helper functions
9. ✅ **Improve Shiny error handling** - Better user feedback
10. ✅ **Add data validation** at function entry points

### Priority 3: Nice to Have

11. ✅ **Create comprehensive vignette**
12. ✅ **Add NEWS.md** for version tracking
13. ✅ **Optimize performance** - Use list accumulation instead of rbind in loops
14. ✅ **Make multinma optional** - Move to Suggests to reduce dependency load
15. ✅ **Add prior selection guidance** with scale-dependent defaults

---

## 9. SPECIFIC CODE CORRECTIONS

### Correction #1: `milestone_nma.R` Zero Event Handling

```r
# BEFORE (around line 50)
events1 <- n1 - summ1$n.risk
events2 <- n2 - summ2$n.risk

pairwise_counts <- rbind(pairwise_counts, data.frame(
  trial = trial_id,
  treat1 = treat1_name,
  treat2 = treat2_name,
  event1 = events1,
  n1 = n1,
  event2 = events2,
  n2 = n2
))

# AFTER
events1 <- n1 - summ1$n.risk
events2 <- n2 - summ2$n.risk

# Apply continuity correction for zero events
needs_correction <- (events1 == 0 || events2 == 0 ||
                     events1 == n1 || events2 == n2)

if (needs_correction) {
  warning(paste("Zero or all events in trial", trial_id, "at time", t,
                "- applying continuity correction"))
  events1 <- events1 + 0.5
  events2 <- events2 + 0.5
  n1 <- n1 + 1
  n2 <- n2 + 1
}

pairwise_counts <- rbind(pairwise_counts, data.frame(
  trial = trial_id,
  treat1 = treat1_name,
  treat2 = treat2_name,
  event1 = events1,
  n1 = n1,
  event2 = events2,
  n2 = n2
))
```

### Correction #2: `rmst_nma.R` Minimum Trial Check

```r
# BEFORE (around line 75)
if(nrow(pairwise) > 0) {
  nma_result <- netmeta::netmeta(...)
}

# AFTER
if(nrow(pairwise) < 2) {
  warning(paste("Insufficient data for NMA at tau =", tau,
                ": only", nrow(pairwise), "comparisons available.",
                "Network meta-analysis requires at least 2 trials."))
  next
}

nma_result <- netmeta::netmeta(...)
results[[paste0("tau_", tau)]] <- nma_result
```

### Correction #3: `bayes_nma.R` Add Convergence Check

```r
# AFTER model fitting
model <- multinma::nma(...)

# Add convergence diagnostics
if (requireNamespace("posterior", quietly = TRUE)) {
  diagnostics <- posterior::summarise_draws(model)

  rhat_issues <- diagnostics$rhat > 1.1
  if (any(rhat_issues, na.rm = TRUE)) {
    problem_params <- diagnostics$variable[rhat_issues]
    warning("Convergence issues detected (Rhat > 1.1) for: ",
            paste(problem_params, collapse = ", "),
            "\nConsider increasing iterations or adjusting priors.")
  }

  # Check effective sample size
  low_ess <- diagnostics$ess_bulk < 100
  if (any(low_ess, na.rm = TRUE)) {
    warning("Low effective sample size (< 100) detected. ",
            "Results may be unreliable.")
  }
}

return(model)
```

---

## 10. DOCUMENTATION IMPROVEMENTS

### Add to Package-Level Documentation

Create `R/cardioTVNMA-package.R`:

```r
#' cardioTVNMA: Time-Varying Network Meta-Analysis for Cardiology
#'
#' @description
#' This package provides methods for conducting time-varying network
#' meta-analyses of clinical trial data, with a focus on cardiovascular outcomes.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{rmst_nma}}: RMST network meta-analysis
#'   \item \code{\link{milestone_nma}}: Milestone survival NMA
#'   \item \code{\link{bayes_nma}}: Bayesian network meta-analysis
#'   \item \code{\link{launch_app}}: Interactive Shiny dashboard
#' }
#'
#' @section Statistical Methods:
#' The package implements two complementary approaches:
#'
#' **RMST Analysis**: Compares restricted mean survival time, representing
#' the average event-free time up to a specified horizon. This provides a
#' clinically interpretable measure of treatment benefit over time.
#'
#' **Milestone Analysis**: Compares the odds of survival at specific time points.
#' This is useful for understanding treatment effects at clinically meaningful
#' moments (e.g., 1-year survival).
#'
#' @references
#' Guyot P, Ades AE, Ouwens MJ, Welton NJ (2012). Enhanced secondary analysis
#' of survival data: reconstructing the data from published Kaplan-Meier
#' survival curves. BMC Med Res Methodol 12:9.
#'
#' Wei Y, Royston P (2017). Reconstructing time-to-event data from published
#' Kaplan-Meier curves. Stata J 17(4):786-802.
#'
#' @docType package
#' @name cardioTVNMA-package
NULL
```

---

## 11. FINAL VERDICT

### Current State
The package demonstrates good understanding of network meta-analysis principles and provides valuable functionality for time-varying treatment effect analysis. The code structure is generally well-organized and follows R package conventions.

### Critical Concerns
1. **IPD reconstruction is oversimplified** and may produce inaccurate results
2. **Zero-event handling is missing**, which will cause failures with sparse data
3. **Limited validation and testing** raises concerns about reliability
4. **Error handling is inconsistent**, potentially leading to silent failures

### Path to Production

**Short-term (v1.0 release):**
1. Fix critical statistical issues (zero events, validation)
2. Add comprehensive error handling
3. Include missing files and proper example data
4. Expand test coverage to at least 80%
5. Add at least one comprehensive vignette

**Medium-term (v1.1-1.2):**
1. Implement proper IPD reconstruction or integrate existing package
2. Add convergence diagnostics for Bayesian models
3. Performance optimization for large networks
4. Additional plotting options
5. Export functions for publication-ready tables

**Long-term (v2.0):**
1. Support for multi-arm trials without splitting
2. Advanced heterogeneity modeling
3. Meta-regression capabilities
4. Integration with other time-to-event NMA methods

---

## 12. POSITIVE ASPECTS WORTH HIGHLIGHTING

Despite the issues identified, several aspects deserve praise:

1. **Comprehensive scope** - Covers multiple important methodologies
2. **User-friendly interface** - Shiny dashboard lowers barrier to entry
3. **Good documentation structure** - Roxygen2 used consistently
4. **Modular design** - Functions are well-separated by purpose
5. **Practical examples** - Example data generation facilitates testing
6. **Modern R practices** - Uses tidyverse, follows current conventions

---

## APPENDIX: Checklist for Authors

- [ ] Fix zero-event handling in milestone_nma
- [ ] Add minimum trial validation in rmst_nma
- [ ] Validate or fix IPD reconstruction algorithm
- [ ] Add convergence diagnostics to bayes_nma
- [ ] Create .Rbuildignore file
- [ ] Add inst/extdata/example_data.csv
- [ ] Expand test coverage (target: 80%+)
- [ ] Write comprehensive vignette
- [ ] Add NEWS.md file
- [ ] Refactor duplicated code
- [ ] Optimize data frame operations
- [ ] Add input validation to all main functions
- [ ] Improve Shiny error handling
- [ ] Review and verify RMST sign convention
- [ ] Add package-level documentation
- [ ] Consider making multinma optional (Suggests)
- [ ] Add statistical references to documentation
- [ ] Create pkgdown website
- [ ] Run R CMD check with no errors, warnings, or notes
- [ ] Submit to CRAN or publish on GitHub with proper README

---

**Review completed: 2025-10-30**
