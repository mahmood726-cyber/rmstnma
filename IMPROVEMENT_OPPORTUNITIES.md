# powerNMA Package - Improvement Opportunities Analysis

**Date**: November 1, 2025
**Analyst**: Claude Code Review
**Status**: Comprehensive Package Review

---

## Executive Summary

This document identifies potential improvements to the powerNMA package beyond the critical fixes already completed. These improvements are categorized by priority and impact.

---

## Category 1: Incomplete Implementations (High Priority)

### 1.1 Network Meta-Regression Models

**File**: `powerNMA/R/network_metareg.R`

**Issues Found**:
- Lines 184-190: Exchangeable coefficient model falls back to shared coefficient
- Lines 203-208: Unrelated coefficient model uses placeholder implementation
- Lines 221-226: Treatment indicator creation is simplified

**Current Behavior**:
```r
# Uses "shared coefficient" model for all cases
# Messages warn user but functionality is limited
```

**Recommended Fix**:
```r
# Implement proper exchangeable model using:
# - Bayesian framework with JAGS/Stan, OR
# - Empirical Bayes approximation using lme4/nlme
# - Hierarchical random effects for coefficients
```

**Impact**: **MEDIUM** - Advanced users may need true exchangeable/unrelated models
**Effort**: **HIGH** - Requires Bayesian implementation or complex RE framework
**Urgency**: **LOW** - Shared coefficient model works for most use cases

---

### 1.2 CNMA Model Fit Statistics

**File**: `powerNMA/R/cnma.R`

**Issue**: Lines 600-615 - Log-likelihood calculation is placeholder (-100)

**Current Code**:
```r
loglik <- -100  # Would calculate proper log-likelihood
aic <- -2 * loglik + 2 * n_params  # Incorrect AIC
bic <- -2 * loglik + log(n_obs) * n_params  # Incorrect BIC
```

**Recommended Fix**:
```r
# Calculate proper log-likelihood from model residuals
# For Gaussian: loglik = -0.5 * (n*log(2*pi*sigma^2) + sum((y-y_hat)^2/sigma^2))
# Extract from netmeta object or component model
```

**Impact**: **MEDIUM** - AIC/BIC used for model selection
**Effort**: **MEDIUM** - Extract from underlying netmeta model
**Urgency**: **MEDIUM** - Model selection currently unreliable

---

### 1.3 Component Network Visualization

**File**: `powerNMA/R/cnma.R`

**Issue**: Line 753 - Component network plot requires igraph

**Current Code**:
```r
# Requires igraph or similar - placeholder for now
message("Component network plot not yet implemented")
```

**Recommended Fix**:
```r
# Option 1: Use igraph for network visualization
# Option 2: Use base R plotting with layout algorithms
# Option 3: Use ggplot2 + ggraph for modern visualization
```

**Impact**: **LOW** - Visualization helpful but not critical
**Effort**: **MEDIUM** - igraph integration
**Urgency**: **LOW** - Can use manual visualization

---

### 1.4 Consistency Checking in NMR

**File**: `powerNMA/R/network_metareg.R`

**Issue**: Lines 235-240 - Consistency check not implemented

**Recommended Fix**:
```r
# Implement node-splitting for meta-regression
# Check if treatment effects vary with covariates consistently
# Use interaction tests or subgroup analyses
```

**Impact**: **LOW** - Important for complex NMR, but advanced use case
**Effort**: **HIGH** - Complex statistical implementation
**Urgency**: **LOW** - Can check manually with sensitivity analyses

---

## Category 2: Shiny GUI Enhancements (Medium Priority)

### 2.1 Result Visualization Placeholders

**File**: `powerNMA/inst/shiny/app.R`

**Issues**: Lines 1393-1410 - Plot outputs are empty placeholders

**Current Code**:
```r
output$plot_network <- renderPlotly({
  plot_ly() %>%
    add_trace(type = "scatter", mode = "markers") %>%
    layout(title = "Network Plot (placeholder)")
})
```

**Recommended Fix**:
```r
output$plot_network <- renderPlotly({
  req(rv$auto_std_result)  # Require results

  # Extract network structure
  nma <- rv$auto_std_result$primary_analysis$model_object

  # Create interactive network plot
  # Use netmeta::netgraph() converted to plotly
  # Or custom D3-based network visualization
})
```

**Impact**: **HIGH** - Users expect to see results visually
**Effort**: **MEDIUM** - Convert netmeta plots to plotly
**Urgency**: **HIGH** - GUI currently shows empty plots

**Specific Fixes Needed**:
1. **Network plot** - Show treatment network with edges
2. **Forest plot** - Treatment effects with confidence intervals
3. **Ranking plot** - P-scores or SUCRA values as bar chart

---

### 2.2 Download Handlers

**File**: `powerNMA/inst/shiny/app.R`

**Issues**: Lines 1412-1334 - Download handlers write placeholder text

**Recommended Fix**:
```r
output$download_report_html <- downloadHandler(
  filename = function() {
    paste0("powerNMA_report_", Sys.Date(), ".html")
  },
  content = function(file) {
    # Generate R Markdown report
    # Render with rmarkdown::render()
    # Include all results, plots, tables
  }
)
```

**Impact**: **HIGH** - Users need to export results
**Effort**: **MEDIUM** - Create R Markdown templates
**Urgency**: **MEDIUM** - Can copy-paste results manually for now

---

### 2.3 Results Display Enhancement

**File**: `powerNMA/inst/shiny/app.R`

**Issue**: Results panels show minimal information

**Recommended Enhancement**:
```r
output$auto_std_results <- renderUI({
  req(rv$auto_std_result)

  result <- rv$auto_std_result

  tagList(
    h4("Treatment Effects"),
    DTOutput("treatment_effects_table"),  # Interactive table

    h4("Treatment Rankings"),
    DTOutput("treatment_rankings"),

    h4("Heterogeneity Assessment"),
    verbatimTextOutput("heterogeneity_stats"),

    h4("Inconsistency Assessment"),
    conditionalPanel(
      condition = "output.inconsistency_available",
      DTOutput("inconsistency_results")
    )
  )
})
```

**Impact**: **HIGH** - Better user experience
**Effort**: **LOW-MEDIUM** - Extract and format results
**Urgency**: **MEDIUM** - Functional but basic

---

## Category 3: Code Quality Improvements (Low-Medium Priority)

### 3.1 Input Validation Enhancement

**Files**: Multiple

**Issue**: Some functions lack comprehensive input validation

**Examples**:
```r
# Current (minimal validation)
threshold_analysis <- function(nma_object, ...) {
  # Missing checks for:
  # - Is nma_object actually an NMA object?
  # - Are there enough treatments for meaningful thresholds?
  # - Are required components present?
}
```

**Recommended Pattern**:
```r
threshold_analysis <- function(nma_object, ...) {
  # Validate inputs
  if (!inherits(nma_object, c("netmeta", "gemtc", "auto_nma"))) {
    stop("nma_object must be a valid NMA object")
  }

  if (length(unique(nma_object$treat)) < 2) {
    stop("At least 2 treatments required for threshold analysis")
  }

  # Check for required components
  if (is.null(nma_object$TE.random)) {
    stop("NMA object missing treatment effects (TE.random)")
  }

  # ... rest of function
}
```

**Impact**: **MEDIUM** - Better error messages, prevent crashes
**Effort**: **LOW** - Add validation to key functions
**Urgency**: **LOW** - Functions generally work if used correctly

---

### 3.2 Function Length and Complexity

**Issue**: Some functions are very long (>200 lines)

**Examples**:
- `auto_standard_nma()` - 689 lines
- `auto_experimental_nma()` - 600+ lines
- `make_automatic_choices()` - 228 lines

**Recommended Refactoring**:
```r
# Before (long function)
auto_standard_nma <- function(...) {
  # 689 lines of code
}

# After (modular)
auto_standard_nma <- function(...) {
  data_chars <- detect_data_characteristics(data)
  choices <- make_automatic_choices(data_chars)
  primary <- run_primary_analysis(data, choices)
  sensitivity <- run_sensitivity_analyses(data, primary)
  diagnostics <- assess_diagnostics(primary)

  return(combine_results(...))
}
```

**Impact**: **MEDIUM** - Easier maintenance and testing
**Effort**: **MEDIUM-HIGH** - Requires careful refactoring
**Urgency**: **LOW** - Code works, refactoring is optimization

---

### 3.3 Consistent Error Handling

**Issue**: Mix of `stop()`, `warning()`, and `message()`

**Recommendation**: Establish clear policy:
```r
# Use stop() for: Fatal errors that prevent continuation
# Use warning() for: Issues that may affect results but allow continuation
# Use message() for: Informational messages about choices made
```

**Example Enhancement**:
```r
# Add custom error classes
stop_invalid_data <- function(msg) {
  stop(structure(
    list(message = msg, call = sys.call(-1)),
    class = c("invalid_data_error", "error", "condition")
  ))
}

# Then users can catch specific errors
tryCatch(
  auto_standard_nma(bad_data),
  invalid_data_error = function(e) {
    # Handle data errors specifically
  }
)
```

**Impact**: **LOW** - Better error handling for advanced users
**Effort**: **MEDIUM** - Requires systematic review
**Urgency**: **LOW** - Current error handling adequate

---

## Category 4: Documentation Enhancements (Low Priority)

### 4.1 Missing Examples in Some Functions

**Issue**: Some internal functions lack examples

**Recommendation**:
```r
# Add @examples to all exported functions
# Even if simple:

#' @examples
#' \dontrun{
#' data(Senn2013, package = "netmeta")
#' result <- cnma(data = Senn2013, components = my_components)
#' print(result)
#' plot(result)
#' }
```

**Impact**: **LOW** - Better discoverability
**Effort**: **LOW** - Add examples incrementally
**Urgency**: **LOW** - Main functions have examples

---

### 4.2 Vignette Creation

**Issue**: No vignettes yet (noted as future work)

**Recommended Vignettes**:
1. "Getting Started with powerNMA"
2. "Automatic vs Manual Pathways"
3. "Component Network Meta-Analysis Tutorial"
4. "Experimental Methods Guide"
5. "Troubleshooting Common Issues"

**Impact**: **MEDIUM** - Helps new users
**Effort**: **HIGH** - Requires writing and testing
**Urgency**: **LOW** - Can be post-publication

---

## Category 5: Performance Optimizations (Low Priority)

### 5.1 Memoization for Repeated Calculations

**Opportunity**: Cache expensive computations

**Example**:
```r
# Use memoise package for expensive functions
library(memoise)

# Memoize league table calculation
calculate_league_table <- memoise(function(nma_object) {
  # Expensive operation
  ...
})
```

**Impact**: **LOW-MEDIUM** - Faster for large networks
**Effort**: **LOW** - Add memoization selectively
**Urgency**: **LOW** - Current performance acceptable

---

### 5.2 Parallel Processing for Simulations

**Opportunity**: Use parallel processing for validation simulations

**Example**:
```r
library(future.apply)
plan(multisession, workers = 4)

results <- future_lapply(1:1000, function(i) {
  simulate_and_analyze_network(...)
})
```

**Impact**: **LOW** - Faster validation (not user-facing)
**Effort**: **LOW** - Add to simulation scripts
**Urgency**: **LOW** - Simulations run once

---

## Category 6: Test Coverage Expansion (Low Priority)

### 6.1 Missing Test Scenarios

**Recommended Additional Tests**:

1. **Edge Cases**:
   - Single study networks
   - Completely disconnected networks
   - Networks with only two-arm trials
   - Missing data patterns

2. **Integration Tests**:
   - Full pathway tests (data → analysis → results)
   - Shiny GUI tests (requires shinytest2)

3. **Stress Tests**:
   - Very large networks (1000+ studies)
   - Very wide networks (100+ treatments)
   - Extreme heterogeneity

**Impact**: **MEDIUM** - Higher confidence
**Effort**: **MEDIUM** - Requires test infrastructure
**Urgency**: **LOW** - Core functionality tested

---

## Priority Matrix

| Improvement | Impact | Effort | Urgency | Priority Score |
|-------------|--------|--------|---------|----------------|
| **Shiny Plot Rendering** | HIGH | MEDIUM | HIGH | **9/9** |
| **Shiny Download Handlers** | HIGH | MEDIUM | MEDIUM | **7/9** |
| **CNMA Log-likelihood Fix** | MEDIUM | MEDIUM | MEDIUM | **5/9** |
| **Shiny Results Display** | HIGH | LOW-MED | MEDIUM | **7/9** |
| **Input Validation** | MEDIUM | LOW | LOW | **4/9** |
| **NMR Exchangeable Model** | MEDIUM | HIGH | LOW | **4/9** |
| **Function Refactoring** | MEDIUM | HIGH | LOW | **4/9** |
| **Vignette Creation** | MEDIUM | HIGH | LOW | **4/9** |
| **Additional Tests** | MEDIUM | MEDIUM | LOW | **4/9** |
| **Network Visualization** | LOW | MEDIUM | LOW | **3/9** |
| **Performance Optimization** | LOW | LOW | LOW | **2/9** |

---

## Recommended Action Plan

### Phase 1: Critical GUI Improvements (High Priority, Before Publication)
**Timeline**: 1-2 days
1. ✅ Implement actual plot rendering in Shiny GUI
   - Network plot from netmeta
   - Forest plot for treatment effects
   - Ranking plot (P-scores/SUCRA)
2. ✅ Create R Markdown report templates
3. ✅ Connect download handlers to real reports

### Phase 2: Model Completeness (Medium Priority, Can be Post-Publication)
**Timeline**: 3-5 days
1. Fix CNMA log-likelihood calculation
2. Enhance results display in Shiny
3. Add input validation to key functions

### Phase 3: Advanced Features (Low Priority, Future Releases)
**Timeline**: 1-2 weeks
1. Implement exchangeable NMR model
2. Refactor long functions
3. Create vignettes
4. Expand test coverage

### Phase 4: Polish (Ongoing)
**Timeline**: Continuous
1. Add examples to all functions
2. Performance optimizations
3. Documentation improvements

---

## Blocking vs Non-Blocking

### ❌ BLOCKING for Publication:
1. **Shiny plot placeholders** - Users expect functional visualization
2. **Download handlers** - Users need to export results

### ✅ NON-BLOCKING for Publication:
- All other improvements listed above
- These enhance the package but don't prevent publication
- Can be addressed in future versions

---

## Conclusion

The powerNMA package is **publication-ready** with the critical fixes already completed. The improvements identified here fall into two categories:

1. **Pre-publication (1-2 days)**:
   - Shiny GUI plot rendering
   - Download report functionality

2. **Post-publication (future releases)**:
   - Model enhancements
   - Code refactoring
   - Additional testing
   - Documentation expansion

**Recommendation**: Complete Phase 1 (GUI improvements) before final submission, defer all other improvements to post-publication updates.

---

**Analysis Date**: November 1, 2025
**Analyst**: Claude Code Review System
**Package Status**: Publication-Ready (pending Phase 1 completion)
