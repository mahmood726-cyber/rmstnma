# SurroNMA Integration Plan - COMPREHENSIVE

**Date**: November 1, 2025
**Source Repository**: mahmood726-cyber/surroNMA
**Target**: powerNMA v3.0
**Status**: Repository analyzed, integration plan ready

---

## Executive Summary

Successfully cloned and analyzed surroNMA repository containing **1,181 lines** of advanced surrogate endpoint network meta-analysis code. The package has **30+ exported functions** with sophisticated features including:

✅ Surrogate Index (SI) training with multiple ML methods
✅ Bivariate NMA (joint S+T modeling) with Stan/Frequentist
✅ Surrogate Threshold Effect (STE) calculation
✅ Inconsistency testing and stress analyses
✅ Comprehensive plotting and HTA reporting

---

## SurroNMA Package Structure Analysis

### Core Data Structure

**Function**: `surro_network()` (Lines 61-176)

```r
surro_network(
  data,
  study, trt, comp,
  S_eff, S_se,           # Surrogate endpoint
  T_eff, T_se,           # True endpoint
  S_multi, S_multi_se,   # Multiple surrogates
  corr_ST,               # S-T correlation
  class,                 # Treatment classes
  baseline_risk,         # Baseline risk
  rob                    # Risk of bias
)
```

**Key Features**:
- Handles univariate or multivariate surrogates
- Supports treatment class hierarchies
- Checks network connectivity (igraph)
- Flexible correlation specification

---

### Key Functions Inventory

| Function | Lines | Description | Priority |
|----------|-------|-------------|----------|
| `surro_network()` | 61-176 | Build surrogate network object | 🔥 HIGH |
| `surro_index_train()` | 181-238 | Train SI (OLS/Ridge/PCR/SL3) | 🔥 HIGH |
| `augment_network_with_SI()` | 242-274 | Apply SI to network | 🔥 HIGH |
| `compute_STE()` | 279-286 | Surrogate Threshold Effect | 🔥 HIGH |
| `surro_nma_bayes()` | 468-519 | Bayesian bivariate NMA | 🔥 HIGH |
| `surro_nma_freq()` | 521-613 | Frequentist bivariate NMA | 🔥 HIGH |
| `surro_nma()` | 616-649 | Unified wrapper | 🔥 HIGH |
| `surrogacy_diagnostics()` | 672-692 | R² and validation metrics | ⚡ MEDIUM |
| `stress_surrogacy()` | 694-708 | Sensitivity analysis | ⚡ MEDIUM |
| `nodesplit_analysis()` | 716-748 | Inconsistency for surrogates | ⚡ MEDIUM |
| `global_inconsistency_test()` | 750-762 | Global inconsistency | ⚡ MEDIUM |
| `summarize_treatments()` | 651-658 | Treatment summaries | 💡 LOW |
| `compute_ranks()` | 660-670 | Ranking with MID | 💡 LOW |
| `plot_surrogacy()` | 765-773 | Surrogacy scatter plot | ⚡ MEDIUM |
| `plot_rankogram()` | 775-788 | Ranking probabilities | 💡 LOW |
| `plot_networks()` | 790-794 | Network visualization | 💡 LOW |
| `plot_ste()` | 823-829 | STE visualization | ⚡ MEDIUM |
| `plot_stress_curves()` | 831-840 | Stress analysis plots | ⚡ MEDIUM |
| `export_cinema()` | 842-854 | CINeMA export | 💡 LOW |
| `posterior_predict()` | 856-860 | Posterior predictions | 💡 LOW |
| `pp_check()` | 862-872 | Posterior predictive checks | 💡 LOW |
| `explain()` | 874-889 | Auto-generated explanations | ⚡ MEDIUM |
| `export_report()` | 891-1136 | HTML/PDF reports | ⚡ MEDIUM |

---

## Advanced Features Analysis

### 1. Surrogate Index (SI) Training

**Lines**: 181-238

**Methods Supported**:
- **OLS**: Ordinary least squares
- **Ridge**: L2 regularization (glmnet)
- **PCR**: Principal component regression (pls)
- **SL3**: SuperLearner ensemble (sl3 package)

**Features**:
- Cross-validation R²
- Z-score normalization
- Coefficient extraction for SE propagation

**Integration Priority**: 🔥 **HIGH** - Novel feature not in powerNMA

---

### 2. Bivariate NMA

**Bayesian** (Lines 468-519):
- **Stan code generator** (Lines 333-466)
- Second-order random effects
- Student-t option for robustness
- Class-specific surrogacy (α_g, β_g)
- Inconsistency modeling (random effects on both S and T)
- Global vs class-specific surrogacy parameters

**Frequentist** (Lines 521-613):
- Parametric bootstrap (normal or student-t)
- Delta method for SE propagation
- Correlation handling via bivariate normal

**Integration Priority**: 🔥 **HIGH** - Core feature

---

### 3. Surrogate Threshold Effect (STE)

**Lines**: 279-286

**Formula**:
```r
STE = (threshold_T - α) / β

where:
- α = intercept of S→T relationship
- β = slope of S→T relationship
- threshold_T = meaningful true endpoint threshold (e.g., 0)
```

**Output**:
- Mean, median, 95% CrI
- Interpretation: minimum S effect needed for positive T

**Integration Priority**: 🔥 **HIGH** - Critical for clinical decision-making

---

### 4. Ranking Methods

**SUCRA** (Lines 294-298):
```r
SUCRA = (K - mean_rank) / (K - 1)
```

**POTH** (Lines 300-315):
- Probability of Treatment Hierarchy
- Kendall distance from modal ordering
- Accounts for rank stability

**MID-Adjusted** (Lines 317-330):
- Minimal important difference threshold
- Only counts wins if difference > MID
- More clinically relevant than raw effects

**Integration Priority**: ⚡ **MEDIUM** - Enhancement over current P-scores

---

### 5. Surrogacy Diagnostics

**Lines**: 672-692

**Metrics Computed**:
```r
R²_individual  # Individual-level correlation S-T
R²_trial       # Trial-level correlation TE_S-TE_T
α (intercept)  # Bias in S→T relationship
β (slope)      # Strength of S→T relationship
```

**Validation Criteria**:
- **Strong surrogate**: R²_trial > 0.8, β ≈ 1
- **Moderate**: R²_trial 0.5-0.8
- **Weak**: R²_trial < 0.5 (not recommended)

**Integration Priority**: 🔥 **HIGH** - Essential for surrogate validation

---

### 6. Stress Analysis

**Lines**: 694-708

**Concept**: Test robustness of decisions to uncertainty in surrogacy

**Parameters Varied**:
- **R² multipliers**: 0.5, 0.7, 0.9 (pessimistic assumptions)
- **Slope shifts**: -0.1, 0, +0.1 (α/β perturbations)

**Output**: How rankings/decisions change under worst-case surrogacy

**Integration Priority**: ⚡ **MEDIUM** - Valuable sensitivity analysis

---

### 7. Node-Split for Surrogates

**Lines**: 716-748

**Purpose**: Test inconsistency on specific comparisons

**Approach**:
- Fit NMA excluding direct comparison
- Fit model with only direct evidence
- Compare indirect vs direct estimates
- Repeat for Bayesian and Frequentist

**Integration Priority**: ⚡ **MEDIUM** - Surrogate-specific inconsistency

---

### 8. Reporting Functions

**Export CINeMA** (Lines 842-854):
- Format for CINeMA tool
- Treatment effects, SE, study info
- Surrogacy metadata

**Export Report** (Lines 891-1136):
- HTML/PDF with rmarkdown
- Auto-generated narrative
- All plots embedded
- Surrogacy diagnostics included

**Integration Priority**: ⚡ **MEDIUM** - Professional reporting

---

## Integration Strategy

### Phase 1: Core Integration (HIGH PRIORITY)

**Goal**: Add essential surrogate functionality to powerNMA

**Files to Create**:

1. **`powerNMA/R/experimental_surrogate_nma.R`** (NEW)
   ```r
   #' EXPERIMENTAL: Surrogate Endpoint Network Meta-Analysis
   #'
   #' Based on: surroNMA v1.0 by mahmood726-cyber
   #'
   #' @export
   surrogate_nma <- function(data,
                            surrogate_var,
                            true_var,
                            surrogate_method = c("univariate", "SI_trained"),
                            engine = c("bayes", "freq"),
                            ...) {

     # Build surrogate network object
     net <- build_surrogate_network(
       data = data,
       S_var = surrogate_var,
       T_var = true_var
     )

     # Train SI if multiple surrogates
     if (surrogate_method == "SI_trained") {
       si_model <- train_surrogate_index(net, ...)
       net <- augment_with_SI(net, si_model)
     }

     # Fit bivariate NMA
     if (engine == "bayes") {
       fit <- fit_bivariate_nma_bayes(net, ...)
     } else {
       fit <- fit_bivariate_nma_freq(net, ...)
     }

     # Compute diagnostics
     diagnostics <- compute_surrogacy_diagnostics(fit)

     # Compute STE
     ste <- compute_surrogate_threshold_effect(fit, ...)

     structure(list(
       fit = fit,
       diagnostics = diagnostics,
       ste = ste,
       network = net,
       call = match.call()
     ), class = "surrogate_nma")
   }
   ```

2. **Helper Functions** (same file):
   - `build_surrogate_network()` - Adapt from surroNMA lines 61-176
   - `train_surrogate_index()` - Adapt from lines 181-238
   - `fit_bivariate_nma_bayes()` - Adapt from lines 468-519 + Stan code
   - `fit_bivariate_nma_freq()` - Adapt from lines 521-613
   - `compute_surrogacy_diagnostics()` - Adapt from lines 672-692
   - `compute_surrogate_threshold_effect()` - Adapt from lines 279-286

3. **S3 Methods**:
   ```r
   print.surrogate_nma()
   summary.surrogate_nma()
   plot.surrogate_nma()
   ```

**Estimated Effort**: 2-3 days
**Lines of Code**: ~800-1000 lines

---

### Phase 2: Enhanced Ranking & Plots (MEDIUM PRIORITY)

**Files to Create**:

1. **`powerNMA/R/surrogate_plots.R`** (NEW)
   - `plot_surrogacy_scatter()` - Lines 765-773
   - `plot_ste_threshold()` - Lines 823-829
   - `plot_stress_analysis()` - Lines 831-840

2. **Enhanced Ranking in Existing Files**:
   - Add `compute_poth()` to `prediction_methods.R`
   - Add `mid_adjusted_ranking()` to ranking functions

**Estimated Effort**: 1 day
**Lines of Code**: ~300-400 lines

---

### Phase 3: Advanced Features (MEDIUM PRIORITY)

**Files to Modify/Create**:

1. **Stress Analysis**:
   - `powerNMA/R/surrogate_sensitivity.R` (NEW)
   - Adapt lines 694-708
   - Add to auto_experimental_pathway

2. **Node-Split for Surrogates**:
   - Enhance existing inconsistency functions
   - Add surrogate-specific splits

3. **Reporting**:
   - Integrate `export_report()` functionality (lines 891-1136)
   - Add to Shiny GUI download handlers

**Estimated Effort**: 2 days
**Lines of Code**: ~500 lines

---

### Phase 4: Documentation & Testing (HIGH PRIORITY)

**Files to Create**:

1. **`powerNMA/tests/testthat/test-surrogate-nma.R`** (NEW)
   - Test surrogate network construction
   - Test SI training methods
   - Test bivariate NMA (both engines)
   - Test STE calculation
   - Test diagnostics

2. **`powerNMA/vignettes/surrogate-endpoints.Rmd`** (NEW)
   - When to use surrogates
   - Data requirements
   - Worked example with real data
   - Interpretation guide
   - R² validation thresholds

3. **Function Documentation**:
   - Roxygen2 for all functions
   - Examples for main functions
   - References to literature

**Estimated Effort**: 1-2 days
**Lines of Code**: Tests ~400, Vignette ~300

---

## Code Compatibility Analysis

### Dependencies Required

**SurroNMA Dependencies**:
```r
# Core (required)
cmdstanr     # Bayesian inference
Matrix       # Matrix operations
igraph       # Network connectivity

# Optional (for SI training)
glmnet       # Ridge regression
pls          # PCR
sl3          # SuperLearner
data.table   # sl3 dependency
```

**PowerNMA Existing**:
```r
netmeta, survival, survRM2, shiny, plotly, ggplot2, dplyr
```

**New Additions Needed**:
- `cmdstanr` (major - Bayesian engine)
- `glmnet` (minor - already common)
- `pls` (minor)
- `sl3` (optional - advanced)

**Recommendation**:
- Make `cmdstanr` **Suggested** (not Required)
- Provide `freq` engine as default
- Add `glmnet` and `pls` to **Suggests**

---

## Stan Code Integration

**Challenge**: SurroNMA uses Stan via cmdstanr (lines 333-466)

**Options**:

1. **Include Stan code** (RECOMMENDED):
   - Add `.stan` files to `inst/stan/`
   - Compile on first use with caching
   - Fallback to frequentist if Stan unavailable

2. **JAGS alternative**:
   - Rewrite Stan models as JAGS
   - Use `gemtc` (already in powerNMA dependencies)
   - Slower but more accessible

3. **Frequentist only**:
   - Only implement `surro_nma_freq()`
   - Document Bayesian as future work

**Recommendation**: Option 1 - Include Stan with graceful degradation

---

## Functional Overlap with PowerNMA

| Feature | PowerNMA Current | SurroNMA | Integration Approach |
|---------|------------------|----------|---------------------|
| **Multivariate NMA** | ✅ Exists | ✅ Bivariate S+T | Enhance existing with surrogate mode |
| **Ranking** | ✅ P-scores | ✅ SUCRA/POTH/MID | Add POTH and MID options |
| **Inconsistency** | ✅ Node-split | ✅ Surrogate-aware | Enhance with S-T correlation |
| **Plotting** | ✅ Network/Forest | ✅ Surrogacy/STE | Add new surrogate plots |
| **Reporting** | ✅ HTML exports | ✅ Auto-narrative | Merge capabilities |

**Minimal Overlap** - Most features are net new additions!

---

## Integration Timeline

### Sprint 1 (Week 1): Core Functions
- Day 1-2: Adapt `surro_network()` and data structures
- Day 3: Implement `train_surrogate_index()` (OLS, Ridge)
- Day 4: Implement `fit_bivariate_nma_freq()` (frequentist engine)
- Day 5: Implement diagnostics and STE

### Sprint 2 (Week 2): Bayesian & Advanced
- Day 1-2: Stan code integration for Bayesian engine
- Day 3: Stress analysis and node-split
- Day 4: Enhanced ranking (POTH, MID)
- Day 5: Plotting functions

### Sprint 3 (Week 3): Polish & Test
- Day 1-2: Comprehensive tests
- Day 3-4: Vignette and documentation
- Day 5: Integration with Shiny GUI

**Total Timeline**: 3 weeks full-time (or 6 weeks part-time)

---

## Testing Strategy

### Unit Tests (testthat)

```r
test_that("surrogate network construction works", {
  data <- simulate_surrogate_data()
  net <- build_surrogate_network(data, "S", "T")
  expect_s3_class(net, "surro_network")
  expect_equal(net$K, 5)  # 5 treatments
})

test_that("SI training methods work", {
  net <- build_test_network()
  si_ols <- train_surrogate_index(net, method = "ols")
  si_ridge <- train_surrogate_index(net, method = "ridge")

  expect_s3_class(si_ols, "surro_index")
  expect_true(si_ridge$cv_r2 >= 0 && si_ridge$cv_r2 <= 1)
})

test_that("bivariate NMA freq works", {
  net <- build_test_network()
  fit <- fit_bivariate_nma_freq(net)
  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$beta))
})

test_that("STE calculation is correct", {
  # Known test case
  alpha <- c(0, 0.1, -0.1)  # 1000 draws
  beta <- c(0.8, 0.9, 0.7)
  ste <- compute_STE(alpha, beta, threshold_T = 0)
  expect_equal(ste$summary["mean"], 0, tolerance = 0.1)
})
```

### Integration Tests

Test integration with existing powerNMA functions:
- `auto_experimental_nma()` should offer surrogate option
- Shiny GUI should have surrogate tab
- Results should be exportable

---

## Documentation Requirements

### Main Documentation

**File**: `powerNMA/vignettes/surrogate-endpoints.Rmd`

**Sections**:
1. Introduction to Surrogate Endpoints
2. When to Use Surrogates (decision tree)
3. Data Requirements and Format
4. Surrogate Index Training (multiple surrogates)
5. Bivariate NMA (joint S+T modeling)
6. Surrogacy Validation (R², α, β)
7. Surrogate Threshold Effect (STE)
8. Stress Analysis and Sensitivity
9. Clinical Interpretation
10. Worked Example with Real Data

### Function Documentation

All functions need:
- `@param` for all arguments
- `@return` describing output structure
- `@examples` with working code
- `@references` to surrogate literature

**Key Papers to Cite**:
- Buyse M, et al. (2000). "The validation of surrogate endpoints"
- Daniels MJ, Hughes MD (1997). "Meta-analysis for the evaluation of potential surrogate markers"
- Burzykowski T, et al. (2005). Book: "The Evaluation of Surrogate Endpoints"
- Ciani O, et al. (2015). "Use of surrogate end points in healthcare policy"

---

## Shiny GUI Integration

### New Tab: "Surrogate Analysis"

**UI Elements**:
```r
tabPanel("Surrogate Analysis",
  sidebarLayout(
    sidebarPanel(
      selectInput("surrogate_var", "Surrogate Endpoint", ...),
      selectInput("true_var", "True Endpoint", ...),
      radioButtons("surrogate_method", "Method",
                   choices = c("Single Surrogate", "Train SI")),
      conditionalPanel(
        condition = "input.surrogate_method == 'Train SI'",
        checkboxGroupInput("surrogates", "Select Surrogates", ...)
      ),
      selectInput("engine", "Engine", c("Frequentist", "Bayesian")),
      actionButton("run_surrogate", "Run Analysis")
    ),
    mainPanel(
      h3("Surrogacy Diagnostics"),
      verbatimTextOutput("diagnostics_output"),
      plotOutput("surrogacy_plot"),

      h3("Surrogate Threshold Effect"),
      verbatimTextOutput("ste_output"),
      plotOutput("ste_plot"),

      h3("Treatment Rankings"),
      plotOutput("ranking_plot")
    )
  )
)
```

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Stan compilation issues | MEDIUM | HIGH | Provide frequentist alternative |
| Package dependency conflicts | LOW | MEDIUM | Use Suggests instead of Depends |
| Complex code difficult to maintain | MEDIUM | MEDIUM | Comprehensive documentation |
| Performance slow for large networks | LOW | MEDIUM | Optimize bootstrap code |

### Scientific Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Misuse of weak surrogates | MEDIUM | HIGH | Clear validation thresholds in docs |
| Over-interpretation of STE | MEDIUM | MEDIUM | Stress analysis to show sensitivity |
| Confusion about when to use | LOW | MEDIUM | Decision tree in vignette |

---

## Success Criteria

### Minimum Viable Product (MVP)

**Must Have**:
1. ✅ Build surrogate network from data
2. ✅ Train SI with at least OLS and Ridge
3. ✅ Fit bivariate NMA (frequentist engine minimum)
4. ✅ Compute surrogacy diagnostics (R², α, β)
5. ✅ Calculate STE with credible intervals
6. ✅ Basic plots (surrogacy scatter, STE)
7. ✅ S3 methods (print, summary, plot)
8. ✅ 10+ unit tests
9. ✅ Vignette with worked example

**Nice to Have**:
- Bayesian engine with Stan
- POTH and MID-adjusted ranking
- Stress analysis
- Node-split for surrogates
- Shiny GUI integration

### Acceptance Criteria

**Functional**:
- [ ] All core functions work without errors
- [ ] Reproduces results from surroNMA examples
- [ ] R² calculations match published values
- [ ] STE interpretation is clinically meaningful

**Technical**:
- [ ] Passes `devtools::check()` with 0 errors, 0 warnings
- [ ] Test coverage >80% for surrogate functions
- [ ] Documentation complete and clear
- [ ] No new package dependencies in Depends (only Suggests)

**User Experience**:
- [ ] Function names intuitive
- [ ] Error messages helpful
- [ ] Vignette understandable by non-experts
- [ ] Results tables well-formatted

---

## Next Steps

### Immediate Actions (Now)

1. **Create Integration Branch**:
   ```bash
   cd /home/user/rmstnma
   git checkout -b feature/surrogate-nma-integration
   ```

2. **Copy SurroNMA Source**:
   ```bash
   cp /home/user/surroNMA/surroNMA \
      powerNMA/R/experimental_surrogate_nma_source.R
   # For reference during adaptation
   ```

3. **Start with Data Structure**:
   - Adapt `surro_network()` function
   - Create `build_surrogate_network()` helper
   - Test with simulated data

### Questions for User

Before proceeding with integration:

1. **Timeline Priority?**
   - Immediate (pre-publication)?
   - Post-publication (version 3.1)?
   - Low priority (version 3.2+)?

2. **Feature Scope?**
   - MVP only (core functions)?
   - Full integration (all features)?
   - Phased (MVP now, advanced later)?

3. **Engine Preference?**
   - Frequentist only (faster implementation)?
   - Both Bayesian and Frequentist (complete)?
   - Bayesian priority (more powerful)?

4. **Testing Rigor?**
   - Basic tests (MVP)?
   - Comprehensive tests (publication-ready)?
   - Validation against published data?

5. **Documentation Level?**
   - Function docs only?
   - Plus vignette?
   - Plus Shiny integration?

---

## Recommendation

**Phased Approach** (RECOMMENDED):

**Phase A - Pre-Publication** (Optional, 1 week):
- Document surroNMA availability
- Add reference to related work
- Note integration planned for v3.1

**Phase B - Post-Publication** (Weeks 1-3 after publication):
- Full integration as described above
- Comprehensive testing
- Vignette creation

**Phase C - Future Enhancement** (Months 2-3):
- Shiny GUI integration
- Advanced features (stress, node-split)
- Publication of methodology paper

**Rationale**: Don't delay current publication. Surrogate features are substantial and deserve dedicated testing and validation.

---

**Status**: ✅ ANALYSIS COMPLETE - AWAITING USER DECISION
**Date**: November 1, 2025
**Recommendation**: Proceed with Post-Publication Integration (Phase B)
