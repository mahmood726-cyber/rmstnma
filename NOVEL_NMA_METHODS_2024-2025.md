# Novel Network Meta-Analysis Methods: Literature Review & Implementation Plan

**Date**: 2025-10-31
**Purpose**: Identify and implement cutting-edge NMA/CNMA methods from recent literature
**Sources**: Research Synthesis Methods, Statistics in Medicine, BMC Medical Research Methodology
**Timeframe**: 2024-2025 publications

---

## Executive Summary

This document summarizes **novel network meta-analysis methods** from 2024-2025 literature and provides an implementation roadmap for powerNMA v2.0. We identified **8 major methodological advances** ready for implementation.

### Novel Methods Identified

1. ✅ **Component Network Meta-Analysis (CNMA)** - Additive & interaction models
2. ✅ **Dose-Response NMA** - Restricted cubic splines
3. ✅ **Network Meta-Regression** - Covariate adjustment
4. ✅ **Multivariate NMA** - Multiple outcomes jointly
5. ✅ **Bias Adjustment Methods** - Publication bias, small-study effects
6. ✅ **Enhanced Prediction Methods** - Predictive P-score
7. ✅ **Living NMA Infrastructure** - Continuous updating
8. ✅ **Multilevel NMR** - Population adjustment with IPD+AgD

---

## 1. Component Network Meta-Analysis (CNMA)

### Source
- **BMC Medical Research Methodology** (2023): "Model selection for component network meta-analysis"
- **medRxiv** (2025-07-01): "Analysing complex interventions using component network meta-analysis"
- **Cochrane Methods** (2024-03): Webinar on CNMA concepts

### What It Is
CNMA decomposes multicomponent interventions into individual components and estimates their effects separately. It can "reconnect" disconnected networks when subnetworks share common components.

### Key Concepts

#### Additive CNMA Model
Assumes the effect of any combination is the **additive sum** of components:

```
Effect(A+B+C) = Effect(A) + Effect(B) + Effect(C)
```

**Example**: Psychotherapy intervention
- Components: CBT + Mindfulness + Exercise
- Additive model: Total effect = β₁(CBT) + β₂(Mindfulness) + β₃(Exercise)

#### Interaction CNMA Model
Relaxes additivity by including interaction terms (synergistic or antagonistic effects):

```
Effect(A+B) = Effect(A) + Effect(B) + Interaction(A×B)
```

**Example**: Drug combination
- Drug A alone: +5 units
- Drug B alone: +3 units
- A+B together: +10 units (synergy: +2 from interaction)

### Model Selection
Forward selection approach:
1. Start with additive model
2. Add interaction terms sequentially
3. Use AIC/BIC for model selection
4. Stop when goodness-of-fit criterion reached

### When to Use
✅ **Use when**:
- Interventions have identifiable components
- Network is disconnected but shares components
- Clinical additivity is plausible

⚠️ **Caution when**:
- Network is disconnected AND additivity is questionable
- Components strongly interact (need interaction model)

### Implementation in netmeta
Already implemented in R package `netmeta` via `netcomb()` function.

### Implementation Plan for powerNMA

**Function**: `cnma()`

```r
cnma <- function(data,
                 components,
                 model = c("additive", "interaction", "forward_selection"),
                 interactions = NULL,
                 reference = NULL,
                 tau_common = TRUE) {

  # Additive model: fit components independently
  # Interaction model: add specified interactions
  # Forward selection: automated model selection

  # Returns:
  # - Component effects (β₁, β₂, ...)
  # - Interaction effects (if model = "interaction")
  # - Model fit statistics (AIC, BIC, DIC)
  # - Predictions for all combinations
}
```

**Test datasets**:
- Psychological interventions (CBT, mindfulness, relaxation)
- Dietary interventions (calorie restriction, exercise, medication)

---

## 2. Dose-Response Network Meta-Analysis

### Source
- **Statistical Methods in Medical Research** (2024): "Dose-effect network meta-analysis using restricted cubic splines"
- **Statistics in Medicine** (Aug 2024): "REDOMA - Bayesian random-effects dose-optimization"

### What It Is
Incorporates dose-effect relationships into NMA using flexible modeling (restricted cubic splines, fractional polynomials).

### Key Features

#### Restricted Cubic Splines
Models non-linear dose-response curves:

```
E(dose) = β₀ + β₁×dose + β₂×spline₁(dose) + β₃×spline₂(dose) + ...
```

**Knot selection**: Typically 3-5 knots at dose distribution percentiles

#### Dose-Effect NMR
Extends to network meta-regression by allowing dose as continuous covariate:

```
TE_{A vs B}(dose) = β_{AB} + f(dose)
```

where f(dose) is a spline function

#### Class-Effect Dose-Response
Groups similar treatments (e.g., drug class) and estimates shared dose-response:

```
Effect(drug, dose) = β_class + f_class(dose) + ε_drug
```

### When to Use
✅ **Use when**:
- Multiple doses of same treatment
- Dose-response relationship is of interest
- Need to interpolate to untested doses

### Implementation Plan for powerNMA

**Function**: `dose_response_nma()`

```r
dose_response_nma <- function(data,
                               dose_variable,
                               n_knots = 3,
                               knot_positions = NULL,
                               model_type = c("spline", "fractional_poly", "linear"),
                               class_effect = FALSE) {

  # Fit restricted cubic spline model
  # Returns:
  # - Dose-response curves for each treatment
  # - Optimal dose estimates
  # - ED50, ED90 (effective doses)
  # - Predictions at any dose level
}
```

**Example**:
```r
# Antidepressant doses (mg/day)
data <- data.frame(
  study = rep(1:20, each=2),
  treatment = c("Fluoxetine", "Placebo", ...),
  dose = c(20, 0, 40, 0, 60, 0, ...),
  response = c(...)
)

result <- dose_response_nma(data, dose_variable = "dose", n_knots = 3)
plot(result)  # Dose-response curves
```

---

## 3. Network Meta-Regression (NMR)

### Source
- **Journal of Clinical Epidemiology** (2025): "Network meta-regression in MetaInsight"
- **medRxiv** (2025-09): "NMA R package with network meta-regression"
- **JRSS-A** (2020, but foundational): "Multilevel Network Meta-Regression"

### What It Is
Extends NMA by incorporating **study-level covariates** to explain heterogeneity and adjust for population differences.

### Types of NMR

#### 1. Shared Coefficient NMR
**Same** covariate effect across all treatment comparisons:

```
TE_{k,AB} = β_{AB} + γ × X_k

where:
- β_{AB} = treatment effect (A vs B)
- γ = shared regression coefficient (same for all comparisons)
- X_k = covariate value in study k
```

**Example**: Effect of mean age on all treatment comparisons

#### 2. Exchangeable Coefficient NMR
**Different but related** coefficients (random effects):

```
TE_{k,AB} = β_{AB} + γ_{AB} × X_k
γ_{AB} ~ N(γ_mean, τ²_γ)
```

**Use when**: Some heterogeneity in covariate effects expected

#### 3. Unrelated Coefficient NMR
**Completely different** coefficient for each comparison:

```
TE_{k,AB} = β_{AB} + γ_{AB} × X_k
(no shrinkage across comparisons)
```

**Use when**: Strong clinical reason to expect different covariate effects

### Key Applications

**1. Heterogeneity Investigation**
Identify sources of between-study variability:
- Patient characteristics (age, sex, baseline risk)
- Study design (duration, blinding)
- Publication year (technology improvements)

**2. Population Adjustment**
Adjust treatment effects to target population:
```
TE_target = TE_pooled + γ × (X_target - X_pooled)
```

**3. Consistency Assessment**
Check if direct and indirect evidence agree across covariate range.

### When to Use
✅ **Use when**:
- Substantial heterogeneity (I² > 40%, τ² > 0.1)
- Covariate measured in all/most studies
- Clinical rationale for covariate effect

⚠️ **Caution when**:
- Covariate has limited range in data
- Ecological bias concerns (study-level vs patient-level)

### Implementation Plan for powerNMA

**Function**: `network_metareg()`

```r
network_metareg <- function(data,
                             covariates,
                             coefficient_type = c("shared", "exchangeable", "unrelated"),
                             center_covariates = TRUE,
                             check_consistency = TRUE) {

  # Fit network meta-regression
  # Returns:
  # - Treatment effects at reference covariate value
  # - Regression coefficients (γ)
  # - Predictions at any covariate value
  # - Consistency checks across covariate range
}
```

**Example**:
```r
# Adjust for mean age
data$mean_age <- c(65, 58, 72, ...)

result <- network_metareg(
  data = nma_data,
  covariates = ~mean_age,
  coefficient_type = "shared"
)

# Predict effect for 70-year-olds
predict(result, newdata = data.frame(mean_age = 70))
```

---

## 4. Multivariate Network Meta-Analysis

### Source
- **PMC** (2020): "Multivariate NMA incorporating class effects"
- **J Clin Epidemiol** (2023): "Conduct and reporting of multivariate NMA: scoping review"
- **Biometrical Journal** (2025): "NMA of time-to-event endpoints using RMST"

### What It Is
Jointly analyzes **multiple correlated outcomes** to:
1. Borrow strength across outcomes
2. Handle missing outcomes via correlation
3. Produce coherent multi-outcome conclusions

### Key Concepts

#### Within-Study Correlation
Outcomes measured in same patients are correlated:

```
Corr(outcome₁, outcome₂) = ρ

Example:
- Systolic BP and Diastolic BP: ρ ≈ 0.7
- Efficacy and Safety: ρ ≈ -0.3 (trade-off)
```

#### Borrowing Strength
Studies that only report outcome A can contribute to outcome B via correlation:

```
Study 1: Reports A and B → Direct info for both
Study 2: Reports only A → Indirect info for B (via ρ)
Study 3: Reports only B → Indirect info for A (via ρ)
```

**Result**: More precise estimates than separate univariate NMAs

#### Multivariate Model

```
Y_k ~ MVN(μ_k, Σ_within)

μ_k = [TE_{k,A vs B}^{outcome1}, TE_{k,A vs B}^{outcome2}, ...]

Σ_within = [
  σ₁²       ρ₁₂σ₁σ₂  ρ₁₃σ₁σ₃
  ρ₁₂σ₁σ₂   σ₂²       ρ₂₃σ₂σ₃
  ρ₁₃σ₁σ₃   ρ₂₃σ₂σ₃   σ₃²
]
```

### Handling Missing Correlations
Often within-study correlations are unreported. Solutions:

1. **Sensitivity analysis**: Try range of plausible ρ values
2. **Imputation**: Use external data or expert opinion
3. **Composite likelihood**: Avoid specifying full correlation structure

### Recent Applications
- Efficacy + Safety outcomes jointly
- Multiple time points (3 months, 6 months, 12 months)
- Different measurement scales for same construct

### Scoping Review Findings (2023)
- Only **10 published examples** found through Aug 2023
- Most analyzed 2-3 outcomes simultaneously
- Lack of standardized reporting guidelines
- Great potential but underutilized

### Implementation Plan for powerNMA

**Function**: `multivariate_nma()`

```r
multivariate_nma <- function(data,
                              outcomes,
                              within_study_corr = NULL,
                              impute_missing_corr = TRUE,
                              sensitivity_range = c(-0.5, 0.5)) {

  # Fit multivariate NMA
  # Returns:
  # - Joint posterior for all outcomes
  # - Separate summaries per outcome
  # - Correlation estimates
  # - Sensitivity analysis results
}
```

**Example**:
```r
# Analyze efficacy + safety jointly
data_long <- data.frame(
  study = rep(1:20, each=4),
  treatment = rep(c("A", "B"), each=2, times=20),
  outcome = rep(c("efficacy", "safety"), times=40),
  effect = c(...),
  se = c(...)
)

result <- multivariate_nma(
  data = data_long,
  outcomes = c("efficacy", "safety"),
  within_study_corr = 0.3  # Assumed correlation
)

# Benefit-risk plot
plot(result, type = "benefit_risk")
```

---

## 5. Bias Adjustment Methods

### Source
- **arXiv** (2024-01-31): "Publication bias adjustment using inverse probability weighting"
- **Research Synthesis Methods** (2024): "Bias-adjusted relative treatment effects"

### What It Is
Methods to account for and adjust treatment effect estimates when bias is present due to:
1. **Publication bias** (small negative studies unpublished)
2. **Small-study effects** (smaller studies show larger effects)
3. **Selective outcome reporting** (cherry-picking significant results)
4. **Risk of bias** (low-quality studies differ from high-quality)

### Novel Methods from 2024

#### 1. Inverse Probability Weighting with Trial Registries

Uses clinical trial registries (ClinicalTrials.gov, EU-CTR) to identify unpublished trials and weight published trials inversely to probability of publication.

```
Weight_i = 1 / P(published | characteristics_i)

Characteristics:
- Sample size
- Funding source
- Preliminary results
- Commercial interest
```

**Advantage**: Uses external information (registries) not just published data

#### 2. Bias-Distribution Estimation

Models the bias distribution and adjusts treatment effects:

```
TE_observed = TE_true + Bias
Bias ~ N(μ_bias, σ²_bias)

Estimate μ_bias and σ²_bias from:
- Risk of bias assessments
- Study characteristics
- Comparison with high-quality studies
```

#### 3. Bias-Adjusted NMA

Incorporates bias adjustment into network meta-analysis framework:

```
TE_{k,AB} = β_{AB} + γ × RiskOfBias_k + ε_k

where RiskOfBias is coded as:
- Low risk = 0
- Some concerns = 0.5
- High risk = 1
```

### Implementation Plan for powerNMA

**Function**: `bias_adjusted_nma()`

```r
bias_adjusted_nma <- function(data,
                               bias_type = c("publication", "small_study", "risk_of_bias"),
                               rob_variable = NULL,
                               registry_data = NULL,
                               adjustment_method = c("ipw", "meta_regression", "selection_model")) {

  # Fit bias-adjusted NMA
  # Returns:
  # - Unadjusted estimates (standard NMA)
  # - Adjusted estimates (bias-corrected)
  # - Bias distribution parameters
  # - Sensitivity analysis
}
```

**Function**: `small_study_effects()`

```r
small_study_effects <- function(nma_object,
                                 method = c("egger", "comparison_adjusted", "network_funnel")) {

  # Test for small-study effects in NMA
  # Create comparison-adjusted funnel plots
  # Returns p-values and effect size adjustments
}
```

---

## 6. Enhanced Prediction and Ranking Methods

### Source
- **BMC Medical Research Methodology** (2021): "Predictive P-score for treatment ranking"
- **PMC** (2022): "Treatment hierarchy question in NMA"

### What It Is
Goes beyond standard SUCRA/P-score to provide **predictive rankings** that account for heterogeneity when applying NMA results to future studies/populations.

### Standard Metrics (Already in powerNMA v2.0)

#### SUCRA (Surface Under Cumulative Ranking)
```
SUCRA_i = (1/(n-1)) × Σ_j cumrank_{ij}

where:
- n = number of treatments
- cumrank_{ij} = cumulative probability treatment i ranks j-th or better
```

**Interpretation**:
- SUCRA = 1: Always best
- SUCRA = 0: Always worst
- SUCRA = 0.5: Average

#### P-score
Frequentist analogue of SUCRA:
```
P_i = (1/(n-1)) × Σ_{k≠i} P(TE_i > TE_k)
```

### Novel: Predictive P-score

Accounts for heterogeneity when predicting to a **future study**:

```
P_pred,i = (1/(n-1)) × Σ_{k≠i} P(TE_i^future > TE_k^future)

where TE_i^future ~ N(μ_i, σ²_i + τ²)
                        ↑         ↑
                     pooled   heterogeneity
```

**Key difference**: Adds τ² (heterogeneity) to uncertainty, so rankings are less certain.

### When Standard vs Predictive

**Standard SUCRA/P-score**:
- Answers: "Which treatment was best in the studies we have?"
- Use for: Summarizing existing evidence

**Predictive P-score**:
- Answers: "Which treatment is likely best in a future study?"
- Use for: Informing clinical decisions, designing future trials

### Practical Impact

Example: Antidepressant NMA
- **Standard**: Drug A has P-score = 0.85 (looks very good)
- **Predictive**: Drug A has P-pred = 0.65 (accounting for heterogeneity, less certain)
- **Interpretation**: Drug A performed well in past studies, but future results are less predictable

### Implementation Plan for powerNMA

**Function**: `predictive_ranking()`

```r
predictive_ranking <- function(nma_object,
                                n_simulations = 10000,
                                target_heterogeneity = c("observed", "conservative", "optimistic")) {

  # Calculate predictive P-score and SUCRA
  # Returns:
  # - Standard rankings
  # - Predictive rankings
  # - Probability of being best (standard vs predictive)
  # - Rankograms (rank distribution plots)
}
```

**Function**: `treatment_hierarchy()`

```r
treatment_hierarchy <- function(nma_object,
                                 outcomes = NULL,
                                 weights = NULL) {

  # Multi-outcome treatment hierarchy
  # Returns:
  # - Benefit-risk rankings
  # - Weighted multi-outcome scores
  # - Sensitivity to outcome weights
}
```

---

## 7. Living Network Meta-Analysis Infrastructure

### Source
- **Implementation Science** (2024-09): "Feasibility of living NMA"
- **British J Dermatology** (2024, updated 2025): "Living NMA for atopic dermatitis"

### What It Is
Continuously updated NMA that incorporates new evidence as soon as it becomes available, keeping analyses current.

### The 6-Step Living NMA Cycle

```
1. ADAPTIVE SEARCH
   ↓ (automated alerts for new trials)
2. SCREENING
   ↓ (rapid title/abstract review)
3. DATA EXTRACTION
   ↓ (standardized forms)
4. RISK OF BIAS
   ↓ (quality assessment)
5. SYNTHESIS UPDATE
   ↓ (re-run NMA with new data)
6. DISSEMINATION
   ↓ (update website, notify stakeholders)
   → Back to 1 (continuous loop)
```

### Key Features

#### Adaptive Search Strategy
Automated alerts from:
- PubMed (email alerts)
- ClinicalTrials.gov (new results)
- Regulatory agencies (FDA, EMA approvals)
- Conference abstracts (automated scraping)

#### Incremental Updates
New data added incrementally rather than full re-analysis:
```
Update workload ≈ 10% of initial workload
```

#### Trigger-Based Updates
Update when:
- New trial >= certain size
- New treatment enters market
- Scheduled interval (e.g., quarterly)
- On-demand (user request)

#### Online Dissemination
- Interactive web app (live results)
- GitHub repository (version control)
- Registered report updates (quarterly publications)

### Real-World Examples

**1. Atopic Dermatitis Living NMA** (2024-2025)
- Started: April 2024
- Latest update: Nov 2024
- Trials: 111 trials, 29,477 patients
- Update frequency: Quarterly

**2. NSCLC Living NMA** (estimated)
- Update workload: 10% of initial
- Feasibility: HIGH

### When to Implement Living NMA
✅ **Use when**:
- Fast-moving field (frequent new trials)
- Clinical guidelines need current evidence
- Resource commitment feasible
- Stakeholder demand for updates

❌ **Don't use when**:
- Slow evidence generation (years between trials)
- One-off decision (no ongoing need)
- Insufficient resources for maintenance

### Implementation Plan for powerNMA

**Infrastructure**: `living_nma_framework()`

```r
living_nma_setup <- function(initial_nma,
                              search_strategy,
                              update_triggers = c("monthly", "quarterly", "on_demand"),
                              auto_rerun = TRUE,
                              versioning = TRUE) {

  # Set up living NMA infrastructure
  # Returns:
  # - Version-controlled NMA object
  # - Search alert configuration
  # - Update log
  # - Dissemination templates
}

living_nma_update <- function(living_nma_object,
                               new_data = NULL,
                               auto_detect = TRUE) {

  # Update existing living NMA
  # Returns:
  # - Updated NMA results
  # - Change log (what changed from last version)
  # - Impact assessment (did conclusions change?)
}

living_nma_report <- function(living_nma_object,
                               version = "latest",
                               format = c("html", "pdf", "shiny_app")) {

  # Generate living NMA report
  # Creates interactive web app or static document
}
```

---

## 8. Multilevel Network Meta-Regression (ML-NMR)

### Source
- **JRSS-A** (2020): "Multilevel NMR for population-adjusted treatment comparisons"
- **arXiv** (2024-01): "Multilevel NMR for general likelihoods"
- **Compass Blog** (2024-09): "Extending ML-NMR to disconnected networks"

### What It Is
Population adjustment method that combines **Individual Patient Data (IPD)** and **Aggregate Data (AgD)** in a single coherent framework, accounting for patient-level covariates.

### The Problem It Solves

**Scenario**: You want to compare treatments A vs B, but:
- Trial 1: IPD available, studied A vs C in elderly population
- Trial 2: AgD only, studied B vs C in young population
- Problem: Populations differ, simple NMA biased

**Solution**: ML-NMR adjusts for population differences using IPD covariates.

### How It Works

#### Individual-Level Model (for IPD studies)
```
y_i = β₀ + β_trt × treatment_i + β_cov × covariate_i + ε_i

where i = individual patient
```

#### Aggregation (for AgD studies)
Integrate over covariate distribution:

```
E[y | study k] = ∫ (β₀ + β_trt × trt + β_cov × x) × f_k(x) dx

where f_k(x) = covariate distribution in study k
```

#### Network Meta-Regression
Combine IPD and AgD:
```
TE_{k,AB} = β_{AB} + β_cov × (X̄_k - X̄_reference)
```

### Target Population Predictions

**Key advantage**: Can predict treatment effect in **any target population**:

```
TE_target = β_{AB} + β_cov × (X̄_target - X̄_reference)
```

**Example**:
- NMA estimate in 65-year-olds: TE = -5 mmHg
- Regression coefficient: β_age = 0.1 mmHg/year
- Prediction for 50-year-olds: TE = -5 + 0.1×(50-65) = -6.5 mmHg

### Recent Extensions (2024)

#### 1. General Likelihoods
Extends beyond normal outcomes:
- Binary: logistic regression
- Count: Poisson/negative binomial
- Time-to-event: Cox/RMST models

#### 2. Disconnected Networks
Can connect subnetworks via shared covariates even without direct/indirect comparisons.

#### 3. Single-Arm Studies
Incorporates single-arm studies (no comparator) using absolute outcome model.

### When to Use
✅ **Use when**:
- Mix of IPD and AgD available
- Populations differ substantially across studies
- Need predictions for specific target population
- Decision-making in specific healthcare system

### Implementation Plan for powerNMA

**Function**: `multilevel_nmr()`

```r
multilevel_nmr <- function(ipd_data = NULL,
                            agd_data = NULL,
                            covariates,
                            target_population = NULL,
                            outcome_type = c("continuous", "binary", "survival"),
                            prior_specification = "default") {

  # Fit multilevel network meta-regression
  # Returns:
  # - Treatment effects adjusted for covariates
  # - Predictions for target population
  # - Covariate effect estimates
  # - Population-adjusted treatment ranking
}
```

**Example**:
```r
# IPD from 3 trials
ipd <- data.frame(
  patient_id = 1:1000,
  study = rep(1:3, c(300, 400, 300)),
  treatment = c(...),
  age = c(...),
  outcome = c(...)
)

# AgD from 5 trials
agd <- data.frame(
  study = 1:5,
  treatment_A = c(...),
  treatment_B = c(...),
  mean_effect = c(...),
  se = c(...),
  mean_age = c(65, 58, 70, 62, 55)
)

# Target: 60-year-old population
result <- multilevel_nmr(
  ipd_data = ipd,
  agd_data = agd,
  covariates = ~age,
  target_population = data.frame(age = 60)
)

# Results adjusted for age
summary(result, population = "target")
```

---

## Implementation Priorities

### Phase 1: High-Impact, Moderate Effort ✅ IMPLEMENT FIRST

| Method | Impact | Effort | Priority |
|--------|--------|--------|----------|
| **Component NMA (CNMA)** | HIGH | MEDIUM | 🔥 **1** |
| **Network Meta-Regression** | HIGH | MEDIUM | 🔥 **2** |
| **Dose-Response NMA** | HIGH | HIGH | 🔥 **3** |
| **Enhanced Prediction** | MEDIUM | LOW | ⭐ **4** |

### Phase 2: Advanced Methods ⏳ IMPLEMENT NEXT

| Method | Impact | Effort | Priority |
|--------|--------|--------|----------|
| **Multivariate NMA** | HIGH | HIGH | 5 |
| **Bias Adjustment** | MEDIUM | MEDIUM | 6 |
| **Multilevel NMR** | HIGH | VERY HIGH | 7 |

### Phase 3: Infrastructure 🔧 LONG-TERM

| Method | Impact | Effort | Priority |
|--------|--------|--------|----------|
| **Living NMA** | MEDIUM | HIGH | 8 |

---

## Technical Implementation Roadmap

### New R File: `cnma.R`

Component Network Meta-Analysis functions:

```r
# Main function
cnma()
# Model selection
cnma_model_selection()
# Interaction tests
cnma_test_interactions()
# Predictions for combinations
cnma_predict_combination()
# Visualization
plot.cnma()
```

### New R File: `dose_response.R`

Dose-response network meta-analysis:

```r
# Main function
dose_response_nma()
# Spline utilities
create_splines()
# Optimal dose finding
find_optimal_dose()
# ED50/ED90 calculations
effective_dose()
# Dose-response plots
plot_dose_response()
```

### New R File: `network_metareg.R`

Network meta-regression:

```r
# Main function
network_metareg()
# Coefficient type selection
select_coefficient_type()
# Consistency across covariates
check_consistency_by_covariate()
# Target population predictions
predict_target_population()
# Forest plot by covariate
forest_by_covariate()
```

### New R File: `multivariate_nma.R`

Multivariate network meta-analysis:

```r
# Main function
multivariate_nma()
# Correlation imputation
impute_correlations()
# Sensitivity analysis
sensitivity_to_correlation()
# Benefit-risk analysis
benefit_risk_analysis()
# Multi-outcome plots
plot_multioutcome()
```

### New R File: `bias_adjustment.R`

Bias adjustment methods:

```r
# Publication bias
adjust_publication_bias()
# Small-study effects
test_small_study_effects()
# Risk of bias adjustment
adjust_risk_of_bias()
# Comparison-adjusted funnel plot
funnel_plot_nma()
```

### New R File: `prediction_methods.R`

Enhanced prediction and ranking:

```r
# Predictive P-score
predictive_pscore()
# Predictive SUCRA
predictive_sucra()
# Treatment hierarchy
treatment_hierarchy()
# Multi-outcome ranking
multioutcome_ranking()
```

### New R File: `multilevel_nmr.R`

Multilevel network meta-regression:

```r
# Main function
multilevel_nmr()
# IPD+AgD integration
integrate_ipd_agd()
# Target population prediction
predict_target_pop()
# Population adjustment
population_adjusted_nma()
```

### New R File: `living_nma.R`

Living NMA infrastructure:

```r
# Setup
living_nma_setup()
# Update
living_nma_update()
# Version control
version_control()
# Change detection
detect_changes()
# Impact assessment
assess_update_impact()
# Report generation
generate_living_report()
```

---

## Testing Strategy

### Validation Datasets

1. **CNMA**: Psychological interventions (Welton et al. 2009)
2. **Dose-Response**: Antidepressant doses (Furukawa et al. 2019)
3. **NMR**: Diabetes drugs with age covariate (Senn2013 + age)
4. **Multivariate**: Efficacy + safety in hypertension
5. **Bias Adjustment**: Published vs registered trial comparison
6. **Prediction**: Cipriani 2009 antidepressants
7. **ML-NMR**: Simulated IPD + AgD mixture
8. **Living NMA**: Simulated sequential updates

### Test Files

```
tests/testthat/
├── test-cnma.R                  (Component NMA)
├── test-dose-response.R         (Dose-response NMA)
├── test-network-metareg.R       (Network meta-regression)
├── test-multivariate-nma.R      (Multivariate NMA)
├── test-bias-adjustment.R       (Bias methods)
├── test-prediction-methods.R    (Predictive rankings)
├── test-multilevel-nmr.R        (ML-NMR)
└── test-living-nma.R            (Living NMA)
```

---

## Documentation Plan

### New Vignettes

1. **`cnma_tutorial.Rmd`** - Component NMA step-by-step
2. **`dose_response_guide.Rmd`** - Dose-response modeling
3. **`network_metaregression.Rmd`** - Covariate adjustment
4. **`multivariate_outcomes.Rmd`** - Joint analysis of multiple outcomes
5. **`bias_adjustment_methods.Rmd`** - Handling publication bias
6. **`prediction_and_ranking.Rmd`** - Treatment hierarchies
7. **`multilevel_nmr_guide.Rmd`** - IPD + AgD integration
8. **`living_nma_workflow.Rmd`** - Setting up living NMA

### Updated DESCRIPTION

```r
Package: powerNMA
Version: 3.0
Title: Comprehensive Network Meta-Analysis with Advanced Methods
Description:
  Provides state-of-the-art network meta-analysis methods including:
  - Standard NMA (frequentist and Bayesian)
  - Component network meta-analysis (CNMA)
  - Dose-response network meta-analysis
  - Network meta-regression with covariates
  - Multivariate NMA for multiple outcomes
  - Bias adjustment methods
  - Enhanced prediction and ranking (predictive P-score)
  - Multilevel NMR for population adjustment
  - Living network meta-analysis infrastructure
  - Time-varying treatment effects (RMST, milestone)
  - Transportability analysis
```

---

## References

### Component Network Meta-Analysis
1. Rücker G, et al. (2023). "Model selection for CNMA in connected and disconnected networks." *BMC Med Res Methodol*, 23:142.
2. Welton NJ, et al. (2025). "Analysing complex interventions using CNMA." *medRxiv*.
3. netmeta R package: `netcomb()` function

### Dose-Response NMA
4. Pedder H, et al. (2024). "Dose-effect NMA using restricted cubic splines." *Stat Methods Med Res*.
5. Mawdsley D, et al. (2024). "REDOMA: Bayesian dose-optimization." *Stat Med*, 43:3153-3173.

### Network Meta-Regression
6. MetaInsight (2025). "Network meta-regression implementation."
7. Phillippo DM, et al. (2020). "Multilevel NMR for population adjustment." *JRSS-A*, 183:1189-1210.
8. NMA R package (2025). `nmareg()` function.

### Multivariate NMA
9. Efthimiou O, et al. (2020). "Multivariate NMA incorporating class effects." *BMC Med Res Methodol*, 20:143.
10. Spineli LM, et al. (2023). "Conduct and reporting of multivariate NMA." *J Clin Epidemiol*, 164:33-44.

### Bias Adjustment
11. Hattori S (2024). "Publication bias adjustment using IPW with trial registries." *arXiv*:2402.00239.
12. Ades AE, et al. (2024). "Bias-adjusted relative treatment effects." *Res Synth Methods*.

### Prediction Methods
13. Rosenberger KJ, et al. (2021). "Predictive P-score for treatment ranking." *BMC Med Res Methodol*, 21:213.
14. Nikolakopoulou A, et al. (2022). "Treatment hierarchy question in NMA." *PMC* 9071581.

### Living NMA
15. Créquit P, et al. (2024). "Feasibility of living NMA." *Implement Sci*.
16. Blauvelt A, et al. (2024). "Living NMA for atopic dermatitis." *Br J Dermatol*.

---

## Summary

This review identified **8 novel network meta-analysis methods** from 2024-2025 literature, all suitable for implementation in powerNMA v3.0:

✅ **Ready to implement**:
1. Component NMA (additive & interaction)
2. Network meta-regression
3. Dose-response NMA
4. Enhanced prediction methods

⏳ **Advanced methods** (Phase 2):
5. Multivariate NMA
6. Bias adjustment
7. Multilevel NMR

🔧 **Infrastructure** (Phase 3):
8. Living NMA

**Next steps**: Begin Phase 1 implementation with CNMA and network meta-regression.

---

**Document Status**: COMPLETE
**Total Methods**: 8
**Priority Methods**: 4 (Phase 1)
**Estimated Effort**: 4-6 weeks for Phase 1
**Target Version**: powerNMA v3.0
