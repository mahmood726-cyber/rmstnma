# powerNMA: Four-Pathway Analysis System

**Version**: 3.0
**Date**: October 31, 2025
**Status**: Complete implementation of 4 distinct NMA pathways

---

## Overview: The Four-Pathway System

powerNMA provides **4 distinct analysis pathways** for network meta-analysis, designed to meet different research needs, expertise levels, and standardization requirements:

```
┌─────────────────────────────────────────────────────────────────┐
│                     powerNMA Analysis Pathways                   │
├─────────────────────────┬────────────────────────────────────────┤
│   STANDARD METHODS      │        EXPERIMENTAL METHODS            │
│   (Phase 1-2, 2019-2023)│     (2024-2025 Literature)            │
├─────────────────────────┼────────────────────────────────────────┤
│ 1. MANUAL STANDARD      │ 2. MANUAL EXPERIMENTAL                 │
│    User controls all    │    User controls experimental          │
│    choices              │    method selection                    │
│    ↓                    │    ↓                                   │
│    Full flexibility     │    Advanced capabilities               │
├─────────────────────────┼────────────────────────────────────────┤
│ 3. AUTO STANDARD        │ 4. AUTO EXPERIMENTAL                   │
│    All choices          │    All experimental choices            │
│    automated            │    automated                           │
│    ↓                    │    ↓                                   │
│    Standardized workflow│    Cutting-edge automation             │
└─────────────────────────┴────────────────────────────────────────┘
```

---

## Pathway 1: Manual Standard (Phase 1-2 Methods)

### Purpose
**Full user control** over all methodological choices using established, validated methods.

### When to Use
- ✅ Need complete control over analysis decisions
- ✅ Specific methodological requirements
- ✅ Regulatory submissions (FDA, EMA)
- ✅ Clinical guidelines (official)
- ✅ Standard systematic reviews
- ✅ Teaching/learning NMA methodology

### Methods Available (7 methods)

#### Phase 1 - Core Methods
1. **Standard NMA** (`netmeta` wrapper)
   - Basic network meta-analysis
   - Fixed and random effects

2. **Component NMA** (`cnma`)
   - Multicomponent interventions
   - Additive, interaction, forward selection models

3. **Network Meta-Regression** (`network_metareg`)
   - Covariate adjustment
   - Shared/exchangeable/unrelated coefficients

4. **Dose-Response NMA** (`dose_response_nma`)
   - Spline-based dose-response
   - ED50/ED90, optimal dose finding

5. **Enhanced Prediction** (`predictive_ranking`, `treatment_hierarchy`)
   - Predictive P-scores accounting for heterogeneity
   - Multi-outcome benefit-risk ranking

#### Phase 2 - Advanced Methods
6. **Multivariate NMA** (`multivariate_nma`)
   - Multiple correlated outcomes
   - Benefit-risk assessment

7. **Missing Data Handling** (`handle_missing_data`)
   - Pattern-mixture models
   - IMOR framework

8. **Cross-Design Synthesis** (`cross_design_synthesis`)
   - Combine RCT + observational data
   - Bias adjustment

### User Decisions Required
```r
# User must specify:
- Which method to use (standard NMA vs CNMA vs meta-regression, etc.)
- Model type (fixed vs random effects)
- Summary measure (MD, SMD, OR, RR)
- Reference treatment
- Heterogeneity model
- Covariates for meta-regression
- Component specification for CNMA
- Sensitivity analyses to run
- Missing data approach
```

### Example Usage

```r
library(powerNMA)

# Example: Manual Component NMA
result <- cnma(
  data = my_data,
  components = component_matrix,
  model = "interaction",  # USER CHOICE
  reference = "Placebo",  # USER CHOICE
  tau_common = TRUE       # USER CHOICE
)

# Example: Manual Network Meta-Regression
result_nmr <- network_metareg(
  data = my_data,
  covariates = c("age", "severity"),  # USER CHOICE
  coefficient_type = "shared",         # USER CHOICE
  center_covariates = TRUE            # USER CHOICE
)

# Example: Manual Multivariate NMA
result_mv <- multivariate_nma(
  data = multi_outcome_data,
  outcomes = c("efficacy", "safety"),  # USER CHOICE
  within_study_corr = 0.3,            # USER CHOICE
  method = "reml"                      # USER CHOICE
)
```

### Advantages
- Complete control
- Transparent decision-making
- Maximum flexibility
- Suitable for publications requiring detailed methodology

### Disadvantages
- Requires methodological expertise
- Time-consuming decision-making
- Potential for suboptimal choices
- Not standardized across users

---

## Pathway 2: Manual Experimental (2024-2025 Methods)

### Purpose
**Full user control** over cutting-edge experimental methods from latest literature.

### When to Use
- ✅ Advanced research questions
- ✅ Methodological innovation
- ✅ Precision medicine research
- ✅ Robustness assessment beyond GRADE
- ✅ Need specific experimental capability
- ❌ Routine clinical guidelines
- ❌ Regulatory submissions (without validation)

### Methods Available (4 experimental methods)

1. **RMST-based NMA** (`rmst_nma`)
   - Restricted Mean Survival Time for time-to-event
   - No proportional hazards assumption
   - Interpretable time units (months gained)
   - Source: Hua et al. (2025) *Biometrical Journal*

2. **Threshold Analysis** (`threshold_analysis`)
   - Quantitative robustness assessment
   - Robustness score (0-100)
   - Tipping point identification
   - Source: Ades et al. (2025) *Research Synthesis Methods*

3. **Individualized Treatment Rules** (`itr_from_nma`)
   - Personalized treatment selection
   - Effect modifier identification
   - Precision medicine from NMA
   - Source: Chen et al. (2025) *PLOS ONE*

4. **Bayesian Model Averaging** (`model_averaging_nma`)
   - Model uncertainty quantification
   - Honest uncertainty intervals
   - Structural uncertainty
   - Source: Phillippo et al. (2024) *Research Synthesis Methods*

### User Decisions Required

```r
# User must specify:
- Which experimental method(s) to use
- Method-specific parameters:
  * RMST: tau (restriction time)
  * Threshold: risk aversion, decision rule
  * ITR: ML algorithm, covariates
  * BMA: models to average, weighting method
- Data preparation for each method
- Interpretation approach
```

### Example Usage

```r
library(powerNMA)

# Example: Manual RMST-based NMA
result_rmst <- rmst_nma(
  data = survival_data,
  tau = 36,                  # USER CHOICE (3 years)
  data_type = "ipd",         # USER CHOICE
  method = "frequentist"     # USER CHOICE
)

# Example: Manual Threshold Analysis
result_threshold <- threshold_analysis(
  nma_object = nma_result,
  outcome_direction = "higher",    # USER CHOICE
  decision_rule = "maximize_benefit",  # USER CHOICE
  risk_aversion = 1.5,             # USER CHOICE (0-2)
  threshold_type = "all"           # USER CHOICE
)

# Example: Manual ITR
result_itr <- itr_from_nma(
  data = ipd_data,
  outcome_var = "depression_score",  # USER CHOICE
  covariate_vars = c("age", "severity", "sex"),  # USER CHOICE
  method = "machine_learning",       # USER CHOICE
  ml_algorithm = "random_forest"     # USER CHOICE
)

# Example: Manual Model Averaging
result_bma <- model_averaging_nma(
  data = pairwise_data,
  models = c("fixed_effect", "random_effects", "inconsistency_re"),  # USER CHOICE
  weighting = "BIC",                 # USER CHOICE (BIC/AIC/DIC/posterior)
  reference = "Placebo"              # USER CHOICE
)
```

### Advantages
- Access to cutting-edge methods
- Specific advanced capabilities
- Innovation in methodology
- Potential for novel insights

### Disadvantages
- Requires expert-level knowledge
- Limited external validation
- Not suitable for regulatory use
- Complex interpretation

---

## Pathway 3: Auto Standard (Automatic + Standard Methods)

### Purpose
**Fully automated workflow** using established methods with **all choices made automatically** based on data characteristics and best practices. Standardizes NMA across users.

### When to Use
- ✅ Want standardized analysis
- ✅ Limited methodological expertise
- ✅ Need reproducible workflow
- ✅ Time constraints
- ✅ Multiple similar analyses
- ✅ Quality control/consistency across projects
- ✅ Teaching standardized approach

### Function

```r
auto_standard_nma(
  data,                      # Only required parameter
  outcome_type = NULL,       # Auto-detected
  reference_treatment = NULL,  # Auto-selected
  components = NULL,         # Optional for CNMA
  verbose = TRUE            # Show automatic choices
)
```

### What Gets Automated

#### 1. **Data Detection (Automatic)**
- Data format (pairwise, arm-based, IPD)
- Outcome type (continuous, binary, time-to-event)
- Network structure (connected, disconnected, sparse)
- Number of studies and treatments
- Presence of covariates, multiple outcomes, missing data

#### 2. **Method Selection (Automatic)**
```r
Automatic Logic:
- Multicomponent interventions detected → Component NMA
- Multiple outcomes detected → Multivariate NMA
- IPD with covariates → Network Meta-Regression
- Standard comparisons → Standard NMA
- Missing data >10% → Pattern-mixture models
- RCT + Observational → Cross-design synthesis
```

#### 3. **Parameter Settings (Automatic)**
- **Model type**: Random effects (for generalizability)
- **Heterogeneity**: Common tau across comparisons
- **Summary measure**:
  - Continuous outcome → MD (Mean Difference)
  - Binary outcome → OR (Odds Ratio)
  - Time-to-event → HR (Hazard Ratio)
- **Reference treatment**: Most studied or alphabetically first
- **Inconsistency check**: If network has loops (≥3 treatments, ≥4 studies)

#### 4. **Sensitivity Analyses (Automatic)**
```r
Automatically runs:
- Fixed vs random effects comparison (always)
- Leave-one-out (if ≥8 studies)
- Missing data sensitivity (if >10% missing)
```

#### 5. **Outputs Generated (Automatic)**
- Primary NMA results
- Predictive rankings (with heterogeneity)
- Inconsistency assessment
- Meta-regression (if covariates available)
- Missing data handling (if needed)
- Comprehensive diagnostics
- Clinical recommendations
- Auto-generated report

### Example Usage

```r
library(powerNMA)

# Minimal usage - everything automatic!
result <- auto_standard_nma(data = my_network_data)

# The function automatically:
# 1. Detects: 15 studies, 5 treatments, binary outcome, pairwise format
# 2. Chooses: Standard random-effects NMA, OR as summary measure
# 3. Selects: "Placebo" as reference (most studied)
# 4. Runs: Primary analysis + fixed/random comparison + inconsistency check
# 5. Generates: Rankings, diagnostics, report, recommendations

print(result)
summary(result)
plot(result)

# With optional component specification
components <- data.frame(
  treatment = c("A", "B", "AB", "ABC"),
  cbt = c(0, 0, 1, 1),
  medication = c(0, 1, 1, 1),
  exercise = c(1, 0, 0, 1)
)

result_cnma <- auto_standard_nma(
  data = my_data,
  components = components  # Automatically detects CNMA and chooses model
)
# Automatically runs Component NMA with forward selection
```

### Output Structure

```r
result$data_characteristics    # What was detected
result$automatic_choices       # All choices made
result$primary_analysis        # Main results
result$sensitivity_analyses    # Sensitivity analyses
result$inconsistency           # Inconsistency assessment
result$prediction              # Predictive rankings
result$meta_regression         # Meta-regression (if applicable)
result$missing_data            # Missing data handling (if applicable)
result$diagnostics             # Comprehensive diagnostics
result$report                  # Auto-generated report
result$recommendations         # Clinical interpretation
```

### Advantages
- **No expertise required**: Makes optimal choices automatically
- **Standardized**: Same data → same analysis (reproducible)
- **Fast**: Minimal user input needed
- **Comprehensive**: Runs full analysis pipeline
- **Transparent**: Reports all choices made
- **Quality-controlled**: Uses best-practice defaults

### Disadvantages
- Less flexibility than manual
- May not match specific methodological preferences
- Automatic choices may differ from traditional approaches

---

## Pathway 4: Auto Experimental (Automatic + Experimental Methods)

### Purpose
**Fully automated workflow** using **cutting-edge experimental methods** (2024-2025) with all choices made automatically based on research question and data.

### When to Use
- ✅ Advanced research with limited time
- ✅ Want cutting-edge methods without deep expertise
- ✅ Exploratory precision medicine analysis
- ✅ Robustness assessment for complex decisions
- ✅ Methodological exploration
- ⚠️ **EXPERIMENTAL** - for research, not routine clinical use

### Function

```r
auto_experimental_nma(
  data,                              # Only required parameter
  outcome_type = NULL,               # Auto-detected
  research_question = "auto",        # Auto-selects methods
  risk_aversion = 1.0,              # Decision-making conservatism
  reference_treatment = NULL,        # Auto-selected
  verbose = TRUE                    # Show choices
)
```

### Research Question → Automatic Method Selection

```r
research_question = "auto" (default)
→ Analyzes data characteristics
→ Selects optimal experimental methods

research_question = "precision_medicine"
→ ITR + Model Averaging

research_question = "decision_making"
→ Threshold Analysis + Model Averaging

research_question = "survival_analysis"
→ RMST-based NMA + Threshold Analysis

research_question = "model_uncertainty"
→ Model Averaging + Threshold Analysis
```

### Automatic Selection Logic

#### If Time-to-Event Data Detected:
```r
Methods: RMST-based NMA + Threshold Analysis
Rationale:
- RMST provides interpretable survival differences
- Threshold analysis assesses robustness
Auto-selects: tau (restriction time) based on follow-up
```

#### If IPD with Covariates Detected:
```r
Methods: ITR + Model Averaging + Threshold Analysis
Rationale:
- ITR identifies effect modifiers for personalization
- Model Averaging handles structural uncertainty
- Threshold Analysis assesses recommendation robustness
Auto-selects: Regression-based ITR (interpretable)
```

#### If Pairwise Data Only:
```r
Methods: Threshold Analysis + Model Averaging
Rationale:
- Threshold Analysis quantifies robustness
- Model Averaging accounts for model uncertainty
Auto-selects: BIC weighting for model averaging
```

### What Gets Automated

#### 1. **Experimental Method Selection**
- Automatically chooses 1-3 experimental methods based on:
  - Data type (IPD vs aggregate)
  - Outcome type (time-to-event, continuous, binary)
  - Available covariates
  - Research question

#### 2. **Method-Specific Parameters**
- **RMST-NMA**: Auto-selects tau (restriction time)
- **Threshold Analysis**: Auto-sets risk aversion, selects all threshold types
- **ITR**: Auto-selects covariates, chooses regression method
- **Model Averaging**: Auto-selects models to average, uses BIC weighting

#### 3. **Comparison with Standard**
- Automatically runs standard NMA for comparison
- Reports agreement between experimental and standard

#### 4. **Risk Aversion Handling**
```r
risk_aversion = 0   # Risk-neutral (point estimates)
risk_aversion = 1   # Moderate caution (default)
risk_aversion = 2   # High caution (lower CI bounds)
```

### Example Usage

```r
library(powerNMA)

# Example 1: Time-to-event data (auto-selects RMST + Threshold)
survival_data <- data.frame(
  study = rep(1:10, each = 100),
  treatment = sample(c("A", "B", "C"), 1000, replace = TRUE),
  time = rexp(1000, 0.05),
  event = rbinom(1000, 1, 0.7)
)

result <- auto_experimental_nma(data = survival_data)

# Automatically:
# 1. Detects: Time-to-event IPD
# 2. Selects: RMST-based NMA + Threshold Analysis
# 3. Auto-selects tau = 36 months (based on follow-up)
# 4. Runs both methods
# 5. Generates: Survival time differences, robustness scores
# 6. Compares with standard HR-based NMA

# Example 2: IPD for precision medicine (auto-selects ITR)
ipd_data <- data.frame(
  study = rep(1:8, each = 150),
  treatment = sample(c("CBT", "Meds", "Combined"), 1200, replace = TRUE),
  outcome = rnorm(1200),
  age = rnorm(1200, 45, 15),
  severity = rnorm(1200, 25, 8)
)

result_prec <- auto_experimental_nma(
  data = ipd_data,
  outcome_var = "outcome",
  research_question = "precision_medicine"
)

# Automatically:
# 1. Detects: IPD with 2 covariates
# 2. Selects: ITR + Model Averaging
# 3. Identifies effect modifiers (age, severity)
# 4. Derives personalized treatment rules
# 5. Provides prediction function for new patients

# Example 3: Pairwise for guideline (auto-selects Threshold)
pairwise_data <- data.frame(
  study = rep(1:12, each = 2),
  treat1 = sample(c("A", "B", "C"), 24, replace = TRUE),
  treat2 = sample(c("A", "B", "C"), 24, replace = TRUE),
  TE = rnorm(24, 0.5, 0.3),
  seTE = runif(24, 0.1, 0.3)
)

result_guide <- auto_experimental_nma(
  data = pairwise_data,
  research_question = "decision_making",
  risk_aversion = 1.5  # More cautious
)

# Automatically:
# 1. Detects: Pairwise format
# 2. Selects: Threshold Analysis + Model Averaging
# 3. Calculates robustness score
# 4. Identifies tipping points
# 5. Risk-adjusted recommendations
```

### Output Structure

```r
result$data_characteristics    # What was detected
result$automatic_choices       # Experimental methods chosen
result$experimental_analyses   # All experimental results
result$threshold_analysis      # Robustness assessment (if run)
result$itr_analysis            # ITR results (if run)
result$model_averaging         # Model averaging (if run)
result$rmst_analysis           # RMST results (if run)
result$standard_comparison     # Comparison with standard NMA
result$diagnostics             # Experimental diagnostics
result$report                  # Auto-generated report
result$recommendations         # Clinical interpretation
```

### Advantages
- **Cutting-edge automation**: Advanced methods without deep expertise
- **Intelligent selection**: Chooses optimal experimental methods
- **Novel insights**: Effect modifiers, robustness, personalization
- **Transparent**: Reports all automatic choices
- **Comparison**: Always compares with standard methods

### Disadvantages
- **EXPERIMENTAL**: Limited validation
- **Not for routine use**: Requires expert interpretation
- **Regulatory limitations**: Not approved for FDA/EMA submissions
- **Complexity**: More complex interpretation than standard

---

## Comparison Matrix: All Four Pathways

| Feature | Manual Standard | Manual Experimental | Auto Standard | Auto Experimental |
|---------|----------------|---------------------|---------------|-------------------|
| **User Control** | Full | Full | Minimal | Minimal |
| **Automation** | None | None | Complete | Complete |
| **Methods** | Phase 1-2 (7) | Experimental (4) | Phase 1-2 (7) | Experimental (4) |
| **Literature** | 2019-2023 | 2024-2025 | 2019-2023 | 2024-2025 |
| **Validation** | Extensive | Limited | Extensive | Limited |
| **Expertise Required** | High | Expert | Minimal | Moderate |
| **Time Required** | Long | Long | Short | Short |
| **Standardization** | Low | Low | High | High |
| **Flexibility** | Maximum | Maximum | Low | Low |
| **Reproducibility** | Medium | Medium | High | High |
| **Clinical Guidelines** | ✅ Yes | ⚠️ Expert review | ✅ Yes | ❌ No |
| **Regulatory Use** | ✅ Yes | ❌ No | ✅ Yes | ❌ No |
| **Research** | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes |
| **Precision Medicine** | ⚠️ Manual ITR | ✅ Yes | ⚠️ If IPD | ✅ Auto ITR |
| **Robustness Assessment** | GRADE | Threshold | GRADE | Auto Threshold |
| **Novel Insights** | Moderate | High | Moderate | High |
| **Learning Curve** | Steep | Very steep | Gentle | Moderate |

---

## Decision Tree: Which Pathway Should I Use?

```
START: What is your goal?
│
├─ Need regulatory approval (FDA/EMA)?
│  └─ YES → Manual Standard (Pathway 1)
│
├─ Developing official clinical guideline?
│  ├─ Standard approach → Auto Standard (Pathway 3)
│  └─ Need robustness assessment → Manual Standard + Threshold Analysis
│
├─ Advanced research question?
│  ├─ Know exact method needed → Manual Experimental (Pathway 2)
│  └─ Want automatic selection → Auto Experimental (Pathway 4)
│
├─ Precision medicine / personalized treatment?
│  ├─ Have IPD with covariates?
│     ├─ Yes, want control → Manual Experimental: ITR (Pathway 2)
│     └─ Yes, want automation → Auto Experimental (Pathway 4)
│  └─ No IPD → Not feasible
│
├─ Time-to-event with non-proportional hazards?
│  ├─ Want control over tau → Manual Experimental: RMST (Pathway 2)
│  └─ Want automation → Auto Experimental (Pathway 4)
│
├─ Want standardized, reproducible analysis?
│  ├─ Standard methods sufficient → Auto Standard (Pathway 3)
│  └─ Need cutting-edge → Auto Experimental (Pathway 4)
│
├─ Limited methodological expertise?
│  ├─ Standard needs → Auto Standard (Pathway 3)
│  └─ Advanced needs + expert available → Auto Experimental (Pathway 4)
│
└─ Learning/teaching NMA?
   └─ Manual Standard (Pathway 1) - Educational
```

---

## Workflow Examples by Use Case

### Use Case 1: Regulatory Submission (FDA/EMA)

**Choose**: Manual Standard (Pathway 1)

```r
# Full control, validated methods, transparent decisions
result <- cnma(
  data = trial_data,
  components = component_specification,
  model = "additive",         # Explicit choice for regulators
  reference = "Placebo",      # Standard comparator
  tau_common = TRUE,          # Justified in protocol
  sm = "OR"
)

# Document all decisions in submission
```

**Why**: Regulators require validated methods with full transparency.

---

### Use Case 2: Clinical Guideline (Standardized)

**Choose**: Auto Standard (Pathway 3)

```r
# Standardized approach for consistency
result <- auto_standard_nma(
  data = guideline_data,
  verbose = TRUE  # Document automatic choices
)

# All choices made using best practices
# Reproducible across guideline developers
```

**Why**: Ensures consistency across different guideline committees.

---

### Use Case 3: Precision Medicine Research

**Choose**: Auto Experimental (Pathway 4)

```r
# Automatic ITR derivation
result <- auto_experimental_nma(
  data = ipd_with_biomarkers,
  outcome_var = "response",
  research_question = "precision_medicine"
)

# Automatically:
# - Identifies effect modifiers
# - Derives treatment rules
# - Provides prediction function

# Predict for new patient
new_patient <- data.frame(age = 65, biomarker = 8.2)
optimal <- predict(result, new_patient)
```

**Why**: Automatic identification of personalization opportunities.

---

### Use Case 4: Methodological Research

**Choose**: Manual Experimental (Pathway 2)

```r
# Compare multiple experimental methods
rmst_result <- rmst_nma(data, tau = 36)
threshold_result <- threshold_analysis(nma_obj)
itr_result <- itr_from_nma(data, method = "machine_learning")
bma_result <- model_averaging_nma(data, weighting = "DIC")

# Full control for methodological comparison
```

**Why**: Need precise control over experimental methods.

---

### Use Case 5: Rapid Analysis (Time-Limited)

**Choose**: Auto Standard or Auto Experimental (Pathways 3 or 4)

```r
# Standard methods
standard_result <- auto_standard_nma(data)

# Or experimental for novel insights
experimental_result <- auto_experimental_nma(data)

# Both complete in minimal time with automatic decisions
```

**Why**: Maximum efficiency with minimal user input.

---

## Migration Between Pathways

Users can **start with one pathway and migrate** to another:

### Common Migration Patterns

#### Pattern 1: Exploration → Production
```r
# Start: Auto Standard (explore)
explore <- auto_standard_nma(data)

# Review automatic choices
print(explore$automatic_choices)

# Migrate: Manual Standard (production with modifications)
production <- cnma(
  data = data,
  model = "interaction",  # Modified from auto choice
  # ... other parameters from explore$automatic_choices
)
```

#### Pattern 2: Standard → Experimental Enhancement
```r
# Start: Manual/Auto Standard
standard <- auto_standard_nma(data)

# Enhance: Add experimental robustness assessment
threshold <- threshold_analysis(
  nma_object = standard$primary_analysis,
  risk_aversion = 1.0
)

# Combined interpretation
```

#### Pattern 3: Automatic → Manual Refinement
```r
# Start: Auto Experimental
auto_exp <- auto_experimental_nma(data, research_question = "precision_medicine")

# Review: Check effect modifiers
print(auto_exp$itr_analysis$effect_modifiers)

# Refine: Manual ITR with specific choices
manual_itr <- itr_from_nma(
  data = data,
  covariate_vars = c("age", "severity"),  # Refined based on auto results
  method = "regression"  # Changed from auto choice
)
```

---

## Best Practices by Pathway

### Manual Standard - Best Practices
1. **Pre-specify choices**: Document in protocol
2. **Justify decisions**: Provide rationale for each choice
3. **Sensitivity analysis**: Test key assumptions
4. **Compare with automatic**: Check against `auto_standard_nma()` choices
5. **Report transparently**: Document all decisions

### Manual Experimental - Best Practices
1. **Know the methods**: Understand experimental method assumptions
2. **Compare with standard**: Always run standard NMA for comparison
3. **Expert review**: Have statistician review experimental results
4. **Acknowledge limitations**: Note experimental status
5. **Avoid routine use**: Only when standard methods insufficient

### Auto Standard - Best Practices
1. **Review automatic choices**: Check `$automatic_choices`
2. **Verify reasonableness**: Ensure auto choices make clinical sense
3. **Document in methods**: Report automatic choices in publication
4. **Consider sensitivity**: If concerned, run manual with alternatives
5. **Use for standardization**: Ideal for consistent workflows

### Auto Experimental - Best Practices
1. **Understand selection logic**: Review which methods were selected
2. **Compare with standard**: Check `$standard_comparison`
3. **Expert interpretation**: Have statistician interpret results
4. **Not for routine guidelines**: Use for research, not clinical guidelines
5. **Acknowledge experimental**: Clearly note experimental methods used

---

## Summary Table: Quick Reference

| If you need... | Use this pathway... | Function |
|---------------|-------------------|----------|
| FDA/EMA submission | Manual Standard | `cnma()`, `network_metareg()`, etc. |
| Official clinical guideline | Auto Standard | `auto_standard_nma()` |
| Standardized workflow | Auto Standard | `auto_standard_nma()` |
| Precision medicine | Auto Experimental | `auto_experimental_nma(research_question="precision_medicine")` |
| Time-to-event (non-proportional) | Auto/Manual Experimental | `auto_experimental_nma()` or `rmst_nma()` |
| Robustness assessment | Auto/Manual Experimental | `threshold_analysis()` |
| Model uncertainty | Manual Experimental | `model_averaging_nma()` |
| Learning NMA | Manual Standard | Individual functions |
| Fast analysis | Auto Standard | `auto_standard_nma()` |
| Methodological research | Manual Experimental | Individual experimental functions |
| Maximum control | Manual Standard | Individual standard functions |
| Cutting-edge automation | Auto Experimental | `auto_experimental_nma()` |

---

## Technical Implementation Details

### Code Organization

```
powerNMA/R/
├── Standard Methods (Manual)
│   ├── cnma.R                          (Component NMA)
│   ├── network_metareg.R               (Meta-regression)
│   ├── dose_response.R                 (Dose-response)
│   ├── prediction_methods.R            (Predictive rankings)
│   ├── multivariate_nma.R              (Multiple outcomes)
│   ├── missing_data.R                  (Missing data)
│   └── cross_design_synthesis.R        (RCT + observational)
│
├── Experimental Methods (Manual)
│   ├── experimental_rmst_nma.R         (RMST-based NMA)
│   ├── experimental_threshold_analysis.R (Threshold analysis)
│   ├── experimental_itr_nma.R          (Individualized treatment rules)
│   └── experimental_model_averaging.R  (Bayesian model averaging)
│
├── Automatic Pathways
│   ├── auto_standard_pathway.R         (Auto Standard)
│   └── auto_experimental_pathway.R     (Auto Experimental)
│
└── NAMESPACE                           (All exports)
```

### NAMESPACE Exports

```r
# Manual Standard (7 main functions + S3 methods)
export(cnma, network_metareg, dose_response_nma, predictive_ranking,
       treatment_hierarchy, multivariate_nma, handle_missing_data,
       recommend_missing_data_strategy, cross_design_synthesis)

# Manual Experimental (4 functions + S3 methods)
export(rmst_nma, threshold_analysis, itr_from_nma, model_averaging_nma)

# Automatic Pathways (2 functions + S3 methods)
export(auto_standard_nma, auto_experimental_nma)
```

---

## Version History

**Version 3.0** (October 31, 2025)
- ✅ Implemented all 4 pathways
- ✅ 7 standard methods (Phase 1-2)
- ✅ 4 experimental methods (2024-2025)
- ✅ 2 automatic pathways
- ✅ Complete documentation

**Statistics**:
- Total methods: 11 (7 standard + 4 experimental)
- Total pathways: 4 (2 manual + 2 automatic)
- Lines of code: 8,100+
- Documentation pages: 140+
- Exports: 43 (functions + S3 methods)

---

## Future Enhancements

Potential additions to the pathway system:

1. **Interactive Pathway Selector**
   ```r
   select_pathway()  # Interactive wizard
   ```

2. **Pathway Comparison Tool**
   ```r
   compare_pathways(data, pathways = c("auto_standard", "auto_experimental"))
   ```

3. **Custom Automatic Pathway**
   ```r
   create_custom_auto_pathway(choices = my_preferences)
   ```

4. **Pathway Validation**
   ```r
   validate_pathway_choice(data, pathway = "auto_experimental")
   ```

---

## Contact and Support

- **Questions**: Open GitHub issue with tag `[PATHWAY]`
- **Pathway selection help**: Tag `[WHICH_PATHWAY]`
- **Bug reports**: Specify pathway used
- **Feature requests**: Suggest pathway enhancements

---

**Document Version**: 1.0
**Last Updated**: October 31, 2025
**Next Review**: January 2026
