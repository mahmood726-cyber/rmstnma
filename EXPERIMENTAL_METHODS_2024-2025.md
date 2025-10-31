# powerNMA: Experimental Methods Pathway (2024-2025)

**Status**: EXPERIMENTAL - Cutting-edge methods from 2024-2025 literature
**Documentation Date**: October 31, 2025
**Package Version**: 3.0 (Experimental Branch)

---

## Overview

This document describes **4 experimental network meta-analysis methods** added to powerNMA based on the most recent peer-reviewed literature (2024-2025). These methods represent the cutting edge of NMA methodology and are marked as **EXPERIMENTAL** because:

1. **Recency**: Published in 2024-2025, limited real-world validation
2. **Complexity**: Advanced statistical methods requiring careful interpretation
3. **Dependencies**: Some methods require additional packages or assumptions
4. **Active Development**: Methods still being refined in the literature

### Experimental vs Standard Methods

| Feature | Standard Methods (Phase 1-2) | Experimental Methods |
|---------|------------------------------|---------------------|
| **Literature** | Established (2019-2023) | Cutting-edge (2024-2025) |
| **Validation** | Extensively validated | Limited validation |
| **Stability** | Production-ready | Under active development |
| **Complexity** | Moderate | High |
| **Use Case** | Routine analyses | Advanced research questions |
| **Risk Level** | Low | Moderate to high |

### When to Use Experimental Methods

✅ **Use when:**
- Researching novel methodological questions
- Need advanced capabilities not in standard methods
- Have expertise to interpret complex results
- Publishing methodological research
- Collaborating with statisticians

❌ **Avoid when:**
- Conducting routine clinical guidelines
- Limited statistical expertise available
- Requiring regulatory approval (e.g., FDA submissions)
- Need maximum reproducibility
- Time-sensitive decisions

---

## Experimental Method 1: RMST-based Network Meta-Analysis

### 📚 Reference
Hua H, et al. (2025). "Network Meta-Analysis of Time-to-Event Endpoints With Individual Participant Data Using Restricted Mean Survival Time Regression." *Biometrical Journal*, 67(1), e70037.

### 🎯 Purpose
Uses **Restricted Mean Survival Time (RMST)** instead of hazard ratios for time-to-event outcomes in network meta-analysis.

### 💡 Key Innovation
**Problem**: Hazard ratios (HRs) require:
- Proportional hazards assumption (often violated)
- Complex interpretation
- Non-collapsibility issues

**Solution**: RMST provides:
- Area under survival curve up to time τ (tau)
- Direct interpretation: "Treatment A gives X additional months of survival"
- No proportional hazards assumption needed
- More clinically meaningful

### 📊 What is RMST?

```
RMST(τ) = ∫₀^τ S(t) dt
```

Where S(t) is the survival function. RMST represents the **mean survival time up to τ**.

**Example**:
- RMST(36 months) = 28.5 months → Patients survive average of 28.5 months in first 3 years
- RMST difference = 2.5 months → Treatment gives 2.5 additional months on average

### 🔧 Implementation

```r
# Function: rmst_nma()
# File: powerNMA/R/experimental_rmst_nma.R

# Example 1: Aggregate data
data_agg <- data.frame(
  study = rep(c("Study1", "Study2", "Study3"), each = 2),
  treatment = rep(c("A", "B"), 3),
  rmst = c(24.5, 26.8, 23.1, 25.9, 22.8, 27.2),  # months
  rmst_se = c(1.2, 1.3, 1.5, 1.4, 1.1, 1.2)
)

result <- rmst_nma(
  data = data_agg,
  tau = 36,  # 3-year RMST
  data_type = "aggregate",
  reference = "A",
  method = "frequentist"
)

print(result)
plot(result, type = "forest")

# Example 2: Individual participant data
data_ipd <- data.frame(
  study = rep(c("Study1", "Study2"), each = 100),
  treatment = sample(c("A", "B", "C"), 200, replace = TRUE),
  time = rexp(200, rate = 0.05),
  event = rbinom(200, 1, 0.7)
)

result_ipd <- rmst_nma(
  data = data_ipd,
  tau = 24,  # 2-year RMST
  data_type = "ipd",
  reference = "A"
)
```

### ✅ Advantages Over Hazard Ratios

| Feature | Hazard Ratio | RMST |
|---------|-------------|------|
| **Interpretation** | Relative instantaneous risk | Mean survival time |
| **Proportional hazards** | Required | Not required |
| **Clinical meaning** | Abstract | Direct (months gained) |
| **Crossing curves** | Problematic | Handles naturally |
| **Patient understanding** | Difficult | Intuitive |

### ⚠️ Limitations

1. **Requires tau specification**: Must choose restriction time
2. **IPD preferred**: Works best with individual participant data
3. **Censoring**: Heavy censoring before tau reduces precision
4. **Package dependencies**: Requires `survival` package

### 📖 Output Interpretation

```r
# Output example:
Treatment Rankings (by RMST at tau=36 months):
  Rank  Treatment  RMST (months)  95% CI
  1     C          31.2          (29.1, 33.3)
  2     B          28.7          (26.9, 30.5)
  3     A          26.5          (24.8, 28.2)

Interpretation: Treatment C provides on average 4.7 additional months
of survival (compared to A) within the first 3 years.
```

### 🎯 Clinical Use Cases

1. **Cancer trials**: Overall survival where proportional hazards violated
2. **Cardiovascular outcomes**: Time to MACE with crossing survival curves
3. **Chronic disease**: Long-term follow-up with non-proportional effects
4. **Regulatory submissions**: When HR interpretation is challenging

---

## Experimental Method 2: Threshold Analysis for Decision-Making

### 📚 Reference
Ades AE, et al. (2025). "Treatment recommendations based on network meta-analysis: Rules for risk-averse decision-makers." *Research Synthesis Methods*.

Phillippo DM, et al. (2019). "Threshold analysis as an alternative to GRADE." *Annals of Internal Medicine*, 170:538-546.

### 🎯 Purpose
Quantifies **how much evidence would need to change** before a treatment recommendation changes, providing a quantitative measure of recommendation robustness.

### 💡 Key Innovation

**Problem**: Traditional certainty assessments (GRADE, CINeMA) provide:
- Qualitative ratings ("low", "moderate", "high" certainty)
- No specific guidance on how much evidence needed
- Difficult to interpret for decision-makers

**Solution**: Threshold analysis provides:
- **Quantitative thresholds**: "Effect would need to change by X"
- **Tipping point**: Exact value where recommendation changes
- **Robustness score**: 0-100 scale of recommendation stability
- **Actionable guidance**: What new study would change recommendation

### 📊 Key Question

> "How much would the evidence have to change (due to bias, new data, or uncertainty) before the treatment recommendation changes?"

### 🔧 Implementation

```r
# Function: threshold_analysis()
# File: powerNMA/R/experimental_threshold_analysis.R

# Example: Threshold analysis for smoking cessation NMA
library(netmeta)
data(smokingcessation)

nma <- netmeta(TE, seTE, treat1, treat2, studlab,
               data = smokingcessation, sm = "OR")

# Perform threshold analysis
threshold_result <- threshold_analysis(
  nma_object = nma,
  outcome_direction = "higher",  # Higher OR = better
  decision_rule = "maximize_benefit",
  risk_aversion = 1.0,  # 0 = risk-neutral, higher = more cautious
  threshold_type = "all"
)

print(threshold_result)
plot(threshold_result, type = "robustness")
plot(threshold_result, type = "thresholds")
```

### ✅ Three Types of Thresholds

#### 1. **Effect Size Threshold**
"How much would the treatment effect need to change?"

```r
# Example output:
Effect Size Threshold: 0.35 (1.8 standard deviations)
Interpretation: Treatment A's effect would need to decrease by 0.35
(or Treatment B's effect increase by 0.35) before recommendation
changes to Treatment B.
```

#### 2. **Bias Threshold**
"How much systematic bias would change the recommendation?"

```r
# Example output:
Bias Threshold: 0.35 (1.8 SDs)
Robust to bias: YES (threshold > 2 SDs indicates robustness)
Interpretation: A systematic bias of 0.35 would change recommendation.
This is ROBUST to plausible bias (typically < 0.5 SD).
```

#### 3. **New Study Threshold**
"What would a new study need to show?"

```r
# Example output:
New Study Threshold: Required effect = -0.70 (2.3 SDs)
Interpretation: A new study would need to show effect of -0.70
(opposite direction) to change recommendation. This is unlikely,
indicating robust recommendation.
```

### 📊 Robustness Score (0-100)

```r
Robustness Score: 85/100
Category: VERY ROBUST
```

| Score | Category | Interpretation |
|-------|----------|----------------|
| **80-100** | Very Robust | Recommendation extremely stable, unlikely to change |
| **50-79** | Moderately Robust | Reasonably stable, substantial new evidence needed |
| **0-49** | Fragile | Sensitive to uncertainty, more evidence strongly recommended |

### 🎯 Decision Rules

```r
decision_rule = "maximize_benefit"       # Choose best treatment
decision_rule = "minimize_risk"          # Risk-averse choice
decision_rule = "cost_effectiveness"     # Consider costs
decision_rule = "multi_criteria"         # Multiple outcomes
```

**Risk aversion parameter**:
```r
risk_aversion = 0    # Risk-neutral: Use point estimates
risk_aversion = 1    # Moderately risk-averse: Penalize uncertainty
risk_aversion = 2    # Highly risk-averse: Use lower confidence bounds
```

### 📖 Output Interpretation

```r
=============================================================
EXPERIMENTAL: Threshold Analysis for Treatment Recommendation
=============================================================

Current Recommendation: Treatment C
Robustness Score: 85 / 100
Category: VERY ROBUST

Thresholds for Recommendation Change:
--------------------------------------

Effect Size Threshold:
Recommended treatment effect would need to change by 0.42
(2.1 SDs) before recommendation changes to Treatment B

Bias Threshold:
A systematic bias of 0.42 (2.1 SDs) would change recommendation.
This is ROBUST to plausible bias.

New Study Threshold:
A new study would need to show effect of -0.84 (assuming SE = 0.20)
to change recommendation. This represents 4.2 SDs.

Alternative Recommendation: Treatment B

Clinical Guidance:
Recommendation is highly stable. Unlikely to change with new
evidence or bias.
```

### 🎯 Clinical Use Cases

1. **Clinical guidelines**: Assess robustness before publishing recommendations
2. **Policy decisions**: Understand stability of recommendations for reimbursement
3. **Research prioritization**: Identify when more research is truly needed
4. **Regulatory submissions**: Demonstrate robustness to bias/uncertainty
5. **Patient communication**: Explain confidence in recommendations

### ⚠️ Limitations

1. **Single outcome focus**: Current implementation focuses on primary outcome
2. **Simplified cost-effectiveness**: Full health-economic modeling needed for costs
3. **Assumes consistency**: Threshold analysis assumes network consistency holds

---

## Experimental Method 3: Individualized Treatment Rules from NMA

### 📚 Reference
Chen et al. (2025). "Developing a multivariable prediction model to support personalized selection among treatments." *PLOS ONE*.

Zhao et al. (2023). "Meta-analysis of individualized treatment rules via sign-coherency." *JASA*.

### 🎯 Purpose
Derives **personalized treatment recommendations** from network meta-analysis using patient characteristics. Identifies which treatments work best for which patients.

### 💡 Key Innovation

**Problem**: Standard NMA provides:
- Average treatment effects across all patients
- Ignores effect modification (treatment × covariate interactions)
- Cannot personalize recommendations

**Solution**: ITR from NMA provides:
- Patient-specific optimal treatment
- Identifies effect modifiers (which characteristics predict best treatment)
- Enables precision medicine from network evidence

### 📊 What are ITRs?

**Individualized Treatment Rule (ITR)**: A function that maps patient characteristics to optimal treatment.

```
ITR: f(patient covariates) → Optimal Treatment

Example:
- If age > 65 and severity > 7 → Treatment C
- If age ≤ 65 and severity ≤ 7 → Treatment A
- Otherwise → Treatment B
```

### 🔧 Implementation

```r
# Function: itr_from_nma()
# File: powerNMA/R/experimental_itr_nma.R

# Example: Derive ITR from depression treatment NMA with IPD
data_ipd <- data.frame(
  study = rep(paste0("Study", 1:5), each = 200),
  treatment = sample(c("CBT", "Medication", "Combined"), 1000, replace = TRUE),
  age = rnorm(1000, mean = 45, sd = 15),
  sex = sample(c("Male", "Female"), 1000, replace = TRUE),
  severity = rnorm(1000, mean = 25, sd = 8),
  outcome = rnorm(1000)  # Depression score improvement
)

# Derive ITR
itr_result <- itr_from_nma(
  data = data_ipd,
  outcome_var = "outcome",
  treatment_var = "treatment",
  covariate_vars = c("age", "sex", "severity"),
  outcome_type = "continuous",
  method = "regression",  # or "machine_learning"
  reference_treatment = "CBT"
)

print(itr_result)
plot(itr_result, type = "effect_modifiers")

# Predict optimal treatment for new patient
new_patient <- data.frame(
  age = 65,
  sex = "Female",
  severity = 30
)

optimal <- predict(itr_result,
                  newdata = new_patient,
                  treatments = c("CBT", "Medication", "Combined"))

print(optimal)
# Output:
# $optimal_treatment
# [1] "Combined"
# $predicted_outcome
# [1] 15.2
```

### ✅ Four ITR Methods

#### 1. **Regression-based ITR**
```r
method = "regression"
```
- Fits linear/logistic model with treatment × covariate interactions
- Interpretable coefficients
- Fast and stable
- **Use for**: Continuous/binary outcomes, limited covariates

#### 2. **Contrast Regression ITR**
```r
method = "contrast_regression"
```
- Models treatment contrasts with covariate interactions
- More robust to multi-arm trials
- Better for network structure
- **Use for**: Complex network structures

#### 3. **Machine Learning ITR**
```r
method = "machine_learning"
ml_algorithm = "random_forest"  # or "gradient_boosting", "lasso"
```
- Uses ML to learn optimal treatment assignment
- Handles non-linear effects
- Automatic interaction detection
- **Use for**: Many covariates, complex interactions

#### 4. **Sign-Coherency Meta-Analysis**
```r
method = "sign_coherency"
```
- Meta-analyzes ITRs across studies
- Assumes sign coherency (e.g., if age predicts better response in Study 1, same in Study 2)
- Borrows strength across studies
- **Use for**: Multiple studies, consistent effect modification

### 📊 Effect Modifiers

The method identifies **which patient characteristics modify treatment effects**:

```r
# Output:
Effect Modifiers:
  covariate  strength  p_value  direction   significant
  severity   0.45      0.001    positive    TRUE
  age        0.32      0.012    negative    TRUE
  sex        0.08      0.456    positive    FALSE

Interpretation:
- Severity is the strongest effect modifier (strength = 0.45)
- Patients with higher severity benefit more from Combined treatment
- Age negatively modifies effect (younger patients respond better)
- Sex does not significantly modify treatment effects
```

### 📖 Treatment Rules

```r
Treatment Rules:
1. For patients with positive severity, treatment effects differ
   significantly. Consider severity when selecting treatment.

2. For patients with negative age, treatment effects differ
   significantly. Consider age when selecting treatment.

Clinical Use:
Use the prediction function to determine optimal treatment for new
patients based on their characteristics. The model considers
interactions between patient covariates and treatment effects.
```

### 🎯 Clinical Example: Depression Treatment

**Scenario**: Meta-analysis of 10 RCTs comparing CBT, Medication, and Combined therapy for depression.

**Effect Modifiers Identified**:
- **Baseline severity**: High severity patients benefit most from Combined
- **Age**: Younger patients respond better to CBT
- **Prior episodes**: Multiple episodes predict better Medication response

**ITR Recommendations**:
```
Patient 1: Age=25, Severity=18, Episodes=1 → Optimal: CBT
Patient 2: Age=55, Severity=32, Episodes=3 → Optimal: Combined
Patient 3: Age=45, Severity=28, Episodes=2 → Optimal: Medication
```

### ⚠️ Limitations

1. **Requires IPD**: Individual participant data needed (aggregate data insufficient)
2. **Sample size**: Need sufficient sample for interaction detection
3. **Consistency assumption**: Assumes effect modification consistent across studies
4. **Overfitting risk**: ML methods may overfit, cross-validation essential
5. **External validity**: ITR derived from trial populations may not generalize

---

## Experimental Method 4: Bayesian Model Averaging for NMA

### 📚 Reference
Phillippo DM, et al. (2024). "Multilevel network meta-regression with model averaging." *Research Synthesis Methods* (in press).

Hoeting JA, et al. (1999). "Bayesian model averaging: a tutorial." *Statistical Science*, 14(4):382-417.

### 🎯 Purpose
Accounts for **structural uncertainty** in NMA by averaging estimates across multiple plausible model specifications, weighted by model fit or posterior probability.

### 💡 Key Innovation

**Problem**: Single-model inference ignores:
- Should we use fixed or random effects?
- Is there inconsistency in the network?
- What heterogeneity prior is appropriate?
- Which covariates to include?

**Solution**: Model averaging:
- Fits multiple plausible models
- Weights models by fit (AIC/BIC/DIC) or posterior probability
- Averages estimates across models
- Provides honest uncertainty intervals

### 📊 Model Averaging Formula

```
θ̂_BMA = Σ w_m × θ̂_m

where:
- θ̂_m = estimate from model m
- w_m = weight for model m
- Σ w_m = 1
```

**Weights** (BIC-based):
```
w_m = exp(-0.5 × ΔBIC_m) / Σ exp(-0.5 × ΔBIC)

where ΔBIC_m = BIC_m - BIC_best
```

### 🔧 Implementation

```r
# Function: model_averaging_nma()
# File: powerNMA/R/experimental_model_averaging.R

# Example: Model averaging across fixed/random effects and inconsistency models
library(netmeta)
data(smokingcessation)

data_pairs <- data.frame(
  study = smokingcessation$studlab,
  treat1 = smokingcessation$treat1,
  treat2 = smokingcessation$treat2,
  TE = smokingcessation$TE,
  seTE = smokingcessation$seTE
)

# Bayesian model averaging
bma_result <- model_averaging_nma(
  data = data_pairs,
  models = c("fixed_effect", "random_effects", "inconsistency_re"),
  weighting = "BIC",  # or "AIC", "DIC", "posterior", "equal"
  reference = "A",
  sm = "OR"
)

print(bma_result)
summary(bma_result)
plot(bma_result, type = "weights")
plot(bma_result, type = "estimates")
```

### ✅ Model Types Available

```r
models = c(
  "fixed_effect",       # Fixed-effect NMA (no heterogeneity)
  "random_effects",     # Random-effects NMA (standard)
  "inconsistency_re",   # Random-effects allowing inconsistency
  "consistency_re",     # Random-effects enforcing consistency
  "uisd"               # Unrelated study effects
)
```

### 📊 Weighting Methods

#### 1. **BIC Weighting** (Bayesian Information Criterion)
```r
weighting = "BIC"
```
- Penalizes model complexity more than AIC
- Preferred for model selection
- Formula: BIC = -2×logLik + k×log(n)

#### 2. **AIC Weighting** (Akaike Information Criterion)
```r
weighting = "AIC"
```
- Less penalty for complexity
- Better for prediction
- Formula: AIC = -2×logLik + 2k

#### 3. **DIC Weighting** (Deviance Information Criterion)
```r
weighting = "DIC"
```
- Bayesian equivalent of AIC
- Requires JAGS/Stan
- Accounts for effective number of parameters

#### 4. **Posterior Probability**
```r
weighting = "posterior"
```
- Full Bayesian model averaging
- Requires JAGS/Stan
- Most theoretically sound

#### 5. **Equal Weights**
```r
weighting = "equal"
```
- Simple averaging across models
- No preference for any model

### 📖 Output Interpretation

```r
=============================================================
EXPERIMENTAL: Bayesian Model Averaging for NMA
=============================================================

Weighting method: BIC
Number of models: 3
Best single model: random_effects

Model Weights:
  fixed_effect        0.123
  random_effects      0.658
  inconsistency_re    0.219

Model Comparison:
  model                AIC      BIC      logLik   n_params
  random_effects      245.3    256.8    -118.7   4
  inconsistency_re    248.1    262.4    -118.1   6
  fixed_effect        251.6    260.2    -122.8   3

Model-Averaged Treatment Effects:
  treatment   effect_avg   se_avg    95% CI
  A           0.00        0.00      -
  B           0.45        0.18      (0.10, 0.80)
  C           0.62        0.21      (0.21, 1.03)
  D           0.38        0.24      (-0.09, 0.85)

Interpretation:
Model random_effects has highest weight (65.8%), but other models
contribute. Model averaging is beneficial.

Model-averaged estimates account for structural uncertainty and
provide more honest confidence intervals.
```

### 📊 Model Dominance

**Strong Dominance** (weight > 90%):
- One model clearly best
- Model uncertainty small
- Consider using single best model

**Moderate Dominance** (weight 60-90%):
- One model preferred
- Other models contribute
- Model averaging recommended

**Weak Dominance** (weight < 60%):
- No clear best model
- Substantial model uncertainty
- Model averaging essential

### 🎯 Why Model-Averaged SEs are Larger

Model-averaged standard errors incorporate **two sources of uncertainty**:

1. **Within-model uncertainty**: Standard SE from each model
2. **Between-model uncertainty**: Variance across models

```
SE²_BMA = Σ w_m × (SE²_m + (θ̂_m - θ̂_BMA)²)
         \_____________/   \_______________/
         Within-model      Between-model
         uncertainty       uncertainty
```

**Result**: Model-averaged SEs are always ≥ single-model SEs (more honest!).

### 🎯 Clinical Use Cases

1. **Guideline development**: Acknowledge structural uncertainty in recommendations
2. **Health technology assessment**: Provide robust estimates for reimbursement decisions
3. **Sensitivity analysis**: Compare single best model to model-averaged estimates
4. **Regulatory submissions**: Demonstrate robustness to model specification
5. **Research synthesis**: When multiple models have similar fit

### ⚠️ Limitations

1. **Computational cost**: Fitting multiple models takes time
2. **Model specification**: Still requires specifying candidate model set
3. **Interpretation**: Model-averaged estimates harder to explain to non-statisticians
4. **Software limitations**: Full Bayesian BMA requires JAGS/Stan

---

## Comparison of All Experimental Methods

| Method | Primary Purpose | Data Required | Complexity | Clinical Impact |
|--------|----------------|---------------|------------|-----------------|
| **RMST-based NMA** | Time-to-event with meaningful units | IPD or aggregate RMST | Moderate | High (intuitive interpretation) |
| **Threshold Analysis** | Assess recommendation robustness | NMA results | Low | Very High (actionable guidance) |
| **Individualized Treatment Rules** | Personalized treatment selection | IPD with covariates | High | Very High (precision medicine) |
| **Bayesian Model Averaging** | Account for model uncertainty | Pairwise data | Moderate | Moderate (honest uncertainty) |

---

## Implementation Statistics

### Code Metrics

```
Total Experimental Methods: 4
Total Lines of Code: 2,400+
Total Functions: 4 main + 30+ helpers
Total S3 Methods: 12

Files:
- powerNMA/R/experimental_rmst_nma.R           (700+ lines)
- powerNMA/R/experimental_threshold_analysis.R (600+ lines)
- powerNMA/R/experimental_itr_nma.R           (650+ lines)
- powerNMA/R/experimental_model_averaging.R    (550+ lines)
```

### Dependencies

**Required Packages**:
- `netmeta` (for network meta-analysis)
- `survival` (for RMST calculations)

**Optional Packages** (for advanced features):
- `gemtc` / `rjags` (for Bayesian methods)
- `randomForest` / `gbm` (for ML-based ITR)
- `mvmeta` (for multivariate extensions)

---

## Usage Workflow

### 1. Start with Standard Methods (Phase 1-2)

```r
# Use established methods first
library(powerNMA)

# Standard NMA
nma_result <- network_metareg(data, ...)

# Standard prediction
pred_result <- predictive_ranking(nma_result, ...)
```

### 2. Consider Experimental Methods When:

```r
# Use RMST-based NMA if:
- Time-to-event outcome
- Proportional hazards violated
- Need intuitive interpretation
→ rmst_nma()

# Use Threshold Analysis if:
- Need to assess recommendation robustness
- Developing clinical guidelines
- Want quantitative sensitivity analysis
→ threshold_analysis()

# Use ITR if:
- Have individual participant data
- Suspect effect modification
- Want personalized recommendations
→ itr_from_nma()

# Use Model Averaging if:
- Uncertain about model specification
- Multiple models have similar fit
- Want honest uncertainty intervals
→ model_averaging_nma()
```

### 3. Report Both Standard and Experimental

```r
# Best practice: Report both
standard_result <- netmeta(...)
experimental_result <- rmst_nma(...)  # or other experimental method

# Compare and contrast
# Discuss alignment or differences
# Explain choice of final method
```

---

## Validation and Testing

### Status of Validation

| Method | Simulation Testing | Real Data Testing | Peer Review | Status |
|--------|-------------------|-------------------|-------------|---------|
| RMST-based NMA | ✅ Published | ✅ Published | ✅ Published 2025 | **Validated** |
| Threshold Analysis | ✅ Published | ✅ Published | ✅ Published 2025 | **Validated** |
| ITR from NMA | ✅ Published | ⚠️ Limited | ✅ Published 2025 | **Partially Validated** |
| Model Averaging | ✅ Published | ⚠️ Limited | ⏳ In Press 2024 | **Under Review** |

### Recommended Testing

Before using experimental methods in production:

1. **Compare with standard methods**: Do results align?
2. **Sensitivity analysis**: Test with different parameters
3. **Simulation validation**: If possible, test on simulated data
4. **Peer review**: Have statistician review results
5. **Cross-validation**: Use CV to assess performance

---

## Future Directions

### Methods Under Development (Not Yet Implemented)

1. **Living Network Meta-Analysis** (2025+)
   - Real-time updating with new studies
   - Automated evidence surveillance

2. **Network Meta-Analysis for Diagnostic Tests** (2024-2025)
   - Sensitivity/specificity networks
   - Test accuracy comparisons

3. **Causal Inference in NMA** (2025+)
   - Addressing confounding in observational networks
   - G-computation for NMA

4. **Deep Learning for NMA** (2025+)
   - Neural networks for treatment effect prediction
   - Automated feature engineering

---

## References - Experimental Methods

### RMST-based NMA
- Hua H, et al. (2025). Network Meta-Analysis of Time-to-Event Endpoints With Individual Participant Data Using Restricted Mean Survival Time Regression. *Biometrical Journal*, 67(1), e70037.
- Royston P, Parmar MKB (2013). Restricted mean survival time: an alternative to the hazard ratio. *BMC Medical Research Methodology*, 13:152.

### Threshold Analysis
- Ades AE, et al. (2025). Treatment recommendations based on network meta-analysis: Rules for risk-averse decision-makers. *Research Synthesis Methods*.
- Phillippo DM, et al. (2019). Threshold analysis as an alternative to GRADE. *Annals of Internal Medicine*, 170:538-546.

### Individualized Treatment Rules
- Chen et al. (2025). Developing a multivariable prediction model to support personalized selection among treatments. *PLOS ONE*.
- Zhao Y, et al. (2023). Meta-analysis of individualized treatment rules via sign-coherency. *JASA*.

### Bayesian Model Averaging
- Phillippo DM, et al. (2024). Multilevel network meta-regression with model averaging. *Research Synthesis Methods* (in press).
- Hoeting JA, et al. (1999). Bayesian model averaging: a tutorial. *Statistical Science*, 14(4):382-417.

---

## Support and Feedback

**Questions?** Open an issue on GitHub with tag `[EXPERIMENTAL]`

**Bug Reports?** Experimental methods are under active development - please report issues

**Collaborations?** Contact for methodological collaborations on experimental methods

---

**Document Version**: 1.0
**Last Updated**: October 31, 2025
**Next Review**: January 2026
