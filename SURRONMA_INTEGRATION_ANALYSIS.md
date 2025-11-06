# SurroNMA Integration Analysis

**Date**: November 1, 2025
**Purpose**: Analyze potential integration of surrogate endpoint features
**Status**: Repository not accessible - providing integration roadmap

---

## Repository Search Results

**Searched for**: surronma, surrogate-related code
**Found**: No surronma repository in current environment
**Available**: Only rmstnma repository with powerNMA package

---

## Likely Features in SurroNMA (Based on Literature)

Surrogate endpoint network meta-analysis typically includes:

### 1. **Surrogate Validation Methods**

#### Individual-Level Surrogacy
- **R² individual**: Correlation between surrogate (S) and true endpoint (T)
- **Adjusted R²**: Accounting for treatment effect
- **Prentice criteria**: Statistical tests for surrogacy

#### Trial-Level Surrogacy
- **R² trial**: Correlation between treatment effects on S vs T across trials
- **Meta-regression approach**: Regress TE_T on TE_S
- **Prediction intervals**: For new trial's true endpoint from surrogate

### 2. **Multivariate Network Meta-Analysis for Surrogates**

- **Bivariate NMA**: Joint modeling of surrogate and true endpoint
- **Borrowing of strength**: Using surrogate data from trials without true endpoint
- **Copula methods**: Flexible correlation structures

### 3. **Surrogate Threshold Effect (STE)**

- Minimum surrogate effect needed to predict true benefit
- Based on prediction intervals
- Clinical decision-making tool

### 4. **Validation Frameworks**

- **Information Gain**: How much surrogate reduces uncertainty about true endpoint
- **Relative Effect**: RE = β_surrogate / β_true endpoint
- **Surrogacy patterns**: Treatment-specific vs class-specific surrogates

---

## Current powerNMA Features That Overlap

### ✅ Already Implemented in powerNMA:

1. **Multivariate NMA** (`multivariate_nma.R`)
   - Lines 1-260: Joint modeling of multiple outcomes
   - Can be used for surrogate + true endpoint
   - Missing: Explicit surrogate validation metrics

2. **Network Meta-Regression** (`network_metareg.R`)
   - Lines 1-250: Covariate-adjusted NMA
   - Could model S→T relationship
   - Missing: Trial-level R² calculation

3. **Correlation Handling**
   - Basic within-study correlation
   - Missing: Surrogate-specific correlation structures

4. **Prediction Methods** (`prediction_methods.R`)
   - Lines 1-200: Predictive distributions
   - Could predict true endpoint from surrogate
   - Missing: Surrogate-specific prediction intervals

### ❌ Missing in powerNMA:

1. **Surrogate Validation Metrics**
   - R² individual
   - R² trial
   - Prentice criteria tests
   - Information gain

2. **Surrogate Threshold Effect**
   - STE calculation
   - Prediction interval bands
   - Clinical interpretation

3. **Surrogacy Evaluation Plots**
   - S vs T scatter (individual level)
   - TE_S vs TE_T scatter (trial level)
   - Surrogate ROC curves
   - Validation plots

4. **Surrogate-Specific Documentation**
   - When to use surrogates
   - Validation requirements
   - Reporting standards

---

## Recommended Integration Plan

### Phase 1: Core Surrogate Functions (HIGH Priority)

**File to Create**: `powerNMA/R/surrogate_nma.R`

```r
#' Surrogate Endpoint Network Meta-Analysis
#'
#' Validates surrogate endpoints and uses them in NMA
#'
#' @param data_surrogate Surrogate endpoint data
#' @param data_true True endpoint data (may be subset)
#' @param surrogate_var Column name for surrogate
#' @param true_var Column name for true endpoint
#' @export
surrogate_nma <- function(data_surrogate, data_true,
                         surrogate_var, true_var, ...) {

  # Step 1: Validate surrogacy
  validation <- validate_surrogate(
    data_surrogate, data_true,
    surrogate_var, true_var
  )

  # Step 2: Bivariate NMA (joint model)
  bivariate_nma <- multivariate_nma(
    data = merge_surrogate_true(data_surrogate, data_true),
    outcomes = c(surrogate_var, true_var),
    method = "bivariate"
  )

  # Step 3: Calculate STE
  ste <- calculate_ste(validation, bivariate_nma)

  # Step 4: Predict true effects from surrogate
  predictions <- predict_true_from_surrogate(
    surrogate_effects = bivariate_nma$effects[[surrogate_var]],
    surrogacy_model = validation
  )

  structure(list(
    validation = validation,
    bivariate_nma = bivariate_nma,
    ste = ste,
    predictions = predictions,
    call = match.call()
  ), class = "surrogate_nma")
}


#' Validate Surrogate Endpoint
#' @keywords internal
validate_surrogate <- function(data_s, data_t, s_var, t_var) {

  # Individual-level validation (if IPD available)
  if (has_ipd(data_s) && has_ipd(data_t)) {
    r2_individual <- calculate_r2_individual(data_s, data_t, s_var, t_var)
  } else {
    r2_individual <- NA
  }

  # Trial-level validation
  trial_level <- calculate_trial_level_surrogacy(data_s, data_t, s_var, t_var)

  # Prentice criteria
  prentice <- test_prentice_criteria(data_s, data_t, s_var, t_var)

  # Information gain
  info_gain <- calculate_information_gain(data_s, data_t, s_var, t_var)

  list(
    r2_individual = r2_individual,
    r2_trial = trial_level$r2_trial,
    regression_coef = trial_level$beta,
    regression_se = trial_level$se,
    prentice = prentice,
    information_gain = info_gain,
    interpretation = interpret_surrogacy(r2_individual, trial_level$r2_trial)
  )
}


#' Calculate Trial-Level R²
#' @keywords internal
calculate_trial_level_surrogacy <- function(data_s, data_t, s_var, t_var) {

  # Get treatment effects for surrogate
  te_surrogate <- extract_treatment_effects(data_s, s_var)

  # Get treatment effects for true endpoint
  te_true <- extract_treatment_effects(data_t, t_var)

  # Merge by study and comparison
  merged <- merge(te_surrogate, te_true,
                  by = c("study", "comparison"))

  # Meta-regression: TE_true ~ TE_surrogate
  # Weight by inverse variance
  weights <- 1 / (merged$se_true^2 + merged$se_surrogate^2)

  fit <- lm(te_true ~ te_surrogate,
            data = merged,
            weights = weights)

  # R² trial
  r2_trial <- summary(fit)$r.squared

  # Prediction interval for new trial
  pred_interval <- predict(fit, interval = "prediction", level = 0.95)

  list(
    r2_trial = r2_trial,
    beta = coef(fit)[2],
    se = summary(fit)$coefficients[2, 2],
    prediction_interval_width = mean(pred_interval[, 3] - pred_interval[, 2]),
    model = fit
  )
}


#' Calculate Surrogate Threshold Effect
#' @keywords internal
calculate_ste <- function(validation, bivariate_nma) {

  # STE: minimum surrogate effect needed to confidently predict true benefit
  # Based on lower bound of 95% prediction interval

  # Get prediction model: True ~ Surrogate
  beta <- validation$regression_coef
  se_beta <- validation$regression_se

  # Desired true effect (e.g., minimal clinically important difference)
  # For now, assume we want positive true effect
  desired_true_effect <- 0  # Could be parameterized

  # Calculate required surrogate effect
  # Using lower 95% PI bound
  z <- qnorm(0.975)

  # Solve: desired_true = beta * TE_surrogate - z * SE_pred
  # Rearrange: TE_surrogate = (desired_true + z * SE_pred) / beta

  ste <- (desired_true_effect + z * se_beta) / beta

  list(
    ste = ste,
    interpretation = paste0(
      "A surrogate effect of at least ", round(ste, 3),
      " is needed to confidently predict a positive true effect."
    )
  )
}


#' Predict True Effects from Surrogate
#' @keywords internal
predict_true_from_surrogate <- function(surrogate_effects, surrogacy_model) {

  beta <- surrogacy_model$regression_coef
  se_beta <- surrogacy_model$regression_se
  r2_trial <- surrogacy_model$r2_trial

  # Point prediction
  predicted_true <- beta * surrogate_effects

  # Prediction SE
  # SE_pred = sqrt(SE_beta^2 * TE_surrogate^2 + residual_var)
  residual_var <- (1 - r2_trial) * var(predicted_true)
  pred_se <- sqrt(se_beta^2 * surrogate_effects^2 + residual_var)

  # 95% Prediction intervals
  pred_lower <- predicted_true - 1.96 * pred_se
  pred_upper <- predicted_true + 1.96 * pred_se

  data.frame(
    surrogate_effect = surrogate_effects,
    predicted_true_effect = predicted_true,
    pred_se = pred_se,
    pred_lower = pred_lower,
    pred_upper = pred_upper
  )
}


#' S3 print method
#' @export
print.surrogate_nma <- function(x, ...) {
  cat("Surrogate Endpoint Network Meta-Analysis\n")
  cat("=========================================\n\n")

  cat("Surrogate Validation:\n")
  cat("  R² individual:", round(x$validation$r2_individual, 3), "\n")
  cat("  R² trial:", round(x$validation$r2_trial, 3), "\n")
  cat("  Interpretation:", x$validation$interpretation, "\n\n")

  cat("Surrogate Threshold Effect:\n")
  cat("  STE:", round(x$ste$ste, 3), "\n")
  cat("  ", x$ste$interpretation, "\n\n")

  cat("Use plot(x) to visualize surrogacy validation\n")
  cat("Use predict(x, new_surrogate_data) for predictions\n")
}


#' S3 plot method
#' @export
plot.surrogate_nma <- function(x, type = c("validation", "ste", "prediction"), ...) {

  type <- match.arg(type)

  if (type == "validation") {
    # Scatter: TE_true vs TE_surrogate at trial level
    # With regression line and prediction interval
  } else if (type == "ste") {
    # Plot showing STE threshold
  } else if (type == "prediction") {
    # Predicted true effects with uncertainty
  }
}
```

### Phase 2: Integration with Existing Functions (MEDIUM Priority)

**Enhance**: `multivariate_nma.R`

Add surrogate-specific correlation structures:
```r
# In multivariate_nma()
if (!is.null(surrogate_pair)) {
  # Use surrogate-validated correlation
  correlation_matrix <- build_surrogate_correlation(
    surrogate_var = surrogate_pair[1],
    true_var = surrogate_pair[2],
    validation = prior_validation
  )
}
```

### Phase 3: Plotting and Diagnostics (MEDIUM Priority)

**File to Create**: `powerNMA/R/surrogate_plots.R`

```r
#' Create surrogate validation plots
#' @export
plot_surrogate_validation <- function(data_s, data_t, ...) {

  # 1. Individual-level scatter (if IPD)
  # 2. Trial-level scatter
  # 3. Galbraith plot for surrogacy
  # 4. Information gain plot
}

#' Create STE plot
#' @export
plot_ste <- function(surrogate_nma_object, ...) {

  # Show threshold with confidence bands
  # Highlight treatments above/below STE
}
```

### Phase 4: Documentation and Examples (HIGH Priority)

**File to Create**: `powerNMA/vignettes/surrogate-endpoints.Rmd`

```markdown
# Surrogate Endpoint Network Meta-Analysis

## When to Use Surrogates

Surrogate endpoints are useful when:
- True endpoint takes long to measure (e.g., survival)
- Surrogate is measured earlier (e.g., tumor response)
- Want to include trials that only measured surrogate

## Validation Requirements

### Strong Surrogate (R² trial > 0.8)
- Can confidently predict true effects
- Use in primary analysis

### Moderate Surrogate (R² trial 0.5-0.8)
- Use in sensitivity analysis
- Report with caution

### Weak Surrogate (R² trial < 0.5)
- Do NOT use for treatment decisions
- Exploratory only

## Example

[Detailed worked example with real data]
```

---

## Integration Priority Matrix

| Feature | Impact | Effort | Priority | Blocking? |
|---------|--------|--------|----------|-----------|
| Surrogate validation (R²) | HIGH | MEDIUM | 🔥 HIGH | No |
| Bivariate NMA integration | HIGH | LOW | 🔥 HIGH | No |
| STE calculation | MEDIUM | LOW | ⚡ MEDIUM | No |
| Surrogate plots | MEDIUM | MEDIUM | ⚡ MEDIUM | No |
| Prentice criteria | LOW | HIGH | 💡 LOW | No |
| Vignette | HIGH | HIGH | ⚡ MEDIUM | No |

**None are blocking** for current publication, but surrogate features would make powerNMA more comprehensive.

---

## Questions for User

To proceed with integration, please clarify:

1. **Do you have access to a surronma repository?**
   - If yes, can you provide the path or URL?
   - If no, should I implement surrogate features from scratch based on literature?

2. **What specific surrogate features are priorities?**
   - R² validation?
   - Bivariate NMA?
   - STE calculation?
   - All of the above?

3. **Are there specific papers/methods to follow?**
   - Buyse et al. (2000) - Information theory approach?
   - Daniels & Hughes (1997) - Meta-analytic validation?
   - Other specific methodology?

4. **Timeline?**
   - Is this needed before journal submission?
   - Or post-publication enhancement?

---

## Alternative: If Repository is Available

If you can provide access to surronma repo:

```bash
# Steps I would take:
cd /home/user
git clone [surronma-repo-url]
cd surronma

# Analyze structure
find . -name "*.R" | head -20
grep -r "surrogate" --include="*.R" | head -50

# Identify key functions
grep -r "function.*surrogate" --include="*.R"

# Check for validation methods
grep -r "r2.*trial|r2.*individual" --include="*.R"

# Integrate useful features into powerNMA
```

---

## Current Status

- ✅ Analyzed current powerNMA for surrogate-related features
- ✅ Identified overlapping functionality (multivariate NMA)
- ✅ Designed integration plan for surrogate features
- ❌ Cannot access surronma repository (not found in environment)
- ⏳ Awaiting user clarification on:
  - Repository location
  - Feature priorities
  - Timeline

---

## Recommendation

**Option 1**: If surronma repo is available
→ Provide access and I'll integrate useful features

**Option 2**: If implementing from literature
→ Start with Phase 1 (Core functions) based on:
- Buyse M, et al. (2000). "The validation of surrogate endpoints"
- Burzykowski T, et al. (2005). "Evaluation of Surrogate Endpoints"

**Option 3**: Post-publication enhancement
→ Add surrogate features in version 3.1 after initial publication

---

**Status**: Awaiting user input on repository access or preferred approach
**Date**: November 1, 2025
