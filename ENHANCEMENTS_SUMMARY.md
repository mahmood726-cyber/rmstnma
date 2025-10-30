# powerNMA v1.1: Enhanced Features Summary

**Date:** October 30, 2025

**Version:** 1.0.0 → 1.1.0

---

## OVERVIEW

This update adds comprehensive features from CNMA v4.0 and systematic review workflow tools based on methodological review findings. The package now supports both advanced statistical methods AND systematic review best practices.

---

## NEW FEATURES

### A. ADVANCED STATISTICAL METHODS (from CNMA v4.0)

#### 1. **Transportability Weighting** ✨

Adjust evidence for target populations using covariate-based weighting.

```r
# Define target population
target <- list(
  age_mean = 70,
  female_pct = 0.50,
  bmi_mean = 29
)

# Compute weights
weights <- compute_transport_weights(
  data,
  target_population = target,
  metric = "mahalanobis",  # or "euclidean"
  kernel = "gaussian",      # or "tricube"
  truncation = 0.02
)

# Get diagnostics
diagnostics <- transportability_diagnostics(weights$weight)
# Effective sample size, weight distribution, etc.
```

**Features:**
- Mahalanobis and Euclidean distance metrics
- Gaussian and tricube kernel functions
- Weight truncation to avoid extreme values
- Comprehensive diagnostics (ESS, distribution stats)

---

#### 2. **PET-PEESE Publication Bias Analysis** 📊

Precision-Effect Test and Precision-Effect Estimate with Standard Error

```r
# Run PET-PEESE by comparison
pet_results <- pet_peese_analysis(data, min_k = 10)

# Returns:
# - PET intercept and p-value
# - PEESE intercept (if PET p < 0.10)
# - Recommendation on which to use
```

**Method:**
- PET: Regress effect on SE
- PEESE: Regress effect on SE²
- Automatic selection based on PET p-value
- Comparison-specific analysis

---

#### 3. **Leave-One-Treatment-Out (LOTO) Sensitivity** 🔍

Assess network robustness to treatment exclusions.

```r
loto_results <- loto_sensitivity(data, ref_treatment, sm = "HR")

# Shows impact on:
# - Heterogeneity (τ²)
# - I² statistic
# - Network structure
```

---

#### 4. **Multiverse Analysis** 🌌

Test robustness across modeling choices.

```r
multiverse_results <- run_multiverse(data, ref_treatment, sm = "HR")

# Tests:
# - Random effects only
# - Fixed effects only
# - Both models
# Compares heterogeneity estimates
```

---

#### 5. **Component Network Meta-Analysis** 🧩

Analyze additive effects of combination treatments.

```r
# For treatments like "DrugA", "DrugB", "DrugA+DrugB"
component_results <- run_component_nma(nma_fit)

# Estimates:
# - Individual component effects
# - Additivity assumption testing
# - Component interactions
```

---

### B. SYSTEMATIC REVIEW WORKFLOW TOOLS 📋

#### 1. **PRISMA-NMA Checklist Generation** ✅

Automated reporting standards compliance.

```r
protocol_info <- list(
  registration = "PROSPERO CRD42025123456",
  date = "2025-01-15"
)

checklist <- generate_prisma_nma_checklist(results, protocol_info)

# Returns complete PRISMA-NMA checklist with:
# - All required items
# - Completion status
# - Location in manuscript
# - Missing items flagged
```

**PRISMA-NMA Items Covered:**
- Study selection and characteristics
- Network geometry description
- Inconsistency assessment methods and results
- Heterogeneity evaluation
- Treatment ranking with uncertainties
- Statistical methods documentation

---

#### 2. **Network Diagram Creation** 🕸️

Publication-quality network visualization.

```r
network_plot <- create_network_diagram(
  data,
  output_file = "network_diagram.pdf"
)

# Features:
# - Node size by connectivity
# - Edge width by number of studies
# - Clean, publication-ready layout
# - Automatic layout optimization
# - Export to PDF/PNG/SVG
```

---

#### 3. **Risk of Bias Integration** ⚖️

Cochrane RoB 2.0 compatible assessment.

```r
# Prepare RoB data
rob_data <- data.frame(
  studlab = c("Study_01", "Study_02"),
  rob_randomization = c("Low", "Some concerns"),
  rob_deviations = c("Low", "Low"),
  rob_missing = c("Low", "High"),
  rob_measurement = c("Low", "Some concerns"),
  rob_selection = c("Low", "Low"),
  rob_overall = c("Low", "Some concerns")
)

# Integrate with data
data_with_rob <- integrate_rob_assessment(data, rob_data)

# Run RoB-restricted analysis
results_lowrob <- rob_sensitivity_analysis(
  data_with_rob,
  ref_treatment,
  rob_threshold = "Low"  # or c("Low", "Some concerns")
)
```

**Features:**
- Full RoB 2.0 domain structure
- Automatic RoB coverage checking
- Warnings for missing assessments
- Sensitivity analysis restricted to low RoB studies
- Comparison of all-studies vs low-RoB results

---

#### 4. **Protocol Specification & Adherence** 📝

Pre-specification and deviation tracking.

```r
# Specify protocol
protocol <- specify_protocol(
  research_question = "Effect of interventions on mortality",
  treatments = c("Placebo", "DrugA", "DrugB", "DrugC"),
  outcome = "All-cause mortality",
  effect_measure = "HR",
  reference = "Placebo",
  subgroups = c("age_group", "sex"),
  sensitivity = c("rob_restriction", "loto")
)

# Run analysis
results <- run_powernma(data, protocol = protocol)

# Check adherence
deviations <- check_protocol_adherence(results, protocol)
# Flags any analyses not pre-specified
```

**Benefits:**
- Reduces selective reporting
- Documents analysis plan
- Tracks deviations automatically
- PROSPERO-compatible structure

---

#### 5. **Methods Text Generation** ✍️

Automated manuscript methods section.

```r
methods_text <- generate_methods_text(results, format = "cochrane")

# Generates formatted text including:
# - Software and version
# - R version
# - Effect measure
# - Model type (random/fixed)
# - Heterogeneity assessment
# - Inconsistency methods
# - Sensitivity analyses performed
```

**Output Example:**
```
"We performed network meta-analysis using the powerNMA package
(version 1.1.0) in R (version 4.3.1). We used a random-effects
model with the hazard ratio as the summary measure. Network
geometry was assessed by calculating network density and
connectivity. Statistical inconsistency was evaluated using
the design-by-treatment interaction test. Heterogeneity was
quantified using tau-squared (τ²) and the I² statistic.
Treatment ranking was performed using P-scores..."
```

---

#### 6. **Evidence Profile Tables** 📊

GRADE-style Summary of Findings tables.

```r
evidence_profile <- generate_evidence_profile(
  results,
  comparisons = c("DrugA vs Placebo", "DrugB vs Placebo")
)

# Returns structured table with:
# - Treatment comparison
# - Effect estimate and 95% CI
# - Number of studies
# - Number of participants
# - Certainty of evidence (placeholder for GRADE)
```

---

#### 7. **Analysis Archiving** 📦

Complete reproducibility package.

```r
archive_analysis(
  results,
  data,
  output_dir = "analysis_archive_v1"
)

# Creates archive with:
# - Complete results (RDS)
# - Input data (CSV)
# - Configuration (JSON)
# - Session info (TXT)
# - README with reproduction instructions
```

**Archiv Structure:**
```
analysis_archive_v1/
├── results.rds           # Complete results object
├── data.csv              # Input data
├── config.json           # Analysis configuration
├── session_info.txt      # R session details
└── README.md             # Reproduction guide
```

**Benefits:**
- Full reproducibility
- Journal/repository submission ready
- Version control friendly
- Independent verification enabled

---

## ENHANCED FEATURES

### Improved Error Handling

All new functions include:
- Input validation
- Informative error messages
- Graceful degradation
- Warnings for edge cases

### Better Documentation

- Roxygen2 documentation for all new functions
- Usage examples
- Parameter descriptions
- Return value specifications
- Cross-references

### Comprehensive Exports

New NAMESPACE exports:
- 6 advanced analysis functions
- 9 systematic review workflow functions
- Additional utility functions

---

## UPDATED DEPENDENCIES

### New Required Packages:
- **igraph** - Network diagram creation
- **jsonlite** - Configuration export

### Additional Imports:
- `stats::coef`, `stats::lm`, `stats::cov`
- `stats::quantile`, `stats::median`, `stats::sd`
- `utils::sessionInfo`, `utils::capture.output`

---

## USAGE EXAMPLES

### Example 1: Complete Systematic Review Workflow

```r
library(powerNMA)

# 1. Specify protocol
protocol <- specify_protocol(
  research_question = "Comparative effectiveness of antihypertensives",
  treatments = c("Placebo", "ACEi", "ARB", "CCB", "Diuretic"),
  outcome = "Cardiovascular mortality",
  effect_measure = "HR",
  reference = "Placebo"
)

# 2. Prepare data with RoB
data <- simulate_nma_data(n_studies = 40)
rob_data <- # ... load from extraction
data_rob <- integrate_rob_assessment(data, rob_data)

# 3. Run main analysis
results <- run_powernma(
  data_rob,
  data_type = "pairwise",
  config = setup_powernma(sm = "HR")
)

# 4. Check protocol adherence
check_protocol_adherence(results, protocol)

# 5. Generate network diagram
create_network_diagram(data, "network_diagram.pdf")

# 6. Run sensitivity analyses
loto <- loto_sensitivity(data, "Placebo", "HR")
rob_sens <- rob_sensitivity_analysis(data_rob, "Placebo", rob_threshold = "Low")
multiverse <- run_multiverse(data, "Placebo", "HR")

# 7. Publication bias
pet_peese <- pet_peese_analysis(data, min_k = 10)

# 8. Generate reporting materials
prisma <- generate_prisma_nma_checklist(results, protocol_info)
methods <- generate_methods_text(results, format = "cochrane")
evidence <- generate_evidence_profile(results)

# 9. Archive complete analysis
archive_analysis(results, data, "final_analysis_2025")
```

### Example 2: Transportability Analysis

```r
# Define target population
target_us_adults <- list(
  age_mean = 65,
  female_pct = 0.52,
  bmi_mean = 28.5,
  charlson = 1.8
)

# Compute transport weights
weights <- compute_transport_weights(
  data,
  target_population = target_us_adults,
  metric = "mahalanobis",
  kernel = "gaussian"
)

# Check diagnostics
diag <- transportability_diagnostics(weights$weight)
print(diag)

# Run weighted analysis (would need to integrate with run_powernma)
```

### Example 3: Publication Bias Assessment

```r
# Run comprehensive bias assessment
pet_peese_results <- pet_peese_analysis(data, min_k = 10)

# Interpret results
# - Check PET p-values for asymmetry
# - Use PEESE intercept if PET significant
# - Compare to unadjusted estimates
```

---

## COMPARISON: v1.0 vs v1.1

| Feature Category | v1.0 | v1.1 |
|-----------------|------|------|
| **Core NMA** | ✓ | ✓ |
| **Time-Varying (RMST/Milestone)** | ✓ | ✓ |
| **Bayesian NMA** | ✓ | ✓ |
| **Basic Sensitivity (LOO)** | ✓ | ✓ |
| **Transportability Weighting** | ✗ | ✓ |
| **PET-PEESE** | ✗ | ✓ |
| **LOTO Sensitivity** | ✗ | ✓ |
| **Multiverse Analysis** | ✗ | ✓ |
| **Component NMA** | ✗ | ✓ |
| **PRISMA-NMA Checklist** | ✗ | ✓ |
| **Network Diagrams** | ✗ | ✓ |
| **RoB Integration** | ✗ | ✓ |
| **Protocol Support** | ✗ | ✓ |
| **Methods Text Generation** | ✗ | ✓ |
| **Evidence Profiles** | ✗ | ✓ |
| **Analysis Archiving** | ✗ | ✓ |
| **Total Functions** | 13 | 28 (+115%) |

---

## REMAINING LIMITATIONS

### Still Missing (from systematic review evaluation):

1. **GRADE Framework** - Certainty of evidence assessment (planned for v1.2)
2. **Meta-Regression** - With splines and CR2 (planned for v1.2)
3. **Selection Models** - weightr integration (planned for v1.2)
4. **Missing Data Handling** - Imputation methods (planned for v1.3)
5. **Multi-Arm Trial Support** - Proper correlation structure (planned for v1.3)
6. **Living SR Features** - Update workflows (planned for v2.0)

### Known Issues (from methodological review):

1. **IPD Reconstruction** - Still simplified (needs full Guyot implementation)
2. **RMST Sign Convention** - Needs mathematical proof
3. **Milestone Extend** - Extrapolation beyond follow-up
4. **Continuity Correction** - Non-standard approach

---

## DEVELOPMENT ROADMAP

### v1.2 (Next Release - 2-3 months)
- [ ] Full GRADE integration
- [ ] Meta-regression with splines
- [ ] Selection models (weightr)
- [ ] Cluster-robust SE (CR2)
- [ ] Improved continuity correction

### v1.3 (3-6 months)
- [ ] Missing data imputation
- [ ] Multi-arm trial methods
- [ ] Enhanced Bayesian features
- [ ] Additional publication bias methods

### v2.0 (6-12 months)
- [ ] Living systematic review support
- [ ] Web-based interface
- [ ] Complete CINeMA integration
- [ ] RevMan compatibility
- [ ] Automated report generation

---

## TESTING

New tests required for v1.1 features:

```r
# Transport weights
test_that("Transport weights computed correctly", {
  # Test Mahalanobis and Euclidean
  # Test kernel functions
  # Test weight normalization
})

# PRISMA checklist
test_that("PRISMA checklist generation works", {
  # Test with minimal results
  # Test with complete results
  # Test completion status
})

# RoB integration
test_that("RoB integration handles edge cases", {
  # Missing RoB for some studies
  # Different RoB thresholds
  # RoB-restricted analysis
})

# Protocol adherence
test_that("Protocol deviation detection works", {
  # Matching protocol
  # Deviating protocol
  # Missing protocol elements
})
```

---

## DOCUMENTATION UPDATES

### New Vignettes Needed:

1. **"Advanced Analysis with powerNMA"**
   - Transportability weighting
   - Publication bias assessment
   - Sensitivity analyses

2. **"Systematic Review Workflow"**
   - Protocol specification
   - PRISMA-NMA compliance
   - RoB integration
   - Report generation

3. **"Interpretation Guide"**
   - Understanding heterogeneity
   - Assessing inconsistency
   - Treatment ranking
   - Publication bias

---

## MIGRATION GUIDE: v1.0 → v1.1

### Backward Compatibility

✅ **All v1.0 code still works** - No breaking changes

### New Recommended Workflow

**Before (v1.0):**
```r
results <- run_powernma(data, data_type = "pairwise")
```

**Now (v1.1):**
```r
# 1. Specify protocol (recommended)
protocol <- specify_protocol(...)

# 2. Integrate RoB (if available)
data <- integrate_rob_assessment(data, rob_data)

# 3. Run main analysis
results <- run_powernma(data, data_type = "pairwise")

# 4. Generate PRISMA checklist
prisma <- generate_prisma_nma_checklist(results)

# 5. Archive analysis
archive_analysis(results, data)
```

---

## ACKNOWLEDGMENTS

This update integrates methods from:
- **CNMA v4.0** - Advanced NMA techniques
- **CINeMA** - GRADE and PRISMA-NMA standards
- **Cochrane Handbook** - Systematic review best practices
- **PRISMA-NMA Statement** - Reporting standards

---

## SUMMARY

**v1.1 Status:**
- ✅ 15 new major features
- ✅ Enhanced systematic review workflow
- ✅ Improved methodological rigor
- ✅ Better reporting compliance
- ✅ Full backward compatibility

**Next Steps:**
- Extensive testing of new features
- Comprehensive vignettes
- GRADE framework (v1.2)
- Meta-regression (v1.2)

**Current Suitability:**
- ✓ Methods research
- ✓ Exploratory analysis
- ✓ Systematic reviews (with caveats)
- △ Clinical guidelines (needs GRADE - v1.2)
- △ Regulatory submissions (needs validation)

---

**Version:** 1.1.0

**Release Date:** October 30, 2025

**Package Status:** Enhanced - Additional features for systematic review workflow

**License:** MIT
