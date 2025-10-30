# METHODOLOGICAL REVIEW: powerNMA R Package

**Review for:** Research Synthesis Methods / Systematic Reviews

**Manuscript Type:** Software/Methods Article

**Reviewer:** Research Synthesis Methods Specialist

**Date:** October 30, 2025

---

## SUMMARY RECOMMENDATION

**Decision:** Major Revisions Required

**Overall Assessment:** The powerNMA package represents an ambitious attempt to unify time-varying network meta-analysis methods with advanced synthesis techniques. While the package demonstrates technical competence and addresses an important gap, several **critical methodological concerns** must be addressed before publication in a methods journal. The software shows promise but requires substantial methodological refinement, validation, and clearer documentation of assumptions and limitations.

---

## MAJOR CONCERNS

### 1. IPD Reconstruction: Fundamental Methodological Flaws

**Location:** `R/data_functions.R`, lines 20-100

**Issue:** The IPD reconstruction algorithm is oversimplified to the point of being **methodologically unsound** for research use.

**Critical Problems:**

```r
# Lines 53-54
if (s_start > 0) {
  n_events <- round(n_start * (s_start - s_end) / s_start)
}
```

**Methodological Critique:**

1. **Lacks Iterative Refinement:** The Guyot et al. (2012) algorithm requires iterative adjustment to match reported summary statistics. This implementation is a single-pass approximation.

2. **Uniform Distribution Assumption:** Assumes events are uniformly distributed within intervals (line 62), which is **rarely valid** in clinical data where hazards typically vary.

3. **No Uncertainty Quantification:** Reconstructed IPD should carry forward uncertainty from the digitization process. This is not captured.

4. **Censoring Time Assumption:** Places all censored observations at interval endpoints (line 69), which **systematically biases** Kaplan-Meier estimates.

5. **No Validation Against Known Ground Truth:** No evidence of validation against datasets with both published KM curves and available IPD.

**Recommendation:**
- Either implement the full Guyot algorithm correctly OR
- Remove this function entirely and recommend validated packages (IPDfromKM, digitize)
- Add **prominent warnings** that reconstructed IPD should be sensitivity-tested
- Provide guidance on when reconstruction is appropriate vs when IPD request is essential

**Severity:** **Critical** - This could lead to biased meta-analytic conclusions

---

### 2. RMST Sign Convention: Inadequately Justified

**Location:** `R/nma_core.R`, line 181

```r
TE = -1 * rmst_result$unadjusted.result[1, "Est."]
```

**Issue:** Manual sign flipping without rigorous statistical justification.

**Concerns:**

1. **Insufficient Documentation:** The code comment says "need to flip" but doesn't explain **why** this is statistically necessary.

2. **survRM2 Output Interpretation:** The `rmst2()` function returns `arm=1 - arm=0`. The pairwise structure defines `treat1` and `treat2`. The relationship between these should be **mathematically proven**, not empirically adjusted.

3. **Risk of Systematic Error:** If the sign convention is incorrect, **all treatment effect estimates will have reversed direction**. This is catastrophic.

4. **No Validation:** No test case comparing manual RMST calculations to package output to verify correctness.

**Recommendation:**
- Provide **mathematical proof** of the sign convention in supplementary methods
- Add validation test: calculate RMST manually for a simple dataset, verify package matches
- Consider restructuring code to avoid sign manipulation entirely
- Add assertion checks that verify treatment effect direction makes sense (e.g., known effective treatment should favor active arm)

**Severity:** **Critical** - Affects validity of all RMST results

---

### 3. Milestone Analysis: Extend Parameter Creates Bias

**Location:** `R/nma_core.R`, line 249

```r
summ1 <- summary(km1, times = t, extend = TRUE)
```

**Issue:** The `extend = TRUE` parameter **extrapolates** survival curves beyond observed follow-up.

**Statistical Problem:**

When `t` (milestone time) exceeds maximum follow-up in a trial:
- Kaplan-Meier is extrapolated using the last observed hazard
- This **assumes constant hazard** beyond observed data
- Systematically **overestimates** survival if hazard is increasing
- Creates **pseudo-data** that appears real but is speculative

**Methodological Concern:**

In network meta-analysis, this can create **spurious heterogeneity** or **spurious consistency** depending on which trials have shorter follow-up. Trials with complete follow-up to time `t` are compared to trials with extrapolated estimates.

**Example:**
- Trial A: Max follow-up = 180 days
- Trial B: Max follow-up = 400 days
- Milestone time = 365 days
- Trial A uses extrapolated estimate (biased)
- Trial B uses observed estimate (valid)
- Comparison is **invalid**

**Recommendation:**
- Set `extend = FALSE` by default
- Check maximum follow-up time for each trial
- Warn or exclude trials where milestone exceeds follow-up
- Document this limitation explicitly
- Provide guidance on selecting appropriate milestone times

**Severity:** **Major** - Affects validity when follow-up varies

---

### 4. Continuity Correction: Non-Standard Implementation

**Location:** `R/nma_core.R`, lines 260-272

**Issue:** Uses 0.5 correction for zero cells, which is **not universally accepted** for odds ratios in meta-analysis.

**Methodological Context:**

Cochrane guidance (Higgins et al., 2023) recommends:
1. **Treatment-arm continuity correction** (add 0.5 to all arms in zero-event studies)
2. **Consider excluding** double-zero studies entirely
3. Use **alternative effect measures** (risk difference, Peto OR) for sparse data

This implementation:
```r
events1 <- events1 + 0.5
events2 <- events2 + 0.5
n1 <- n1 + 1
n2 <- n2 + 1
```

Adds 0.5 events AND increases denominator by 1, which is **unconventional**.

**Statistical Implications:**
- Standard correction: Add 0.5 to numerator only
- This approach: Adds 0.5 to numerator AND 1 to denominator
- Changes the effective rate differently than standard methods

**Concern:** Results may not be comparable to meta-analyses using standard methods.

**Recommendation:**
- Use standard 0.5 continuity correction (add to numerator only)
- Make correction method a configurable parameter
- Document choice and provide sensitivity analysis option
- Cite methodological literature supporting chosen approach
- Consider Peto OR method for sparse data as alternative

**Severity:** **Major** - Affects comparability and reproducibility

---

### 5. Transportability Weighting: Incomplete Implementation

**Location:** Reference to methods but not fully implemented in reviewed code

**Missing Components:**

1. **No Covariate Balance Assessment:** Should report standardized mean differences (SMDs) for target population vs trial populations

2. **No Positivity Checks:** Transportability requires **positivity assumption** - that target population characteristics are within the convex hull of trial populations. No diagnostic provided.

3. **No Effective Sample Size:** Should report effective N after weighting (N_eff = (Σw)² / Σw²)

4. **No Sensitivity Analysis:** Should explore impact of different distance metrics and kernel choices

5. **Limited Guidance:** No decision rule for when transportability is appropriate vs when it introduces more bias

**Methodological Literature:**

Recent work (Dahabreh et al. 2020, Stuart et al. 2011) emphasizes:
- **Covariate overlap** is essential
- **Strong confounding assumptions** required
- Results can be **highly sensitive** to weighting choices

**Current Implementation Risk:**

Without diagnostics, users may apply transportability inappropriately, creating **false precision** about target population effects.

**Recommendation:**
- Implement balance diagnostics (before/after weighting SMDs)
- Add positivity checks with warnings
- Report effective sample size
- Provide sensitivity analysis across weighting schemes
- Add extensive documentation about assumptions and limitations
- Include example showing when transportability fails

**Severity:** **Major** - Core methodological feature needs completion

---

### 6. Multi-Arm Trials: Incorrectly Handled

**Location:** `R/nma_core.R`, RMST and Milestone functions

**Issue:** Code iterates through pairwise comparisons within multi-arm trials:

```r
for (trial_id in unique(ipd$trial)) {
  trial_data <- ipd[ipd$trial == trial_id, ]
  treatments <- sort(unique(trial_data$treatment))

  if (length(treatments) < 2) next

  # Selects only ONE comparison even if >2 arms
  treat1 <- ...
  treat2 <- setdiff(treatments, reference)[1]
```

**Problem:** For a 3-arm trial (A vs B vs C), this only analyzes **ONE comparison** and discards information from other arms.

**Correct Approach (Standard in NMA):**
1. Account for correlation between comparisons from same trial
2. Either:
   - Use multi-arm extensions of meta-analysis models
   - Split into all pairwise comparisons with variance adjustment
   - Use contrast-based approach with covariance structure

**Impact:** Loss of precision and potentially biased treatment rankings.

**Recommendation:**
- Implement proper multi-arm trial handling
- Document clearly if multi-arm trials are not supported
- Provide guidance to users on how to handle multi-arm trials
- Consider using netmeta's built-in multi-arm capabilities

**Severity:** **Major** - Common data structure handled incorrectly

---

## MODERATE CONCERNS

### 7. Reference Treatment Selection

**Issue:** Auto-selection uses most connected treatment (line 67 in powernma.R).

**Problem:** Most connected ≠ clinically meaningful comparator. For time-varying analyses, should be treatment with most complete follow-up.

**Recommendation:** Allow user specification with clear warnings; provide guidance on appropriate reference choice.

### 8. Fixed Seed for Reproducibility

**Issue:** Default seed = 42 (line 32 in nma_core.R)

**Problem:** Users may not realize results depend on this seed; IPD reconstruction includes randomness.

**Recommendation:**
- Document that IPD reconstruction is stochastic
- Recommend users run multiple seeds and assess variability
- Provide function to assess reconstruction uncertainty

### 9. Heterogeneity Assessment

**Issue:** Reports τ and I² but no prediction intervals, between-study variance CIs, or heterogeneity exploration beyond ML approach.

**Recommendation:**
- Add 95% prediction intervals for treatment effects
- Provide tau² confidence intervals
- Implement Cochran's Q test formally
- Guide users on interpreting heterogeneity magnitudes

### 10. Missing Reporting Standards

**Issue:** No integration with reporting guidelines (PRISMA-NMA, GRADE).

**Recommendation:**
- Provide checklist generators for PRISMA-NMA
- Include GRADE assessment framework integration
- Output structured tables meeting journal requirements

---

## MINOR CONCERNS

### 11. Documentation of Assumptions

**Issue:** Insufficient documentation of assumptions for:
- Consistency/transitivity assumption in NMA
- Proportional hazards for time-varying methods
- Similarity assumption for transportability
- Missing data mechanisms

**Recommendation:** Add detailed assumptions section to package documentation.

### 12. Publication Bias Methods

**Issue:** PET-PEESE, selection models, etc. mentioned but not implemented in reviewed code.

**Recommendation:** Either implement fully or remove from documentation.

### 13. Model Diagnostics

**Issue:** Limited model diagnostics (residuals, influential studies, leverage).

**Recommendation:** Add diagnostic plots and influence statistics.

### 14. Effect Measure Guidance

**Issue:** No guidance on choosing sm (HR vs OR vs RR) for time-to-event data.

**Recommendation:** Provide decision tree for effect measure selection.

---

## STATISTICAL VALIDITY ASSESSMENT

### Frequentist NMA (Standard Component)
**Assessment:** ✅ **Valid** - Uses established netmeta package correctly

**Evidence:**
- Proper random-effects model implementation
- Correct handling of study-level effects
- Appropriate inference framework

### RMST NMA
**Assessment:** ⚠️ **Potentially Valid** - But requires verification

**Concerns:**
- Sign convention needs mathematical proof
- Multi-arm trial handling incomplete
- No validation against known results

**Recommendation:** Provide worked example with hand-calculated results

### Milestone NMA
**Assessment:** ⚠️ **Conditionally Valid** - If follow-up adequate

**Concerns:**
- `extend = TRUE` can introduce bias
- Continuity correction non-standard
- Need follow-up adequacy checks

**Recommendation:** Add follow-up diagnostics; validate continuity correction choice

### IPD Reconstruction
**Assessment:** ❌ **Not Valid** for research use

**Verdict:** Oversimplified algorithm; recommend removal or complete rewrite

### Transportability
**Assessment:** ⚠️ **Incomplete** - Cannot assess without full implementation

**Concerns:**
- Missing key diagnostics
- No sensitivity analyses
- Insufficient user guidance

---

## CODE QUALITY ASSESSMENT

### Strengths
✅ Clear function structure and naming
✅ Comprehensive error handling
✅ Good use of existing packages
✅ Unit tests present
✅ roxygen2 documentation

### Weaknesses
❌ Insufficient statistical validation
❌ Limited integration tests
❌ No benchmark comparisons to existing methods
❌ Sparse code comments explaining statistical logic

---

## COMPARISON TO EXISTING METHODS

### Advantages Over Existing Software
1. **Unified interface** - Combines methods typically requiring multiple packages
2. **Time-varying capability** - Novel integration of RMST with NMA
3. **User-friendly** - Lower barrier to entry than individual packages

### Disadvantages vs Existing Software
1. **Less validated** - netmeta, gemtc have extensive validation
2. **Simplified methods** - Some implementations sacrifice methodological rigor
3. **Limited documentation** - Statistical assumptions less explicit than specialist packages

---

## RECOMMENDATIONS FOR AUTHORS

### Essential (Must Address for Acceptance)

1. **Fix or Remove IPD Reconstruction**
   - Current implementation is methodologically flawed
   - Either implement Guyot algorithm correctly or remove entirely
   - If kept, add **bold warnings** about limitations

2. **Verify RMST Sign Convention**
   - Provide mathematical proof
   - Add validation tests
   - Document thoroughly

3. **Fix Milestone Extend Parameter**
   - Default to `extend = FALSE`
   - Check follow-up adequacy
   - Document limitation

4. **Standardize Continuity Correction**
   - Use established methods
   - Make configurable
   - Document choice with citations

5. **Implement Multi-Arm Trials**
   - Current approach discards information
   - Use proper correlation structure

6. **Complete Transportability Methods**
   - Add diagnostics
   - Provide sensitivity analyses
   - Document assumptions

### Important (Should Address)

7. Improve heterogeneity assessment (prediction intervals, τ² CIs)
8. Add model diagnostics and influence statistics
9. Provide reporting checklist generators
10. Add worked examples with validation

### Recommended (Enhance Quality)

11. Benchmark against existing packages
12. Add more comprehensive tests
13. Improve statistical documentation
14. Create detailed vignettes with real data examples

---

## VALIDATION REQUIREMENTS

Before publication, authors must demonstrate:

1. **Mathematical Correctness**
   - Prove sign conventions
   - Verify variance calculations
   - Document all transformations

2. **Numerical Accuracy**
   - Compare to hand calculations
   - Benchmark against established packages
   - Validate continuity correction effects

3. **Reproducibility**
   - Provide complete worked examples
   - Show sensitivity to random seed
   - Document stochastic components

4. **Clinical Validity**
   - Apply to published datasets
   - Compare results to published meta-analyses
   - Show results are consistent with known findings

---

## SPECIFIC RECOMMENDATIONS BY FUNCTION

### `rmst_nma()`
- [ ] Add validation test against manual calculation
- [ ] Document sign convention with proof
- [ ] Implement multi-arm trial support
- [ ] Add sensitivity analysis for tau selection
- [ ] Check proportional hazards assumption

### `milestone_nma()`
- [ ] Change extend = FALSE default
- [ ] Add follow-up adequacy checks
- [ ] Use standard continuity correction
- [ ] Provide Peto OR alternative for sparse data
- [ ] Add sensitivity analysis for milestone time selection

### `reconstruct_ipd()`
- [ ] Implement full Guyot algorithm OR remove
- [ ] Add uncertainty quantification
- [ ] Provide validation examples
- [ ] Add strong warnings
- [ ] Recommend validated alternatives

### `robust_netmeta()`
- [ ] Add diagnostic plots
- [ ] Provide influence statistics
- [ ] Check consistency assumption
- [ ] Add prediction intervals

---

## SUITABILITY FOR PUBLICATION

**Current Status:** Not suitable for publication without major revisions

**Rationale:**
- Several **methodologically critical** issues
- Incomplete implementation of advertised features
- Insufficient validation
- Limited documentation of assumptions

**Path to Acceptance:**

1. Address all ESSENTIAL recommendations
2. Provide comprehensive validation
3. Add worked examples with known results
4. Document all statistical assumptions clearly
5. Provide comparison to existing methods

**Estimated Revision Scope:** 3-6 months of methodological work

---

## ADDITIONAL COMMENTS

### For a Software Journal (e.g., Journal of Statistical Software)

The package would need:
- More extensive vignettes with real data
- Complete implementation of all advertised features
- Benchmark comparisons
- User guide with interpretation guidance

### For a Methods Journal (e.g., Research Synthesis Methods)

The package would need:
- Clear methodological contribution (what's novel?)
- Validation against existing methods
- Simulation study demonstrating performance
- Discussion of when methods are appropriate

### For Clinical Journal (e.g., BMC Medical Research Methodology)

The package would need:
- Clinical examples with interpretation
- Comparison to published meta-analyses
- Guidance for applied researchers
- Discussion of practical considerations

---

## CONCLUSION

The powerNMA package represents an ambitious and potentially valuable contribution to evidence synthesis methodology. However, several **critical methodological issues** must be addressed before the software is suitable for research use. The authors demonstrate programming competence, but need to strengthen the statistical rigor and validation of their methods.

**Primary Concerns:**
1. IPD reconstruction algorithm is oversimplified
2. RMST sign convention inadequately justified
3. Milestone analysis may extrapolate inappropriately
4. Transportability methods incomplete
5. Multi-arm trials handled incorrectly

**Recommendation:** **Major Revisions Required**

The reviewers encourage resubmission after addressing the methodological concerns outlined above. With appropriate revision, this could become a valuable tool for the research synthesis community.

---

**Reviewer Expertise:** Network meta-analysis, survival analysis, evidence synthesis methods

**Conflict of Interest:** None

**Review Completed:** October 30, 2025
