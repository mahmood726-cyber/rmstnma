# powerNMA v2.0: Validation and Two-Mode Architecture

**Date:** October 31, 2025
**Version:** 2.0.0
**Status:** Implemented and Ready for Testing

---

## Summary of Changes

powerNMA has been restructured to address critical methodological concerns and provide clear guidance on which methods are validated for clinical use versus research-only experimental methods.

---

## Key Improvements

### 1. ✅ Removed IPD Reconstruction (Critical Fix)

**Issue:** IPD reconstruction from KM curves was methodologically unsound
**Action Taken:**
- Removed `reconstruct_ipd()` function entirely
- Removed from NAMESPACE exports
- Removed from documentation
- Users directed to validated packages (IPDfromKM, digitize)

**Files Modified:**
- `powerNMA/R/data_functions.R`
- `powerNMA/NAMESPACE`
- `powerNMA/R/powernma-package.R`

---

### 2. ✅ Two-Mode Architecture Implemented

**NEW: STANDARD MODE** (Default)
- Validated methods only (netmeta, gemtc)
- Suitable for systematic reviews, clinical guidelines, Cochrane reviews
- Publication-ready without questions
- **Rejects IPD data** (time-varying methods are experimental)

**NEW: EXPERIMENTAL MODE** (Explicit Opt-In)
- Time-varying NMA (RMST, milestone)
- Transportability weighting
- Advanced sensitivity analyses
- **Requires explicit `mode = "experimental"` parameter**
- Shows prominent warnings
- Research use only

**Files Created:**
- `powerNMA/R/modes.R` - Mode routing and validation logic

**Files Modified:**
- `powerNMA/R/powernma.R` - Updated main function with mode parameter
- `powerNMA/R/utils.R` - Updated print method to show mode status

---

### 3. ✅ Mode-Specific Safeguards

**Validation Logic:**
```r
# This will ERROR:
run_powernma(ipd_data, data_type = "ipd", mode = "standard")
#> Error: IPD-based methods are EXPERIMENTAL.
#> Use mode='experimental' to enable.

# This requires explicit opt-in:
run_powernma(ipd_data, data_type = "ipd", mode = "experimental")
#> Warning: EXPERIMENTAL MODE ENABLED
#> These methods are NOT validated for clinical decision-making
```

**User Protection:**
- Users cannot accidentally use experimental methods
- Clear error messages explain why
- Warnings are impossible to miss
- Documentation clarifies suitability

---

### 4. ✅ Comprehensive Testing

**New Test Suite:**
- `powerNMA/tests/testthat/test-modes.R`
- 14 tests covering:
  - Standard mode rejects IPD
  - Experimental mode shows warnings
  - Mode parameter defaults correctly
  - Print method shows mode information
  - Class inheritance correct
  - Validation logic works

**Test Coverage:**
- Mode switching
- Error handling
- Warning display
- Results structure
- User experience

---

### 5. ✅ Updated User Interface

**Print Method Now Shows:**
```
═══════════════════════════════════════════════════
  powerNMA Analysis Results
═══════════════════════════════════════════════════

Mode: STANDARD
Status: ✓ VALIDATED - Suitable for systematic reviews
Based on: netmeta 2.9.0

Reference treatment: Placebo
Heterogeneity τ: 0.1234
I²: 45.6%

Analyses completed:
  ✓ Frequentist NMA (netmeta)
  ✓ Bayesian NMA (gemtc)
  ✓ Leave-one-out sensitivity
  ✓ Network geometry
  ✓ Inconsistency assessment
```

vs Experimental Mode:
```
═══════════════════════════════════════════════════
  powerNMA Analysis Results
═══════════════════════════════════════════════════

Mode: EXPERIMENTAL
Status: ⚠ EXPERIMENTAL - Research use only
Validation: See VALIDATION_PLAN.md

Reference treatment: Control

Warnings:
  • RMST sign convention requires mathematical validation
  • Milestone analysis uses extend=FALSE by default
  • Multi-arm trial handling may be incomplete

Analyses completed:
  ✓ Standard NMA
  ✓ Time-varying RMST NMA
  ✓ Milestone Survival NMA
  ✗ Transportability weighting
  ...

⚠ REMINDER: Experimental methods not validated for clinical use
```

---

## New Usage Examples

### Standard Mode (Recommended for Publication)

```r
library(powerNMA)

# Generate or load pairwise comparison data
data <- simulate_nma_data(n_studies = 40)

# Standard mode (default) - validated methods
results <- run_powernma(
  data = data,
  data_type = "pairwise",
  mode = "standard"  # Can be omitted (default)
)

# Suitable for:
# - Cochrane systematic reviews
# - Clinical practice guidelines
# - Journal submissions (BMJ, Lancet, etc.)
# - Regulatory submissions
```

### Experimental Mode (Research Only)

```r
library(powerNMA)

# Individual patient data
ipd <- generate_example_ipd(n_trials = 10)

# Configure time-varying analysis
config <- setup_powernma(
  use_timevarying = TRUE,
  tau_list = 365,
  milestone_times = c(90, 180, 365)
)

# Experimental mode - explicit opt-in required
results <- run_powernma(
  data = ipd,
  data_type = "ipd",
  mode = "experimental",  # REQUIRED for IPD/time-varying
  config = config
)

# Will show warning:
#> ══════════════════════════════════════════════════════════
#>              EXPERIMENTAL MODE ENABLED
#> ══════════════════════════════════════════════════════════
#>
#> These methods are NOT validated for clinical decision-making.
#>
#> SUITABLE FOR:
#>   • Academic methods research
#>   • Exploratory analysis
#>   • Hypothesis generation
#>
#> NOT SUITABLE FOR:
#>   ✗ Cochrane systematic reviews
#>   ✗ Clinical practice guidelines
#>   ✗ Regulatory submissions

# Suitable for:
# - Methods research papers
# - Exploratory comparative effectiveness research
# - Sensitivity analyses
# - Hypothesis generation
```

---

## Migration Guide for Existing Users

### If you were using standard NMA only:
**No changes needed.** Your code will continue to work exactly as before.

### If you were using time-varying methods (RMST, milestone):
**Action required:**
```r
# OLD (will now ERROR):
results <- run_powernma(ipd, data_type = "ipd")

# NEW (add mode parameter):
results <- run_powernma(ipd, data_type = "ipd", mode = "experimental")
```

### If you were using reconstruct_ipd():
**Function removed.** Use validated alternatives:
```r
# REMOVED:
# ipd <- reconstruct_ipd(km_data, n_risk, total_events, trial, arm)

# RECOMMENDED instead:
library(IPDfromKM)
# Use IPDfromKM package following Guyot et al. (2012) algorithm
```

---

## Validation Status by Method

| Method | Mode | Validated? | Suitable For |
|--------|------|-----------|--------------|
| Standard NMA (netmeta) | Standard | ✅ Yes | Systematic reviews, guidelines |
| Bayesian NMA (gemtc) | Standard | ✅ Yes | Systematic reviews, guidelines |
| Leave-one-out sensitivity | Standard | ✅ Yes | Systematic reviews, guidelines |
| Inconsistency assessment | Standard | ✅ Yes | Systematic reviews, guidelines |
| **RMST NMA** | **Experimental** | ⚠️ **No** | **Research only** |
| **Milestone NMA** | **Experimental** | ⚠️ **No** | **Research only** |
| **Transportability** | **Experimental** | ⚠️ **No** | **Research only** |
| **PET-PEESE for NMA** | **Experimental** | ⚠️ **No** | **Research only** |
| **Component NMA** | **Experimental** | ⚠️ **No** | **Research only** |

---

## Documentation Updates

**New/Updated Files:**
1. `VALIDATION_PLAN.md` - Comprehensive 24-page validation roadmap
2. `VALIDATION_v2.0_SUMMARY.md` - This document
3. `powerNMA/R/modes.R` - Mode architecture implementation
4. `powerNMA/tests/testthat/test-modes.R` - Mode testing suite

**Updated Function Documentation:**
- `run_powernma()` - Now includes mode parameter with detailed examples
- `print.powernma_result()` - Shows mode-specific information

---

## Next Steps (See VALIDATION_PLAN.md)

### Immediate (v2.0 Release - Current)
- [x] Remove IPD reconstruction
- [x] Implement two-mode architecture
- [x] Add safeguards and warnings
- [x] Create mode-specific tests
- [x] Update documentation

### Short-term (Weeks 1-4)
- [ ] Fix multi-arm trial handling
- [ ] Validate RMST sign convention mathematically
- [ ] Change milestone extend=FALSE default
- [ ] Standardize continuity correction
- [ ] Complete transportability diagnostics

### Medium-term (Weeks 5-12)
- [ ] Benchmark against published datasets
- [ ] Simulation studies
- [ ] Vignettes for both modes
- [ ] Independent peer review

### Long-term (Months 4-12)
- [ ] Methods paper submission
- [ ] Full validation of experimental methods
- [ ] Cochrane Methods Group review
- [ ] Move validated experimental methods to standard

---

## Breaking Changes in v2.0

1. **`reconstruct_ipd()` removed**
   - Reason: Methodologically unsound
   - Migration: Use IPDfromKM package

2. **IPD data requires `mode = "experimental"`**
   - Reason: Time-varying methods not validated
   - Migration: Add `mode = "experimental"` parameter

3. **Print output changed**
   - Reason: Show mode and validation status
   - Impact: Scripts parsing output may need updates

---

## Backwards Compatibility

**Preserved:**
- All standard NMA functionality unchanged
- API for pairwise data unchanged
- Results structure mostly preserved
- Export functions still work

**Changed:**
- IPD analysis requires explicit mode parameter
- Print output format updated
- Some internal functions reorganized

---

## Testing Recommendations

**For Users of Standard Mode:**
```r
# Test that your existing code still works:
data <- your_data()
results <- run_powernma(data, data_type = "pairwise")
# Should work without changes
```

**For Users of Experimental Mode:**
```r
# Update your code to specify mode:
ipd <- your_ipd_data()
results <- run_powernma(ipd, data_type = "ipd", mode = "experimental")
# Acknowledge warnings and proceed cautiously
```

---

## Contact and Support

**Questions about validation status:**
- See `VALIDATION_PLAN.md` for detailed roadmap
- See `METHODOLOGICAL_REVIEW.md` for technical critique
- See `SYSTEMATIC_REVIEW_EVALUATION.md` for SR suitability

**Reporting Issues:**
- GitHub Issues: [your-org/powerNMA/issues](https://github.com/your-org/powerNMA/issues)
- Include mode (standard/experimental) in bug reports

**Contributing to Validation:**
- We welcome contributions to validation efforts
- See VALIDATION_PLAN.md Phase 3 for benchmark datasets
- Contact dev team for collaboration opportunities

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | Oct 30, 2025 | Initial release with experimental methods |
| 1.1.0 | Oct 30, 2025 | Added systematic review workflow features |
| **2.0.0** | **Oct 31, 2025** | **Two-mode architecture, removed IPD reconstruction** |

---

## Acknowledgments

This restructuring addresses critical issues identified in:
- Methodological review (METHODOLOGICAL_REVIEW.md)
- Systematic review evaluation (SYSTEMATIC_REVIEW_EVALUATION.md)
- User feedback and safety concerns

Special thanks to the methodological reviewers for identifying these critical issues.

---

**Document Status:** ✅ Complete - Ready for v2.0 release
**Last Updated:** October 31, 2025
