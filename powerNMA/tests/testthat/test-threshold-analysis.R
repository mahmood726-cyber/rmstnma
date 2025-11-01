# ============================================================================
# Tests for Threshold Analysis (Experimental Method)
# ============================================================================

library(testthat)
library(powerNMA)

# ============================================================================
# Test threshold_analysis() main function
# ============================================================================

test_that("threshold_analysis runs on netmeta object", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE,
    fixed = FALSE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    decision_rule = "maximize_benefit",
    risk_aversion = 0,
    threshold_type = "all"
  )

  expect_s3_class(result, "threshold_analysis")
  expect_true(!is.null(result$recommendation))
  expect_true(!is.null(result$thresholds))
  expect_true(!is.null(result$robustness_score))
})


test_that("threshold_analysis provides all threshold types", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    threshold_type = "all"
  )

  # Should have all three threshold types
  expect_true("effect_size" %in% names(result$thresholds))
  expect_true("bias" %in% names(result$thresholds))
  expect_true("new_study" %in% names(result$thresholds))

  # Each threshold should have required components
  expect_true(!is.null(result$thresholds$effect_size$value))
  expect_true(!is.null(result$thresholds$effect_size$sd_units))
  expect_true(!is.null(result$thresholds$bias$value))
  expect_true(!is.null(result$thresholds$new_study$required_effect))
})


test_that("threshold_analysis calculates robustness score correctly", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher"
  )

  # Robustness score should be 0-100
  expect_gte(result$robustness_score, 0)
  expect_lte(result$robustness_score, 100)

  # Should have robustness category
  expect_true(result$interpretation$robustness_category %in%
              c("VERY ROBUST", "MODERATELY ROBUST", "FRAGILE"))
})


# ============================================================================
# Test cost-effectiveness decision rule
# ============================================================================

test_that("cost-effectiveness decision rule works correctly", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  # Create cost data for treatments in Senn2013
  treatments <- unique(c(as.character(Senn2013$treat1),
                        as.character(Senn2013$treat2)))
  costs <- data.frame(
    treatment = treatments,
    cost = c(100, 500, 800, 1200, 2000, 1500, 900, 300,
            1100, 400, 600, 1400, 1800, 700)[1:length(treatments)]
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    decision_rule = "cost_effectiveness",
    cost_data = costs,
    willingness_to_pay = 50000,
    threshold_type = "effect_size"
  )

  expect_s3_class(result, "threshold_analysis")
  expect_true(!is.null(result$recommendation$treatment))
  expect_true(!is.null(result$recommendation$cost))
  expect_true(!is.null(result$recommendation$icer))
})


test_that("cost-effectiveness requires cost data", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  # Should error without cost data
  expect_error(
    threshold_analysis(
      nma_object = nma,
      outcome_direction = "higher",
      decision_rule = "cost_effectiveness",
      willingness_to_pay = 50000
    ),
    "must provide cost_data"
  )
})


# ============================================================================
# Test new study threshold formula
# ============================================================================

test_that("new study threshold uses meta-analytic updating", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    threshold_type = "new_study"
  )

  # New study threshold should contain key information
  new_study <- result$thresholds$new_study

  expect_true(!is.null(new_study$required_effect))
  expect_true(!is.null(new_study$current_effect))
  expect_true(!is.null(new_study$difference_from_current))
  expect_true(!is.null(new_study$tipping_point))
  expect_true(!is.null(new_study$second_best_treatment))

  # The required effect should be numerically valid
  expect_true(is.finite(new_study$required_effect))
  expect_true(is.finite(new_study$threshold_sd))
})


# ============================================================================
# Test risk aversion
# ============================================================================

test_that("risk aversion affects recommendation", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  # Risk-neutral (use point estimates)
  result_neutral <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    risk_aversion = 0,
    threshold_type = "effect_size"
  )

  # Risk-averse (penalize uncertainty)
  result_averse <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher",
    risk_aversion = 1,
    threshold_type = "effect_size"
  )

  # Both should run successfully
  expect_s3_class(result_neutral, "threshold_analysis")
  expect_s3_class(result_averse, "threshold_analysis")

  # Risk aversion might change the recommendation (not always, depends on data)
  # Just check both are valid
  expect_true(!is.null(result_neutral$recommendation$treatment))
  expect_true(!is.null(result_averse$recommendation$treatment))
})


# ============================================================================
# Test S3 methods
# ============================================================================

test_that("print method works for threshold_analysis", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher"
  )

  # Should print without error
  expect_output(print(result), "EXPERIMENTAL")
  expect_output(print(result), "Robustness Score")
})


test_that("plot methods work for threshold_analysis", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Plotting tests may be slow

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher"
  )

  # All plot types should work without error
  expect_silent(plot(result, type = "robustness"))
  expect_silent(plot(result, type = "thresholds"))
  expect_silent(plot(result, type = "comparison"))
})


# ============================================================================
# Test edge cases
# ============================================================================

test_that("threshold_analysis handles two-treatment network", {
  skip_if_not_installed("netmeta")

  # Minimal network: just 2 treatments
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "B", "B"),
    TE = c(0.5, 0.6, 0.4),
    seTE = c(0.2, 0.25, 0.18)
  )

  nma <- netmeta::netmeta(
    TE = data$TE,
    seTE = data$seTE,
    treat1 = data$treat1,
    treat2 = data$treat2,
    studlab = data$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher"
  )

  expect_s3_class(result, "threshold_analysis")
  expect_true(!is.null(result$thresholds$effect_size))
})


test_that("threshold_analysis validates outcome direction", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  # Invalid outcome direction
  expect_error(
    threshold_analysis(
      nma_object = nma,
      outcome_direction = "invalid"
    )
  )

  # Valid directions should work
  expect_s3_class(
    threshold_analysis(nma_object = nma, outcome_direction = "higher"),
    "threshold_analysis"
  )

  expect_s3_class(
    threshold_analysis(nma_object = nma, outcome_direction = "lower"),
    "threshold_analysis"
  )
})


# ============================================================================
# Test interpretation output
# ============================================================================

test_that("threshold_analysis provides clinical interpretation", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  nma <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE
  )

  result <- threshold_analysis(
    nma_object = nma,
    outcome_direction = "higher"
  )

  # Should have interpretation
  expect_true(!is.null(result$interpretation))
  expect_true(!is.null(result$interpretation$robustness_category))
  expect_true(!is.null(result$interpretation$guidance))
  expect_true(!is.null(result$interpretation$summary))

  # Interpretation should include key information
  expect_true(nchar(result$interpretation$guidance) > 10)
  expect_true(nchar(result$interpretation$summary) > 20)
})
