# ============================================================================
# Tests for Automatic Pathways
# ============================================================================

library(testthat)
library(powerNMA)

# ============================================================================
# Test auto_standard_nma()
# ============================================================================

test_that("auto_standard_nma runs on pairwise data", {
  skip_if_not_installed("netmeta")

  # Use netmeta example dataset
  data(Senn2013, package = "netmeta")

  result <- auto_standard_nma(
    data = Senn2013,
    verbose = FALSE
  )

  expect_s3_class(result, "auto_standard_nma")
  expect_equal(result$primary_analysis$status, "completed")
  expect_true(!is.null(result$primary_analysis$model_object))
})


test_that("auto_standard_nma matches netmeta results numerically", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  # Run automatic pathway
  auto_result <- auto_standard_nma(
    data = Senn2013,
    verbose = FALSE
  )

  # Run manual netmeta
  manual_result <- netmeta::netmeta(
    TE = Senn2013$TE,
    seTE = Senn2013$seTE,
    treat1 = Senn2013$treat1,
    treat2 = Senn2013$treat2,
    studlab = Senn2013$studlab,
    sm = "MD",
    random = TRUE,
    fixed = FALSE
  )

  # Treatment effects should match within tolerance
  expect_equal(
    auto_result$primary_analysis$model_object$TE.random,
    manual_result$TE.random,
    tolerance = 0.001
  )

  # Heterogeneity should match
  expect_equal(
    auto_result$diagnostics$tau2,
    manual_result$tau^2,
    tolerance = 0.001
  )
})


test_that("auto_standard_nma provides all expected components", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Check all required components exist
  expect_true(!is.null(result$data_characteristics))
  expect_true(!is.null(result$automatic_choices))
  expect_true(!is.null(result$primary_analysis))
  expect_true(!is.null(result$sensitivity_analyses))
  expect_true(!is.null(result$diagnostics))
  expect_true(!is.null(result$report))
  expect_true(!is.null(result$recommendations))

  # Check data characteristics
  expect_true(result$data_characteristics$n_studies > 0)
  expect_true(result$data_characteristics$n_treatments > 0)

  # Check automatic choices
  expect_true(!is.null(result$automatic_choices$primary_method))
  expect_true(!is.null(result$automatic_choices$reference))
})


test_that("auto_standard_nma sensitivity analyses run", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should have at least fixed_vs_random
  expect_true(length(result$sensitivity_analyses) > 0)
  expect_true("fixed_vs_random" %in% names(result$sensitivity_analyses))

  # Fixed vs random should complete
  expect_equal(
    result$sensitivity_analyses$fixed_vs_random$status,
    "completed"
  )
})


test_that("auto_standard_nma inconsistency assessment works", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should have inconsistency assessment
  expect_true(!is.null(result$inconsistency))
  expect_true(result$inconsistency$status %in% c("completed", "consistent", "unclear"))
})


test_that("auto_standard_nma rankings are generated", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should have prediction/ranking
  expect_true(!is.null(result$prediction))
  expect_equal(result$prediction$status, "completed")
})


test_that("auto_standard_nma handles arm-based data", {
  skip_if_not_installed("netmeta")

  # Create arm-based binary data
  data <- data.frame(
    study = rep(1:5, each = 3),
    treatment = rep(c("A", "B", "C"), 5),
    events = c(10, 15, 12, 8, 14, 11, 9, 16, 13, 11, 17, 14, 10, 15, 12),
    n = rep(50, 15)
  )

  result <- auto_standard_nma(
    data = data,
    outcome_type = "binary",
    verbose = FALSE
  )

  expect_s3_class(result, "auto_standard_nma")
  expect_equal(result$primary_analysis$status, "completed")
})


# ============================================================================
# Test auto_experimental_nma()
# ============================================================================

test_that("auto_experimental_nma runs on pairwise data", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_experimental_nma(
    data = Senn2013,
    research_question = "decision_making",
    verbose = FALSE
  )

  expect_s3_class(result, "auto_experimental_nma")
  expect_true(!is.null(result$experimental_analyses))
})


test_that("auto_experimental_nma selects appropriate methods", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  # Test decision_making question
  result_dm <- auto_experimental_nma(
    data = Senn2013,
    research_question = "decision_making",
    verbose = FALSE
  )

  # Should include threshold analysis
  expect_true("threshold_analysis" %in% result_dm$automatic_choices$methods_to_run)

  # Test model_uncertainty question
  result_mu <- auto_experimental_nma(
    data = Senn2013,
    research_question = "model_uncertainty",
    verbose = FALSE
  )

  # Should include model averaging
  expect_true("model_averaging" %in% result_mu$automatic_choices$methods_to_run)
})


test_that("auto_experimental_nma provides standard comparison", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_experimental_nma(
    data = Senn2013,
    research_question = "decision_making",
    verbose = FALSE
  )

  # Should have standard comparison
  expect_true(!is.null(result$standard_comparison))
  expect_true(!is.null(result$standard_comparison$best_treatment))
})


test_that("auto_experimental_nma handles time-to-event data", {
  skip_on_cran()  # Time-to-event test may be slow

  # Create simulated survival data
  data <- data.frame(
    study = rep(1:3, each = 100),
    treatment = rep(c("A", "B", "C"), each = 100),
    time = rexp(300, rate = 0.05),
    event = rbinom(300, 1, 0.7)
  )

  result <- auto_experimental_nma(
    data = data,
    outcome_type = "time_to_event",
    research_question = "survival_analysis",
    verbose = FALSE
  )

  expect_s3_class(result, "auto_experimental_nma")

  # Should select RMST method
  expect_true("rmst_nma" %in% result$automatic_choices$methods_to_run)
})


# ============================================================================
# Test S3 methods
# ============================================================================

test_that("print method works for auto_standard_nma", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")
  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should print without error
  expect_output(print(result), "AUTOMATIC STANDARD NMA")
})


test_that("summary method works for auto_standard_nma", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")
  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should summarize without error
  expect_output(summary(result), "Detailed Summary")
})


test_that("print method works for auto_experimental_nma", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")
  result <- auto_experimental_nma(
    data = Senn2013,
    research_question = "decision_making",
    verbose = FALSE
  )

  # Should print without error
  expect_output(print(result), "EXPERIMENTAL")
})


# ============================================================================
# Test edge cases
# ============================================================================

test_that("auto_standard_nma handles small networks", {
  skip_if_not_installed("netmeta")

  # Small network: 3 studies, 3 treatments
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "B"),
    treat2 = c("B", "C", "C"),
    TE = c(0.5, 0.8, 0.3),
    seTE = c(0.2, 0.25, 0.18)
  )

  result <- auto_standard_nma(data = data, verbose = FALSE)

  expect_s3_class(result, "auto_standard_nma")
  expect_equal(result$primary_analysis$status, "completed")
})


test_that("auto_standard_nma detects missing data", {
  skip_if_not_installed("netmeta")

  # Data with missing values
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "B", "A"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.5, NA, 0.3, 0.7),
    seTE = c(0.2, 0.25, 0.18, 0.22)
  )

  result <- auto_standard_nma(data = data, verbose = FALSE)

  # Should detect missing data
  expect_true(result$data_characteristics$missing_pct > 0)
})


test_that("auto_standard_nma gracefully handles errors", {
  # Malformed data (negative SE)
  data <- data.frame(
    studlab = c("S1"),
    treat1 = c("A"),
    treat2 = c("B"),
    TE = c(0.5),
    seTE = c(-0.1)  # Invalid negative SE
  )

  # Should fail gracefully
  expect_error(auto_standard_nma(data = data, verbose = FALSE))
})


# ============================================================================
# Test recommendations
# ============================================================================

test_that("auto_standard_nma provides interpretable recommendations", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")
  result <- auto_standard_nma(data = Senn2013, verbose = FALSE)

  # Should have recommendations
  expect_true(!is.null(result$recommendations))
  expect_true(!is.null(result$recommendations$best_treatment))
  expect_true(!is.null(result$recommendations$confidence))
  expect_true(!is.null(result$recommendations$certainty))
  expect_true(!is.null(result$recommendations$clinical_interpretation))
})


test_that("auto_experimental_nma provides robustness assessment", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  result <- auto_experimental_nma(
    data = Senn2013,
    research_question = "decision_making",
    verbose = FALSE
  )

  # Should have threshold analysis results if it ran
  if (!is.null(result$threshold_analysis) &&
      result$threshold_analysis$status == "completed") {
    expect_true(!is.null(result$recommendations$robustness))
  }
})
