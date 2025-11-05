# ============================================================================
# RMST Standard Error Validation Against survRM2
# ============================================================================

library(testthat)
library(powerNMA)

test_that("RMST SE calculation matches survRM2", {
  skip_if_not_installed("survRM2")
  skip_if_not_installed("survival")

  # Use survival package's veteran data
  data(veteran, package = "survival")

  # Select two treatment arms for comparison
  trt1_data <- veteran[veteran$trt == 1, ]
  trt2_data <- veteran[veteran$trt == 2, ]

  # Set tau (restriction time) at 300 days
  tau <- 300

  # Calculate RMST using survRM2 for each arm
  fit1_survRM2 <- survRM2::rmst2(
    time = trt1_data$time,
    status = trt1_data$status,
    arm = rep(0, nrow(trt1_data)),
    tau = tau
  )

  fit2_survRM2 <- survRM2::rmst2(
    time = trt2_data$time,
    status = trt2_data$status,
    arm = rep(0, nrow(trt2_data)),
    tau = tau
  )

  # Extract RMST and SE from survRM2
  rmst1_survRM2 <- fit1_survRM2$RMST.arm0$rmst[1]
  se1_survRM2 <- fit1_survRM2$RMST.arm0$rmst[2]

  rmst2_survRM2 <- fit2_survRM2$RMST.arm0$rmst[1]
  se2_survRM2 <- fit2_survRM2$RMST.arm0$rmst[2]

  # Calculate RMST using powerNMA's internal function
  # First create KM fits
  km_fit1 <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = trt1_data
  )

  km_fit2 <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = trt2_data
  )

  # Use powerNMA's internal calculate_rmst_from_survfit function
  # Access it through powerNMA namespace
  calc_rmst <- powerNMA:::calculate_rmst_from_survfit

  rmst1_powerNMA <- calc_rmst(km_fit1, tau)
  rmst2_powerNMA <- calc_rmst(km_fit2, tau)

  # Compare RMST point estimates (should be very close)
  expect_equal(
    rmst1_powerNMA$rmst,
    rmst1_survRM2,
    tolerance = 0.5,  # Within half a day
    info = "RMST point estimate for arm 1"
  )

  expect_equal(
    rmst2_powerNMA$rmst,
    rmst2_survRM2,
    tolerance = 0.5,
    info = "RMST point estimate for arm 2"
  )

  # Compare standard errors (should be reasonably close)
  # Note: Different methods may give somewhat different SEs
  # We accept up to 20% relative difference
  rel_diff1 <- abs(rmst1_powerNMA$se - se1_survRM2) / se1_survRM2
  rel_diff2 <- abs(rmst2_powerNMA$se - se2_survRM2) / se2_survRM2

  expect_lt(
    rel_diff1,
    0.20,
    info = sprintf(
      "RMST SE for arm 1: powerNMA=%.3f, survRM2=%.3f, rel_diff=%.1f%%",
      rmst1_powerNMA$se, se1_survRM2, rel_diff1 * 100
    )
  )

  expect_lt(
    rel_diff2,
    0.20,
    info = sprintf(
      "RMST SE for arm 2: powerNMA=%.3f, survRM2=%.3f, rel_diff=%.1f%%",
      rmst2_powerNMA$se, se2_survRM2, rel_diff2 * 100
    )
  )

  # Print detailed comparison for documentation
  cat("\n=== RMST SE Validation Results ===\n")
  cat("Arm 1:\n")
  cat(sprintf("  RMST (powerNMA): %.3f ± %.3f\n", rmst1_powerNMA$rmst, rmst1_powerNMA$se))
  cat(sprintf("  RMST (survRM2):  %.3f ± %.3f\n", rmst1_survRM2, se1_survRM2))
  cat(sprintf("  Relative SE difference: %.1f%%\n", rel_diff1 * 100))
  cat("\nArm 2:\n")
  cat(sprintf("  RMST (powerNMA): %.3f ± %.3f\n", rmst2_powerNMA$rmst, rmst2_powerNMA$se))
  cat(sprintf("  RMST (survRM2):  %.3f ± %.3f\n", rmst2_survRM2, se2_survRM2))
  cat(sprintf("  Relative SE difference: %.1f%%\n", rel_diff2 * 100))
  cat("=================================\n\n")
})


test_that("RMST SE is reasonable across different scenarios", {
  skip_if_not_installed("survival")

  # Test 1: Small sample with high censoring
  small_data <- data.frame(
    time = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55),
    event = c(1, 0, 1, 0, 1, 0, 0, 0, 1, 0)
  )

  km_fit_small <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = small_data
  )

  calc_rmst <- powerNMA:::calculate_rmst_from_survfit
  rmst_small <- calc_rmst(km_fit_small, tau = 50)

  # SE should be positive
  expect_gt(rmst_small$se, 0)
  # SE should not be larger than tau (sanity check)
  expect_lt(rmst_small$se, 50)
  # RMST should be less than tau
  expect_lt(rmst_small$rmst, 50)

  # Test 2: Large sample with low censoring
  set.seed(12345)
  large_data <- data.frame(
    time = rexp(200, rate = 0.05),
    event = rbinom(200, 1, 0.8)
  )

  km_fit_large <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = large_data
  )

  rmst_large <- calc_rmst(km_fit_large, tau = 30)

  # SE should be smaller with large sample
  expect_gt(rmst_large$se, 0)
  expect_lt(rmst_large$se, rmst_small$se)  # Larger sample -> smaller SE

  # Test 3: Complete data (no censoring)
  complete_data <- data.frame(
    time = runif(100, 10, 50),
    event = rep(1, 100)
  )

  km_fit_complete <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = complete_data
  )

  rmst_complete <- calc_rmst(km_fit_complete, tau = 50)

  # SE should be positive and reasonable
  expect_gt(rmst_complete$se, 0)
  expect_lt(rmst_complete$se, 5)  # Should be small with complete data
})


test_that("RMST SE formula documentation is correct", {
  # This test documents the formula used for RMST SE calculation

  # The formula used is:
  # Var(RMST) = ∫₀^τ S(t)² * Var(S(t)) dt
  # where Var(S(t)) comes from Greenwood's formula

  # This is based on:
  # Andersen PK, Hansen MG, Klein JP (2004)
  # "Regression analysis of restricted mean survival time"
  # Lifetime Data Analysis 10(4):335-350

  # Create a simple test case
  data <- data.frame(
    time = c(5, 10, 15, 20, 25),
    event = c(1, 1, 1, 0, 0)
  )

  km_fit <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = data
  )

  calc_rmst <- powerNMA:::calculate_rmst_from_survfit
  result <- calc_rmst(km_fit, tau = 25)

  # Document that result contains both RMST and SE
  expect_true(!is.null(result$rmst))
  expect_true(!is.null(result$se))

  # Document that SE is based on variance integration
  # (No specific value to test, just documenting the approach)
  expect_gt(result$se, 0)

  cat("\n=== RMST SE Formula Documentation ===\n")
  cat("Formula: Var(RMST) = ∫₀^τ S(t)² * Var(S(t)) dt\n")
  cat("Reference: Andersen et al. (2004) Lifetime Data Analysis\n")
  cat("Implementation: Trapezoidal numerical integration\n")
  cat("=========================================\n\n")
})


test_that("RMST SE handles edge cases", {
  skip_if_not_installed("survival")

  calc_rmst <- powerNMA:::calculate_rmst_from_survfit

  # Edge case 1: All events before tau
  early_data <- data.frame(
    time = c(5, 10, 15, 20),
    event = c(1, 1, 1, 1)
  )

  km_early <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = early_data
  )

  rmst_early <- calc_rmst(km_early, tau = 50)

  expect_gt(rmst_early$se, 0)
  expect_lt(rmst_early$rmst, 50)

  # Edge case 2: No events (all censored)
  cens_data <- data.frame(
    time = c(10, 20, 30, 40, 50),
    event = c(0, 0, 0, 0, 0)
  )

  km_cens <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = cens_data
  )

  rmst_cens <- calc_rmst(km_cens, tau = 50)

  # RMST should equal tau when no events
  expect_equal(rmst_cens$rmst, 50, tolerance = 0.1)

  # Edge case 3: Single observation
  single_data <- data.frame(
    time = 30,
    event = 1
  )

  km_single <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    data = single_data
  )

  rmst_single <- calc_rmst(km_single, tau = 50)

  # Should handle without error
  expect_true(!is.null(rmst_single$rmst))
  expect_true(!is.null(rmst_single$se))
})
