# RMST Sign Convention Validation Tests
#
# MATHEMATICAL PROOF OF SIGN CONVENTION:
#
# survRM2::rmst2() calculates: RMST(arm=1) - RMST(arm=0)
# where arm=0 is the first group in the data, arm=1 is the second group
#
# powerNMA pairwise structure uses:
#   treat1 (typically reference/control)
#   treat2 (typically active treatment)
#
# For the comparison "Treatment vs Control":
#   - We want: RMST(Treatment) - RMST(Control)
#   - survRM2 coding: arm=0 for Control, arm=1 for Treatment
#   - survRM2 output: RMST(arm=1) - RMST(arm=0) = RMST(Treatment) - RMST(Control)
#
# Therefore, NO sign flip is needed when:
#   - treat1 = Control (arm=0)
#   - treat2 = Treatment (arm=1)
#   - We want TE = treat2 - treat1
#
# HOWEVER, the current code uses: TE = -1 * rmst_result
# This suggests the coding may be reversed. Let's validate.

test_that("RMST sign convention is correct: better treatment has positive RMST", {
  skip_if_not_installed("survRM2")

  set.seed(42)

  # Create IPD where Treatment is BETTER than Control
  # Lower hazard = longer survival = higher RMST
  ipd <- data.frame(
    trial = "ValidationTrial",
    treatment = rep(c("Control", "Treatment"), each = 200),
    time = c(
      rexp(200, rate = 0.05),  # Control: hazard = 0.05
      rexp(200, rate = 0.02)   # Treatment: hazard = 0.02 (BETTER)
    ),
    status = 1  # All events
  )

  # Manual RMST calculation at tau = 365
  tau <- 365

  control_data <- ipd[ipd$treatment == "Control", ]
  treatment_data <- ipd[ipd$treatment == "Treatment", ]

  km_control <- survival::survfit(survival::Surv(time, status) ~ 1, data = control_data)
  km_treatment <- survival::survfit(survival::Surv(time, status) ~ 1, data = treatment_data)

  # RMST is area under KM curve up to tau
  # Approximate using summary
  times <- seq(0, tau, by = 1)
  surv_control <- summary(km_control, times = times, extend = TRUE)$surv
  surv_treatment <- summary(km_treatment, times = times, extend = TRUE)$surv

  rmst_control_manual <- sum(surv_control)
  rmst_treatment_manual <- sum(surv_treatment)

  # Treatment effect (Treatment - Control)
  manual_TE <- rmst_treatment_manual - rmst_control_manual

  # Treatment is better, so RMST(Treatment) > RMST(Control)
  # Therefore: Treatment - Control should be POSITIVE
  expect_true(manual_TE > 0,
    info = sprintf("Manual TE = %.2f should be positive (Treatment better)",
                   manual_TE))

  # Now test package function
  result <- rmst_nma(ipd, tau_list = tau, reference = "Control")

  # Extract the treatment effect for Treatment vs Control
  nma_data <- result$tau_365$data
  te_row <- nma_data[nma_data$treat1 == "Treatment" & nma_data$treat2 == "Control", ]

  if (nrow(te_row) == 0) {
    # Try reverse
    te_row <- nma_data[nma_data$treat1 == "Control" & nma_data$treat2 == "Treatment", ]
  }

  package_TE <- te_row$TE[1]

  # The package TE should have the SAME SIGN as manual TE
  expect_true(sign(package_TE) == sign(manual_TE),
    info = sprintf(
      "Package TE = %.2f should have same sign as manual TE = %.2f\nBetter treatment should have positive RMST difference",
      package_TE, manual_TE
    ))

  # They should also be numerically close (within 5%)
  expect_equal(abs(package_TE), abs(manual_TE), tolerance = 0.05 * abs(manual_TE))
})

test_that("RMST sign convention validated with survRM2 directly", {
  skip_if_not_installed("survRM2")

  set.seed(123)

  # Simple test: Treatment has longer survival
  data <- data.frame(
    treatment = rep(c("Control", "Treatment"), each = 100),
    time = c(rexp(100, 0.1), rexp(100, 0.05)),  # Treatment: half the hazard
    status = 1
  )

  # Use survRM2 directly
  combined_data <- data
  combined_data$arm <- ifelse(combined_data$treatment == "Control", 0, 1)

  rmst_direct <- survRM2::rmst2(
    time = combined_data$time,
    status = combined_data$status,
    arm = combined_data$arm,
    tau = 100
  )

  # survRM2 output: arm=1 - arm=0 = Treatment - Control
  survRM2_estimate <- rmst_direct$unadjusted.result[1, "Est."]

  # Treatment has better survival, so Treatment - Control should be POSITIVE
  expect_true(survRM2_estimate > 0,
    info = sprintf("survRM2 estimate = %.3f should be positive", survRM2_estimate))

  # Now check what powerNMA does
  ipd <- data.frame(
    trial = "DirectTest",
    treatment = data$treatment,
    time = data$time,
    status = data$status
  )

  result <- rmst_nma(ipd, tau_list = 100, reference = "Control")

  # Extract TE from result
  pw_data <- result$tau_100$data
  te_package <- pw_data$TE[1]

  # Document the relationship
  cat("\n=== RMST Sign Convention Validation ===\n")
  cat(sprintf("survRM2 direct estimate: %.3f (Treatment - Control)\n", survRM2_estimate))
  cat(sprintf("powerNMA TE:             %.3f\n", te_package))
  cat(sprintf("Sign flip in code:       %s\n", ifelse(sign(te_package) != sign(survRM2_estimate), "YES", "NO")))

  # The key test: does the sign convention preserve the direction?
  # Better treatment should have positive TE in NMA context
  # (where positive = favors treatment over reference)
  expect_true(te_package > 0,
    info = "powerNMA should show positive TE when Treatment is better than Control")
})

test_that("RMST manual calculation matches survRM2", {
  skip_if_not_installed("survRM2")

  set.seed(789)

  # Create simple IPD
  ipd <- data.frame(
    trial = "ManualTest",
    treatment = rep(c("A", "B"), each = 50),
    time = c(rexp(50, 0.08), rexp(50, 0.04)),
    status = 1
  )

  tau <- 50

  # Method 1: survRM2
  combined <- ipd
  combined$arm <- ifelse(combined$treatment == "A", 0, 1)

  rmst_survRM2 <- survRM2::rmst2(
    time = combined$time,
    status = combined$status,
    arm = combined$arm,
    tau = tau
  )

  # Method 2: Manual integration under KM curve
  km_A <- survival::survfit(survival::Surv(time, status) ~ 1,
                            data = ipd[ipd$treatment == "A", ])
  km_B <- survival::survfit(survival::Surv(time, status) ~ 1,
                            data = ipd[ipd$treatment == "B", ])

  # RMST = integral of S(t) from 0 to tau
  # Approximate with sum
  t_grid <- seq(0, tau, length.out = 1000)
  surv_A <- summary(km_A, times = t_grid, extend = FALSE)$surv
  surv_B <- summary(km_B, times = t_grid, extend = FALSE)$surv

  # Handle NAs at end
  surv_A[is.na(surv_A)] <- tail(na.omit(surv_A), 1)
  surv_B[is.na(surv_B)] <- tail(na.omit(surv_B), 1)

  rmst_A <- sum(surv_A) * (tau / length(t_grid))
  rmst_B <- sum(surv_B) * (tau / length(t_grid))

  manual_diff <- rmst_B - rmst_A  # B - A (arm=1 - arm=0)
  survRM2_diff <- rmst_survRM2$unadjusted.result[1, "Est."]

  # These should be very close
  expect_equal(manual_diff, survRM2_diff, tolerance = 1,
    info = sprintf("Manual = %.2f, survRM2 = %.2f", manual_diff, survRM2_diff))
})

test_that("RMST NMA sign convention documented and correct", {
  skip_if_not_installed("survRM2")

  # DOCUMENTATION TEST
  # This test serves as documentation of the correct sign convention

  set.seed(999)

  # Scenario: Drug A is standard, Drug B is better
  ipd <- data.frame(
    trial = "SignTest",
    treatment = rep(c("DrugA", "DrugB"), each = 100),
    time = c(
      rexp(100, 0.06),  # DrugA: median survival = ln(2)/0.06 ≈ 11.6
      rexp(100, 0.03)   # DrugB: median survival = ln(2)/0.03 ≈ 23.1 (BETTER)
    ),
    status = 1
  )

  # Calculate medians to confirm
  median_A <- median(ipd$time[ipd$treatment == "DrugA"])
  median_B <- median(ipd$time[ipd$treatment == "DrugB"])

  cat("\n=== Sign Convention Documentation ===\n")
  cat(sprintf("DrugA median survival: %.2f\n", median_A))
  cat(sprintf("DrugB median survival: %.2f\n", median_B))
  cat("Expected: DrugB > DrugA (DrugB is better)\n\n")

  # Run RMST NMA with DrugA as reference
  result <- rmst_nma(ipd, tau_list = 50, reference = "DrugA")

  pw_data <- result$tau_50$data
  te <- pw_data$TE[1]

  cat(sprintf("RMST Treatment Effect (TE): %.3f\n", te))
  cat("Interpretation:\n")
  cat("  - treat1 = DrugA (reference)\n")
  cat("  - treat2 = DrugB (intervention)\n")
  cat("  - TE = RMST difference\n")
  cat("  - Positive TE = DrugB has longer RMST = DrugB is better\n")
  cat("  - Negative TE = DrugB has shorter RMST = DrugB is worse\n\n")

  # The test: TE should be POSITIVE (DrugB is better)
  expect_true(te > 0,
    label = sprintf("TE = %.3f should be positive (DrugB better than DrugA)", te))

  cat("CONCLUSION: Sign convention VALIDATED ✓\n")
  cat("=========================================\n\n")
})
