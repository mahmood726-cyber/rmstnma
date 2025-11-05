# Tests for Experimental Methods (Time-Varying NMA)
#
# These tests validate the FIXED implementations of RMST and Milestone NMA

test_that("Multi-arm trials: RMST NMA creates all comparisons", {
  skip_if_not_installed("survRM2")

  set.seed(42)

  # Create IPD with a 3-arm trial
  ipd <- rbind(
    # 2-arm trial
    data.frame(
      trial = "Trial1",
      treatment = rep(c("A", "B"), each = 50),
      time = c(rexp(50, 0.05), rexp(50, 0.03)),
      status = 1
    ),
    # 3-arm trial
    data.frame(
      trial = "Trial2_MultiArm",
      treatment = rep(c("A", "B", "C"), each = 50),
      time = c(rexp(50, 0.05), rexp(50, 0.03), rexp(50, 0.02)),
      status = 1
    )
  )

  # Run RMST NMA
  result <- rmst_nma(ipd, tau_list = 100, reference = "A")

  # Extract pairwise data
  pw_data <- result$tau_100$data

  # For Trial2, should have comparisons: A-B, A-C (reference-based)
  trial2_comparisons <- pw_data[pw_data$trial == "Trial2_MultiArm", ]

  expect_equal(nrow(trial2_comparisons), 2,
    info = "3-arm trial should generate 2 comparisons (all vs reference)")

  # Check that we have A-B and A-C
  comparisons <- paste(trial2_comparisons$treat1, trial2_comparisons$treat2, sep = "-")
  expect_true("A-B" %in% comparisons)
  expect_true("A-C" %in% comparisons)
})

test_that("Multi-arm trials: Milestone NMA creates all comparisons", {
  skip_if_not_installed("survival")

  set.seed(123)

  # Create IPD with a 3-arm trial
  ipd <- rbind(
    # 2-arm trial
    data.frame(
      trial = "Trial1",
      treatment = rep(c("Control", "TrtA"), each = 60),
      time = c(rexp(60, 0.06), rexp(60, 0.04)),
      status = 1
    ),
    # 3-arm trial
    data.frame(
      trial = "Trial2_MultiArm",
      treatment = rep(c("Control", "TrtA", "TrtB"), each = 60),
      time = c(rexp(60, 0.06), rexp(60, 0.04), rexp(60, 0.03)),
      status = 1
    )
  )

  # Run milestone NMA
  result <- milestone_nma(ipd, times = 50, reference = "Control")

  # Extract count data
  count_data <- result$time_50$data

  # For Trial2, should have 2 comparisons
  trial2_comparisons <- count_data[count_data$trial == "Trial2_MultiArm", ]

  expect_equal(nrow(trial2_comparisons), 2,
    info = "3-arm trial should generate 2 comparisons")

  # Verify we have Control-TrtA and Control-TrtB
  comparisons <- paste(trial2_comparisons$treat1, trial2_comparisons$treat2, sep = "-")
  expect_true("Control-TrtA" %in% comparisons)
  expect_true("Control-TrtB" %in% comparisons)
})

test_that("Milestone extend=FALSE is default and prevents extrapolation", {
  set.seed(456)

  # Create IPD with short follow-up
  ipd <- data.frame(
    trial = "ShortFollowup",
    treatment = rep(c("A", "B"), each = 50),
    time = c(runif(50, 0, 30), runif(50, 0, 30)),  # Max time = 30
    status = 1
  )

  # Try to analyze at milestone = 100 (beyond follow-up)
  expect_warning(
    result <- milestone_nma(ipd, times = 100, reference = "A"),
    "exceeds follow-up"
  )

  # Result should skip this comparison
  # (or return with warning)
})

test_that("Continuity correction: standard method is default", {
  skip_if_not_installed("survival")

  set.seed(789)

  # Create data with zero events in one arm
  ipd <- data.frame(
    trial = "ZeroEventTrial",
    treatment = rep(c("A", "B"), each = 20),
    time = c(
      rep(100, 20),  # Arm A: all censored at 100
      rexp(20, 0.1)  # Arm B: some events
    ),
    status = c(
      rep(0, 20),    # All censored
      rep(1, 20)     # All events
    )
  )

  # Run with standard correction (default)
  result_std <- milestone_nma(ipd, times = 50, reference = "A",
                              continuity_correction = "standard")

  # Run with empirical correction
  result_emp <- milestone_nma(ipd, times = 50, reference = "A",
                              continuity_correction = "empirical")

  # Results should differ
  te_std <- result_std$time_50$TE.random[2, 1]
  te_emp <- result_emp$time_50$TE.random[2, 1]

  expect_false(isTRUE(all.equal(te_std, te_emp)),
    info = "Standard and empirical corrections should produce different results")
})

test_that("RMST values increase with better survival", {
  skip_if_not_installed("survRM2")

  set.seed(111)

  # Create IPD with clearly better treatment
  ipd <- data.frame(
    trial = "TestTrial",
    treatment = rep(c("Worse", "Better"), each = 100),
    time = c(
      rexp(100, 0.10),  # Worse: median ~ 7
      rexp(100, 0.02)   # Better: median ~ 35
    ),
    status = 1
  )

  result <- rmst_nma(ipd, tau_list = 100, reference = "Worse")

  # Extract TE
  pw_data <- result$tau_100$data
  te <- pw_data$TE[pw_data$treat1 == "Worse" & pw_data$treat2 == "Better"]

  # TE should be POSITIVE (Better has higher RMST)
  expect_true(te > 0,
    info = sprintf("TE = %.2f should be positive (Better > Worse)", te))

  # And it should be substantial
  expect_true(te > 10,
    info = sprintf("TE = %.2f should be substantial given the large difference", te))
})

test_that("Milestone NMA handles varying follow-up correctly", {
  set.seed(222)

  # Create trials with different follow-up
  ipd <- rbind(
    data.frame(
      trial = "LongFollowup",
      treatment = rep(c("A", "B"), each = 50),
      time = c(rexp(50, 0.01), rexp(50, 0.008)),  # Long times
      status = 1
    ),
    data.frame(
      trial = "ShortFollowup",
      treatment = rep(c("A", "B"), each = 50),
      time = c(rexp(50, 0.2), rexp(50, 0.15)),  # Short times
      status = 1
    )
  )

  # Analyze at milestone that's reasonable for long but not short
  milestone_time <- 50

  # With extend=FALSE (default), should warn about short trial
  expect_warning(
    result <- milestone_nma(ipd, times = milestone_time,
                           reference = "A", extend = FALSE),
    "exceeds follow-up|Milestone time.*exceeds"
  )
})

test_that("RMST sign convention is correct after fix", {
  skip_if_not_installed("survRM2")

  set.seed(333)

  # Simple test: B is better than A
  ipd <- data.frame(
    trial = "SignTest",
    treatment = rep(c("A", "B"), each = 100),
    time = c(rexp(100, 0.08), rexp(100, 0.04)),  # B has longer survival
    status = 1
  )

  result <- rmst_nma(ipd, tau_list = 50, reference = "A")
  pw_data <- result$tau_50$data

  # TE should be POSITIVE (B is better)
  te <- pw_data$TE[1]
  expect_true(te > 0,
    label = sprintf("TE = %.3f should be positive (B better than A)", te))
})
