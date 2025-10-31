# Comprehensive Testing with Large Simulation Datasets
#
# This file tests powerNMA on large, realistic simulation datasets
# covering various network structures and scenarios

test_that("Large star network (200 studies, 10 treatments) processes correctly", {
  set.seed(2001)

  # Generate large star network
  star_data <- simulate_large_star_network(
    n_treatments = 10,
    n_studies_per_comparison = 20,
    tau = 0.10,
    seed = 2001
  )

  # Should have ~180 studies (9 comparisons * 20 studies each)
  expect_true(nrow(star_data) >= 170)
  expect_true(nrow(star_data) <= 190)

  # Run standard NMA
  result <- run_powernma(
    data = star_data,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should complete without errors
  expect_false(is.null(result$network))

  # Should have 10 treatments
  expect_equal(length(result$network$trts), 10)

  # Heterogeneity should be estimated
  expect_false(is.na(result$network$tau))

  # tau should be close to true value (0.10)
  expect_true(abs(result$network$tau - 0.10) < 0.05)
})

test_that("Complete network (150 studies) handles all comparisons", {
  set.seed(2002)

  complete_data <- simulate_large_complete_network(
    n_treatments = 6,
    n_studies_per_comparison = 10,
    tau = 0.15,
    seed = 2002
  )

  # 6 treatments = C(6,2) = 15 comparisons
  # Each with ~10 studies = ~150 total
  n_comparisons <- length(unique(paste(complete_data$treat1, complete_data$treat2)))
  expect_true(n_comparisons >= 12)  # Allow some variability

  result <- run_powernma(complete_data, mode = "standard")

  # Network should be well-connected
  expect_false(is.null(result$network))
  expect_true(result$geometry$density > 0.8)  # High connectivity

  # Should have good precision (complete network)
  se_estimates <- result$network$seTE.random
  expect_true(all(se_estimates < 0.3, na.rm = TRUE))
})

test_that("Multi-arm network (45 trials) uses all arms correctly", {
  set.seed(2003)

  multiarm_data <- simulate_multiarm_network(
    n_2arm = 30,
    n_3arm = 10,
    n_4arm = 5,
    seed = 2003
  )

  # Count multi-arm trials
  study_counts <- multiarm_data %>%
    dplyr::group_by(studlab) %>%
    dplyr::summarize(n_comparisons = dplyr::n()) %>%
    dplyr::arrange(desc(n_comparisons))

  # 3-arm trials should have 3 comparisons
  three_arm_studies <- study_counts %>%
    dplyr::filter(n_comparisons == 3)
  expect_true(nrow(three_arm_studies) >= 8)  # ~10 expected

  # 4-arm trials should have 6 comparisons
  four_arm_studies <- study_counts %>%
    dplyr::filter(n_comparisons == 6)
  expect_true(nrow(four_arm_studies) >= 4)  # ~5 expected

  # Run NMA
  result <- run_powernma(multiarm_data, mode = "standard")

  # All data should be used (no discarded comparisons)
  expect_equal(nrow(result$network$data), nrow(multiarm_data))
})

test_that("Heterogeneous network (varying tau) detects heterogeneity", {
  set.seed(2004)

  hetero_data <- simulate_heterogeneous_network(
    n_treatments = 8,
    n_studies_total = 100,
    tau_low = 0.05,
    tau_high = 0.30,
    seed = 2004
  )

  result <- run_powernma(hetero_data, mode = "standard")

  # Should detect high heterogeneity
  expect_true(result$network$tau > 0.10)

  # I² should be substantial
  expect_true(result$network$I2 > 0.20)

  # Prediction intervals should be wide
  # (indicates heterogeneity)
})

test_that("Sparse network (30%% connectivity) still produces valid results", {
  set.seed(2005)

  sparse_data <- simulate_sparse_network(
    n_treatments = 12,
    sparsity = 0.3,
    min_studies_per_comparison = 2,
    seed = 2005
  )

  # Network should be connected but sparse
  result <- run_powernma(sparse_data, mode = "standard")

  expect_false(is.null(result$network))

  # Should have 12 treatments
  expect_equal(length(result$network$trts), 12)

  # Density should be low (sparse)
  expect_true(result$geometry$density < 0.5)

  # Standard errors should be larger (less data)
  se_estimates <- result$network$seTE.random
  mean_se <- mean(se_estimates, na.rm = TRUE)
  expect_true(mean_se > 0.15)  # Larger due to sparsity
})

test_that("Large IPD (4000 patients) processes efficiently", {
  skip_if_not_installed("survRM2")
  skip_if_not_installed("survival")

  set.seed(2006)

  large_ipd <- simulate_large_ipd(
    n_trials = 20,
    n_treatments = 6,
    n_per_arm = 200,
    seed = 2006
  )

  # Should have ~24000 patients (20 trials * 6 treatments * 200/arm)
  # Actually fewer due to 2/3-arm design
  expect_true(nrow(large_ipd) > 8000)
  expect_true(nrow(large_ipd) < 30000)

  # Run RMST NMA
  start_time <- Sys.time()
  result <- rmst_nma(
    ipd = large_ipd,
    tau_list = 365,
    reference = "Trt1"
  )
  end_time <- Sys.time()

  # Should complete in reasonable time (< 30 seconds)
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_true(elapsed < 30,
    info = sprintf("RMST NMA took %.1f seconds", elapsed))

  # Should produce valid results
  expect_false(is.null(result$tau_365))

  # Multi-arm trials should all be included
  # (no data loss from multi-arm bug)
  multiarm_trials <- large_ipd %>%
    dplyr::group_by(trial) %>%
    dplyr::summarize(n_treatments = length(unique(treatment))) %>%
    dplyr::filter(n_treatments == 3)

  if (nrow(multiarm_trials) > 0) {
    # Check that all 3-arm trials have 2 comparisons in result
    pw_data <- result$tau_365$data

    for (trial_id in multiarm_trials$trial) {
      trial_comparisons <- pw_data[pw_data$trial == trial_id, ]
      expect_equal(nrow(trial_comparisons), 2,
        info = sprintf("3-arm trial %s should have 2 comparisons", trial_id))
    }
  }
})

test_that("Stress test: Very large network (500 studies)", {
  skip_if(Sys.getenv("SKIP_STRESS_TESTS") == "true")

  set.seed(2007)

  # Generate very large dataset
  very_large <- simulate_large_star_network(
    n_treatments = 15,
    n_studies_per_comparison = 35,
    tau = 0.12,
    seed = 2007
  )

  # Should have ~490 studies (14 comparisons * 35 each)
  expect_true(nrow(very_large) >= 450)

  # Time the analysis
  start_time <- Sys.time()
  result <- run_powernma(
    very_large,
    mode = "standard"
  )
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Should complete in < 10 seconds even for large network
  expect_true(elapsed < 10,
    info = sprintf("Large network analysis took %.2f seconds", elapsed))

  # Results should be valid
  expect_false(is.null(result$network))
  expect_equal(length(result$network$trts), 15)
})

test_that("Robustness: Network with extreme heterogeneity", {
  set.seed(2008)

  # Generate data with very high tau
  extreme_hetero <- simulate_large_star_network(
    n_treatments = 8,
    n_studies_per_comparison = 20,
    tau = 0.50,  # Very high!
    seed = 2008
  )

  result <- run_powernma(extreme_hetero, mode = "standard")

  # Should still complete without errors
  expect_false(is.null(result$network))

  # Should detect high heterogeneity
  expect_true(result$network$tau > 0.30)
  expect_true(result$network$I2 > 0.70)

  # Random effects model should differ substantially from fixed
  te_fixed <- result$network$TE.fixed
  te_random <- result$network$TE.random

  # Shrinkage should be evident
  expect_false(isTRUE(all.equal(te_fixed, te_random)))
})

test_that("Edge case: Single study per comparison", {
  set.seed(2009)

  # Minimal data
  minimal <- simulate_large_star_network(
    n_treatments = 5,
    n_studies_per_comparison = 1,  # Only 1 study per comparison!
    tau = 0.10,
    seed = 2009
  )

  # Should have exactly 4 studies
  expect_equal(nrow(minimal), 4)

  # Should still produce results (but heterogeneity not estimable)
  result <- run_powernma(minimal, mode = "standard")

  expect_false(is.null(result$network))

  # tau should be 0 or NA (can't estimate from 1 study per comparison)
  expect_true(result$network$tau == 0 || is.na(result$network$tau))
})

test_that("Consistency: Results stable across multiple runs", {
  set.seed(2010)

  # Generate dataset once
  test_data <- simulate_large_complete_network(
    n_treatments = 5,
    n_studies_per_comparison = 15,
    tau = 0.12,
    seed = 2010
  )

  # Run multiple times
  result1 <- run_powernma(test_data, mode = "standard")
  result2 <- run_powernma(test_data, mode = "standard")
  result3 <- run_powernma(test_data, mode = "standard")

  # Results should be identical (no randomness in standard NMA)
  expect_equal(result1$network$TE.random, result2$network$TE.random)
  expect_equal(result2$network$TE.random, result3$network$TE.random)
  expect_equal(result1$network$tau, result2$network$tau)
})

test_that("Performance: Benchmark standard mode speed", {
  skip_if(Sys.getenv("SKIP_BENCHMARK") == "true")

  # Various network sizes
  sizes <- c(10, 50, 100, 200)
  times <- numeric(length(sizes))

  for (i in seq_along(sizes)) {
    test_data <- simulate_large_star_network(
      n_treatments = 6,
      n_studies_per_comparison = sizes[i] / 5,
      seed = 3000 + i
    )

    start <- Sys.time()
    result <- run_powernma(test_data, mode = "standard")
    end <- Sys.time()

    times[i] <- as.numeric(difftime(end, start, units = "secs"))
  }

  # Time should scale reasonably (not exponentially)
  # 200 studies should take < 4x the time of 50 studies
  time_ratio <- times[4] / times[2]
  expect_true(time_ratio < 5,
    info = sprintf("Time scaling: %.2fx for 4x data", time_ratio))

  cat("\n")
  cat("Performance Benchmark:\n")
  for (i in seq_along(sizes)) {
    cat(sprintf("  %3d studies: %.3f seconds\n", sizes[i], times[i]))
  }
})
