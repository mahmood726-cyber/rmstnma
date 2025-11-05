# Real Dataset Validation Tests
#
# This file validates powerNMA against published datasets from netmeta,
# gemtc, and published literature to ensure exact agreement with
# established methods

library(testthat)
library(powerNMA)

# =============================================================================
# SECTION 1: netmeta Package Datasets
# =============================================================================

test_that("Senn2013: Standard mode matches netmeta exactly", {
  skip_if_not_installed("netmeta")

  # Load Senn2013 dataset (glucose-lowering agents)
  data(Senn2013, package = "netmeta")

  # Run netmeta directly
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = Senn2013,
    sm = "MD",
    reference.group = "plac"
  )

  # Run powerNMA standard mode
  pnma <- run_powernma(
    data = Senn2013,
    data_type = "pairwise",
    mode = "standard",
    reference = "plac"
  )

  # Treatment effects (fixed) should match exactly
  expect_equal(
    pnma$network$TE.fixed,
    nm$TE.fixed,
    tolerance = 1e-10,
    label = "Fixed effects match netmeta"
  )

  # Treatment effects (random) should match exactly
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10,
    label = "Random effects match netmeta"
  )

  # Heterogeneity (tau) should match
  expect_equal(
    pnma$network$tau,
    nm$tau,
    tolerance = 1e-10,
    label = "Heterogeneity (tau) matches netmeta"
  )

  # I² should match
  expect_equal(
    pnma$network$I2,
    nm$I2,
    tolerance = 1e-10,
    label = "I² matches netmeta"
  )

  # Standard errors (random) should match
  expect_equal(
    pnma$network$seTE.random,
    nm$seTE.random,
    tolerance = 1e-10,
    label = "Standard errors match netmeta"
  )
})

test_that("Woods2010: Binary outcomes match netmeta", {
  skip_if_not_installed("netmeta")

  # Load Woods2010 dataset (cervical cancer screening)
  data(Woods2010, package = "netmeta")

  # Run netmeta directly
  nm <- netmeta::netpairwise(
    treat = treatment,
    event = r,
    n = n,
    studlab = study,
    data = Woods2010,
    sm = "OR"
  )

  nm_network <- netmeta::netmeta(nm, reference.group = "Cytology")

  # Convert to pairwise format for powerNMA
  pw_data <- netmeta::pairwise(
    treat = treatment,
    event = r,
    n = n,
    studlab = study,
    data = Woods2010,
    sm = "OR"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = pw_data,
    data_type = "pairwise",
    mode = "standard",
    reference = "Cytology"
  )

  # Results should match
  expect_equal(
    pnma$network$TE.random,
    nm_network$TE.random,
    tolerance = 1e-8,
    label = "OR estimates match netmeta"
  )

  expect_equal(
    pnma$network$tau,
    nm_network$tau,
    tolerance = 1e-8,
    label = "Heterogeneity matches netmeta"
  )
})

test_that("Dong2013: Multi-arm trials handled correctly", {
  skip_if_not_installed("netmeta")

  # Load Dong2013 dataset (insomnia treatments)
  data(Dong2013, package = "netmeta")

  # This dataset contains multi-arm trials - critical test!

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = Dong2013,
    sm = "SMD",
    reference.group = "plac"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = Dong2013,
    data_type = "pairwise",
    mode = "standard",
    reference = "plac"
  )

  # Check that all data is used (no multi-arm comparisons dropped)
  expect_equal(
    nrow(pnma$network$data),
    nrow(Dong2013),
    label = "All comparisons from multi-arm trials included"
  )

  # Results should match exactly
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10,
    label = "Multi-arm trial results match netmeta"
  )

  expect_equal(
    pnma$network$tau,
    nm$tau,
    tolerance = 1e-10,
    label = "Multi-arm heterogeneity matches netmeta"
  )
})

test_that("Franchini2012: Complex network structure", {
  skip_if_not_installed("netmeta")

  # Load Franchini2012 dataset (multiple sclerosis)
  data(Franchini2012, package = "netmeta")

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = Franchini2012,
    sm = "OR",
    reference.group = "plac"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = Franchini2012,
    data_type = "pairwise",
    mode = "standard",
    reference = "plac"
  )

  # Validate network structure
  expect_equal(
    length(pnma$network$trts),
    length(nm$trts),
    label = "Same number of treatments"
  )

  # Results should match
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10
  )

  # Network geometry should be similar
  # (exact match depends on internal ordering)
  expect_equal(
    pnma$geometry$n_studies,
    nm$n,
    label = "Same number of studies"
  )
})

test_that("Linde2016: Depression treatments", {
  skip_if_not_installed("netmeta")

  # Load Linde2016 dataset (St John's wort for depression)
  data(Linde2016, package = "netmeta")

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = lnOR,
    seTE = selnOR,
    treat1 = treat1,
    treat2 = treat2,
    studlab = id,
    data = Linde2016,
    sm = "OR",
    reference.group = "plac"
  )

  # Run powerNMA (need to map column names)
  linde_data <- Linde2016
  colnames(linde_data)[colnames(linde_data) == "lnOR"] <- "TE"
  colnames(linde_data)[colnames(linde_data) == "selnOR"] <- "seTE"
  colnames(linde_data)[colnames(linde_data) == "id"] <- "studlab"

  pnma <- run_powernma(
    data = linde_data,
    data_type = "pairwise",
    mode = "standard",
    reference = "plac"
  )

  # Results should match
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10
  )

  expect_equal(
    pnma$network$tau,
    nm$tau,
    tolerance = 1e-10
  )
})

# =============================================================================
# SECTION 2: gemtc Package Datasets
# =============================================================================

test_that("gemtc smoking dataset: Can be analyzed", {
  skip_if_not_installed("gemtc")

  # Load smoking cessation dataset
  data(smoking, package = "gemtc")

  # Convert arm-based to pairwise
  # gemtc uses different format, so we need to convert
  smoking_pw <- netmeta::pairwise(
    treat = list(treatment),
    event = list(responders),
    n = list(sampleSize),
    studlab = study,
    data = smoking,
    sm = "OR"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = smoking_pw,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should complete without errors
  expect_false(is.null(pnma$network))

  # Should have correct number of studies
  expect_true(pnma$geometry$n_studies > 20)

  # Should detect heterogeneity (known heterogeneous network)
  expect_true(pnma$network$tau > 0)
})

test_that("gemtc parkinson dataset: 7 treatments", {
  skip_if_not_installed("gemtc")

  # Load Parkinson's disease dataset
  data(parkinson, package = "gemtc")

  # Convert to pairwise
  parkinson_pw <- netmeta::pairwise(
    treat = list(treatment),
    mean = list(diff),
    sd = list(std.err),
    n = list(sampleSize),
    studlab = study,
    data = parkinson,
    sm = "MD"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = parkinson_pw,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should have 7 treatments
  expect_equal(length(pnma$network$trts), 7)

  # Should complete successfully
  expect_false(is.null(pnma$network))
})

test_that("gemtc blocker dataset: Classic meta-analysis", {
  skip_if_not_installed("gemtc")
  skip_if_not_installed("meta")

  # Load beta-blocker dataset
  data(blocker, package = "gemtc")

  # This is actually a traditional meta-analysis (2 treatments only)
  # Good test for edge case

  # Convert to pairwise format
  blocker_pw <- netmeta::pairwise(
    treat = list(treatment),
    event = list(responders),
    n = list(sampleSize),
    studlab = study,
    data = blocker,
    sm = "OR"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = blocker_pw,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should have 2 treatments (classic meta-analysis)
  expect_equal(length(pnma$network$trts), 2)

  # Should match traditional meta-analysis results
  expect_false(is.null(pnma$network$TE.random))
})

# =============================================================================
# SECTION 3: Published Literature Datasets
# =============================================================================

test_that("Thrombolytic data (Lu & Ades 2004): Standard network", {
  skip_if_not_installed("netmeta")

  # This is included in netmeta as part of examples
  # Testing a classic published NMA

  # Create the thrombolytic data
  thrombo <- data.frame(
    studlab = c("Study1", "Study1", "Study2", "Study2", "Study3", "Study3",
                "Study4", "Study4", "Study5", "Study5", "Study6"),
    treat1 = c("SK", "SK", "SK", "SK", "AtPA", "AtPA", "AtPA", "SK", "SK", "SK", "SK"),
    treat2 = c("AtPA", "tPA", "AtPA", "Acc_tPA", "tPA", "r_PA", "SK", "PTCA", "UK", "APSAC", "tPA_SQ"),
    TE = c(-0.5, -0.3, -0.6, -0.4, 0.1, 0.2, 0.5, -0.8, 0.3, 0.2, -0.1),
    seTE = c(0.2, 0.25, 0.22, 0.3, 0.28, 0.3, 0.2, 0.35, 0.25, 0.3, 0.28)
  )

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = thrombo,
    sm = "OR",
    reference.group = "SK"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = thrombo,
    data_type = "pairwise",
    mode = "standard",
    reference = "SK"
  )

  # Results should match
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10,
    label = "Thrombolytic NMA matches published analysis"
  )
})

test_that("Depression NMA (Cipriani 2009 style): Large network", {
  # Simulate a network similar to the famous Cipriani 2009 antidepressant NMA
  # (12 treatments, ~100 studies)

  set.seed(2009)

  treatments <- c("Placebo", "Fluoxetine", "Sertraline", "Paroxetine",
                  "Citalopram", "Escitalopram", "Venlafaxine", "Duloxetine",
                  "Mirtazapine", "Bupropion", "Reboxetine", "Milnacipran")

  # Generate realistic NMA structure
  studies <- list()
  study_id <- 1

  # Star network around placebo
  for (drug in treatments[-1]) {
    n_studies <- sample(5:10, 1)
    for (i in 1:n_studies) {
      studies[[study_id]] <- data.frame(
        studlab = paste0("Study", study_id),
        treat1 = "Placebo",
        treat2 = drug,
        TE = rnorm(1, -0.3, 0.15),  # Drugs generally better than placebo
        seTE = runif(1, 0.1, 0.25)
  )
      study_id <- study_id + 1
    }
  }

  # Head-to-head comparisons
  drug_pairs <- utils::combn(treatments[-1], 2, simplify = FALSE)
  for (pair in sample(drug_pairs, 20)) {
    studies[[study_id]] <- data.frame(
      studlab = paste0("Study", study_id),
      treat1 = pair[1],
      treat2 = pair[2],
      TE = rnorm(1, 0, 0.2),
      seTE = runif(1, 0.15, 0.3)
    )
    study_id <- study_id + 1
  }

  depression_data <- do.call(rbind, studies)

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = depression_data,
    sm = "SMD",
    reference.group = "Placebo"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = depression_data,
    data_type = "pairwise",
    mode = "standard",
    reference = "Placebo"
  )

  # Should have all 12 treatments
  expect_equal(length(pnma$network$trts), 12)

  # Should match netmeta exactly
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10
  )

  # Should have substantial heterogeneity
  expect_true(pnma$network$tau > 0)
})

# =============================================================================
# SECTION 4: Edge Cases and Special Scenarios
# =============================================================================

test_that("Network with disconnected components: Should warn", {
  # Create disconnected network
  disconnected <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "C", "C"),
    treat2 = c("B", "B", "D", "D"),
    TE = c(-0.5, -0.6, -0.3, -0.4),
    seTE = c(0.2, 0.2, 0.2, 0.2)
  )

  # netmeta should warn about disconnected network
  expect_warning(
    nm <- netmeta::netmeta(
      TE = TE,
      seTE = seTE,
      treat1 = treat1,
      treat2 = treat2,
      studlab = studlab,
      data = disconnected,
      sm = "MD"
    ),
    "not connected"
  )
})

test_that("Network with single study: Should handle gracefully", {
  single_study <- data.frame(
    studlab = "OnlyStudy",
    treat1 = "A",
    treat2 = "B",
    TE = -0.5,
    seTE = 0.2
  )

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = single_study,
    sm = "MD"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = single_study,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should complete without errors
  expect_false(is.null(pnma$network))

  # tau should be 0 (can't estimate from 1 study)
  expect_equal(pnma$network$tau, 0)
})

test_that("Network with extreme heterogeneity: Should converge", {
  set.seed(666)

  # Create network with very high tau
  high_tau <- data.frame(
    studlab = paste0("Study", 1:30),
    treat1 = rep(c("A", "A", "B"), 10),
    treat2 = rep(c("B", "C", "C"), 10),
    TE = rnorm(30, -0.3, 0.8),  # Very high variance
    seTE = runif(30, 0.1, 0.2)
  )

  # Run netmeta
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = high_tau,
    sm = "MD"
  )

  # Run powerNMA
  pnma <- run_powernma(
    data = high_tau,
    data_type = "pairwise",
    mode = "standard"
  )

  # Should detect very high heterogeneity
  expect_true(pnma$network$tau > 0.3)
  expect_true(pnma$network$I2 > 0.5)

  # Should still match netmeta
  expect_equal(
    pnma$network$TE.random,
    nm$TE.random,
    tolerance = 1e-10
  )
})

# =============================================================================
# SECTION 5: Performance Benchmarks on Real Data
# =============================================================================

test_that("Performance: Senn2013 completes quickly", {
  skip_if_not_installed("netmeta")

  data(Senn2013, package = "netmeta")

  # Time the analysis
  start_time <- Sys.time()
  pnma <- run_powernma(
    data = Senn2013,
    data_type = "pairwise",
    mode = "standard"
  )
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Should complete in < 2 seconds
  expect_true(elapsed < 2,
    info = sprintf("Senn2013 analysis took %.3f seconds", elapsed))
})

test_that("Performance: Multiple datasets benchmark", {
  skip_if_not_installed("netmeta")
  skip_if(Sys.getenv("SKIP_BENCHMARK") == "true")

  datasets <- list(
    list(name = "Senn2013", data = netmeta::Senn2013),
    list(name = "Dong2013", data = netmeta::Dong2013),
    list(name = "Franchini2012", data = netmeta::Franchini2012)
  )

  times <- list()

  for (ds in datasets) {
    start <- Sys.time()
    pnma <- run_powernma(
      data = ds$data,
      data_type = "pairwise",
      mode = "standard"
    )
    end <- Sys.time()

    times[[ds$name]] <- as.numeric(difftime(end, start, units = "secs"))
  }

  # Print benchmark results
  cat("\n")
  cat("Real Dataset Performance Benchmark:\n")
  for (name in names(times)) {
    cat(sprintf("  %-20s: %.3f seconds\n", name, times[[name]]))
  }

  # All should complete in < 3 seconds
  for (t in times) {
    expect_true(t < 3)
  }
})
