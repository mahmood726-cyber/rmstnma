test_that("Data generation works", {
  # IPD generation using CI helper for faster tests
  ipd <- create_small_test_data("ipd")
  expect_true(validate_ipd(ipd))
  expect_equal(nrow(ipd), 2 * 20 * 2)

  # NMA data generation using CI helper for faster tests
  nma_data <- create_small_test_data("nma")
  expect_true(validate_nma_input(nma_data))
  expect_true("studlab" %in% names(nma_data))
  expect_true("TE" %in% names(nma_data))
})

test_that("Data validation catches errors", {
  bad_ipd <- data.frame(
    trial = c("T1", "T1"),
    treatment = c("A", "B"),
    time = c(10, -5),  # Negative time
    status = c(1, 1)
  )

  expect_warning(validate_ipd(bad_ipd))

  bad_nma <- data.frame(
    studlab = c("S1"),
    treat1 = c("A"),
    treat2 = c("B"),
    TE = c(0.5),
    seTE = c(-0.1)  # Negative SE
  )

  expect_error(validate_nma_input(bad_nma))
})

test_that("Configuration setup works", {
  config <- setup_powernma(sm = "OR", use_bayesian = FALSE)

  expect_equal(config$sm, "OR")
  expect_false(config$use_bayesian)
  expect_true("tau_list" %in% names(config))
})

test_that("Network geometry calculation works", {
  # Use small test data for faster CI execution
  data <- create_small_test_data("nma")
  geom <- network_geometry(data)

  expect_true("n_studies" %in% names(geom))
  expect_true("n_treatments" %in% names(geom))
  expect_gt(geom$n_studies, 0)
  expect_gt(geom$n_treatments, 0)
})

test_that("RMST NMA runs without error", {
  skip_if_not_installed("survRM2")

  # Use small test data for faster CI execution
  ipd <- create_small_test_data("ipd")
  results <- rmst_nma(ipd, tau_list = c(180, 365), reference = "Control")

  expect_s3_class(results, "rmst_nma")
  expect_true("tau_180" %in% names(results) || "tau_365" %in% names(results))
})

test_that("Milestone NMA runs correctly", {
  # Use small test data for faster CI execution
  ipd <- create_small_test_data("ipd")
  results <- milestone_nma(ipd, times = c(180), reference = "Control")

  expect_s3_class(results, "milestone_nma")
  if (length(results) > 0) {
    expect_true("day_180" %in% names(results))
  }
})

test_that("Standard NMA runs", {
  ci_safe_require("netmeta")  # Quiet package loading in CI

  # Use small test data for faster CI execution
  data <- create_small_test_data("nma")
  results <- run_powernma(
    data,
    data_type = "pairwise",
    config = setup_powernma(
      use_bayesian = FALSE,
      run_sensitivity = FALSE,
      export_results = FALSE
    )
  )

  expect_s3_class(results, "powernma_result")
  expect_true(!is.null(results$results$main_nma))
})
