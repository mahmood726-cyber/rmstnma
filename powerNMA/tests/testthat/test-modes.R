# Tests for Standard vs Experimental Mode Architecture

test_that("Standard mode rejects IPD data", {
  ipd <- generate_example_ipd(n_trials = 3, seed = 123)

  expect_error(
    run_powernma(ipd, data_type = "ipd", mode = "standard"),
    "IPD-based.*EXPERIMENTAL"
  )
})

test_that("Experimental mode shows warnings", {
  ipd <- generate_example_ipd(n_trials = 3, seed = 123)

  expect_warning(
    run_powernma(ipd, data_type = "ipd", mode = "experimental"),
    "EXPERIMENTAL MODE"
  )
})

test_that("Standard mode accepts pairwise data", {
  data <- simulate_nma_data(n_studies = 20, seed = 123)

  expect_silent({
    result <- run_powernma(data, data_type = "pairwise", mode = "standard")
  })

  result <- run_powernma(data, data_type = "pairwise", mode = "standard")

  expect_equal(result$mode, "standard")
  expect_true(result$validated)
  expect_false(is.null(result$network))
})

test_that("Experimental mode runs time-varying methods", {
  ipd <- generate_example_ipd(n_trials = 5, seed = 123)
  config <- setup_powernma(use_timevarying = TRUE, tau_list = 365, milestone_times = c(180, 365))

  expect_warning({
    result <- run_powernma(ipd, data_type = "ipd", mode = "experimental", config = config)
  }, "EXPERIMENTAL MODE")

  expect_equal(result$mode, "experimental")
  expect_false(result$validated)

  # Should have RMST and milestone results
  expect_false(is.null(result$rmst))
  expect_false(is.null(result$milestone))
})

test_that("Mode parameter defaults to standard", {
  data <- simulate_nma_data(n_studies = 10, seed = 123)

  # Not specifying mode should default to standard
  result <- run_powernma(data, data_type = "pairwise")

  expect_equal(result$mode, "standard")
})

test_that("Print method shows mode information", {
  data <- simulate_nma_data(n_studies = 10, seed = 123)

  result_std <- run_powernma(data, data_type = "pairwise", mode = "standard")

  output <- capture.output(print(result_std))
  expect_true(any(grepl("STANDARD", output)))
  expect_true(any(grepl("VALIDATED", output)))
})

test_that("Experimental mode has warnings field", {
  ipd <- generate_example_ipd(n_trials = 3, seed = 123)
  config <- setup_powernma(use_timevarying = TRUE, tau_list = 365)

  suppressWarnings({
    result <- run_powernma(ipd, data_type = "ipd", mode = "experimental", config = config)
  })

  expect_false(is.null(result$warnings))
  expect_true(length(result$warnings) > 0)
  expect_true(any(grepl("RMST", result$warnings)))
})

test_that("Standard mode runs netmeta correctly", {
  data <- simulate_nma_data(n_studies = 20, seed = 123)

  result <- run_powernma(data, data_type = "pairwise", mode = "standard")

  # Should have network results
  expect_false(is.null(result$network))
  expect_true(inherits(result$network, "netmeta"))

  # Should have heterogeneity stats
  expect_false(is.na(result$network$tau))
  expect_false(is.na(result$network$I2))
})

test_that("Validate mode-datatype combination function works", {
  # Standard + pairwise: should pass
  expect_silent(validate_mode_datatype("standard", "pairwise"))

  # Standard + ipd: should error
  expect_error(
    validate_mode_datatype("standard", "ipd"),
    "IPD-based.*EXPERIMENTAL"
  )

  # Experimental + both: should pass
  expect_silent(validate_mode_datatype("experimental", "pairwise"))
  expect_silent(validate_mode_datatype("experimental", "ipd"))
})

test_that("Mode description function works", {
  expect_match(get_mode_description("standard"), "Validated")
  expect_match(get_mode_description("experimental"), "research use")
})

test_that("Results have correct class based on mode", {
  data <- simulate_nma_data(n_studies = 10, seed = 123)

  result_std <- run_powernma(data, data_type = "pairwise", mode = "standard")
  expect_true(inherits(result_std, "powernma_standard"))

  ipd <- generate_example_ipd(n_trials = 3, seed = 123)
  suppressWarnings({
    result_exp <- run_powernma(ipd, data_type = "ipd", mode = "experimental")
  })
  expect_true(inherits(result_exp, "powernma_experimental"))
})
