context("Phase 8: Network Analysis, Simulation & Integration")

library(testthat)
library(powerNMA)

# Generate test data
set.seed(456)
test_data <- data.frame(
  studlab = rep(paste0("Study", 1:10), each = 3),
  treat1 = rep(c("Placebo", "Placebo", "A"), 10),
  treat2 = rep(c("A", "B", "B"), 10),
  TE = rnorm(30, 0.5, 0.3),
  seTE = runif(30, 0.1, 0.3),
  stringsAsFactors = FALSE
)

# ==============================================================================
# Network Geometry Tests
# ==============================================================================

test_that("Network geometry analysis works", {
  geometry <- analyze_network_geometry(test_data)

  expect_is(geometry, "network_geometry")
  expect_true("connectivity" %in% names(geometry))
  expect_true("characteristics" %in% names(geometry))
  expect_true("graph_metrics" %in% names(geometry))
  expect_is(geometry$adjacency_matrix, "matrix")
})

test_that("Network connectivity is assessed correctly", {
  geometry <- analyze_network_geometry(test_data)

  expect_is(geometry$connectivity$is_connected, "logical")
  expect_is(geometry$connectivity$density, "numeric")
  expect_true(geometry$connectivity$density >= 0 && geometry$connectivity$density <= 1)
})

test_that("Network robustness is calculated", {
  geometry <- analyze_network_geometry(test_data)

  expect_true("robustness" %in% names(geometry))
  expect_is(geometry$robustness$robustness_score, "numeric")
  expect_true(geometry$robustness$robustness_score >= 0 &&
             geometry$robustness$robustness_score <= 1)
})

# ==============================================================================
# Simulation & Power Analysis Tests
# ==============================================================================

test_that("NMA data simulation works", {
  sim_data <- simulate_nma_data(
    n_treatments = 4,
    n_studies = 15,
    network_structure = "star",
    effect_size = 0.5
  )

  expect_is(sim_data, "data.frame")
  expect_true("studlab" %in% names(sim_data))
  expect_true("treat1" %in% names(sim_data))
  expect_true("treat2" %in% names(sim_data))
  expect_true("TE" %in% names(sim_data))
  expect_true("seTE" %in% names(sim_data))
})

test_that("Different network structures can be simulated", {
  star_data <- simulate_nma_data(n_treatments = 4, n_studies = 10, network_structure = "star")
  line_data <- simulate_nma_data(n_treatments = 4, n_studies = 10, network_structure = "line")
  complete_data <- simulate_nma_data(n_treatments = 4, n_studies = 10, network_structure = "complete")

  expect_is(star_data, "data.frame")
  expect_is(line_data, "data.frame")
  expect_is(complete_data, "data.frame")
})

test_that("Power analysis runs", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Power analysis takes time

  power_result <- nma_power_analysis(
    n_treatments = 4,
    n_studies_per_comparison = 3,
    effect_size = 0.5,
    n_simulations = 100  # Reduced for testing
  )

  expect_is(power_result, "nma_power")
  expect_is(power_result$power, "numeric")
  expect_true(power_result$power >= 0 && power_result$power <= 1)
})

test_that("Sample size calculation works", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Takes time

  ss_result <- calculate_nma_sample_size(
    n_treatments = 4,
    effect_size = 0.5,
    desired_power = 0.80
  )

  expect_is(ss_result, "nma_sample_size")
  expect_true("required_studies_per_comparison" %in% names(ss_result))
})

test_that("Study allocation optimization works", {
  allocation <- optimize_study_allocation(
    n_treatments = 5,
    total_studies = 30
  )

  expect_is(allocation, "data.frame")
  expect_equal(sum(allocation$Allocated_Studies), 30)
})

# ==============================================================================
# Export & Interoperability Tests
# ==============================================================================

test_that("NMA results can be exported", {
  skip_if_not_installed("netmeta")

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data
  )

  temp_dir <- tempdir()

  expect_error({
    files <- export_nma_results(
      nma_result = nma,
      output_dir = file.path(temp_dir, "export_test"),
      formats = c("csv", "json", "html"),
      prefix = "test_nma"
    )
  }, NA)

  # Cleanup
  unlink(file.path(temp_dir, "export_test"), recursive = TRUE)
})

test_that("GRADE template can be created", {
  skip_if_not_installed("netmeta")

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data
  )

  temp_file <- tempfile(fileext = ".csv")

  grade_df <- export_for_grade(nma, temp_file)

  expect_is(grade_df, "data.frame")
  expect_true("Risk_of_Bias" %in% names(grade_df))
  expect_true("Overall_Quality" %in% names(grade_df))
  expect_true(file.exists(temp_file))

  # Cleanup
  unlink(temp_file)
})

test_that("BUGSnet conversion works", {
  bugsnet_data <- convert_to_bugsnet(test_data)

  expect_is(bugsnet_data, "data.frame")
  expect_true("studyID" %in% names(bugsnet_data))
  expect_true("T" %in% names(bugsnet_data))
  expect_true("C" %in% names(bugsnet_data))
})

# ==============================================================================
# Comprehensive Wrapper Tests
# ==============================================================================

test_that("Ultimate comprehensive NMA wrapper works", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Comprehensive analysis takes time

  expect_error({
    result <- run_ultimate_nma(
      data = test_data,
      sm = "MD",
      reference = "Placebo",
      assess_geometry = TRUE,
      calculate_rankings = TRUE,
      assess_inconsistency = FALSE,  # Skip to save time
      assess_publication_bias = FALSE,
      generate_visualizations = FALSE,
      export_results = FALSE,
      perform_simulation = FALSE
    )
  }, NA)

  expect_is(result, "comprehensive_nma")
  expect_true("core_nma" %in% names(result))
  expect_true("geometry" %in% names(result))
})

test_that("Quick NMA analysis works", {
  skip_if_not_installed("netmeta")

  expect_error({
    result <- quick_nma(test_data, sm = "MD", reference = "Placebo")
  }, NA)

  expect_is(result, "quick_nma")
  expect_true("nma" %in% names(result))
  expect_true("sucra" %in% names(result))
  expect_true("heterogeneity" %in% names(result))
})

# ==============================================================================
# Print/Plot Methods Tests
# ==============================================================================

test_that("Network geometry has print method", {
  geometry <- analyze_network_geometry(test_data)

  expect_error(print(geometry), NA)
  expect_output(print(geometry), "Network Geometry Analysis")
})

test_that("Power analysis has print method", {
  skip_if_not_installed("netmeta")
  skip_on_cran()

  power_result <- nma_power_analysis(
    n_treatments = 4,
    n_studies_per_comparison = 3,
    n_simulations = 50
  )

  expect_error(print(power_result), NA)
  expect_output(print(power_result), "Power Analysis")
})

test_that("Comprehensive NMA has print method", {
  skip_if_not_installed("netmeta")
  skip_on_cran()

  result <- run_ultimate_nma(
    data = test_data,
    assess_geometry = FALSE,
    calculate_rankings = FALSE,
    assess_inconsistency = FALSE,
    assess_publication_bias = FALSE,
    generate_visualizations = FALSE
  )

  expect_error(print(result), NA)
  expect_output(print(result), "COMPREHENSIVE")
})

# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("Phase 8 methods integrate with Phase 7", {
  skip_if_not_installed("netmeta")
  skip_on_cran()

  # Should be able to use geometry with other features
  geometry <- analyze_network_geometry(test_data)

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data
  )

  # Should work together without errors
  expect_true(geometry$n_treatments == length(nma$trts))
})

# ==============================================================================
# Edge Cases
# ==============================================================================

test_that("Handles disconnected networks", {
  # Create disconnected network
  disconnected_data <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "C"),
    treat2 = c("B", "D"),
    TE = c(0.5, 0.5),
    seTE = c(0.2, 0.2),
    stringsAsFactors = FALSE
  )

  geometry <- analyze_network_geometry(disconnected_data)

  expect_false(geometry$connectivity$is_connected)
  expect_equal(geometry$connectivity$n_components, 2)
})

test_that("Handles small networks", {
  small_data <- test_data[1:3, ]

  expect_error(analyze_network_geometry(small_data), NA)
})

test_that("Simulation handles binary outcomes", {
  binary_sim <- simulate_nma_data(
    n_treatments = 3,
    n_studies = 10,
    outcome_type = "binary"
  )

  expect_is(binary_sim, "data.frame")
  expect_true(all(is.finite(binary_sim$TE)))
})
