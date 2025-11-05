# ============================================================================
# Integration Tests for powerNMA
# ============================================================================
#
# These tests verify that different components work together correctly
# in end-to-end workflows.
#

test_that("Complete standard NMA workflow executes successfully", {
  skip_if_not_installed("netmeta")

  # Generate data
  data <- simulate_nma_data(n_studies = 15, seed = 123)

  # Validate data
  expect_true(validate_nma_input(data))

  # Run diagnostics
  diag <- diagnose_nma_data(data, return_details = TRUE)
  expect_equal(diag$status, "PASS")

  # Run analysis
  config <- setup_powernma(
    use_bayesian = FALSE,
    run_sensitivity = FALSE,
    export_results = FALSE
  )

  results <- run_powernma(
    data = data,
    data_type = "pairwise",
    mode = "standard",
    config = config
  )

  # Check results structure
  expect_s3_class(results, "powernma_result")
  expect_equal(results$mode, "standard")
  expect_false(is.null(results$network))
})

test_that("Data diagnostics detect and report issues correctly", {
  # Create problematic data
  bad_data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "B", "C"),
    TE = c(0.5, Inf, 0.3),  # One Inf value
    seTE = c(0.1, 0.1, -0.1)  # One negative value
  )

  # Diagnose should detect issues
  diag <- diagnose_nma_data(bad_data, return_details = TRUE)

  expect_true(diag$has_inf_TE)
  expect_true(diag$has_negative_seTE)
  expect_gt(diag$n_issues, 0)
  expect_equal(diag$status, "FAILED")
})

test_that("Export functionality works for all modes", {
  skip_if_not_installed("netmeta")

  data <- simulate_nma_data(n_studies = 10, seed = 456)

  # Test standard mode export
  config_standard <- setup_powernma(
    use_bayesian = FALSE,
    run_sensitivity = TRUE,
    export_results = TRUE,
    output_dir = tempdir()
  )

  results_standard <- run_powernma(
    data = data,
    data_type = "pairwise",
    mode = "standard",
    config = config_standard
  )

  # Check that files were created
  expect_true(file.exists(file.path(tempdir(), "powernma_results.rds")))

  # Clean up
  unlink(file.path(tempdir(), "*.rds"))
  unlink(file.path(tempdir(), "*.csv"))
})

test_that("Data transformation helpers work correctly", {
  data <- simulate_nma_data(n_studies = 20, seed = 789)

  # Test filtering
  filtered <- filter_treatments(data, exclude_treatments = "DrugD")
  expect_false("DrugD" %in% c(filtered$treat1, filtered$treat2))

  # Test standardization
  data_with_spaces <- data
  data_with_spaces$treat1 <- paste0(" ", data_with_spaces$treat1, " ")
  standardized <- standardize_treatment_names(data_with_spaces)
  expect_false(any(grepl("^\\s|\\s$", standardized$treat1)))

  # Test comparison matrix
  matrix <- create_comparison_matrix(data)
  expect_true(is.matrix(matrix))
  expect_equal(nrow(matrix), ncol(matrix))

  # Test treatment pairs
  pairs <- get_treatment_pairs(data)
  expect_true(all(c("treat1", "treat2", "comparison_type") %in% names(pairs)))
})

test_that("Progress and performance utilities work", {
  # Test memory reporting
  expect_silent(report_memory_usage())

  # Test memory estimation
  mem_estimate <- estimate_memory_needs(
    n_studies = 50,
    n_treatments = 10,
    use_bayesian = FALSE
  )
  expect_true(is.numeric(mem_estimate))
  expect_gt(mem_estimate, 0)

  # Test timing
  result <- time_operation({
    Sys.sleep(0.1)
    42
  }, label = "Test operation")

  expect_equal(result, 42)
})

test_that("Validation helpers catch all edge cases", {
  # NULL checks
  expect_error(
    assert_not_null(NULL, "test_arg", "test_func"),
    "\\[test_func\\].*test_arg.*cannot be NULL"
  )

  # Data frame checks
  expect_error(
    assert_data_frame(list(a = 1), "test_arg", "test_func"),
    "\\[test_func\\].*must be a data.frame"
  )

  # Numeric checks
  expect_error(
    assert_numeric("text", "test_arg", "test_func"),
    "\\[test_func\\].*must be numeric"
  )

  expect_error(
    assert_numeric(Inf, "test_arg", "test_func"),
    "\\[test_func\\].*non-finite values"
  )

  expect_error(
    assert_numeric(5, "test_arg", "test_func", min_value = 10),
    "\\[test_func\\].*must be in range"
  )

  # Positive integer checks
  expect_error(
    assert_positive_integer(0, "test_arg", "test_func"),
    "\\[test_func\\].*must be a positive integer"
  )

  expect_error(
    assert_positive_integer(2.5, "test_arg", "test_func"),
    "\\[test_func\\].*must be a positive integer"
  )

  # Column checks
  test_df <- data.frame(a = 1, b = 2)
  expect_error(
    assert_columns_exist(test_df, c("a", "c"), "test_func"),
    "\\[test_func\\].*Missing required columns: c"
  )
})

test_that("Subgroup analysis workflow works", {
  data <- simulate_nma_data(n_studies = 30, seed = 321)

  # Split by study design
  subgroups <- split_by_subgroup(data, "study_design")

  expect_true(is.list(subgroups))
  expect_gt(length(subgroups), 0)

  # Each subgroup should be analyzable
  for (sg_name in names(subgroups)) {
    sg_data <- subgroups[[sg_name]]

    if (nrow(sg_data) >= .MIN_COMPARISONS_FOR_NMA) {
      expect_true(validate_nma_input(sg_data))
    }
  }
})

test_that("Diagnostic functions provide comprehensive information", {
  data <- simulate_nma_data(n_studies = 25, seed = 654)

  # Full diagnostics
  diag <- diagnose_nma_data(data, return_details = TRUE)

  # Check all expected fields are present
  expected_fields <- c(
    "has_required_columns", "n_studies", "n_comparisons",
    "n_treatments", "TE_range", "TE_mean", "TE_sd",
    "seTE_range", "status", "issues", "warnings",
    "recommendations", "suggested_reference"
  )

  for (field in expected_fields) {
    expect_true(field %in% names(diag),
               info = sprintf("Missing field: %s", field))
  }

  # Test quick summary
  summary <- quick_summary(data)
  expect_equal(summary$n_studies, diag$n_studies)
})

test_that("Helper function generate_pairwise_comparisons handles all cases", {
  skip_if_not_installed("netmeta")

  # Case 1: With reference treatment
  treatments1 <- c("A", "B", "C", "D")
  comps1 <- powerNMA:::generate_pairwise_comparisons(treatments1, reference = "A")
  expect_equal(length(comps1), 3)  # A-B, A-C, A-D
  expect_true(all(sapply(comps1, function(x) "A" %in% x)))

  # Case 2: Without reference (all pairwise)
  treatments2 <- c("A", "B", "C")
  comps2 <- powerNMA:::generate_pairwise_comparisons(treatments2)
  expect_equal(length(comps2), 3)  # A-B, A-C, B-C

  # Case 3: Single treatment
  treatments3 <- c("A")
  comps3 <- powerNMA:::generate_pairwise_comparisons(treatments3)
  expect_equal(length(comps3), 0)

  # Case 4: Empty vector
  treatments4 <- character(0)
  comps4 <- powerNMA:::generate_pairwise_comparisons(treatments4)
  expect_equal(length(comps4), 0)

  # Case 5: Reference not in treatment list
  treatments5 <- c("B", "C", "D")
  comps5 <- powerNMA:::generate_pairwise_comparisons(treatments5, reference = "A")
  expect_equal(length(comps5), 3)  # Falls back to all pairwise
})

test_that("Constants are properly defined and accessible", {
  # Check key constants exist
  expect_true(exists(".DEFAULT_SE_MIN"))
  expect_true(exists(".DEFAULT_SE_MAX"))
  expect_true(exists(".MIN_STUDIES_FOR_LOO"))
  expect_true(exists(".MIN_COMPARISONS_FOR_NMA"))

  # Check helper functions work
  sim_defaults <- get_simulation_defaults()
  expect_true(is.list(sim_defaults))
  expect_true("se_min" %in% names(sim_defaults))

  event_rates <- get_event_rates()
  expect_true(is.numeric(event_rates))
  expect_true("Control" %in% names(event_rates))

  true_effects <- get_true_effects()
  expect_true(is.numeric(true_effects))
})

test_that("Export function handles errors gracefully", {
  skip_if_not_installed("netmeta")

  data <- simulate_nma_data(n_studies = 10, seed = 999)
  config <- setup_powernma(
    use_bayesian = FALSE,
    run_sensitivity = FALSE,
    export_results = FALSE
  )

  results <- run_powernma(
    data = data,
    data_type = "pairwise",
    mode = "standard",
    config = config
  )

  # Test export to invalid directory (should create it)
  test_dir <- file.path(tempdir(), "nonexistent_dir")
  if (dir.exists(test_dir)) unlink(test_dir, recursive = TRUE)

  exported <- export_powernma_tables_unified(results, test_dir, verbose = FALSE)

  expect_true(dir.exists(test_dir))
  expect_true(is.character(exported))

  # Clean up
  unlink(test_dir, recursive = TRUE)
})

test_that("Complete workflow with all features", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Comprehensive test, skip on CRAN

  # 1. Generate data
  data <- simulate_nma_data(n_studies = 30, seed = 111)

  # 2. Run diagnostics
  diag <- diagnose_nma_data(data, return_details = TRUE)
  expect_equal(diag$status, "PASS")

  # 3. Standardize treatment names
  data_clean <- standardize_treatment_names(data, trim_whitespace = TRUE)

  # 4. Create comparison matrix
  matrix <- create_comparison_matrix(data_clean)
  expect_true(is.matrix(matrix))

  # 5. Quick summary
  summary <- quick_summary(data_clean)
  expect_gt(summary$n_studies, 0)

  # 6. Run analysis with timing
  result <- time_operation({
    run_powernma(
      data = data_clean,
      data_type = "pairwise",
      mode = "standard",
      config = setup_powernma(
        use_bayesian = FALSE,
        run_sensitivity = TRUE,
        export_results = FALSE
      )
    )
  }, label = "Complete NMA")

  expect_s3_class(result, "powernma_result")

  # 7. Verify all components
  expect_false(is.null(result$network))
  expect_false(is.null(result$geometry))
  expect_equal(result$mode, "standard")
  expect_equal(result$validated, TRUE)
})
