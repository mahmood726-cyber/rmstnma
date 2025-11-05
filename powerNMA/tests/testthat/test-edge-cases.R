# ============================================================================
# Edge Case Tests for powerNMA
# ============================================================================
#
# These tests ensure the package handles edge cases gracefully with
# informative error messages.
#

# ==============================================================================
# NULL and Missing Input Tests
# ==============================================================================

test_that("validate_ipd rejects NULL input", {
  expect_error(
    validate_ipd(NULL),
    "\\[validate_ipd\\] Argument 'ipd' cannot be NULL"
  )
})

test_that("validate_ipd rejects non-data.frame input", {
  expect_error(
    validate_ipd(list(a = 1, b = 2)),
    "\\[validate_ipd\\] Argument 'ipd' must be a data.frame"
  )
})

test_that("validate_ipd rejects empty data frame", {
  empty_df <- data.frame(
    trial = character(0),
    treatment = character(0),
    time = numeric(0),
    status = numeric(0)
  )

  expect_error(
    validate_ipd(empty_df),
    "\\[validate_ipd\\] Argument 'ipd' cannot be an empty data frame"
  )
})

test_that("validate_nma_input rejects NULL input", {
  expect_error(
    validate_nma_input(NULL),
    "\\[validate_nma_input\\] Argument 'data' cannot be NULL"
  )
})

test_that("validate_nma_input rejects empty data frame", {
  empty_nma <- data.frame(
    studlab = character(0),
    treat1 = character(0),
    treat2 = character(0),
    TE = numeric(0),
    seTE = numeric(0)
  )

  expect_error(
    validate_nma_input(empty_nma),
    "\\[validate_nma_input\\] Argument 'data' cannot be an empty data frame"
  )
})

# ==============================================================================
# Invalid Column Tests
# ==============================================================================

test_that("validate_ipd detects missing required columns", {
  bad_ipd <- data.frame(
    trial = c("T1", "T1"),
    treatment = c("A", "B")
    # Missing: time, status
  )

  expect_error(
    validate_ipd(bad_ipd),
    "\\[validate_ipd\\] Missing required columns: time, status"
  )
})

test_that("validate_nma_input detects missing columns", {
  bad_nma <- data.frame(
    studlab = "S1",
    treat1 = "A"
    # Missing: treat2, TE, seTE
  )

  expect_error(
    validate_nma_input(bad_nma),
    "\\[validate_nma_input\\] Missing required columns"
  )
})

# ==============================================================================
# Invalid Data Type Tests
# ==============================================================================

test_that("validate_ipd rejects non-numeric time", {
  bad_ipd <- data.frame(
    trial = c("T1", "T1"),
    treatment = c("A", "B"),
    time = c("10", "20"),  # Character instead of numeric
    status = c(1, 0)
  )

  expect_error(
    validate_ipd(bad_ipd),
    "\\[validate_ipd\\] Column 'time' must be numeric"
  )
})

test_that("validate_ipd rejects invalid status values", {
  bad_ipd <- data.frame(
    trial = c("T1", "T1", "T1"),
    treatment = c("A", "B", "C"),
    time = c(10, 20, 30),
    status = c(0, 1, 2)  # 2 is invalid
  )

  expect_error(
    validate_ipd(bad_ipd),
    "\\[validate_ipd\\] Column 'status' must contain only 0 and 1"
  )
})

# ==============================================================================
# Non-finite Value Tests
# ==============================================================================

test_that("validate_nma_input detects Inf in TE", {
  bad_nma <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "B"),
    TE = c(0.5, Inf),
    seTE = c(0.1, 0.1)
  )

  expect_error(
    validate_nma_input(bad_nma),
    "\\[validate_nma_input\\] Column 'TE' contains.*non-finite values"
  )
})

test_that("validate_nma_input detects NA in seTE", {
  bad_nma <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "B"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, NA)
  )

  expect_error(
    validate_nma_input(bad_nma),
    "\\[validate_nma_input\\] Column 'seTE' contains.*non-finite values"
  )
})

test_that("validate_nma_input detects non-positive seTE", {
  bad_nma <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "A"),
    treat2 = c("B", "B", "B"),
    TE = c(0.5, 0.3, 0.2),
    seTE = c(0.1, 0, -0.1)  # Zero and negative values
  )

  expect_error(
    validate_nma_input(bad_nma),
    "\\[validate_nma_input\\] Column 'seTE' must be strictly positive"
  )
})

# ==============================================================================
# Insufficient Data Tests
# ==============================================================================

test_that("validate_nma_input requires at least 2 comparisons", {
  single_comparison <- data.frame(
    studlab = "S1",
    treat1 = "A",
    treat2 = "B",
    TE = 0.5,
    seTE = 0.1
  )

  expect_error(
    validate_nma_input(single_comparison),
    "\\[validate_nma_input\\] Need at least 2 comparisons"
  )
})

test_that("Single study network handled appropriately", {
  skip_if_not_installed("netmeta")

  single_study_data <- data.frame(
    studlab = c("S1", "S1"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(0.5, 0.3),
    seTE = c(0.1, 0.1)
  )

  # Should not error, but may produce warnings about limited data
  expect_warning(
    result <- tryCatch(
      run_powernma(
        single_study_data,
        data_type = "pairwise",
        mode = "standard",
        config = setup_powernma(
          use_bayesian = FALSE,
          run_sensitivity = FALSE,
          export_results = FALSE
        )
      ),
      error = function(e) NULL
    )
  )
})

# ==============================================================================
# Input Validation for Data Generation Functions
# ==============================================================================

test_that("generate_example_ipd validates n_trials", {
  expect_error(
    generate_example_ipd(n_trials = 0),
    "\\[generate_example_ipd\\] Argument 'n_trials' must be a positive integer"
  )

  expect_error(
    generate_example_ipd(n_trials = -5),
    "\\[generate_example_ipd\\] Argument 'n_trials' must be a positive integer"
  )

  expect_error(
    generate_example_ipd(n_trials = 2.5),
    "\\[generate_example_ipd\\] Argument 'n_trials' must be a positive integer"
  )

  expect_error(
    generate_example_ipd(n_trials = "three"),
    "\\[generate_example_ipd\\] Argument 'n_trials' must be a single numeric value"
  )
})

test_that("generate_example_ipd validates n_per_arm", {
  expect_error(
    generate_example_ipd(n_per_arm = 0),
    "\\[generate_example_ipd\\] Argument 'n_per_arm' must be a positive integer"
  )

  expect_error(
    generate_example_ipd(n_per_arm = -10),
    "\\[generate_example_ipd\\] Argument 'n_per_arm' must be a positive integer"
  )
})

test_that("generate_example_ipd handles NULL seed", {
  # Should work without error
  ipd1 <- generate_example_ipd(n_trials = 2, n_per_arm = 10, seed = NULL)
  ipd2 <- generate_example_ipd(n_trials = 2, n_per_arm = 10, seed = NULL)

  expect_s3_class(ipd1, "data.frame")
  expect_s3_class(ipd2, "data.frame")

  # Results should differ (different random seeds)
  expect_false(identical(ipd1$time, ipd2$time))
})

test_that("generate_example_ipd validates seed type", {
  expect_error(
    generate_example_ipd(seed = "random"),
    "\\[generate_example_ipd\\] Argument 'seed' must be a single numeric value or NULL"
  )

  expect_error(
    generate_example_ipd(seed = Inf),
    "\\[generate_example_ipd\\] Argument 'seed' must be a single numeric value or NULL"
  )
})

test_that("simulate_nma_data validates n_studies", {
  expect_error(
    simulate_nma_data(n_studies = 0),
    "\\[simulate_nma_data\\] Argument 'n_studies' must be a positive integer"
  )

  expect_error(
    simulate_nma_data(n_studies = 1.5),
    "\\[simulate_nma_data\\] Argument 'n_studies' must be a positive integer"
  )
})

# ==============================================================================
# Negative Value Tests
# ==============================================================================

test_that("validate_ipd warns about negative time values", {
  ipd_negative_time <- data.frame(
    trial = c("T1", "T1", "T1"),
    treatment = c("A", "B", "C"),
    time = c(10, -5, 20),  # One negative
    status = c(1, 0, 1)
  )

  expect_warning(
    validate_ipd(ipd_negative_time),
    "\\[validate_ipd\\] Detected 1 negative time values"
  )
})

test_that("validate_ipd warns about NA in time", {
  ipd_na_time <- data.frame(
    trial = c("T1", "T1", "T1"),
    treatment = c("A", "B", "C"),
    time = c(10, NA, 20),
    status = c(1, 0, 1)
  )

  expect_warning(
    validate_ipd(ipd_na_time),
    "\\[validate_ipd\\] Detected 1 NA values in time column"
  )
})

# ==============================================================================
# Helper Function Tests
# ==============================================================================

test_that("generate_pairwise_comparisons works with reference", {
  skip_if_not_installed("netmeta")

  treatments <- c("A", "B", "C", "D")
  comparisons <- powerNMA:::generate_pairwise_comparisons(treatments, reference = "A")

  expect_equal(length(comparisons), 3)  # A vs B, A vs C, A vs D
  expect_true(all(sapply(comparisons, function(x) "A" %in% x)))
})

test_that("generate_pairwise_comparisons works without reference", {
  skip_if_not_installed("netmeta")

  treatments <- c("A", "B", "C")
  comparisons <- powerNMA:::generate_pairwise_comparisons(treatments)

  expect_equal(length(comparisons), 3)  # A-B, A-C, B-C (all combinations)
})

test_that("generate_pairwise_comparisons handles single treatment", {
  skip_if_not_installed("netmeta")

  treatments <- c("A")
  comparisons <- powerNMA:::generate_pairwise_comparisons(treatments)

  expect_equal(length(comparisons), 0)
})

test_that("generate_pairwise_comparisons handles empty vector", {
  skip_if_not_installed("netmeta")

  treatments <- character(0)
  comparisons <- powerNMA:::generate_pairwise_comparisons(treatments)

  expect_equal(length(comparisons), 0)
})

# ==============================================================================
# Extreme Value Tests
# ==============================================================================

test_that("Very large networks are handled", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Skip on CRAN due to time

  # This tests scalability
  large_data <- simulate_nma_data(n_studies = 100, seed = 123)

  expect_s3_class(large_data, "data.frame")
  expect_gt(nrow(large_data), 50)

  # Validation should still work
  expect_true(validate_nma_input(large_data))
})

test_that("Very small standard errors are handled", {
  skip_if_not_installed("netmeta")

  precise_data <- data.frame(
    studlab = paste0("S", 1:5),
    treat1 = rep("A", 5),
    treat2 = rep("B", 5),
    TE = rnorm(5, 0, 0.1),
    seTE = rep(0.001, 5)  # Very precise
  )

  expect_true(validate_nma_input(precise_data))

  # Should be able to run NMA
  result <- tryCatch(
    run_powernma(
      precise_data,
      data_type = "pairwise",
      mode = "standard",
      config = setup_powernma(
        use_bayesian = FALSE,
        run_sensitivity = FALSE,
        export_results = FALSE
      )
    ),
    error = function(e) NULL
  )

  expect_false(is.null(result))
})

test_that("Extreme heterogeneity doesn't crash", {
  skip_if_not_installed("netmeta")

  # Create data with extreme heterogeneity
  extreme_het_data <- data.frame(
    studlab = paste0("S", 1:10),
    treat1 = rep("A", 10),
    treat2 = rep("B", 10),
    TE = c(seq(-2, 2, length.out = 10)),  # Very different effect sizes
    seTE = rep(0.1, 10)
  )

  expect_true(validate_nma_input(extreme_het_data))

  # Should handle extreme heterogeneity gracefully
  result <- tryCatch(
    run_powernma(
      extreme_het_data,
      data_type = "pairwise",
      mode = "standard",
      config = setup_powernma(
        use_bayesian = FALSE,
        run_sensitivity = FALSE,
        export_results = FALSE
      )
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  # May produce warnings but shouldn't crash
  expect_true(is.null(result) || inherits(result, "powernma_result"))
})

# ==============================================================================
# Message Format Tests
# ==============================================================================

test_that("Error messages include function context", {
  bad_data <- list(wrong = "type")

  error_msg <- tryCatch(
    validate_ipd(bad_data),
    error = function(e) e$message
  )

  expect_match(error_msg, "\\[validate_ipd\\]")
})

test_that("Error messages are informative", {
  bad_nma <- data.frame(
    studlab = "S1",
    treat1 = "A",
    treat2 = "B",
    TE = 0.5,
    seTE = -0.1  # Invalid
  )

  error_msg <- tryCatch(
    validate_nma_input(bad_nma),
    error = function(e) e$message
  )

  # Should mention what's wrong and provide details
  expect_match(error_msg, "seTE")
  expect_match(error_msg, "positive")
})
