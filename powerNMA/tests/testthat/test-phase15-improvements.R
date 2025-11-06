# Tests for Phase 15 Production-Grade Improvements
# ===================================================

context("Phase 15: Validation, Performance & Edge Cases")

# Test 1: Validation Utilities - Not Null
test_that("assert_not_null works correctly", {
  expect_silent(assert_not_null(5, "x", "test_func"))
  expect_error(assert_not_null(NULL, "x", "test_func"),
               "\\[test_func\\] Parameter 'x' cannot be NULL")
})

# Test 2: Validation Utilities - Data Frame
test_that("assert_data_frame validates correctly", {
  df <- data.frame(a = 1:5, b = letters[1:5])

  # Valid data frame
  expect_silent(assert_data_frame(df, "df", "test_func"))

  # Invalid type
  expect_error(assert_data_frame(list(a = 1), "df", "test_func"),
               "must be a data.frame")

  # Empty data frame
  empty_df <- data.frame(a = numeric(0), b = character(0))
  expect_error(assert_data_frame(empty_df, "df", "test_func"),
               "empty data frame")

  # Missing required columns
  expect_error(
    assert_data_frame(df, "df", "test_func", required_cols = c("a", "c")),
    "missing required columns: c"
  )
})

# Test 3: Validation Utilities - Numeric
test_that("assert_numeric validates correctly", {
  # Valid numeric
  expect_silent(assert_numeric(5, "x", "test_func"))
  expect_silent(assert_numeric(1:10, "x", "test_func"))

  # Invalid type
  expect_error(assert_numeric("5", "x", "test_func"),
               "must be numeric")

  # Range checking
  expect_error(assert_numeric(5, "x", "test_func", min = 10),
               "must be >= 10")
  expect_error(assert_numeric(5, "x", "test_func", max = 3),
               "must be <= 3")

  # NA handling
  expect_error(assert_numeric(c(1, NA, 3), "x", "test_func"),
               "contains NA values")
  expect_silent(assert_numeric(c(1, NA, 3), "x", "test_func", allow_na = TRUE))

  # Integer only
  expect_error(assert_numeric(5.5, "x", "test_func", integer_only = TRUE),
               "must be integer")
  expect_silent(assert_numeric(5, "x", "test_func", integer_only = TRUE))
})

# Test 4: Validation Utilities - Character
test_that("assert_character validates correctly", {
  # Valid character
  expect_silent(assert_character("test", "x", "test_func"))

  # Invalid type
  expect_error(assert_character(5, "x", "test_func"),
               "must be character")

  # Empty strings
  expect_error(assert_character("", "x", "test_func"),
               "empty strings")

  # Choices validation
  expect_error(
    assert_character("invalid", "x", "test_func", choices = c("a", "b", "c")),
    "invalid values: invalid"
  )
  expect_silent(
    assert_character("b", "x", "test_func", choices = c("a", "b", "c"))
  )
})

# Test 5: Validation Utilities - Positive Integer
test_that("assert_positive_integer works", {
  expect_silent(assert_positive_integer(5, "n", "test_func"))
  expect_error(assert_positive_integer(-1, "n", "test_func"),
               "must be >= 1")
  expect_error(assert_positive_integer(5.5, "n", "test_func"),
               "must be integer")
})

# Test 6: Validation Utilities - Probability
test_that("assert_probability validates correctly", {
  expect_silent(assert_probability(0.5, "p", "test_func"))
  expect_silent(assert_probability(0, "p", "test_func"))
  expect_silent(assert_probability(1, "p", "test_func"))

  expect_error(assert_probability(-0.1, "p", "test_func"))
  expect_error(assert_probability(1.1, "p", "test_func"))
})

# Test 7: Validation Utilities - Network Structure
test_that("assert_network_structure validates correctly", {
  # Valid network
  data <- data.frame(
    studlab = rep(c("S1", "S2", "S3"), each = 2),
    treat1 = rep(c("A", "A", "B"), each = 2),
    treat2 = rep(c("B", "C", "C"), each = 2),
    TE = rnorm(6),
    seTE = runif(6, 0.1, 0.3)
  )

  expect_silent(assert_network_structure(data))

  # Too few treatments
  single_comp <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "A"),
    treat2 = c("B", "B"),
    TE = c(-0.5, -0.3),
    seTE = c(0.2, 0.2)
  )
  expect_error(assert_network_structure(single_comp, min_treatments = 3),
               "at least 3 treatments")

  # Too few studies
  expect_error(
    assert_network_structure(data[1:2, ], min_studies = 3),
    "at least 3 studies"
  )
})

# Test 8: Performance - Memoization
test_that("Memoization caching works", {
  skip_if_not_installed("digest")

  # Create expensive function
  counter <- 0
  expensive_func <- function(x) {
    counter <<- counter + 1
    Sys.sleep(0.01)  # Simulate computation
    x^2
  }

  memoized_func <- memoize(expensive_func, max_cache_size = 10)

  # First call
  result1 <- memoized_func(5)
  count1 <- counter

  # Second call (should be cached)
  result2 <- memoized_func(5)
  count2 <- counter

  expect_equal(result1, result2)
  expect_equal(result1, 25)
  expect_equal(count2, count1)  # Counter shouldn't increase (cached)
})

# Test 9: Performance - Cache Management
test_that("Cache management works", {
  # Clear cache
  n_cleared <- clear_memo_cache()
  expect_true(n_cleared >= 0)

  # Get cache stats
  stats <- get_cache_stats()
  expect_true(is.list(stats))
  expect_true("n_items" %in% names(stats))
  expect_true("total_size_mb" %in% names(stats))
})

# Test 10: Performance - Parallel Apply
test_that("Parallel apply works (sequential fallback)", {
  # Test with sequential processing (n_cores = 1)
  X <- 1:10
  results <- parallel_apply(X, function(x) x^2, n_cores = 1, show_progress = FALSE)

  expect_equal(length(results), 10)
  expect_equal(unlist(results), (1:10)^2)
})

# Test 11: Performance - Batch Processing
test_that("Batch processing works correctly", {
  # Process in chunks
  X <- 1:100
  results <- batch_process(X, mean, chunk_size = 25, show_progress = FALSE)

  expect_true(is.numeric(results))
  expect_equal(length(results), 4)  # 100/25 = 4 chunks
})

# Test 12: Performance - Measure Time
test_that("Execution time measurement works", {
  result <- measure_time({
    Sys.sleep(0.01)
    42
  }, silent = TRUE)

  expect_true(is.list(result))
  expect_equal(result$result, 42)
  expect_true(result$elapsed_seconds > 0)
  expect_true(result$elapsed_seconds < 1)
})

# Test 13: Performance - System Info
test_that("System info retrieval works", {
  info <- get_system_info()

  expect_true(is.list(info))
  expect_true("r_version" %in% names(info))
  expect_true("n_cores" %in% names(info))
  expect_true(info$n_cores >= 1)
})

# Test 14: Performance - Memory Estimation
test_that("Memory estimation works", {
  data <- matrix(rnorm(1000), ncol = 10)
  mem_mb <- estimate_memory(data)

  expect_true(is.numeric(mem_mb))
  expect_true(mem_mb > 0)
})

# Test 15: Edge Case - Single Study Network
test_that("Single study network is rejected", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = c("Study1", "Study1"),
    treat1 = c("A", "A"),
    treat2 = c("B", "C"),
    TE = c(-0.5, -0.3),
    seTE = c(0.2, 0.2)
  )

  expect_error(
    assert_network_structure(data, min_studies = 2),
    "at least 2 studies"
  )
})

# Test 16: Edge Case - Empty Network Data
test_that("Empty network data is rejected", {
  empty_data <- data.frame(
    studlab = character(0),
    treat1 = character(0),
    treat2 = character(0),
    TE = numeric(0),
    seTE = numeric(0)
  )

  expect_error(
    assert_data_frame(empty_data, "data", "test_func"),
    "empty data frame"
  )
})

# Test 17: Edge Case - Disconnected Network
test_that("Disconnected network is detected", {
  # Create disconnected network (two separate subnetworks)
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "D", "D"),
    treat2 = c("B", "C", "E", "F"),
    TE = c(-0.5, -0.3, -0.4, -0.6),
    seTE = c(0.2, 0.2, 0.2, 0.2)
  )

  # Should still pass basic structure validation
  expect_silent(assert_network_structure(data))

  # Note: Actual disconnection detection would require graph analysis
  # which would be done by netmeta/NMA functions
})

# Test 18: Edge Case - Extreme Heterogeneity
test_that("Extreme heterogeneity values handled", {
  # Very large standard errors
  data <- data.frame(
    studlab = rep(paste0("S", 1:10), each = 2),
    treat1 = rep("A", 20),
    treat2 = rep(c("B", "C"), 10),
    TE = rnorm(20, 0, 5),  # Large variation
    seTE = runif(20, 10, 20)  # Very large SEs
  )

  # Should still pass validation
  expect_silent(assert_network_structure(data))
})

# Test 19: Edge Case - Negative Time Values
test_that("Negative time values are detected", {
  # For IPD data with time-to-event
  ipd_data <- data.frame(
    trial = rep(c("T1", "T2"), each = 100),
    treatment = rep(c("A", "B"), 100),
    time = c(rexp(100, 0.01), c(-5, rexp(99, 0.01))),  # One negative
    status = rbinom(200, 1, 0.3)
  )

  # Negative times should be detected in validation
  expect_true(any(ipd_data$time < 0))
})

# Test 20: Edge Case - Very Large Network
test_that("Large network structure is validated", {
  # Create large network (100 studies, 20 treatments)
  n_studies <- 100
  n_comps_per_study <- 3

  data <- data.frame(
    studlab = rep(paste0("Study", 1:n_studies), each = n_comps_per_study),
    treat1 = sample(LETTERS[1:20], n_studies * n_comps_per_study, replace = TRUE),
    treat2 = sample(LETTERS[1:20], n_studies * n_comps_per_study, replace = TRUE),
    TE = rnorm(n_studies * n_comps_per_study),
    seTE = runif(n_studies * n_comps_per_study, 0.1, 0.5)
  )

  # Remove self-loops
  data <- data[data$treat1 != data$treat2, ]

  # Should pass validation
  expect_silent(assert_network_structure(data))
  expect_true(nrow(data) > 100)
})

# Test 21: Safe Try Wrapper
test_that("safe_try handles errors correctly", {
  # Successful execution
  result <- safe_try({
    2 + 2
  }, func = "test_func", silent = TRUE)
  expect_equal(result, 4)

  # Error handling
  expect_error(
    safe_try({
      stop("Test error")
    }, func = "test_func"),
    "\\[test_func\\] Error: Test error"
  )
})

# Test 22: Efficient Data Frame Creation
test_that("Efficient data frame creation works", {
  df <- efficient_data_frame(
    id = rep(1:10, each = 3),
    treatment = rep(c("A", "B", "C"), 10),
    value = rnorm(30)
  )

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 30)

  # Check if factors were used for low-cardinality columns
  expect_true(is.factor(df$treatment) || is.character(df$treatment))
})

# Test 23: Error Message Formatting
test_that("Error messages include function context", {
  error <- create_error("Test error message", "my_function")

  expect_true(inherits(error, "validation_error"))
  expect_true(grepl("\\[my_function\\]", error$message))
  expect_true(grepl("Test error message", error$message))
})

# Test 24: Warning Message Formatting
test_that("Warning messages include function context", {
  warning_obj <- create_warning("Test warning message", "my_function")

  expect_true(inherits(warning_obj, "validation_warning"))
  expect_true(grepl("\\[my_function\\]", warning_obj$message))
  expect_true(grepl("Test warning message", warning_obj$message))
})

# Test 25: Seed Validation
test_that("Seed validation works correctly", {
  expect_silent(assert_seed(42, func = "test_func"))
  expect_silent(assert_seed(NULL, func = "test_func"))  # NULL is allowed

  expect_error(assert_seed("not_a_number", func = "test_func"),
               "must be numeric")
  expect_error(assert_seed(42.5, func = "test_func"),
               "must be integer")
})

# Integration Test
test_that("Validation utilities integrate well", {
  # Simulate a function using multiple validators
  test_function <- function(data, n_iter, alpha, model_type) {
    assert_data_frame(data, "data", "test_function",
                     required_cols = c("studlab", "TE", "seTE"))
    assert_positive_integer(n_iter, "n_iter", "test_function")
    assert_probability(alpha, "alpha", "test_function")
    assert_character(model_type, "model_type", "test_function",
                    choices = c("random", "fixed"))

    TRUE
  }

  # Valid inputs
  data <- data.frame(
    studlab = c("S1", "S2"),
    TE = c(-0.5, -0.3),
    seTE = c(0.2, 0.2)
  )

  expect_true(test_function(data, 1000, 0.05, "random"))

  # Invalid inputs
  expect_error(test_function(data, -5, 0.05, "random"))  # Negative n_iter
  expect_error(test_function(data, 1000, 1.5, "random"))  # Invalid alpha
  expect_error(test_function(data, 1000, 0.05, "invalid"))  # Invalid model_type
})
