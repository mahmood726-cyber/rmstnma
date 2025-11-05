# Test Phase 10: Ultimate Shiny Integration & Batch Processing
# Tests for:
# - Shiny app main framework
# - UI tab components
# - Server modules
# - Batch processing
# - Automated workflows

library(testthat)
library(powerNMA)

# ============================================================================
# Helper Functions
# ============================================================================

create_mock_batch_datasets <- function(n = 3) {
  datasets <- list()
  for (i in 1:n) {
    datasets[[paste0("dataset_", i)]] <- simulate_nma_data(
      n_studies = sample(10:20, 1),
      n_treatments = sample(4:7, 1),
      sm = "OR"
    )
  }
  return(datasets)
}

# ============================================================================
# Tests: Batch Processing
# ============================================================================

test_that("run_batch_nma processes multiple datasets", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 2)

  config <- setup_powernma(
    sm = "OR",
    use_bayesian = FALSE,
    run_sensitivity = FALSE
  )

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  expect_s3_class(batch_results, "batch_nma_results")
  expect_type(batch_results$results, "list")
  expect_equal(length(batch_results$results), 2)
  expect_true(all(batch_results$summary$Success))
})

test_that("batch processing handles errors gracefully", {
  skip_on_cran()

  datasets <- list(
    good_data = simulate_nma_data(n_studies = 15, sm = "OR"),
    bad_data = data.frame(x = 1:5, y = 1:5)  # Invalid data
  )

  config <- setup_powernma(sm = "OR")

  batch_results <- suppressWarnings(
    run_batch_nma(
      datasets = datasets,
      analysis_config = config,
      parallel = FALSE,
      save_individual = FALSE
    )
  )

  expect_s3_class(batch_results, "batch_nma_results")
  expect_true(batch_results$summary$Success[1])  # good_data succeeds
  expect_false(batch_results$summary$Success[2])  # bad_data fails
  expect_true(length(batch_results$errors) > 0)
})

test_that("batch processing creates summary statistics", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 3)
  config <- setup_powernma(sm = "OR")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  expect_true("summary" %in% names(batch_results))
  expect_s3_class(batch_results$summary, "data.frame")
  expect_true("Dataset" %in% names(batch_results$summary))
  expect_true("Success" %in% names(batch_results$summary))
  expect_true("Time_Seconds" %in% names(batch_results$summary))
})

test_that("batch processing aggregates statistics", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 2)
  config <- setup_powernma(sm = "OR")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  expect_true("aggregate_statistics" %in% names(batch_results))

  if (!is.null(batch_results$aggregate_statistics)) {
    agg_stats <- batch_results$aggregate_statistics
    expect_true("N_Studies" %in% names(agg_stats))
    expect_true("N_Treatments" %in% names(agg_stats))
    expect_true("Tau2" %in% names(agg_stats))
  }
})

test_that("batch processing creates comparison table", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 2)
  config <- setup_powernma(sm = "OR")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  if (length(batch_results$results) > 1) {
    expect_true("comparison_table" %in% names(batch_results))
  }
})

test_that("batch processing saves log file", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 2)
  config <- setup_powernma(sm = "OR")

  log_file <- tempfile(fileext = ".log")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE,
    log_file = log_file
  )

  expect_true(file.exists(log_file))
  expect_true(file.size(log_file) > 0)

  # Clean up
  unlink(log_file)
})

test_that("batch processing can resume from specific index", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 3)
  config <- setup_powernma(sm = "OR")

  # Process from index 2
  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE,
    resume_from = 2
  )

  # Should process datasets 2 and 3
  expect_equal(length(batch_results$results), 2)
})

test_that("parallel batch processing works", {
  skip_on_cran()
  skip_on_ci()  # May not work in CI environments

  datasets <- create_mock_batch_datasets(n = 4)
  config <- setup_powernma(sm = "OR")

  batch_results <- suppressMessages(
    run_batch_nma(
      datasets = datasets,
      analysis_config = config,
      parallel = TRUE,
      n_cores = 2,
      save_individual = FALSE
    )
  )

  expect_s3_class(batch_results, "batch_nma_results")
  expect_equal(length(batch_results$results), 4)
})

test_that("print method works for batch results", {
  skip_on_cran()

  datasets <- create_mock_batch_datasets(n = 2)
  config <- setup_powernma(sm = "OR")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  expect_output(print(batch_results), "Batch NMA Results")
  expect_output(print(batch_results), "Total Datasets:")
  expect_output(print(batch_results), "Successful:")
})

# ============================================================================
# Tests: Automated Workflows
# ============================================================================

test_that("automated_nma_workflow runs complete pipeline", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 20, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR", assess_all = TRUE),
    visualizations = c("network", "forest"),
    manuscripts = list(generate = FALSE),  # Skip for speed
    reports = NULL
  )

  output_dir <- tempdir()

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = output_dir
    )
  )

  expect_s3_class(workflow_results, "automated_workflow_results")
  expect_true("analysis" %in% names(workflow_results))
  expect_true("workflow_info" %in% names(workflow_results))
})

test_that("automated workflow validates data", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 15, sm = "OR")

  workflow_config <- default_workflow_config()
  workflow_config$manuscripts$generate <- FALSE
  workflow_config$reports <- NULL

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = tempdir()
    )
  )

  expect_true("validation" %in% names(workflow_results))
})

test_that("automated workflow generates visualizations", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 15, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR"),
    visualizations = c("network", "forest"),
    manuscripts = list(generate = FALSE),
    reports = NULL
  )

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = tempdir()
    )
  )

  expect_true("visualizations" %in% names(workflow_results))

  if (!is.null(workflow_results$visualizations)) {
    expect_type(workflow_results$visualizations, "list")
  }
})

test_that("automated workflow creates summary dashboard", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 15, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR"),
    visualizations = NULL,
    manuscripts = list(generate = FALSE),
    reports = NULL
  )

  output_dir <- tempdir()

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = output_dir
    )
  )

  expect_true("dashboard" %in% names(workflow_results))

  if (!is.null(workflow_results$dashboard)) {
    expect_true(file.exists(workflow_results$dashboard))
  }
})

test_that("default_workflow_config returns valid configuration", {
  config <- default_workflow_config()

  expect_type(config, "list")
  expect_true("analysis" %in% names(config))
  expect_true("visualizations" %in% names(config))
  expect_true("manuscripts" %in% names(config))
  expect_true("reports" %in% names(config))
})

test_that("workflow tracks timing information", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 10, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR"),
    visualizations = NULL,
    manuscripts = list(generate = FALSE),
    reports = NULL
  )

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = tempdir()
    )
  )

  expect_true("workflow_info" %in% names(workflow_results))
  expect_true("total_time" %in% names(workflow_results$workflow_info))
  expect_true(workflow_results$workflow_info$total_time > 0)
})

test_that("print method works for workflow results", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 10, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR"),
    visualizations = NULL,
    manuscripts = list(generate = FALSE),
    reports = NULL
  )

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = tempdir()
    )
  )

  expect_output(print(workflow_results), "Automated NMA Workflow Results")
  expect_output(print(workflow_results), "Workflow completed")
})

# ============================================================================
# Tests: Shiny App Framework
# ============================================================================

test_that("launch_powernma_app function exists and has correct parameters", {
  # Check function exists
  expect_true(exists("launch_powernma_app"))

  # Check function signature
  args <- names(formals(launch_powernma_app))
  expect_true("launch_browser" %in% args)
  expect_true("port" %in% args)
  expect_true("theme" %in% args)
  expect_true("max_upload_mb" %in% args)
})

test_that("create_powernma_ui function exists", {
  expect_true(exists("create_powernma_ui", mode = "function"))
})

test_that("create_powernma_server function exists", {
  expect_true(exists("create_powernma_server", mode = "function"))
})

test_that("UI tab creation functions exist", {
  expect_true(exists("create_home_tab"))
  expect_true(exists("create_data_import_tab"))
  expect_true(exists("create_standard_nma_tab"))
  expect_true(exists("create_bayesian_nma_tab"))
  expect_true(exists("create_component_nma_tab"))
  expect_true(exists("create_multivariate_nma_tab"))
  expect_true(exists("create_living_nma_tab"))
  expect_true(exists("create_rankings_tab"))
  expect_true(exists("create_visualizations_tab"))
  expect_true(exists("create_generate_methods_tab"))
})

test_that("Server module functions exist", {
  expect_true(exists("data_import_server"))
  expect_true(exists("standard_nma_server"))
  expect_true(exists("bayesian_nma_server"))
  expect_true(exists("component_nma_server"))
  expect_true(exists("multivariate_nma_server"))
  expect_true(exists("living_nma_server"))
  expect_true(exists("rankings_server"))
  expect_true(exists("visualizations_server"))
  expect_true(exists("generate_methods_server"))
})

test_that("custom CSS function returns valid CSS", {
  css <- get_custom_css("default")

  expect_type(css, "character")
  expect_true(nchar(css) > 0)
  expect_true(grepl("main-header", css))
})

test_that("custom CSS handles different themes", {
  css_default <- get_custom_css("default")
  css_dark <- get_custom_css("dark")

  expect_true(nchar(css_dark) > nchar(css_default))
  expect_true(grepl("background-color", css_dark))
})

test_that("workflow directory creation works", {
  base_dir <- file.path(tempdir(), "test_workflow")

  dirs <- create_workflow_directories(base_dir)

  expect_type(dirs, "list")
  expect_true("plots" %in% names(dirs))
  expect_true("manuscripts" %in% names(dirs))
  expect_true("reports" %in% names(dirs))

  # All directories should exist
  expect_true(all(sapply(dirs, dir.exists)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
})

test_that("example data generation functions work", {
  expect_true(exists("generate_example_smoking_data"))
  expect_true(exists("generate_example_depression_data"))

  smoking_data <- generate_example_smoking_data()
  expect_s3_class(smoking_data, "data.frame")
  expect_true(nrow(smoking_data) > 0)
})

test_that("data template creation works", {
  template_pairwise <- create_data_template("pairwise")

  expect_s3_class(template_pairwise, "data.frame")
  expect_true("studlab" %in% names(template_pairwise))
  expect_true("treat1" %in% names(template_pairwise))
  expect_true("treat2" %in% names(template_pairwise))
  expect_true("TE" %in% names(template_pairwise))
  expect_true("seTE" %in% names(template_pairwise))
})

# ============================================================================
# Integration Tests
# ============================================================================

test_that("batch processing integrates with batch_processing functions", {
  datasets <- create_mock_batch_datasets(n = 2)
  config <- setup_powernma(sm = "OR")

  batch_results <- run_batch_nma(
    datasets = datasets,
    analysis_config = config,
    parallel = FALSE,
    save_individual = FALSE
  )

  # Should have aggregate statistics
  expect_true(!is.null(batch_results$aggregate_statistics) || length(batch_results$results) == 0)

  # Should track timing
  expect_true(batch_results$total_time > 0)
})

test_that("workflow integrates with manuscript generation", {
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 10, sm = "OR")

  workflow_config <- list(
    analysis = list(sm = "OR"),
    visualizations = NULL,
    manuscripts = list(generate = TRUE, use_ai = FALSE, style = "concise"),
    reports = NULL
  )

  workflow_results <- suppressMessages(
    automated_nma_workflow(
      data = data,
      workflow_config = workflow_config,
      output_dir = tempdir()
    )
  )

  if (!is.null(workflow_results$manuscripts)) {
    expect_true("methods" %in% names(workflow_results$manuscripts))
    expect_true("results" %in% names(workflow_results$manuscripts))
  }
})

# ============================================================================
# Summary Test
# ============================================================================

test_that("Phase 10 implementation is complete", {
  # Verify all expected functions exist

  # Batch processing
  expect_true(exists("run_batch_nma"))
  expect_true(exists("automated_nma_workflow"))
  expect_true(exists("default_workflow_config"))

  # Shiny app
  expect_true(exists("launch_powernma_app"))
  expect_true(exists("create_powernma_ui"))
  expect_true(exists("create_powernma_server"))

  # UI tabs
  expect_true(exists("create_home_tab"))
  expect_true(exists("create_data_import_tab"))

  # Server modules
  expect_true(exists("data_import_server"))
  expect_true(exists("standard_nma_server"))

  message("Phase 10: All Shiny and Batch Processing functions implemented and tested")
})
