context("Phase 7: Breakthrough 2024-2025 Methods")

library(testthat)
library(powerNMA)

# Generate test data
set.seed(123)
test_data <- data.frame(
  studlab = rep(paste0("Study", 1:10), each = 3),
  treat1 = rep(c("Placebo", "Placebo", "A"), 10),
  treat2 = rep(c("A", "B", "B"), 10),
  TE = rnorm(30, 0.5, 0.3),
  seTE = runif(30, 0.1, 0.3),
  stringsAsFactors = FALSE
)

# ==============================================================================
# Component Network Meta-Analysis Tests
# ==============================================================================

test_that("Component matrix creation works", {
  treatments <- c("Placebo", "A", "B", "A+B")
  components <- c("Component_A", "Component_B")

  comp_matrix <- create_component_matrix(treatments, components, interactive = FALSE)

  expect_is(comp_matrix, "matrix")
  expect_equal(nrow(comp_matrix), 4)
  expect_equal(ncol(comp_matrix), 2)
  expect_equal(rownames(comp_matrix), treatments)
  expect_equal(colnames(comp_matrix), components)
})

test_that("Additive CNMA runs successfully", {
  skip_if_not_installed("netmeta")

  # Create simple component structure
  comp_matrix <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1
  ), nrow = 4, byrow = TRUE)

  rownames(comp_matrix) <- c("Placebo", "A", "B", "A+B")
  colnames(comp_matrix) <- c("Comp_A", "Comp_B")

  # Modify test data for component structure
  cnma_data <- test_data
  cnma_data$treat1 <- factor(cnma_data$treat1, levels = rownames(comp_matrix))
  cnma_data$treat2 <- factor(cnma_data$treat2, levels = rownames(comp_matrix))

  expect_error({
    result <- additive_cnma(
      data = cnma_data,
      comp_matrix = comp_matrix,
      reference = "Placebo"
    )
  }, NA)  # Expect no error
})

test_that("Interaction CNMA handles interaction terms", {
  skip_if_not_installed("netmeta")

  comp_matrix <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1
  ), nrow = 4, byrow = TRUE)

  rownames(comp_matrix) <- c("Placebo", "A", "B", "A+B")
  colnames(comp_matrix) <- c("Comp_A", "Comp_B")

  cnma_data <- test_data
  cnma_data$treat1 <- factor(cnma_data$treat1, levels = rownames(comp_matrix))
  cnma_data$treat2 <- factor(cnma_data$treat2, levels = rownames(comp_matrix))

  expect_error({
    result <- interaction_cnma(
      data = cnma_data,
      comp_matrix = comp_matrix,
      interactions = c("Comp_A:Comp_B"),
      reference = "Placebo"
    )
  }, NA)
})

# ==============================================================================
# Multivariate NMA Tests
# ==============================================================================

test_that("Multivariate NMA data preparation works", {
  # Create two outcome datasets
  outcome1 <- test_data
  outcome2 <- test_data
  outcome2$TE <- outcome2$TE + rnorm(nrow(outcome2), 0, 0.1)

  mvnma_data <- prepare_mvnma_data(
    data_list = list(outcome1, outcome2),
    outcome_names = c("Efficacy", "Safety")
  )

  expect_is(mvnma_data, "mvnma_data")
  expect_equal(mvnma_data$n_outcomes, 2)
  expect_equal(mvnma_data$outcome_names, c("Efficacy", "Safety"))
  expect_is(mvnma_data$correlation, "matrix")
})

test_that("Multivariate NMA fitting works", {
  skip_if_not_installed("netmeta")

  outcome1 <- test_data
  outcome2 <- test_data
  outcome2$TE <- outcome2$TE + rnorm(nrow(outcome2), 0, 0.1)

  mvnma_data <- prepare_mvnma_data(
    data_list = list(outcome1, outcome2),
    outcome_names = c("Efficacy", "Safety")
  )

  expect_error({
    result <- fit_mvnma(
      mvnma_data = mvnma_data,
      reference = "Placebo",
      sm = c("MD", "MD"),
      method = "separate"
    )
  }, NA)

  expect_is(result, "mvnma_result")
  expect_equal(length(result$separate_nmas), 2)
})

test_that("Net benefit calculation works", {
  skip_if_not_installed("netmeta")

  outcome1 <- test_data
  outcome2 <- test_data

  mvnma_data <- prepare_mvnma_data(
    data_list = list(outcome1, outcome2),
    outcome_names = c("Efficacy", "Safety")
  )

  result <- fit_mvnma(mvnma_data, method = "separate")

  net_benefit <- calculate_net_benefit(
    mvnma_result = result,
    weights = c(Efficacy = 0.6, Safety = 0.4),
    directions = c(Efficacy = "higher", Safety = "lower")
  )

  expect_is(net_benefit, "data.frame")
  expect_true("Net_Benefit" %in% names(net_benefit))
  expect_true("Treatment" %in% names(net_benefit))
})

# ==============================================================================
# Living NMA Tests
# ==============================================================================

test_that("Living NMA initialization works", {
  skip_if_not_installed("netmeta")

  temp_dir <- tempdir()
  project_name <- "test_living_nma"

  expect_error({
    project <- initialize_living_nma(
      baseline_data = test_data,
      project_name = project_name,
      update_frequency = "quarterly",
      output_dir = file.path(temp_dir, project_name)
    )
  }, NA)

  expect_is(project, "living_nma_project")
  expect_equal(project$current_version, 1)
  expect_true(dir.exists(project$output_dir))

  # Cleanup
  unlink(file.path(temp_dir, project_name), recursive = TRUE)
})

test_that("Living NMA update works", {
  skip_if_not_installed("netmeta")

  temp_dir <- tempdir()
  project_name <- "test_living_update"

  project <- initialize_living_nma(
    baseline_data = test_data,
    project_name = project_name,
    output_dir = file.path(temp_dir, project_name)
  )

  # Create new studies
  new_data <- data.frame(
    studlab = rep(paste0("NewStudy", 1:2), each = 3),
    treat1 = rep(c("Placebo", "Placebo", "A"), 2),
    treat2 = rep(c("A", "B", "B"), 2),
    TE = rnorm(6, 0.5, 0.3),
    seTE = runif(6, 0.1, 0.3),
    stringsAsFactors = FALSE
  )

  expect_error({
    updated_project <- update_living_nma(
      project = project,
      new_data = new_data,
      force_update = TRUE,
      update_note = "Test update"
    )
  }, NA)

  expect_equal(updated_project$current_version, 2)

  # Cleanup
  unlink(file.path(temp_dir, project_name), recursive = TRUE)
})

test_that("Version history retrieval works", {
  skip_if_not_installed("netmeta")

  temp_dir <- tempdir()
  project_name <- "test_version_history"

  project <- initialize_living_nma(
    baseline_data = test_data,
    project_name = project_name,
    output_dir = file.path(temp_dir, project_name)
  )

  history <- get_version_history(project)

  expect_is(history, "data.frame")
  expect_true("Version" %in% names(history))
  expect_true("N_Studies" %in% names(history))
  expect_equal(nrow(history), 1)

  # Cleanup
  unlink(file.path(temp_dir, project_name), recursive = TRUE)
})

# ==============================================================================
# Personalized Prediction Tests
# ==============================================================================

test_that("IPD meta-regression fitting works", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("netmeta")

  # Create IPD data
  ipd <- data.frame(
    study = rep(1:5, each = 100),
    patient_id = 1:500,
    treatment = sample(c("A", "B", "C"), 500, replace = TRUE),
    outcome = rnorm(500),
    age = rnorm(500, 50, 10),
    baseline_severity = rnorm(500, 0, 1),
    stringsAsFactors = FALSE
  )

  expect_error({
    model <- fit_ipd_metaregression(
      ipd_data = ipd,
      outcome_var = "outcome",
      treatment_var = "treatment",
      study_var = "study",
      effect_modifiers = c("age", "baseline_severity"),
      reference = "A",
      family = "gaussian"
    )
  }, NA)

  expect_is(model, "ipd_nma")
  expect_equal(model$reference, "A")
})

test_that("Personalized predictions work", {
  skip_if_not_installed("lme4")

  ipd <- data.frame(
    study = rep(1:5, each = 100),
    treatment = sample(c("A", "B", "C"), 500, replace = TRUE),
    outcome = rnorm(500),
    age = rnorm(500, 50, 10),
    baseline_severity = rnorm(500, 0, 1),
    stringsAsFactors = FALSE
  )

  model <- fit_ipd_metaregression(
    ipd_data = ipd,
    outcome_var = "outcome",
    treatment_var = "treatment",
    study_var = "study",
    effect_modifiers = c("age", "baseline_severity"),
    family = "gaussian"
  )

  new_patient <- data.frame(
    age = 65,
    baseline_severity = 1.5
  )

  expect_error({
    predictions <- predict_personalized_effects(
      ipd_model = model,
      newdata = new_patient
    )
  }, NA)

  expect_is(predictions, "data.frame")
  expect_true("Treatment" %in% names(predictions))
  expect_true("Effect_vs_Reference" %in% names(predictions))
})

test_that("Best treatment identification works", {
  skip_if_not_installed("lme4")

  ipd <- data.frame(
    study = rep(1:5, each = 100),
    treatment = sample(c("A", "B", "C"), 500, replace = TRUE),
    outcome = rnorm(500),
    age = rnorm(500, 50, 10),
    baseline_severity = rnorm(500, 0, 1),
    stringsAsFactors = FALSE
  )

  model <- fit_ipd_metaregression(
    ipd_data = ipd,
    outcome_var = "outcome",
    treatment_var = "treatment",
    study_var = "study",
    effect_modifiers = c("age", "baseline_severity"),
    family = "gaussian"
  )

  new_patient <- data.frame(
    age = 65,
    baseline_severity = 1.5
  )

  expect_error({
    best <- identify_best_treatment(
      ipd_model = model,
      newdata = new_patient,
      outcome_direction = "maximize"
    )
  }, NA)

  expect_is(best, "data.frame")
  expect_true("Recommended_Treatment" %in% names(best))
  expect_true("Expected_Effect" %in% names(best))
})

# ==============================================================================
# Value of Information Tests
# ==============================================================================

test_that("Threshold analysis works", {
  skip_if_not_installed("netmeta")

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data,
    reference.group = "Placebo"
  )

  costs <- data.frame(
    treatment = c("Placebo", "A", "B"),
    cost = c(0, 1000, 1500)
  )

  expect_error({
    threshold_result <- threshold_analysis(
      nma_result = nma,
      cost_data = costs,
      threshold_range = seq(0, 10000, by = 1000)
    )
  }, NA)

  expect_is(threshold_result, "threshold_analysis")
  expect_true("threshold_results" %in% names(threshold_result))
  expect_true("switch_points" %in% names(threshold_result))
})

test_that("EVPI calculation works", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("MASS")

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data,
    reference.group = "Placebo"
  )

  costs <- data.frame(
    treatment = c("Placebo", "A", "B"),
    cost = c(0, 1000, 1500)
  )

  expect_error({
    evpi_result <- calculate_evpi(
      nma_result = nma,
      cost_data = costs,
      threshold = 5000,
      n_sim = 100  # Small for testing
    )
  }, NA)

  expect_is(evpi_result, "evpi")
  expect_true("evpi" %in% names(evpi_result))
  expect_is(evpi_result$evpi, "numeric")
})

test_that("Population EVPI calculation works", {
  evpi_result <- list(
    evpi = 100,
    threshold = 5000,
    n_sim = 100
  )
  class(evpi_result) <- "evpi"

  expect_error({
    pop_evpi <- calculate_population_evpi(
      evpi_result = evpi_result,
      population_size = 10000,
      time_horizon = 5,
      discount_rate = 0.035
    )
  }, NA)

  expect_is(pop_evpi, "numeric")
  expect_true(pop_evpi > evpi_result$evpi)  # Should be larger than per-patient
})

# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("Phase 7 methods integrate with each other", {
  skip_if_not_installed("netmeta")

  # This test ensures Phase 7 methods can work together
  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data,
    reference.group = "Placebo"
  )

  # Should be able to run multiple Phase 7 analyses on same data
  expect_error({
    # Multivariate setup
    mvnma_data <- prepare_mvnma_data(
      data_list = list(test_data, test_data),
      outcome_names = c("Primary", "Secondary")
    )

    # Living NMA setup
    temp_dir <- tempdir()
    living_project <- initialize_living_nma(
      baseline_data = test_data,
      project_name = "integration_test",
      output_dir = file.path(temp_dir, "integration_test")
    )

    # Cleanup
    unlink(file.path(temp_dir, "integration_test"), recursive = TRUE)
  }, NA)
})

test_that("All Phase 7 classes have print methods", {
  # Check that all major classes have print methods
  expect_true(exists("print.cnma_additive"))
  expect_true(exists("print.cnma_interaction"))
  expect_true(exists("print.mvnma_result"))
  expect_true(exists("print.living_nma_project"))
  expect_true(exists("print.ipd_nma"))
  expect_true(exists("print.threshold_analysis"))
  expect_true(exists("print.evpi"))
})

test_that("All Phase 7 classes have plot methods where appropriate", {
  # Check that plotting classes have plot methods
  expect_true(exists("plot.cnma_additive"))
  expect_true(exists("plot.cnma_interaction"))
  expect_true(exists("plot.mvnma_result"))
  expect_true(exists("plot.threshold_analysis"))
})

# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("Phase 7 functions handle errors gracefully", {
  # Component NMA with invalid matrix
  expect_error(
    create_component_matrix(treatments = c(), components = c("A")),
    NA  # Should handle gracefully
  )

  # Multivariate NMA with mismatched data
  expect_error(
    prepare_mvnma_data(data_list = "not_a_list", outcome_names = c("A")),
    "data_list must be a list"
  )

  # Living NMA with invalid directory
  expect_error(
    load_living_nma("/nonexistent/path"),
    "not found"
  )
})

# ==============================================================================
# Performance Tests
# ==============================================================================

test_that("Phase 7 methods complete in reasonable time", {
  skip_if_not_installed("netmeta")
  skip_on_cran()  # Performance tests may take longer

  start_time <- Sys.time()

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = test_data,
    reference.group = "Placebo"
  )

  costs <- data.frame(
    treatment = c("Placebo", "A", "B"),
    cost = c(0, 1000, 1500)
  )

  threshold_result <- threshold_analysis(
    nma_result = nma,
    cost_data = costs,
    threshold_range = seq(0, 10000, by = 1000)
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  expect_true(elapsed < 30)  # Should complete in under 30 seconds
})
