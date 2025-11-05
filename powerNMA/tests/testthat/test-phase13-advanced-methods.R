# Tests for Phase 13: Advanced Statistical Methods from Latest Journals

context("Phase 13: Cutting-Edge Statistical Methods")

# Test 1: Multivariate NMA
test_that("Multivariate NMA works", {
  skip_if_not_installed("netmeta")
  
  mv_data <- simulate_multivariate_nma_data(
    n_studies = 10,
    n_outcomes = 3,
    correlation = 0.5
  )
  
  expect_true(is.list(mv_data))
  expect_equal(length(mv_data$datasets), 3)
  
  mv_result <- run_multivariate_nma(
    data_list = mv_data$datasets,
    outcome_names = c("Efficacy", "Safety_AE", "Safety_SAE"),
    outcome_types = c("benefit", "harm", "harm"),
    method = "frequentist"
  )
  
  expect_s3_class(mv_result, "multivariate_nma")
  expect_equal(mv_result$n_outcomes, 3)
  expect_true(!is.null(mv_result$benefit_risk))
})

# Test 2: ROPE Analysis
test_that("ROPE analysis works", {
  data <- simulate_nma_data(n_studies = 15)
  nma_result <- netmeta::netmeta(
    TE = data$TE, seTE = data$seTE,
    treat1 = data$treat1, treat2 = data$treat2,
    studlab = data$studlab
  )
  
  rope_result <- analyze_rope(
    nma_result = nma_result,
    mcid = 0.5,
    method = "frequentist"
  )
  
  expect_s3_class(rope_result, "rope_analysis")
  expect_true(!is.null(rope_result$decisions))
  expect_output(print(rope_result), "ROPE")
})

# Test 3: Non-inferiority Testing
test_that("Non-inferiority testing works", {
  data <- simulate_nma_data(n_studies = 10)
  nma_result <- netmeta::netmeta(
    TE = data$TE, seTE = data$seTE,
    treat1 = data$treat1, treat2 = data$treat2,
    studlab = data$studlab
  )
  
  test_result <- test_noninferiority(
    nma_result = nma_result,
    test_type = "non_inferiority",
    margin = 0.3,
    test_treatment = "T2",
    reference_treatment = "T1",
    method = "frequentist"
  )
  
  expect_s3_class(test_result, "noninferiority_test")
  expect_true(!is.null(test_result$conclusion))
})

# Test 4: MCID Assessment
test_that("MCID assessment works", {
  data <- simulate_nma_data(n_studies = 12)
  nma_result <- netmeta::netmeta(
    TE = data$TE, seTE = data$seTE,
    treat1 = data$treat1, treat2 = data$treat2,
    studlab = data$studlab
  )
  
  mcid_result <- assess_clinical_importance(
    nma_result = nma_result,
    mcid = 0.5,
    mcid_source = "literature"
  )
  
  expect_s3_class(mcid_result, "mcid_assessment")
  expect_true(!is.null(mcid_result$summary))
})

# Test 5: Transitivity Assessment
test_that("Transitivity assessment works", {
  data <- simulate_nma_data(n_studies = 20)
  data$age_mean <- rnorm(nrow(data), 60, 10)
  data$comparison <- paste(data$treat1, "vs", data$treat2)
  
  transit_result <- assess_transitivity(
    nma_data = data,
    effect_modifiers = c("age_mean"),
    statistical_test = "kruskal_wallis"
  )
  
  expect_s3_class(transit_result, "transitivity_assessment")
  expect_true(!is.null(transit_result$transitivity_score))
})

# Test 6: Treatment Selection Model
test_that("Treatment selection model works", {
  data <- simulate_nma_data(n_studies = 15)
  nma_result <- netmeta::netmeta(
    TE = data$TE, seTE = data$seTE,
    treat1 = data$treat1, treat2 = data$treat2,
    studlab = data$studlab
  )
  
  model <- build_treatment_selection_model(
    nma_result = nma_result,
    patient_characteristics = c("age", "sex"),
    method = "regression"
  )
  
  expect_s3_class(model, "treatment_selection_model")
  
  prediction <- predict_optimal_treatment(
    model = model,
    new_patient = data.frame(age = 65, sex = "M")
  )
  
  expect_true(!is.null(prediction$optimal_treatment))
})

# Test 7: RCT + RWE Integration
test_that("RCT-RWE integration works", {
  rct_data <- simulate_nma_data(n_studies = 10)
  rwe_data <- simulate_nma_data(n_studies = 5)
  
  integrated <- integrate_rct_rwe(
    rct_data = rct_data,
    rwe_data = rwe_data,
    bias_adjustment = TRUE,
    method = "hierarchical"
  )
  
  expect_s3_class(integrated, "integrated_rct_rwe")
  expect_true(!is.null(integrated$integrated_model))
})

# Test 8: Sequential NMA
test_that("Sequential NMA works", {
  # Create sequential data updates
  data_updates <- list(
    simulate_nma_data(n_studies = 5),
    simulate_nma_data(n_studies = 10),
    simulate_nma_data(n_studies = 15)
  )
  
  seq_result <- run_sequential_nma(
    data_updates = data_updates,
    alpha = 0.05,
    spending_function = "obrien_fleming"
  )
  
  expect_s3_class(seq_result, "sequential_nma")
  expect_equal(length(seq_result$sequential_results), 3)
})

# Test 9: Benefit-Risk MCDA
test_that("Benefit-risk MCDA works", {
  mv_data <- simulate_multivariate_nma_data(n_studies = 10, n_outcomes = 2)
  
  mv_result <- run_multivariate_nma(
    data_list = mv_data$datasets,
    outcome_names = c("Benefit", "Harm"),
    outcome_types = c("benefit", "harm"),
    method = "frequentist"
  )
  
  mcda_result <- benefit_risk_mcda(
    mv_nma_result = mv_result,
    method = "weighted_sum"
  )
  
  expect_s3_class(mcda_result, "benefit_risk_mcda")
  expect_true(!is.null(mcda_result$optimal_treatment))
})

message("Phase 13 Advanced Methods Tests Complete")
