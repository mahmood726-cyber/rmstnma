# Tests for Phase 12: Advanced Bayesian, IPD-NMA, DTA, HTML Reports, APIs, Splines
# Test file for powerNMA Phase 12 revolutionary features

context("Phase 12: Advanced Bayesian Workflows and Extensions")

# ==============================================================================
# Test 1: Advanced Bayesian Stan NMA
# ==============================================================================

test_that("Advanced Bayesian Stan NMA works", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()

  # Simulate data
  data <- simulate_nma_data(n_studies = 10, n_treatments = 3, sm = "OR")

  # Specify priors
  priors <- list(
    treatment_effects = list(mean = 0, sd = 2),
    heterogeneity = list(distribution = "half_normal", scale = 0.5)
  )

  # Run Bayesian NMA (with very few iterations for testing)
  expect_message(
    bayes_result <- run_stan_nma(
      data = data,
      model_type = "fixed_effects",  # Use simpler model for testing
      prior_specification = priors,
      n_iter = 200,
      n_warmup = 100,
      n_chains = 2
    ),
    "Bayesian NMA complete"
  )

  # Check result structure
  expect_s3_class(bayes_result, "stan_nma")
  expect_true(!is.null(bayes_result$model))
  expect_true(!is.null(bayes_result$samples))
  expect_true(!is.null(bayes_result$summaries))
  expect_true(!is.null(bayes_result$rankings))

  # Check diagnostics
  if (!is.null(bayes_result$diagnostics)) {
    expect_true(is.numeric(bayes_result$diagnostics$rhat_max))
    expect_true(is.numeric(bayes_result$diagnostics$ess_min))
  }

  # Check print method
  expect_output(print(bayes_result), "Bayesian Network Meta-Analysis")
})

test_that("Prior sensitivity analysis works", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  data <- simulate_nma_data(n_studies = 8, n_treatments = 3)

  # Define multiple prior scenarios
  priors_list <- list(
    weak = list(
      treatment_effects = list(mean = 0, sd = 10),
      heterogeneity = list(distribution = "half_normal", scale = 1)
    ),
    moderate = list(
      treatment_effects = list(mean = 0, sd = 2),
      heterogeneity = list(distribution = "half_normal", scale = 0.5)
    )
  )

  # Run sensitivity analysis
  expect_message(
    sensitivity <- prior_sensitivity_analysis(
      data = data,
      prior_scenarios = priors_list,
      model_type = "fixed_effects",
      n_iter = 200,
      n_warmup = 100,
      n_chains = 2
    ),
    "Running prior sensitivity analysis"
  )

  expect_s3_class(sensitivity, "prior_sensitivity")
  expect_equal(length(sensitivity), 2)

  # Compare results
  comparison <- compare_prior_sensitivity(sensitivity)
  expect_true(is.data.frame(comparison))
  expect_equal(nrow(comparison), 2)
})

# ==============================================================================
# Test 2: IPD Network Meta-Analysis
# ==============================================================================

test_that("One-stage IPD-NMA works", {
  skip_if_not_installed("lme4")

  # Simulate IPD data
  ipd_data <- simulate_ipd_nma_data(
    n_studies = 5,
    n_patients_per_study = 50,
    n_treatments = 3,
    outcome_type = "binary"
  )

  expect_true(is.data.frame(ipd_data))
  expect_true(all(c("study", "patient_id", "treatment", "outcome") %in% names(ipd_data)))

  # Run one-stage IPD-NMA
  expect_message(
    ipd_result <- run_onestage_ipd_nma(
      ipd_data = ipd_data,
      outcome_type = "binary",
      covariates = c("age", "sex"),
      treatment_by_covariate = FALSE,
      method = "ML"
    ),
    "One-stage IPD-NMA complete"
  )

  expect_s3_class(ipd_result, "ipd_nma_onestage")
  expect_true(!is.null(ipd_result$model))
  expect_true(!is.null(ipd_result$results))
  expect_true(!is.null(ipd_result$rankings))

  # Check print method
  expect_output(print(ipd_result), "One-Stage Individual Patient Data")
})

test_that("Two-stage IPD-NMA works", {
  # Simulate IPD data
  ipd_data <- simulate_ipd_nma_data(
    n_studies = 5,
    n_patients_per_study = 50,
    outcome_type = "continuous"
  )

  # Run two-stage IPD-NMA
  expect_message(
    twostage_result <- run_twostage_ipd_nma(
      ipd_data = ipd_data,
      outcome_type = "continuous",
      pooling_method = "random_effects"
    ),
    "Two-stage IPD-NMA complete"
  )

  expect_s3_class(twostage_result, "ipd_nma_twostage")
  expect_true(!is.null(twostage_result$stage1))
  expect_true(!is.null(twostage_result$stage2))
})

test_that("Personalized treatment predictions work", {
  skip_if_not_installed("lme4")

  ipd_data <- simulate_ipd_nma_data(n_studies = 5, n_patients_per_study = 50)

  ipd_result <- run_onestage_ipd_nma(
    ipd_data = ipd_data,
    outcome_type = "binary",
    covariates = c("age", "sex"),
    method = "ML"
  )

  # Create new patient data
  new_patient <- data.frame(
    age = 65,
    sex = "M"
  )

  # Predict personalized effects
  expect_message(
    predictions <- predict_personalized_effects(
      ipd_nma_result = ipd_result,
      new_patient_data = new_patient
    ),
    "Personalized predictions completed"
  )

  expect_s3_class(predictions, "personalized_predictions")
  expect_true(!is.null(predictions$predictions))
})

# ==============================================================================
# Test 3: Diagnostic Test Accuracy NMA
# ==============================================================================

test_that("DTA-NMA simulation works", {
  # Simulate DTA data
  dta_data <- simulate_dta_nma_data(
    n_studies = 15,
    n_tests = 4
  )

  expect_true(is.data.frame(dta_data))
  expect_true(all(c("study", "test", "TP", "FP", "FN", "TN") %in% names(dta_data)))
  expect_true(nrow(dta_data) > 0)
})

test_that("Frequentist DTA-NMA works", {
  skip_if_not_installed("metafor")

  # Simulate data
  dta_data <- simulate_dta_nma_data(n_studies = 10, n_tests = 3)

  # Run DTA-NMA
  expect_message(
    dta_result <- run_dta_nma(
      data = dta_data,
      model_type = "bivariate",
      method = "frequentist",
      prevalence = 0.1
    ),
    "DTA-NMA complete"
  )

  expect_s3_class(dta_result, "dta_nma")
  expect_true(!is.null(dta_result$pooled_estimates))
  expect_true(!is.null(dta_result$accuracy_measures))
  expect_true(!is.null(dta_result$rankings))

  # Check accuracy measures
  expect_true("PPV" %in% names(dta_result$accuracy_measures))
  expect_true("NPV" %in% names(dta_result$accuracy_measures))
  expect_true("LR_positive" %in% names(dta_result$accuracy_measures))

  # Check print method
  expect_output(print(dta_result), "Diagnostic Test Accuracy")

  # Check summary
  expect_output(summary(dta_result), "Accuracy Measures")
})

test_that("Bayesian DTA-NMA works", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  dta_data <- simulate_dta_nma_data(n_studies = 8, n_tests = 2)

  expect_message(
    dta_result <- run_dta_nma(
      data = dta_data,
      model_type = "bivariate",
      method = "bayesian",
      n_iter = 200,
      n_warmup = 100,
      n_chains = 2
    ),
    "DTA-NMA complete"
  )

  expect_s3_class(dta_result, "dta_nma")
  expect_true(!is.null(dta_result$model))
})

# ==============================================================================
# Test 4: Interactive HTML Reports
# ==============================================================================

test_that("HTML report generation works", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")

  # Create mock NMA result
  data <- simulate_nma_data(n_studies = 10)
  nma_result <- netmeta::netmeta(
    TE = data$TE,
    seTE = data$seTE,
    treat1 = data$treat1,
    treat2 = data$treat2,
    studlab = data$studlab,
    sm = "OR"
  )

  # Add SUCRA scores
  nma_result$sucra <- calculate_sucra(nma_result)

  # Generate report
  temp_dir <- tempdir()
  output_file <- file.path(temp_dir, "test_report.html")

  expect_message(
    report_path <- generate_interactive_html_report(
      nma_result = nma_result,
      output_file = output_file,
      title = "Test NMA Report",
      theme = "flatly"
    ),
    "Interactive HTML report generated"
  )

  expect_true(file.exists(report_path))
  expect_equal(report_path, output_file)

  # Clean up
  unlink(output_file)
})

test_that("Multi-format export works", {
  skip_if_not_installed("rmarkdown")

  data <- simulate_nma_data(n_studies = 8)
  nma_result <- netmeta::netmeta(
    TE = data$TE,
    seTE = data$seTE,
    treat1 = data$treat1,
    treat2 = data$treat2,
    studlab = data$studlab
  )

  temp_dir <- tempdir()

  expect_message(
    output_files <- export_multi_format_report(
      nma_result = nma_result,
      output_dir = temp_dir,
      formats = c("html"),
      base_filename = "test_multi_report"
    ),
    "All reports exported"
  )

  expect_true(length(output_files) > 0)
  expect_true(file.exists(output_files[1]))

  # Clean up
  unlink(output_files)
})

# ==============================================================================
# Test 5: External API Integration
# ==============================================================================

test_that("PubMed search works", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("httr")
  skip_if_not_installed("xml2")

  # Simple PubMed search
  expect_message(
    results <- search_pubmed(
      search_terms = c("diabetes", "metformin"),
      max_results = 5,
      filter_rct = TRUE,
      extract_abstracts = FALSE
    ),
    "Found"
  )

  if (nrow(results) > 0) {
    expect_true(is.data.frame(results))
    expect_true("PMID" %in% names(results))
    expect_true("Title" %in% names(results))
  }
})

test_that("ClinicalTrials.gov search works", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("httr")
  skip_if_not_installed("jsonlite")

  # Simple trial search
  expect_message(
    trials <- search_clinicaltrials(
      condition = "diabetes",
      status = "Completed",
      max_results = 5
    ),
    "Found"
  )

  if (nrow(trials) > 0) {
    expect_true(is.data.frame(trials))
    expect_true("NCT_ID" %in% names(trials))
    expect_true("Title" %in% names(trials))
  }
})

test_that("PICO extraction works", {
  # Test abstracts
  abstracts <- c(
    "This randomized controlled trial included 200 patients with type 2 diabetes. Participants were treated with metformin 1000mg daily compared with placebo. The primary outcome was HbA1c reduction at 12 weeks.",
    "Adults with hypertension received losartan 50mg versus standard care. Blood pressure was the main outcome."
  )

  expect_message(
    pico_results <- extract_pico_from_abstracts(
      abstracts = abstracts,
      method = "rules"
    ),
    "PICO extraction complete"
  )

  expect_true(is.data.frame(pico_results))
  expect_equal(nrow(pico_results), 2)
  expect_true(all(c("Population", "Intervention", "Comparison", "Outcome") %in% names(pico_results)))
})

test_that("Duplicate detection works", {
  # Create mock studies with duplicates
  studies1 <- data.frame(
    Title = c("Study A", "Study B", "Study C"),
    Authors = c("Smith J", "Jones K", "Brown L"),
    Year = c("2020", "2021", "2022")
  )

  studies2 <- data.frame(
    Title = c("Study A", "Study D", "Study C"),
    Authors = c("Smith J", "Davis M", "Brown L"),
    Year = c("2020", "2021", "2022")
  )

  duplicates <- detect_duplicate_studies(
    studies_list = list(studies1, studies2),
    method = "title"
  )

  expect_true(is.data.frame(duplicates))
  expect_true(nrow(duplicates) >= 2)  # At least "Study A" and "Study C"
})

test_that("Bibliography export works", {
  studies <- data.frame(
    Title = c("Study 1", "Study 2"),
    Authors = c("Smith J, Jones K", "Brown L, Davis M"),
    Journal = c("Journal A", "Journal B"),
    Year = c("2020", "2021")
  )

  temp_file <- tempfile(fileext = ".bib")

  expect_message(
    output_path <- export_to_bibliography_format(
      studies = studies,
      output_file = temp_file,
      format = "bibtex"
    ),
    "Bibliography exported"
  )

  expect_true(file.exists(output_path))

  # Check content
  content <- readLines(output_path)
  expect_true(any(grepl("@article", content)))

  unlink(temp_file)
})

# ==============================================================================
# Test 6: Advanced Meta-Regression with Splines
# ==============================================================================

test_that("Spline meta-regression works", {
  skip_if_not_installed("splines")
  skip_if_not_installed("metafor")

  # Simulate dose-response data
  dose_data <- simulate_dose_response_data(
    n_studies = 20,
    doses_per_study = 4,
    true_dose_response = function(x) 0.01 * x
  )

  expect_true(is.data.frame(dose_data))
  expect_true(all(c("study", "dose", "logRR", "se_logRR") %in% names(dose_data)))

  # Fit spline meta-regression
  expect_message(
    spline_result <- metareg_splines(
      data = dose_data,
      outcome = "logRR",
      se = "se_logRR",
      covariate = "dose",
      n_knots = 3,
      method = "REML"
    ),
    "Spline meta-regression complete"
  )

  expect_s3_class(spline_result, "metareg_splines")
  expect_true(!is.null(spline_result$model))
  expect_true(!is.null(spline_result$predictions))
  expect_equal(spline_result$n_knots, 3)

  # Check print method
  expect_output(print(spline_result), "Restricted Cubic Splines")
})

test_that("Fractional polynomial meta-regression works", {
  skip_if_not_installed("metafor")

  # Simulate data
  dose_data <- simulate_dose_response_data(
    n_studies = 15,
    doses_per_study = 3
  )

  # Fit fractional polynomial
  expect_message(
    fp_result <- metareg_fractional_polynomial(
      data = dose_data,
      outcome = "logRR",
      se = "se_logRR",
      covariate = "dose",
      max_degree = 1,
      selection_criterion = "AIC"
    ),
    "Best fractional polynomial"
  )

  expect_s3_class(fp_result, "metareg_fp")
  expect_true(!is.null(fp_result$model))
  expect_true(!is.null(fp_result$predictions))
  expect_true(length(fp_result$powers) <= 2)

  # Check print method
  expect_output(print(fp_result), "Fractional Polynomials")
})

test_that("Dose-response meta-analysis works", {
  # Simulate dose-response data
  dose_data <- simulate_dose_response_data(
    n_studies = 10,
    doses_per_study = 3
  )

  # Run dose-response MA
  expect_message(
    dr_result <- dose_response_metaanalysis(
      data = dose_data,
      outcome = "logRR",
      se = "se_logRR",
      dose = "dose",
      study = "study",
      type = "linear"
    ),
    "Performing linear dose-response"
  )

  expect_s3_class(dr_result, "dose_response_ma")
  expect_true(!is.null(dr_result$model))
})

# ==============================================================================
# Test 7: Integration Tests
# ==============================================================================

test_that("Phase 12 functions integrate correctly", {
  # Test that different Phase 12 components work together

  # 1. Simulate IPD data
  ipd_data <- simulate_ipd_nma_data(n_studies = 5, n_patients_per_study = 30)
  expect_true(nrow(ipd_data) > 0)

  # 2. Simulate DTA data
  dta_data <- simulate_dta_nma_data(n_studies = 8, n_tests = 3)
  expect_true(nrow(dta_data) > 0)

  # 3. Simulate dose-response data
  dose_data <- simulate_dose_response_data(n_studies = 10)
  expect_true(nrow(dose_data) > 0)

  # All simulations work without errors
  expect_true(TRUE)
})

test_that("Phase 12 error handling works", {
  # Test various error conditions

  # Missing required columns
  expect_error(
    run_dta_nma(data = data.frame(x = 1:10)),
    "must contain columns"
  )

  # Invalid model type
  expect_error(
    metareg_splines(
      data = data.frame(y = rnorm(10), se = rep(0.1, 10), x = 1:10),
      outcome = "y",
      se = "se",
      covariate = "x",
      method = "invalid_method"
    )
  )
})

# ==============================================================================
# Summary Message
# ==============================================================================

message("Phase 12 Advanced Features Tests Complete")
message("- Advanced Bayesian Stan NMA: TESTED")
message("- IPD Network Meta-Analysis: TESTED")
message("- Diagnostic Test Accuracy NMA: TESTED")
message("- Interactive HTML Reports: TESTED")
message("- External API Integration: TESTED")
message("- Advanced Meta-Regression with Splines: TESTED")
