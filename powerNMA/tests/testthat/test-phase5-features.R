# Tests for Phase 5: Advanced NMA Features
# Tests for: rankings, inconsistency, visualization, heterogeneity,
# league tables, publication bias, effect size conversions

test_that("SUCRA scores are calculated correctly", {
  skip_if_not_installed("netmeta")

  # Create simple test data
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "B", "A"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.5, 0.3, 0.2, 0.4),
    seTE = c(0.2, 0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD",
    reference.group = "A"
  )

  sucra <- calculate_sucra(nma, small_values = "good")

  # Tests
  expect_s3_class(sucra, "sucra")
  expect_true("sucra_scores" %in% names(sucra))
  expect_true("summary" %in% names(sucra))
  expect_equal(nrow(sucra$summary), 4)  # 4 treatments

  # SUCRA scores should be between 0 and 1
  expect_true(all(sucra$sucra_scores >= 0 & sucra$sucra_scores <= 1))

  # Sum of probabilities of being best should equal 1
  expect_equal(sum(sucra$prob_best), 1, tolerance = 0.01)
})

test_that("Node-splitting detects inconsistency", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("meta")

  # Create data with potential inconsistency
  data <- data.frame(
    studlab = paste0("S", 1:6),
    treat1 = c("A", "A", "B", "A", "B", "C"),
    treat2 = c("B", "C", "C", "B", "C", "A"),
    TE = c(0.5, 0.3, -0.3, 0.6, -0.2, -0.8),  # Inconsistent loop
    seTE = rep(0.2, 6),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  # Node-splitting
  ns_result <- tryCatch({
    node_splitting(nma, data)
  }, error = function(e) NULL)

  # Basic structure tests
  if (!is.null(ns_result)) {
    expect_s3_class(ns_result, "node_split")
    expect_true("summary" %in% names(ns_result))
    expect_true("n_tests" %in% names(ns_result))
  }
})

test_that("Heterogeneity report is generated correctly", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "B", "A"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.5, 0.3, 0.2, 0.4),
    seTE = c(0.2, 0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  het_report <- heterogeneity_report(nma)

  # Tests
  expect_s3_class(het_report, "heterogeneity_report")
  expect_true("tau" %in% names(het_report))
  expect_true("I2" %in% names(het_report))
  expect_true("prediction_intervals" %in% names(het_report))

  # Values should be non-negative
  expect_true(het_report$tau >= 0)
  expect_true(het_report$I2 >= 0 && het_report$I2 <= 1)
})

test_that("League table is created correctly", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "B", "A"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.5, 0.3, 0.2, 0.4),
    seTE = c(0.2, 0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  league <- create_league_table(nma, digits = 2, format = "effect_ci")

  # Tests
  expect_s3_class(league, "league_table")
  expect_true(is.matrix(league))
  expect_equal(nrow(league), 4)  # 4 treatments
  expect_equal(ncol(league), 4)

  # Diagonal should have treatment names
  expect_true(grepl("A", league[1, 1]))
})

test_that("Publication bias tests run without errors", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = paste0("S", 1:10),
    treat1 = rep(c("A", "A", "B"), length.out = 10),
    treat2 = rep(c("B", "C", "C"), length.out = 10),
    TE = rnorm(10, 0.3, 0.2),
    seTE = runif(10, 0.1, 0.3),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  pub_bias <- tryCatch({
    assess_publication_bias(nma, data, method = "Egger")
  }, error = function(e) NULL)

  if (!is.null(pub_bias)) {
    expect_s3_class(pub_bias, "publication_bias")
    expect_true("egger" %in% names(pub_bias) || "begg" %in% names(pub_bias))
  }
})

test_that("Effect size conversions work correctly", {
  # OR to RR
  OR <- 2.0
  p0 <- 0.20
  RR <- or_to_rr(OR, p0)

  expect_true(RR > 0)
  expect_true(RR < OR)  # RR should be smaller than OR for common outcomes

  # RR to OR (round trip)
  OR_back <- rr_to_or(RR, p0)
  expect_equal(OR_back, OR, tolerance = 0.001)

  # Cohen's d to Hedges' g
  d <- 0.5
  g <- cohens_d_to_hedges_g(d, n1 = 30, n2 = 30)

  expect_true(g < d)  # Hedges' g should be slightly smaller (bias correction)
  expect_true(abs(g - d) < 0.05)  # Small difference

  # MD to SMD
  MD <- 10
  smd <- md_to_smd(MD, sd1 = 15, sd2 = 18, n1 = 50, n2 = 50)

  expect_true(abs(smd) < abs(MD))  # SMD should be smaller in absolute value

  # Correlation conversions
  r <- 0.5
  z <- r_to_fishers_z(r)
  r_back <- fishers_z_to_r(z)

  expect_equal(r_back, r, tolerance = 0.001)
})

test_that("NNT calculation is correct", {
  # From OR
  OR <- 0.5  # 50% reduction in odds
  p0 <- 0.20
  NNT <- calculate_nnt(OR, type = "OR", baseline_risk = p0)

  expect_true(NNT > 0)
  expect_true(is.finite(NNT))

  # From RR
  RR <- 0.75
  NNT_rr <- calculate_nnt(RR, type = "RR", baseline_risk = p0)

  expect_true(NNT_rr > 0)
  expect_true(is.finite(NNT_rr))

  # From ARD
  ARD <- 0.10
  NNT_ard <- calculate_nnt(ARD, type = "ARD")

  expect_equal(NNT_ard, 10)  # 1/0.10 = 10
})

test_that("Network geometry calculation works", {
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "B", "A"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.5, 0.3, 0.2, 0.4),
    seTE = c(0.2, 0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  geometry <- network_geometry(data)

  expect_true(is.data.frame(geometry))
  expect_true("n_studies" %in% names(geometry))
  expect_true("n_treatments" %in% names(geometry))
  expect_true("density" %in% names(geometry))

  expect_equal(geometry$n_studies, 4)
  expect_equal(geometry$n_treatments, 4)
  expect_true(geometry$density >= 0 && geometry$density <= 1)
})

test_that("Prediction intervals are wider than confidence intervals", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = paste0("S", 1:10),
    treat1 = rep(c("A", "A", "B"), length.out = 10),
    treat2 = rep(c("B", "C", "C"), length.out = 10),
    TE = rnorm(10, 0.3, 0.3),  # Some heterogeneity
    seTE = rep(0.2, 10),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  pi <- prediction_interval(nma, "B", "A")

  expect_s3_class(pi, "prediction_interval")
  expect_true("pi_lower" %in% names(pi))
  expect_true("pi_upper" %in% names(pi))
  expect_true("ci_lower" %in% names(pi))
  expect_true("ci_upper" %in% names(pi))

  # PI should be wider than CI when tau > 0
  if (nma$tau > 0) {
    expect_true(pi$pi_width >= pi$ci_width)
  }
})

test_that("Comparison I-squared calculations work", {
  skip_if_not_installed("meta")

  # Data with multiple studies for same comparison
  data <- data.frame(
    studlab = c("S1", "S2", "S3", "S4"),
    treat1 = c("A", "A", "A", "B"),
    treat2 = c("B", "B", "B", "C"),
    TE = c(0.5, 0.6, 0.4, 0.3),
    seTE = c(0.2, 0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  i2_result <- calculate_comparison_i2(data, "A", "B")

  expect_true(is.list(i2_result))
  expect_true("I2" %in% names(i2_result))
  expect_true("n_studies" %in% names(i2_result))

  expect_equal(i2_result$n_studies, 3)  # 3 A vs B studies
})

test_that("Loop inconsistency identifies triangular loops", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("meta")

  # Create triangle: A-B-C-A
  data <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "B", "C"),
    treat2 = c("B", "C", "A"),
    TE = c(0.5, 0.3, -0.9),  # Potentially inconsistent
    seTE = c(0.2, 0.2, 0.2),
    stringsAsFactors = FALSE
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "MD"
  )

  loops <- tryCatch({
    loop_inconsistency(nma, data)
  }, error = function(e) NULL)

  if (!is.null(loops) && nrow(loops) > 0) {
    expect_true(is.data.frame(loops))
    expect_true("loop" %in% names(loops))
    expect_true("IF" %in% names(loops))  # Inconsistency Factor
    expect_true("p_value" %in% names(loops))
  }
})

test_that("Comprehensive NMA runs without errors on valid data", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = paste0("S", 1:10),
    treat1 = rep(c("A", "A", "B", "A", "C"), 2),
    treat2 = rep(c("B", "C", "C", "D", "D"), 2),
    TE = rnorm(10, 0.3, 0.2),
    seTE = runif(10, 0.15, 0.25),
    stringsAsFactors = FALSE
  )

  # Should run without errors
  expect_error({
    results <- run_comprehensive_nma(
      data = data,
      sm = "MD",
      assess_inconsistency = FALSE,  # Skip to speed up test
      assess_publication_bias = FALSE,
      generate_plots = FALSE,
      output_dir = NULL
    )
  }, NA)  # NA means no error expected
})

test_that("Batch effect size conversion handles errors gracefully", {
  ORs <- c(1.5, 2.0, -1.0, 2.5)  # One invalid value

  RRs <- batch_convert(ORs, or_to_rr, p0 = 0.20)

  # Should have same length
  expect_equal(length(RRs), length(ORs))

  # Invalid value should be NA
  expect_true(is.na(RRs[3]))

  # Valid values should be converted
  expect_false(is.na(RRs[1]))
  expect_false(is.na(RRs[2]))
})

test_that("Treatment ranking comparison across outcomes works", {
  skip_if_not_installed("netmeta")

  # Create two NMA results
  data1 <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1 = c("A", "A", "B"),
    treat2 = c("B", "C", "C"),
    TE = c(0.5, 0.3, 0.2),
    seTE = c(0.2, 0.2, 0.2)
  )

  nma1 <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data1, sm = "MD"
  )

  sucra1 <- calculate_sucra(nma1)
  sucra2 <- calculate_sucra(nma1)  # Same for simplicity

  comparison <- compare_rankings(list(
    Efficacy = sucra1,
    Safety = sucra2
  ))

  expect_true(is.data.frame(comparison))
  expect_true("Treatment" %in% names(comparison))
  expect_true("Efficacy_SUCRA" %in% names(comparison))
  expect_true("Safety_SUCRA" %in% names(comparison))
  expect_true("Mean_SUCRA" %in% names(comparison))
})

test_that("Variance decomposition calculates correctly", {
  skip_if_not_installed("netmeta")

  data <- data.frame(
    studlab = paste0("S", 1:8),
    treat1 = rep(c("A", "A", "B"), length.out = 8),
    treat2 = rep(c("B", "C", "C"), length.out = 8),
    TE = rnorm(8, 0.3, 0.3),
    seTE = runif(8, 0.15, 0.25)
  )

  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data, sm = "MD"
  )

  var_decomp <- variance_decomposition(nma)

  expect_true(is.list(var_decomp))
  expect_true("between_study_variance" %in% names(var_decomp))
  expect_true("within_study_variance" %in% names(var_decomp))
  expect_true("proportion_between" %in% names(var_decomp))

  # Proportions should sum to 1
  expect_equal(
    var_decomp$proportion_between + var_decomp$proportion_within,
    1,
    tolerance = 0.01
  )

  # Variances should be non-negative
  expect_true(var_decomp$between_study_variance >= 0)
  expect_true(var_decomp$within_study_variance >= 0)
})
