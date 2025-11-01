# ==============================================================================
# Tests for Surrogate Endpoint Network Meta-Analysis
# ==============================================================================

library(testthat)

# Helper function to create test data
create_test_surrogate_data <- function(n_studies = 5, n_treatments = 4,
                                       missing_true_pct = 0.3) {

  # Create treatment pairs
  comparisons <- expand.grid(
    treat1 = 1:n_treatments,
    treat2 = 1:n_treatments
  )
  comparisons <- comparisons[comparisons$treat1 < comparisons$treat2, ]
  comparisons <- comparisons[sample(nrow(comparisons), n_studies), ]

  # Treatment labels
  trt_labels <- LETTERS[1:n_treatments]

  data.frame(
    study = paste0("Study", 1:n_studies),
    treat1 = trt_labels[comparisons$treat1],
    treat2 = trt_labels[comparisons$treat2],
    # Surrogate endpoint (always observed)
    pfs_effect = rnorm(n_studies, mean = 0.5, sd = 0.2),
    pfs_se = runif(n_studies, 0.1, 0.2),
    # True endpoint (some missing)
    os_effect = ifelse(
      runif(n_studies) > missing_true_pct,
      rnorm(n_studies, mean = 0.3, sd = 0.15),
      NA
    ),
    os_se = ifelse(
      !is.na(rnorm(n_studies)),
      runif(n_studies, 0.15, 0.25),
      NA
    )
  )
}


# ==============================================================================
# Test: build_surrogate_network()
# ==============================================================================

test_that("build_surrogate_network creates valid surrogate_network object", {

  data <- create_test_surrogate_data(n_studies = 6, n_treatments = 4)

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  # Check class
  expect_s3_class(net, "surrogate_network")
  expect_s3_class(net, "powernma_data")

  # Check structure
  expect_true("S_eff" %in% names(net))
  expect_true("T_eff" %in% names(net))
  expect_true("K" %in% names(net))
  expect_true("J" %in% names(net))

  # Check dimensions
  expect_equal(net$K, 4)  # 4 treatments
  expect_equal(length(net$S_eff), 6)  # 6 comparisons

  # Check treatment coding
  expect_type(net$treat1, "integer")
  expect_type(net$treat2, "integer")
  expect_true(all(net$treat1 >= 1 & net$treat1 <= net$K))
  expect_true(all(net$treat2 >= 1 & net$treat2 <= net$K))

  # Check flags
  expect_true(net$has_true_endpoint)
  expect_true(net$n_true_observed >= 0)
  expect_true(net$n_true_observed <= net$n_comparisons)
})


test_that("build_surrogate_network handles missing true endpoint", {

  data <- create_test_surrogate_data(n_studies = 5, missing_true_pct = 1.0)
  data$os_effect <- NULL
  data$os_se <- NULL

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se
  )

  expect_s3_class(net, "surrogate_network")
  expect_false(net$has_true_endpoint)
  expect_equal(net$n_true_observed, 0)
  expect_null(net$T_eff)
})


test_that("build_surrogate_network validates inputs", {

  data <- data.frame(
    study = c("S1", "S2"),
    treat1 = c("A", "B"),
    treat2 = c("B", "C"),
    pfs = c(0.5, 0.6),
    pfs_se = c(0.1, 0.1)
  )

  # Missing required columns
  expect_error(
    build_surrogate_network(
      data = data,
      study = study,
      treat1 = treat1,
      treat2 = treat2,
      S_eff = missing_col,
      S_se = pfs_se
    ),
    "object 'missing_col' not found|cannot find"
  )

  # Missing S_se
  expect_error(
    build_surrogate_network(
      data = data,
      study = study,
      treat1 = treat1,
      treat2 = treat2,
      S_eff = pfs
      # S_se missing
    ),
    "S_eff and S_se are required"
  )
})


test_that("build_surrogate_network handles multiple surrogates", {

  data <- data.frame(
    study = c("S1", "S2", "S3"),
    treat1 = c("A", "B", "A"),
    treat2 = c("B", "C", "C"),
    S1 = rnorm(3, 0.5, 0.1),
    S2 = rnorm(3, 0.4, 0.1),
    S3 = rnorm(3, 0.6, 0.1),
    se = rep(0.1, 3),
    T_eff = c(0.3, 0.2, 0.4),
    T_se = rep(0.15, 3)
  )

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = S1,
    S_se = se,
    T_eff = T_eff,
    T_se = T_se,
    S_multi = c("S1", "S2", "S3")
  )

  expect_s3_class(net, "surrogate_network")
  expect_false(is.null(net$S_multi))
  expect_equal(ncol(net$S_multi), 3)
  expect_equal(nrow(net$S_multi), 3)
})


# ==============================================================================
# Test: train_surrogate_index()
# ==============================================================================

test_that("train_surrogate_index with OLS method", {

  # Create data with multiple surrogates
  data <- data.frame(
    study = paste0("S", 1:20),
    treat1 = sample(LETTERS[1:4], 20, replace = TRUE),
    treat2 = sample(LETTERS[1:4], 20, replace = TRUE),
    S1 = rnorm(20, 0.5, 0.2),
    S2 = rnorm(20, 0.4, 0.2),
    S3 = rnorm(20, 0.6, 0.2),
    se_surr = rep(0.1, 20),
    T_eff = rnorm(20, 0.3, 0.15),  # True endpoint (all observed)
    T_se = rep(0.15, 20)
  )

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = S1,
    S_se = se_surr,
    T_eff = T_eff,
    T_se = T_se,
    S_multi = c("S1", "S2", "S3")
  )

  si_model <- train_surrogate_index(
    net = net,
    method = "ols",
    standardize = TRUE,
    seed = 123
  )

  # Check structure
  expect_s3_class(si_model, "surrogate_index")
  expect_equal(si_model$method, "ols")
  expect_equal(si_model$n_surrogates, 3)
  expect_true(!is.null(si_model$fit$coef))
  expect_true(si_model$r_squared >= 0 && si_model$r_squared <= 1)
})


test_that("train_surrogate_index with Ridge method", {

  skip_if_not_installed("glmnet")

  data <- data.frame(
    study = paste0("S", 1:15),
    treat1 = sample(LETTERS[1:3], 15, replace = TRUE),
    treat2 = sample(LETTERS[1:3], 15, replace = TRUE),
    S1 = rnorm(15, 0.5, 0.2),
    S2 = rnorm(15, 0.4, 0.2),
    se_surr = rep(0.1, 15),
    T_eff = rnorm(15, 0.3, 0.15),
    T_se = rep(0.15, 15)
  )

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = S1,
    S_se = se_surr,
    T_eff = T_eff,
    T_se = T_se,
    S_multi = c("S1", "S2")
  )

  si_model <- train_surrogate_index(
    net = net,
    method = "ridge",
    standardize = TRUE
  )

  expect_s3_class(si_model, "surrogate_index")
  expect_equal(si_model$method, "ridge")
  expect_true(!is.null(si_model$fit$lambda))
  expect_true(si_model$fit$lambda > 0)
})


test_that("train_surrogate_index validates inputs", {

  data <- create_test_surrogate_data(n_studies = 10, missing_true_pct = 0.2)

  # Network without S_multi
  net_no_multi <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  expect_error(
    train_surrogate_index(net_no_multi, method = "ols"),
    "S_multi required"
  )

  # Network without true endpoint
  data_no_true <- data
  data_no_true$os_effect <- NULL
  data_no_true$os_se <- NULL

  net_no_true <- build_surrogate_network(
    data = data_no_true,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    S_multi = c("pfs_effect")
  )

  expect_error(
    train_surrogate_index(net_no_true, method = "ols"),
    "True endpoint.*required"
  )
})


# ==============================================================================
# Test: apply_surrogate_index()
# ==============================================================================

test_that("apply_surrogate_index updates network correctly", {

  skip_if_not_installed("glmnet")

  data <- data.frame(
    study = paste0("S", 1:15),
    treat1 = sample(LETTERS[1:3], 15, replace = TRUE),
    treat2 = sample(LETTERS[1:3], 15, replace = TRUE),
    S1 = rnorm(15, 0.5, 0.1),
    S2 = rnorm(15, 0.4, 0.1),
    se_s = rep(0.1, 15),
    T_eff = rnorm(15, 0.3, 0.1),
    T_se = rep(0.15, 15)
  )

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = S1,
    S_se = se_s,
    T_eff = T_eff,
    T_se = T_se,
    S_multi = c("S1", "S2")
  )

  # Train SI
  si_model <- train_surrogate_index(net, method = "ols")

  # Apply SI
  net_augmented <- apply_surrogate_index(net, si_model)

  # Check that S_eff was updated
  expect_false(identical(net$S_eff, net_augmented$S_eff))

  # Check attributes
  expect_true(attr(net_augmented, "has_SI"))
  expect_equal(attr(net_augmented, "SI_method"), "ols")
  expect_true(!is.null(attr(net_augmented, "SI_r_squared")))
})


# ==============================================================================
# Test: compute_surrogate_threshold_effect()
# ==============================================================================

test_that("compute_surrogate_threshold_effect calculates STE correctly", {

  # Simulate posterior draws
  set.seed(456)
  n_draws <- 500

  # α ~ N(0.1, 0.05), β ~ N(0.8, 0.1)
  alpha_draws <- rnorm(n_draws, mean = 0.1, sd = 0.05)
  beta_draws <- rnorm(n_draws, mean = 0.8, sd = 0.1)

  ste_result <- compute_surrogate_threshold_effect(
    alpha_draws = alpha_draws,
    beta_draws = beta_draws,
    threshold_T = 0.0,
    conf_level = 0.95
  )

  # Check structure
  expect_s3_class(ste_result, "surrogate_threshold_effect")
  expect_true(!is.null(ste_result$median))
  expect_true(!is.null(ste_result$mean))
  expect_true(!is.null(ste_result$ci_lower))
  expect_true(!is.null(ste_result$ci_upper))

  # Check formula: STE = (threshold_T - alpha) / beta
  # Expected STE ≈ (0 - 0.1) / 0.8 = -0.125
  expect_true(abs(ste_result$mean - (-0.125)) < 0.05)

  # Check CI ordering
  expect_true(ste_result$ci_lower < ste_result$median)
  expect_true(ste_result$median < ste_result$ci_upper)
})


test_that("compute_surrogate_threshold_effect handles edge cases", {

  # Case 1: Few draws
  expect_error(
    compute_surrogate_threshold_effect(
      alpha_draws = c(0.1, 0.2),
      beta_draws = c(0.8, 0.9),
      threshold_T = 0
    ),
    "at least 10 draws"
  )

  # Case 2: Mismatched lengths
  expect_error(
    compute_surrogate_threshold_effect(
      alpha_draws = rnorm(100),
      beta_draws = rnorm(50),
      threshold_T = 0
    ),
    "same length"
  )

  # Case 3: Beta near zero (weak surrogacy)
  alpha_draws <- rnorm(100, 0, 0.1)
  beta_draws <- rnorm(100, 0.01, 0.01)  # Very small beta

  expect_warning(
    compute_surrogate_threshold_effect(alpha_draws, beta_draws, 0),
    "non-finite STE|weak surrogacy"
  )
})


# ==============================================================================
# Test: fit_bivariate_nma_freq()
# ==============================================================================

test_that("fit_bivariate_nma_freq produces valid results", {

  skip_if_not_installed("MASS")

  # Create larger dataset for stable estimation
  set.seed(789)
  data <- create_test_surrogate_data(
    n_studies = 20,
    n_treatments = 5,
    missing_true_pct = 0.3
  )

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  fit <- fit_bivariate_nma_freq(
    net = net,
    n_boot = 50,  # Small for speed
    boot_method = "normal",
    seed = 123
  )

  # Check structure
  expect_s3_class(fit, "bivariate_nma_fit")
  expect_s3_class(fit, "powernma_result")
  expect_equal(fit$engine, "frequentist")

  # Check treatment effects
  expect_equal(length(fit$dS), net$K)
  expect_equal(length(fit$dT), net$K)

  # Reference treatment should be 0
  expect_equal(fit$dS[1], 0)
  expect_equal(fit$dT[1], 0)

  # Check surrogacy parameters
  expect_true(!is.null(fit$surrogacy$alpha))
  expect_true(!is.null(fit$surrogacy$beta))
  expect_equal(length(fit$surrogacy$alpha_draws), 50)
  expect_equal(length(fit$surrogacy$beta_draws), 50)

  # Check rankings
  expect_true(!is.null(fit$ranks$sucra))
  expect_equal(length(fit$ranks$sucra), net$K)
  expect_true(all(fit$ranks$sucra >= 0 & fit$ranks$sucra <= 1))

  # Check bootstrap draws
  expect_equal(dim(fit$draws_T), c(50, net$K))
})


test_that("fit_bivariate_nma_freq handles missing true endpoint gracefully", {

  # All true endpoint missing
  data <- create_test_surrogate_data(n_studies = 10, missing_true_pct = 1.0)
  data$os_effect <- NA
  data$os_se <- NA

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  expect_warning(
    fit <- fit_bivariate_nma_freq(net, n_boot = 20),
    "No true endpoint data"
  )

  # Should still fit surrogate NMA
  expect_s3_class(fit, "bivariate_nma_fit")
  expect_true(!is.null(fit$dS))
})


test_that("fit_bivariate_nma_freq with student-t bootstrap", {

  data <- create_test_surrogate_data(n_studies = 15, missing_true_pct = 0.2)

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  fit <- fit_bivariate_nma_freq(
    net = net,
    n_boot = 30,
    boot_method = "student",
    df = 5
  )

  expect_s3_class(fit, "bivariate_nma_fit")
  expect_equal(fit$uncertainty$boot_method, "student")
})


# ==============================================================================
# Test: Print Methods
# ==============================================================================

test_that("print methods work without errors", {

  data <- create_test_surrogate_data(n_studies = 8, missing_true_pct = 0.25)

  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  # Print surrogate_network
  expect_output(print(net), "Surrogate Network")
  expect_output(print(net), "Studies")
  expect_output(print(net), "Treatments")

  # Fit and print bivariate_nma_fit
  fit <- fit_bivariate_nma_freq(net, n_boot = 20, seed = 999)
  expect_output(print(fit), "Bivariate Network Meta-Analysis")
  expect_output(print(fit), "Treatment Effects")
  expect_output(print(fit), "Surrogacy Parameters")

  # Print STE
  ste <- compute_surrogate_threshold_effect(
    alpha_draws = rnorm(100, 0.1, 0.05),
    beta_draws = rnorm(100, 0.8, 0.1),
    threshold_T = 0
  )
  expect_output(print(ste), "Surrogate Threshold Effect")
  expect_output(print(ste), "STE estimate")
})


# ==============================================================================
# Integration Test: Full Workflow
# ==============================================================================

test_that("Full surrogate NMA workflow", {

  skip_if_not_installed("MASS")

  # 1. Create realistic dataset
  set.seed(12345)
  data <- create_test_surrogate_data(
    n_studies = 25,
    n_treatments = 6,
    missing_true_pct = 0.4
  )

  # 2. Build network
  net <- build_surrogate_network(
    data = data,
    study = study,
    treat1 = treat1,
    treat2 = treat2,
    S_eff = pfs_effect,
    S_se = pfs_se,
    T_eff = os_effect,
    T_se = os_se
  )

  expect_s3_class(net, "surrogate_network")

  # 3. Fit bivariate NMA
  fit <- fit_bivariate_nma_freq(
    net = net,
    n_boot = 100,
    boot_method = "normal",
    seed = 456
  )

  expect_s3_class(fit, "bivariate_nma_fit")

  # 4. Compute STE
  ste <- compute_surrogate_threshold_effect(
    alpha_draws = fit$surrogacy$alpha_draws,
    beta_draws = fit$surrogacy$beta_draws,
    threshold_T = 0.0,
    conf_level = 0.95
  )

  expect_s3_class(ste, "surrogate_threshold_effect")
  expect_true(is.finite(ste$median))

  # 5. Check consistency
  # Best treatment on surrogate should correlate with best on true
  best_S <- which.max(fit$dS)
  best_T <- which.max(fit$dT)

  # They should be similar (allowing for some variation)
  expect_true(abs(best_S - best_T) <= 2)  # Within 2 ranks
})
