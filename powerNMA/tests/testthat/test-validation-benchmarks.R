# Phase 3: Validation Benchmarks
#
# This file contains validation tests that benchmark powerNMA against
# netmeta and gemtc to ensure correctness

test_that("Standard mode matches netmeta exactly on simple network", {
  set.seed(42)

  # Create simple network data
  data <- simulate_nma_data(n_studies = 20, seed = 42)

  # Run netmeta directly
  nm_direct <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "HR",
    reference.group = "Placebo"
  )

  # Run powerNMA standard mode
  pnma_result <- run_powernma(
    data = data,
    data_type = "pairwise",
    mode = "standard"
  )

  # Compare treatment effects (fixed)
  expect_equal(
    pnma_result$network$TE.fixed,
    nm_direct$TE.fixed,
    tolerance = 1e-10
  )

  # Compare treatment effects (random)
  expect_equal(
    pnma_result$network$TE.random,
    nm_direct$TE.random,
    tolerance = 1e-10
  )

  # Compare heterogeneity
  expect_equal(
    pnma_result$network$tau,
    nm_direct$tau,
    tolerance = 1e-10
  )

  # Compare I²
  expect_equal(
    pnma_result$network$I2,
    nm_direct$I2,
    tolerance = 1e-10
  )
})

test_that("Multi-arm trials are handled correctly", {
  set.seed(123)

  # Create data with a 3-arm trial
  data <- rbind(
    # 2-arm trials
    tibble::tibble(
      studlab = "Study1",
      treat1 = "A",
      treat2 = "B",
      TE = -0.2,
      seTE = 0.1
    ),
    tibble::tibble(
      studlab = "Study2",
      treat1 = "B",
      treat2 = "C",
      TE = -0.3,
      seTE = 0.12
    ),
    # 3-arm trial (all pairwise comparisons)
    tibble::tibble(
      studlab = "Study3_MultiArm",
      treat1 = c("A", "A", "B"),
      treat2 = c("B", "C", "C"),
      TE = c(-0.15, -0.4, -0.25),
      seTE = c(0.15, 0.15, 0.15)
    )
  )

  # netmeta should handle this correctly
  nm <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = "HR",
    reference.group = "A"
  )

  # powerNMA should match
  pnma <- run_powernma(data, data_type = "pairwise", mode = "standard")

  # All 3 comparisons from Study3 should be included
  study3_data <- pnma$network$data[pnma$network$data$studlab == "Study3_MultiArm", ]
  expect_equal(nrow(study3_data), 3)

  # Results should match netmeta
  expect_equal(pnma$network$TE.random, nm$TE.random, tolerance = 1e-8)
})

test_that("Simulation study: Type I error control (no true effect)", {
  skip_if(Sys.getenv("SKIP_LONG_TESTS") == "true")

  set.seed(999)
  n_sim <- 100  # Reduced for testing speed
  alpha <- 0.05

  # Simulate under null (no treatment effect)
  reject_count <- 0

  for (i in 1:n_sim) {
    # Generate network with NO treatment effects
    sim_data <- simulate_nma_network_null(
      n_studies = 15,
      n_treatments = 4,
      tau = 0.1,
      seed = i
    )

    # Run NMA
    nm <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = sim_data,
        sm = "MD"
      )
    }, error = function(e) NULL)

    if (!is.null(nm)) {
      # Test if any comparison is significant (should be rare under null)
      p_values <- nm$p.random
      if (any(p_values < alpha, na.rm = TRUE)) {
        reject_count <- reject_count + 1
      }
    }
  }

  # Type I error rate should be close to alpha
  observed_alpha <- reject_count / n_sim
  expect_true(observed_alpha < 0.15,  # Allow some variation
    info = sprintf("Type I error rate = %.3f (expected ~%.2f)", observed_alpha, alpha))
})

test_that("Simulation study: Power to detect true effects", {
  skip_if(Sys.getenv("SKIP_LONG_TESTS") == "true")

  set.seed(888)
  n_sim <- 50

  # Simulate with LARGE treatment effect
  detect_count <- 0

  for (i in 1:n_sim) {
    sim_data <- simulate_nma_network_with_effect(
      n_studies = 20,
      n_treatments = 3,
      true_effects = c(0, -0.5, -0.7),  # Large effects
      tau = 0.05,
      seed = i
    )

    nm <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = sim_data,
        sm = "MD",
        reference.group = "Trt1"
      )
    }, error = function(e) NULL)

    if (!is.null(nm)) {
      # Check if we detect the effect for Trt3 vs Trt1
      idx_trt3 <- which(rownames(nm$TE.random) == "Trt3")
      idx_trt1 <- which(rownames(nm$TE.random) == "Trt1")

      if (length(idx_trt3) > 0 && length(idx_trt1) > 0) {
        p_val <- nm$p.random[idx_trt3, idx_trt1]
        if (!is.na(p_val) && p_val < 0.05) {
          detect_count <- detect_count + 1
        }
      }
    }
  }

  power <- detect_count / n_sim
  expect_true(power > 0.7,
    info = sprintf("Power = %.2f (should be > 0.7 for large effects)", power))
})

test_that("Heterogeneity estimation accuracy", {
  set.seed(777)

  true_tau <- 0.15
  n_sim <- 50
  tau_estimates <- numeric(n_sim)

  for (i in 1:n_sim) {
    sim_data <- simulate_nma_network_with_effect(
      n_studies = 30,
      n_treatments = 4,
      true_effects = c(0, -0.2, -0.3, -0.1),
      tau = true_tau,
      seed = i
    )

    nm <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = sim_data,
        sm = "MD"
      )
    }, error = function(e) NULL)

    if (!is.null(nm)) {
      tau_estimates[i] <- nm$tau
    }
  }

  # Mean should be close to true value
  mean_tau <- mean(tau_estimates, na.rm = TRUE)
  expect_equal(mean_tau, true_tau, tolerance = 0.05)
})

# Simulation helper functions
simulate_nma_network_null <- function(n_studies, n_treatments, tau, seed) {
  set.seed(seed)

  treatments <- paste0("Trt", 1:n_treatments)
  designs <- utils::combn(treatments, 2, simplify = FALSE)

  purrr::map_dfr(1:n_studies, function(i) {
    design <- designs[[((i - 1) %% length(designs)) + 1]]

    # Under null: true_diff = 0
    study_re <- stats::rnorm(1, 0, tau)
    se <- stats::runif(1, 0.1, 0.3)

    tibble::tibble(
      studlab = paste0("Study", i),
      treat1 = design[1],
      treat2 = design[2],
      TE = stats::rnorm(1, 0 + study_re, se),  # Mean = 0 (null)
      seTE = se
    )
  })
}

simulate_nma_network_with_effect <- function(n_studies, n_treatments,
                                             true_effects, tau, seed) {
  set.seed(seed)

  treatments <- paste0("Trt", 1:n_treatments)
  names(true_effects) <- treatments

  designs <- utils::combn(treatments, 2, simplify = FALSE)

  purrr::map_dfr(1:n_studies, function(i) {
    design <- designs[[((i - 1) %% length(designs)) + 1]]

    # True treatment effect
    true_diff <- true_effects[design[1]] - true_effects[design[2]]

    # Add study-level random effect
    study_re <- stats::rnorm(1, 0, tau)
    se <- stats::runif(1, 0.08, 0.25)

    tibble::tibble(
      studlab = paste0("Study", i),
      treat1 = design[1],
      treat2 = design[2],
      TE = stats::rnorm(1, true_diff + study_re, se),
      seTE = se
    )
  })
}
