#' Generate Comprehensive Simulation Datasets for Testing
#'
#' This file contains functions to generate large, realistic simulation datasets
#' for comprehensive validation of powerNMA across various scenarios.

#' Generate Large Star Network
#'
#' Creates a star-shaped network where all treatments are compared to one reference
#'
#' @param n_treatments Number of treatments (including reference)
#' @param n_studies_per_comparison Average number of studies per comparison
#' @param true_effects Vector of true treatment effects (length = n_treatments)
#' @param tau Between-study heterogeneity (SD)
#' @param avg_study_size Average sample size per study arm
#' @param effect_measure "MD", "OR", "HR", etc.
#' @param seed Random seed
#' @return Data frame with pairwise comparisons
#' @export
simulate_large_star_network <- function(n_treatments = 10,
                                       n_studies_per_comparison = 20,
                                       true_effects = NULL,
                                       tau = 0.10,
                                       avg_study_size = 200,
                                       effect_measure = "MD",
                                       seed = 42) {
  set.seed(seed)

  if (is.null(true_effects)) {
    # Generate realistic treatment effects
    true_effects <- c(0, sort(stats::rnorm(n_treatments - 1, -0.3, 0.15)))
  }

  treatments <- paste0("Trt", seq_len(n_treatments))
  names(true_effects) <- treatments

  # Star network: all compare to Trt1 (reference)
  reference <- treatments[1]

  pairwise_data <- list()
  study_counter <- 1

  for (i in 2:n_treatments) {
    active_trt <- treatments[i]

    # Generate studies for this comparison
    for (j in 1:n_studies_per_comparison) {
      # True treatment effect
      true_diff <- true_effects[active_trt] - true_effects[reference]

      # Add study-level random effect
      study_re <- stats::rnorm(1, 0, tau)

      # Sample size variability
      n_per_arm <- round(stats::runif(1, avg_study_size * 0.5, avg_study_size * 1.5))

      # Standard error depends on sample size and effect measure
      if (effect_measure == "MD") {
        # Continuous outcome
        within_study_sd <- 1.0
        se <- within_study_sd * sqrt(2 / n_per_arm)
      } else {
        # Binary outcome (OR, RR, HR)
        baseline_risk <- 0.30
        se <- stats::runif(1, 0.10, 0.30)
      }

      # Observed effect
      observed_TE <- stats::rnorm(1, true_diff + study_re, se)

      pairwise_data[[study_counter]] <- tibble::tibble(
        studlab = paste0("Study_", study_counter),
        treat1 = reference,
        treat2 = active_trt,
        TE = observed_TE,
        seTE = se,
        n1 = n_per_arm,
        n2 = n_per_arm,
        year = sample(2000:2023, 1),
        quality = sample(c("High", "Moderate", "Low"), 1)
      )

      study_counter <- study_counter + 1
    }
  }

  dplyr::bind_rows(pairwise_data)
}

#' Generate Large Complete Network
#'
#' Creates a complete network where all treatments are compared to each other
#'
#' @param n_treatments Number of treatments
#' @param n_studies_per_comparison Average studies per comparison
#' @param true_effects Vector of true effects
#' @param tau Heterogeneity
#' @param seed Random seed
#' @return Data frame
#' @export
simulate_large_complete_network <- function(n_treatments = 6,
                                           n_studies_per_comparison = 10,
                                           true_effects = NULL,
                                           tau = 0.15,
                                           seed = 42) {
  set.seed(seed)

  if (is.null(true_effects)) {
    true_effects <- c(0, sort(stats::rnorm(n_treatments - 1, -0.3, 0.2)))
  }

  treatments <- paste0("Trt", seq_len(n_treatments))
  names(true_effects) <- treatments

  # All pairwise combinations
  comparisons <- utils::combn(treatments, 2, simplify = FALSE)

  pairwise_data <- list()
  study_counter <- 1

  for (comp in comparisons) {
    trt1 <- comp[1]
    trt2 <- comp[2]

    # Random number of studies (Poisson distributed)
    n_studies <- max(1, stats::rpois(1, lambda = n_studies_per_comparison))

    for (j in 1:n_studies) {
      true_diff <- true_effects[trt1] - true_effects[trt2]
      study_re <- stats::rnorm(1, 0, tau)
      se <- stats::runif(1, 0.08, 0.25)

      pairwise_data[[study_counter]] <- tibble::tibble(
        studlab = paste0("Study_", study_counter),
        treat1 = trt1,
        treat2 = trt2,
        TE = stats::rnorm(1, true_diff + study_re, se),
        seTE = se
      )

      study_counter <- study_counter + 1
    }
  }

  dplyr::bind_rows(pairwise_data)
}

#' Generate Network with Multi-Arm Trials
#'
#' Creates a network specifically including multi-arm trials
#'
#' @param n_2arm Number of 2-arm trials
#' @param n_3arm Number of 3-arm trials
#' @param n_4arm Number of 4-arm trials
#' @param n_treatments Total number of treatments in network
#' @param tau Heterogeneity
#' @param seed Random seed
#' @return Data frame
#' @export
simulate_multiarm_network <- function(n_2arm = 30,
                                     n_3arm = 10,
                                     n_4arm = 5,
                                     n_treatments = 8,
                                     tau = 0.12,
                                     seed = 42) {
  set.seed(seed)

  treatments <- paste0("Trt", seq_len(n_treatments))
  true_effects <- c(0, sort(stats::rnorm(n_treatments - 1, -0.25, 0.15)))
  names(true_effects) <- treatments

  pairwise_data <- list()
  study_counter <- 1

  # Helper function to generate one multi-arm trial
  generate_multiarm_trial <- function(study_id, n_arms) {
    # Randomly select treatments for this trial
    trial_treatments <- sample(treatments, n_arms)

    # Generate all pairwise comparisons
    pairs <- utils::combn(trial_treatments, 2, simplify = FALSE)

    study_re <- stats::rnorm(1, 0, tau)

    lapply(pairs, function(pair) {
      true_diff <- true_effects[pair[1]] - true_effects[pair[2]]
      se <- stats::runif(1, 0.10, 0.25)

      tibble::tibble(
        studlab = study_id,
        treat1 = pair[1],
        treat2 = pair[2],
        TE = stats::rnorm(1, true_diff + study_re, se),
        seTE = se,
        n_arms = n_arms
      )
    })
  }

  # Generate 2-arm trials
  for (i in 1:n_2arm) {
    pairs <- generate_multiarm_trial(paste0("Study2arm_", i), 2)
    pairwise_data <- c(pairwise_data, pairs)
  }

  # Generate 3-arm trials
  for (i in 1:n_3arm) {
    pairs <- generate_multiarm_trial(paste0("Study3arm_", i), 3)
    pairwise_data <- c(pairwise_data, pairs)
  }

  # Generate 4-arm trials
  for (i in 1:n_4arm) {
    pairs <- generate_multiarm_trial(paste0("Study4arm_", i), 4)
    pairwise_data <- c(pairwise_data, pairs)
  }

  dplyr::bind_rows(pairwise_data)
}

#' Generate Network with Varying Heterogeneity
#'
#' Creates a network where different comparisons have different heterogeneity levels
#'
#' @param n_treatments Number of treatments
#' @param n_studies_total Total number of studies
#' @param tau_low Low heterogeneity (for some comparisons)
#' @param tau_high High heterogeneity (for other comparisons)
#' @param seed Random seed
#' @return Data frame
#' @export
simulate_heterogeneous_network <- function(n_treatments = 8,
                                          n_studies_total = 100,
                                          tau_low = 0.05,
                                          tau_high = 0.30,
                                          seed = 42) {
  set.seed(seed)

  treatments <- paste0("Trt", seq_len(n_treatments))
  true_effects <- c(0, sort(stats::rnorm(n_treatments - 1, -0.3, 0.2)))
  names(true_effects) <- treatments

  comparisons <- utils::combn(treatments, 2, simplify = FALSE)

  # Assign heterogeneity level to each comparison
  comparison_tau <- stats::runif(length(comparisons), tau_low, tau_high)

  pairwise_data <- list()

  for (i in 1:n_studies_total) {
    # Random comparison
    comp_idx <- sample(length(comparisons), 1)
    comp <- comparisons[[comp_idx]]
    tau <- comparison_tau[comp_idx]

    true_diff <- true_effects[comp[1]] - true_effects[comp[2]]
    study_re <- stats::rnorm(1, 0, tau)
    se <- stats::runif(1, 0.08, 0.30)

    pairwise_data[[i]] <- tibble::tibble(
      studlab = paste0("Study_", i),
      treat1 = comp[1],
      treat2 = comp[2],
      TE = stats::rnorm(1, true_diff + study_re, se),
      seTE = se,
      comparison_tau = tau  # Store true tau for validation
    )
  }

  dplyr::bind_rows(pairwise_data)
}

#' Generate Sparse Network
#'
#' Creates a network with some disconnected comparisons and sparse data
#'
#' @param n_treatments Number of treatments
#' @param sparsity Proportion of possible comparisons to include (0-1)
#' @param min_studies_per_comparison Minimum studies per comparison
#' @param seed Random seed
#' @return Data frame
#' @export
simulate_sparse_network <- function(n_treatments = 12,
                                   sparsity = 0.3,
                                   min_studies_per_comparison = 2,
                                   seed = 42) {
  set.seed(seed)

  treatments <- paste0("Trt", seq_len(n_treatments))
  true_effects <- c(0, sort(stats::rnorm(n_treatments - 1, -0.25, 0.15)))
  names(true_effects) <- treatments

  # All possible comparisons
  all_comparisons <- utils::combn(treatments, 2, simplify = FALSE)

  # Randomly select subset based on sparsity
  n_comparisons <- ceiling(length(all_comparisons) * sparsity)
  selected_comparisons <- sample(all_comparisons, n_comparisons)

  pairwise_data <- list()
  study_counter <- 1

  for (comp in selected_comparisons) {
    # Very few studies per comparison (sparse)
    n_studies <- sample(min_studies_per_comparison:(min_studies_per_comparison + 3), 1)

    for (j in 1:n_studies) {
      true_diff <- true_effects[comp[1]] - true_effects[comp[2]]
      study_re <- stats::rnorm(1, 0, 0.15)
      se <- stats::runif(1, 0.15, 0.40)  # Large SE (small studies)

      pairwise_data[[study_counter]] <- tibble::tibble(
        studlab = paste0("Study_", study_counter),
        treat1 = comp[1],
        treat2 = comp[2],
        TE = stats::rnorm(1, true_diff + study_re, se),
        seTE = se
      )

      study_counter <- study_counter + 1
    }
  }

  dplyr::bind_rows(pairwise_data)
}

#' Generate Large IPD for Time-Varying NMA
#'
#' Creates large individual patient data for RMST/milestone analysis
#'
#' @param n_trials Number of trials
#' @param n_treatments Number of treatments
#' @param n_per_arm Sample size per arm
#' @param true_hazards Vector of true hazard rates
#' @param censoring_rate Censoring proportion
#' @param seed Random seed
#' @return IPD data frame
#' @export
simulate_large_ipd <- function(n_trials = 20,
                              n_treatments = 6,
                              n_per_arm = 200,
                              true_hazards = NULL,
                              censoring_rate = 0.20,
                              seed = 42) {
  set.seed(seed)

  if (is.null(true_hazards)) {
    # Realistic hazard rates (lower = better survival)
    true_hazards <- seq(0.010, 0.002, length.out = n_treatments)
  }

  treatments <- paste0("Trt", seq_len(n_treatments))
  names(true_hazards) <- treatments

  # Create trial designs (mix of 2-arm and 3-arm)
  ipd_list <- list()

  for (trial_id in 1:n_trials) {
    # Randomly choose 2 or 3 treatments for this trial
    n_arms <- sample(c(2, 3), 1, prob = c(0.7, 0.3))
    trial_treatments <- sample(treatments, n_arms)

    # Generate IPD for each arm
    for (trt in trial_treatments) {
      # Event times from exponential distribution
      event_times <- stats::rexp(n_per_arm, rate = true_hazards[trt])

      # Censoring times
      censor_times <- stats::runif(n_per_arm, 100, 500)

      # Observed times and status
      observed_times <- pmin(event_times, censor_times)
      status <- as.integer(event_times <= censor_times)

      # Add some random censoring
      random_censor <- stats::rbinom(n_per_arm, 1, censoring_rate)
      status[random_censor == 1] <- 0

      ipd_list[[length(ipd_list) + 1]] <- tibble::tibble(
        trial = paste0("Trial", sprintf("%02d", trial_id)),
        treatment = trt,
        time = observed_times,
        status = status,
        patient_id = paste0("T", trial_id, "_", trt, "_P", 1:n_per_arm)
      )
    }
  }

  dplyr::bind_rows(ipd_list)
}

#' Generate Complete Test Suite
#'
#' Creates all simulation datasets and saves them
#'
#' @param output_dir Directory to save datasets
#' @export
generate_all_test_datasets <- function(output_dir = "inst/testdata") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("Generating test datasets...\n\n")

  # 1. Large star network
  cat("1. Large star network (200 studies, 10 treatments)...\n")
  star <- simulate_large_star_network(
    n_treatments = 10,
    n_studies_per_comparison = 20,
    seed = 1001
  )
  readr::write_csv(star, file.path(output_dir, "large_star_network.csv"))
  cat(sprintf("   Created: %d studies, %d comparisons\n\n",
              length(unique(star$studlab)), nrow(star)))

  # 2. Complete network
  cat("2. Complete network (150 studies, 6 treatments)...\n")
  complete <- simulate_large_complete_network(
    n_treatments = 6,
    n_studies_per_comparison = 10,
    seed = 1002
  )
  readr::write_csv(complete, file.path(output_dir, "complete_network.csv"))
  cat(sprintf("   Created: %d studies\n\n", length(unique(complete$studlab))))

  # 3. Multi-arm network
  cat("3. Multi-arm network (45 trials with 2/3/4 arms)...\n")
  multiarm <- simulate_multiarm_network(
    n_2arm = 30,
    n_3arm = 10,
    n_4arm = 5,
    seed = 1003
  )
  readr::write_csv(multiarm, file.path(output_dir, "multiarm_network.csv"))
  cat(sprintf("   Created: %d unique trials, %d comparisons\n\n",
              length(unique(multiarm$studlab)), nrow(multiarm)))

  # 4. Heterogeneous network
  cat("4. Heterogeneous network (varying tau)...\n")
  hetero <- simulate_heterogeneous_network(
    n_treatments = 8,
    n_studies_total = 100,
    seed = 1004
  )
  readr::write_csv(hetero, file.path(output_dir, "heterogeneous_network.csv"))
  cat(sprintf("   Created: %d studies\n\n", nrow(hetero)))

  # 5. Sparse network
  cat("5. Sparse network (30%% connectivity)...\n")
  sparse <- simulate_sparse_network(
    n_treatments = 12,
    sparsity = 0.3,
    seed = 1005
  )
  readr::write_csv(sparse, file.path(output_dir, "sparse_network.csv"))
  cat(sprintf("   Created: %d studies\n\n", length(unique(sparse$studlab))))

  # 6. Large IPD
  cat("6. Large IPD dataset (20 trials, 6 treatments, 200/arm)...\n")
  large_ipd <- simulate_large_ipd(
    n_trials = 20,
    n_treatments = 6,
    n_per_arm = 200,
    seed = 1006
  )
  readr::write_csv(large_ipd, file.path(output_dir, "large_ipd.csv"))
  cat(sprintf("   Created: %d patients across %d trials\n\n",
              nrow(large_ipd), length(unique(large_ipd$trial))))

  cat("All datasets generated successfully!\n")
  cat(sprintf("Location: %s/\n", output_dir))

  invisible(list(
    star = star,
    complete = complete,
    multiarm = multiarm,
    heterogeneous = hetero,
    sparse = sparse,
    large_ipd = large_ipd
  ))
}
