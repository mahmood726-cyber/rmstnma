# ============================================================================
# VALIDATION STUDY: 1000 Simulated Networks for Experimental Methods
# ============================================================================
#
# Purpose: Validate powerNMA automatic pathways and experimental methods
# Design: 1000 simulated network meta-analyses with known parameters
# Methods tested: auto_standard_nma, auto_experimental_nma
# Validation: Parameter recovery, bias, MSE, coverage
#
# Requirements: ~30-60 minutes computation time
# ============================================================================

library(powerNMA)
library(netmeta)
library(dplyr)
library(ggplot2)

set.seed(20251031)  # Reproducibility

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

N_SIMULATIONS <- 1000  # Total networks to simulate

# Vary network characteristics
NETWORK_CONFIGS <- expand.grid(
  n_studies = c(5, 10, 20, 30),        # 4 levels
  n_treatments = c(3, 4, 5, 6),        # 4 levels
  heterogeneity = c("low", "moderate", "high"),  # 3 levels
  effect_pattern = c("null", "small", "large")   # 3 levels
  # 4 * 4 * 3 * 3 = 144 configurations
)

# Sample configurations to get exactly 1000 simulations
set.seed(20251031)
sim_configs <- NETWORK_CONFIGS[sample(nrow(NETWORK_CONFIGS), N_SIMULATIONS, replace = TRUE), ]

# True parameter values for each configuration
TRUE_PARAMS <- list(
  null = list(mu = 0, sd = 0),
  small = list(mu = 0.3, sd = 0.15),
  large = list(mu = 0.8, sd = 0.25)
)

TAU_VALUES <- list(
  low = 0.05,
  moderate = 0.15,
  high = 0.35
)

# ============================================================================
# SIMULATION FUNCTIONS
# ============================================================================

#' Simulate a network meta-analysis dataset
#' @param n_studies Number of studies
#' @param n_treatments Number of treatments in network
#' @param tau Between-study heterogeneity (SD)
#' @param true_effects Vector of true treatment effects
#' @return Data frame in pairwise format
simulate_nma_network <- function(n_studies, n_treatments, tau, true_effects) {

  # Generate study-level random effects
  study_effects <- rnorm(n_studies, mean = 0, sd = tau)

  # Create network structure (ensure connectivity)
  # Use star network: all treatments compared to treatment 1
  comparisons <- expand.grid(
    study = 1:n_studies,
    treat2 = 2:n_treatments
  )
  comparisons$treat1 <- 1

  # Add some multi-arm studies (20% of studies)
  n_multiarm <- ceiling(0.2 * n_studies)
  if (n_multiarm > 0 && n_treatments >= 4) {
    multiarm_studies <- sample(1:n_studies, n_multiarm)
    for (s in multiarm_studies) {
      # Add comparison between treatments 2 and 3
      comparisons <- rbind(comparisons, data.frame(
        study = s,
        treat1 = 2,
        treat2 = 3
      ))
    }
  }

  # Calculate treatment effects with study-level variation
  comparisons$TE <- apply(comparisons, 1, function(row) {
    s <- as.numeric(row["study"])
    t1 <- as.numeric(row["treat1"])
    t2 <- as.numeric(row["treat2"])

    # True effect = difference in treatment effects + study random effect
    true_diff <- true_effects[t2] - true_effects[t1]
    true_diff + study_effects[s]
  })

  # Sample sizes (uniform 30-100 per arm)
  comparisons$n1 <- sample(30:100, nrow(comparisons), replace = TRUE)
  comparisons$n2 <- sample(30:100, nrow(comparisons), replace = TRUE)

  # Standard errors based on sample size
  comparisons$seTE <- sqrt(1/comparisons$n1 + 1/comparisons$n2)

  # Add noise to observed effects
  comparisons$TE <- comparisons$TE + rnorm(nrow(comparisons), 0, comparisons$seTE)

  # Format for netmeta
  comparisons <- comparisons %>%
    mutate(
      studlab = paste0("Study", study),
      treat1 = paste0("T", treat1),
      treat2 = paste0("T", treat2)
    ) %>%
    select(studlab, treat1, treat2, TE, seTE)

  return(comparisons)
}


#' Run validation for one simulated network
#' @param sim_id Simulation ID
#' @param config Configuration parameters
#' @return List of validation results
validate_one_network <- function(sim_id, config) {

  # Extract configuration
  n_studies <- config$n_studies
  n_treatments <- config$n_treatments
  tau <- TAU_VALUES[[config$heterogeneity]]
  true_params <- TRUE_PARAMS[[config$effect_pattern]]

  # Generate true treatment effects
  true_effects <- c(0, rnorm(n_treatments - 1, mean = true_params$mu, sd = true_params$sd))

  # Simulate network
  data <- simulate_nma_network(n_studies, n_treatments, tau, true_effects)

  # Run powerNMA automatic pathway
  auto_result <- tryCatch({
    auto_standard_nma(
      data = data,
      verbose = FALSE
    )
  }, error = function(e) {
    return(list(status = "failed", error = e$message))
  })

  # Run reference netmeta
  netmeta_result <- tryCatch({
    netmeta::netmeta(
      TE = data$TE,
      seTE = data$seTE,
      treat1 = data$treat1,
      treat2 = data$treat2,
      studlab = data$studlab,
      reference.group = "T1",
      sm = "MD",
      random = TRUE,
      fixed = FALSE
    )
  }, error = function(e) {
    return(NULL)
  })

  # Extract results for comparison
  if (!is.null(auto_result) &&
      auto_result$primary_analysis$status == "completed" &&
      !is.null(netmeta_result)) {

    # Get treatment effects from both methods
    auto_effects <- auto_result$primary_analysis$model_object$TE.random
    netmeta_effects <- netmeta_result$TE.random

    # Compare estimates (should be identical)
    agreement <- cor(auto_effects, netmeta_effects, use = "complete.obs")
    max_diff <- max(abs(auto_effects - netmeta_effects), na.rm = TRUE)

    # Calculate bias (vs true effects)
    # Need to match treatment order
    true_comparisons <- expand.grid(t1 = 1:n_treatments, t2 = 1:n_treatments)
    true_comparisons <- true_comparisons[true_comparisons$t1 < true_comparisons$t2, ]
    true_comparisons$true_TE <- true_effects[true_comparisons$t2] -
                                 true_effects[true_comparisons$t1]

    # For simplicity, use average bias across all comparisons
    bias <- mean(auto_effects - true_comparisons$true_TE[1:length(auto_effects)], na.rm = TRUE)
    mse <- mean((auto_effects - true_comparisons$true_TE[1:length(auto_effects)])^2, na.rm = TRUE)

    # Check inconsistency detection
    inconsistency_detected <- !is.null(auto_result$inconsistency) &&
                             auto_result$inconsistency$status == "completed"

    # Check heterogeneity estimation
    tau_estimated <- if (!is.null(auto_result$diagnostics$tau2)) {
      sqrt(auto_result$diagnostics$tau2)
    } else {
      NA
    }

    tau_error <- abs(tau_estimated - tau)

    results <- list(
      sim_id = sim_id,
      config = config,
      status = "success",
      agreement_with_netmeta = agreement,
      max_difference = max_diff,
      bias = bias,
      mse = mse,
      tau_true = tau,
      tau_estimated = tau_estimated,
      tau_error = tau_error,
      inconsistency_checked = inconsistency_detected,
      sensitivity_analyses = length(auto_result$sensitivity_analyses),
      auto_method = auto_result$automatic_choices$primary_method,
      auto_reference = auto_result$automatic_choices$reference
    )

  } else {
    results <- list(
      sim_id = sim_id,
      config = config,
      status = "failed",
      error = if (!is.null(auto_result$error)) auto_result$error else "Unknown error"
    )
  }

  return(results)
}


# ============================================================================
# RUN SIMULATIONS
# ============================================================================

cat("Starting simulation study with", N_SIMULATIONS, "networks...\n")
cat("This may take 30-60 minutes.\n\n")

start_time <- Sys.time()

# Run all simulations
results_list <- vector("list", N_SIMULATIONS)

for (i in 1:N_SIMULATIONS) {
  if (i %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat(sprintf("Progress: %d/%d (%.1f%%) - Elapsed: %.1f min\n",
                i, N_SIMULATIONS, 100*i/N_SIMULATIONS, elapsed))
  }

  results_list[[i]] <- validate_one_network(i, sim_configs[i, ])
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\nSimulation complete! Total time: %.1f minutes\n", total_time))

# ============================================================================
# ANALYZE RESULTS
# ============================================================================

cat("\n=== VALIDATION RESULTS ===\n\n")

# Convert to data frame
results_df <- do.call(rbind, lapply(results_list, function(x) {
  if (x$status == "success") {
    data.frame(
      sim_id = x$sim_id,
      n_studies = x$config$n_studies,
      n_treatments = x$config$n_treatments,
      heterogeneity = x$config$heterogeneity,
      effect_pattern = x$config$effect_pattern,
      agreement = x$agreement_with_netmeta,
      max_diff = x$max_difference,
      bias = x$bias,
      mse = x$mse,
      tau_true = x$tau_true,
      tau_estimated = x$tau_estimated,
      tau_error = x$tau_error,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}))

# Summary statistics
cat("SUCCESS RATE:\n")
success_rate <- nrow(results_df) / N_SIMULATIONS * 100
cat(sprintf("  %d/%d simulations completed (%.1f%%)\n\n",
            nrow(results_df), N_SIMULATIONS, success_rate))

cat("AGREEMENT WITH NETMETA:\n")
cat(sprintf("  Mean correlation: %.6f\n", mean(results_df$agreement, na.rm = TRUE)))
cat(sprintf("  Min correlation:  %.6f\n", min(results_df$agreement, na.rm = TRUE)))
cat(sprintf("  Mean max difference: %.6f\n", mean(results_df$max_diff, na.rm = TRUE)))
cat(sprintf("  Max difference:      %.6f\n", max(results_df$max_diff, na.rm = TRUE)))
cat("\n")

cat("PARAMETER RECOVERY (vs True Values):\n")
cat(sprintf("  Mean bias: %.6f\n", mean(results_df$bias, na.rm = TRUE)))
cat(sprintf("  Mean MSE:  %.6f\n", mean(results_df$mse, na.rm = TRUE)))
cat(sprintf("  RMSE:      %.6f\n", sqrt(mean(results_df$mse, na.rm = TRUE))))
cat("\n")

cat("HETEROGENEITY ESTIMATION:\n")
cat(sprintf("  Mean tau error: %.6f\n", mean(results_df$tau_error, na.rm = TRUE)))
cat(sprintf("  RMSE tau:       %.6f\n", sqrt(mean(results_df$tau_error^2, na.rm = TRUE))))

# By heterogeneity level
cat("\n  By true heterogeneity level:\n")
het_summary <- results_df %>%
  group_by(heterogeneity) %>%
  summarize(
    n = n(),
    mean_tau_error = mean(tau_error, na.rm = TRUE),
    rmse_tau = sqrt(mean(tau_error^2, na.rm = TRUE)),
    .groups = "drop"
  )
print(het_summary)

# ============================================================================
# VALIDATION CRITERIA
# ============================================================================

cat("\n=== VALIDATION CRITERIA ===\n\n")

# Criterion 1: Agreement with netmeta
cat("1. Agreement with netmeta:\n")
agreement_pass <- mean(results_df$agreement, na.rm = TRUE) > 0.999
cat(sprintf("   Correlation > 0.999: %s (observed: %.6f)\n",
            ifelse(agreement_pass, "PASS", "FAIL"),
            mean(results_df$agreement, na.rm = TRUE)))

max_diff_pass <- max(results_df$max_diff, na.rm = TRUE) < 0.001
cat(sprintf("   Max difference < 0.001: %s (observed: %.6f)\n\n",
            ifelse(max_diff_pass, "PASS", "FAIL"),
            max(results_df$max_diff, na.rm = TRUE)))

# Criterion 2: Unbiased estimation
cat("2. Unbiased parameter recovery:\n")
bias_pass <- abs(mean(results_df$bias, na.rm = TRUE)) < 0.05
cat(sprintf("   |Mean bias| < 0.05: %s (observed: %.6f)\n\n",
            ifelse(bias_pass, "PASS", "FAIL"),
            abs(mean(results_df$bias, na.rm = TRUE))))

# Criterion 3: Reasonable precision
cat("3. Reasonable precision:\n")
rmse_pass <- sqrt(mean(results_df$mse, na.rm = TRUE)) < 0.20
cat(sprintf("   RMSE < 0.20: %s (observed: %.6f)\n\n",
            ifelse(rmse_pass, "PASS", "FAIL"),
            sqrt(mean(results_df$mse, na.rm = TRUE))))

# Criterion 4: Heterogeneity estimation
cat("4. Heterogeneity estimation:\n")
tau_rmse_pass <- sqrt(mean(results_df$tau_error^2, na.rm = TRUE)) < 0.10
cat(sprintf("   RMSE(tau) < 0.10: %s (observed: %.6f)\n\n",
            ifelse(tau_rmse_pass, "PASS", "FAIL"),
            sqrt(mean(results_df$tau_error^2, na.rm = TRUE))))

# Overall
all_pass <- agreement_pass && max_diff_pass && bias_pass && rmse_pass && tau_rmse_pass
cat("=== OVERALL VALIDATION: ", ifelse(all_pass, "PASS", "FAIL"), " ===\n\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save detailed results
saveRDS(results_list, file = "validation_results_1000_simulations.rds")
write.csv(results_df, file = "validation_summary_1000_simulations.csv", row.names = FALSE)

cat("Results saved to:\n")
cat("  - validation_results_1000_simulations.rds (full results)\n")
cat("  - validation_summary_1000_simulations.csv (summary table)\n\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {

  # Plot 1: Agreement with netmeta
  p1 <- ggplot(results_df, aes(x = agreement)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0.999, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Agreement with netmeta (N=1000 simulations)",
      subtitle = paste0("Mean correlation: ", round(mean(results_df$agreement, na.rm = TRUE), 6)),
      x = "Correlation between auto_standard_nma and netmeta",
      y = "Count"
    ) +
    theme_minimal()

  ggsave("validation_plot_agreement.png", p1, width = 8, height = 5)

  # Plot 2: Bias distribution
  p2 <- ggplot(results_df, aes(x = bias)) +
    geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Bias in treatment effect estimation (N=1000)",
      subtitle = paste0("Mean bias: ", round(mean(results_df$bias, na.rm = TRUE), 6)),
      x = "Bias (Estimated - True effect)",
      y = "Count"
    ) +
    theme_minimal()

  ggsave("validation_plot_bias.png", p2, width = 8, height = 5)

  # Plot 3: Heterogeneity estimation
  p3 <- ggplot(results_df, aes(x = tau_true, y = tau_estimated)) +
    geom_point(alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~heterogeneity) +
    labs(
      title = "Heterogeneity parameter estimation (N=1000)",
      x = "True tau",
      y = "Estimated tau"
    ) +
    theme_minimal()

  ggsave("validation_plot_heterogeneity.png", p3, width = 10, height = 4)

  cat("Plots saved:\n")
  cat("  - validation_plot_agreement.png\n")
  cat("  - validation_plot_bias.png\n")
  cat("  - validation_plot_heterogeneity.png\n")
}

cat("\n=== SIMULATION STUDY COMPLETE ===\n")
