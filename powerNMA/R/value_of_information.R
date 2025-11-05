#' Threshold Analysis and Value of Information
#'
#' Tools for threshold analysis and value of information (VOI) analysis to
#' inform treatment decisions and research prioritization.
#'
#' @name value_of_information
#' @references
#' Claxton K, et al. (2005). Probabilistic sensitivity analysis for NICE
#' technology assessment: not an optional extra. Health Economics, 14(4):339-347.
#'
#' Strong M, et al. (2014). Estimating multiparameter partial expected value
#' of perfect information from a probabilistic sensitivity analysis sample.
#' Medical Decision Making, 34(3):311-326.
NULL

#' Threshold Analysis for Treatment Selection
#'
#' Determine at what threshold (willingness-to-pay, risk tolerance) the
#' optimal treatment changes.
#'
#' @param nma_result NMA result object
#' @param cost_data Data frame with treatment costs
#' @param benefit_measure Which outcome to use as benefit ("effect" or custom)
#' @param threshold_range Range of thresholds to test
#' @param reference Reference treatment
#' @return threshold_analysis object
#' @export
#' @examples
#' \dontrun{
#' # Treatment costs
#' costs <- data.frame(
#'   treatment = c("Placebo", "Drug_A", "Drug_B", "Drug_C"),
#'   cost = c(0, 1000, 1500, 2000)
#' )
#'
#' threshold <- threshold_analysis(
#'   nma_result = nma,
#'   cost_data = costs,
#'   threshold_range = seq(0, 50000, by = 1000)
#' )
#' plot(threshold)
#' }
threshold_analysis <- function(nma_result, cost_data, benefit_measure = "effect",
                              threshold_range = seq(0, 50000, by = 1000),
                              reference = NULL) {

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  treatments <- rownames(nma_result$TE.random)

  # Extract treatment effects (benefits)
  benefits <- nma_result$TE.random[, reference]
  benefits_se <- nma_result$seTE.random[, reference]

  # Merge with costs
  treatment_data <- data.frame(
    Treatment = treatments,
    Benefit = benefits,
    Benefit_SE = benefits_se,
    stringsAsFactors = FALSE
  )

  treatment_data <- merge(treatment_data, cost_data,
                         by.x = "Treatment", by.y = "treatment",
                         all.x = TRUE)

  # Incremental cost and benefit vs reference
  ref_cost <- treatment_data$cost[treatment_data$Treatment == reference]
  ref_benefit <- treatment_data$Benefit[treatment_data$Treatment == reference]

  treatment_data$Incremental_Cost <- treatment_data$cost - ref_cost
  treatment_data$Incremental_Benefit <- treatment_data$Benefit - ref_benefit

  message(sprintf("Running threshold analysis for %d thresholds...", length(threshold_range)))

  # For each threshold, calculate net benefit and find optimal treatment
  threshold_results <- lapply(threshold_range, function(threshold) {

    # Net benefit = Benefit * Threshold - Cost
    # Or: Net Benefit = Incremental_Benefit * Threshold - Incremental_Cost
    treatment_data$Net_Benefit <- treatment_data$Incremental_Benefit * threshold -
                                  treatment_data$Incremental_Cost

    # Find optimal treatment
    optimal_idx <- which.max(treatment_data$Net_Benefit)
    optimal_treatment <- treatment_data$Treatment[optimal_idx]
    optimal_nb <- treatment_data$Net_Benefit[optimal_idx]

    data.frame(
      Threshold = threshold,
      Optimal_Treatment = optimal_treatment,
      Net_Benefit = optimal_nb,
      stringsAsFactors = FALSE
    )
  })

  threshold_df <- do.call(rbind, threshold_results)

  # Identify threshold switches (where optimal treatment changes)
  switches <- which(threshold_df$Optimal_Treatment[-1] != threshold_df$Optimal_Treatment[-nrow(threshold_df)])

  switch_points <- data.frame(
    Threshold = threshold_df$Threshold[switches],
    From_Treatment = threshold_df$Optimal_Treatment[switches],
    To_Treatment = threshold_df$Optimal_Treatment[switches + 1],
    stringsAsFactors = FALSE
  )

  result <- list(
    threshold_results = threshold_df,
    switch_points = switch_points,
    treatment_data = treatment_data,
    threshold_range = threshold_range,
    reference = reference
  )

  class(result) <- c("threshold_analysis", "list")
  result
}

#' Expected Value of Perfect Information (EVPI)
#'
#' Calculate the value of eliminating all uncertainty about treatment effects.
#'
#' @param nma_result NMA result object
#' @param cost_data Treatment costs
#' @param threshold Willingness-to-pay threshold
#' @param n_sim Number of simulations
#' @param reference Reference treatment
#' @return evpi object
#' @export
#' @examples
#' \dontrun{
#' evpi_result <- calculate_evpi(
#'   nma_result = nma,
#'   cost_data = costs,
#'   threshold = 20000,
#'   n_sim = 10000
#' )
#' print(evpi_result)
#' }
calculate_evpi <- function(nma_result, cost_data, threshold, n_sim = 10000,
                          reference = NULL) {

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  message(sprintf("Calculating EVPI with %d simulations...", n_sim))

  treatments <- rownames(nma_result$TE.random)

  # Extract treatment effects and covariance
  effects <- nma_result$TE.random[, reference]
  vcov_matrix <- nma_result$Cov.random

  # Merge costs
  treatment_data <- data.frame(
    Treatment = treatments,
    Effect = effects,
    stringsAsFactors = FALSE
  )

  treatment_data <- merge(treatment_data, cost_data,
                         by.x = "Treatment", by.y = "treatment",
                         all.x = TRUE)

  ref_cost <- treatment_data$cost[treatment_data$Treatment == reference]

  treatment_data$Incremental_Cost <- treatment_data$cost - ref_cost
  treatment_data$Incremental_Effect <- treatment_data$Effect

  # Simulate from multivariate normal
  simulated_effects <- MASS::mvrnorm(n = n_sim, mu = effects, Sigma = vcov_matrix)

  # For each simulation, calculate net benefits and find optimal treatment
  simulation_nb <- matrix(NA, nrow = n_sim, ncol = length(treatments))
  colnames(simulation_nb) <- treatments

  for (i in seq_len(n_sim)) {
    sim_effects <- simulated_effects[i, ]

    # Net benefit for each treatment in this simulation
    nb <- (sim_effects - sim_effects[reference]) * threshold -
          (treatment_data$cost - ref_cost)

    simulation_nb[i, ] <- nb
  }

  # Find optimal treatment in each simulation
  optimal_per_sim <- apply(simulation_nb, 1, function(x) {
    treatments[which.max(x)]
  })

  max_nb_per_sim <- apply(simulation_nb, 1, max)

  # Expected net benefit with perfect information
  # = average of maximum NB across simulations
  expected_nb_perfect_info <- mean(max_nb_per_sim)

  # Expected net benefit with current information
  # = NB of treatment that is optimal on average (using mean effects)
  mean_nb <- colMeans(simulation_nb)
  expected_nb_current_info <- max(mean_nb)
  optimal_treatment_current <- treatments[which.max(mean_nb)]

  # EVPI = difference
  evpi <- expected_nb_perfect_info - expected_nb_current_info

  # Probability each treatment is optimal
  prob_optimal <- table(factor(optimal_per_sim, levels = treatments)) / n_sim

  result <- list(
    evpi = evpi,
    expected_nb_perfect_info = expected_nb_perfect_info,
    expected_nb_current_info = expected_nb_current_info,
    optimal_treatment = optimal_treatment_current,
    prob_optimal = as.vector(prob_optimal),
    treatment_names = names(prob_optimal),
    threshold = threshold,
    n_sim = n_sim
  )

  class(result) <- c("evpi", "list")

  message(sprintf("✓ EVPI = %.2f per patient", evpi))
  message(sprintf("  Optimal treatment: %s (%.1f%% probability)",
                 optimal_treatment_current,
                 max(prob_optimal) * 100))

  result
}

#' EVPI for Population
#'
#' Calculate total EVPI for a population over time.
#'
#' @param evpi_result evpi object
#' @param population_size Number of patients per year
#' @param time_horizon Years
#' @param discount_rate Annual discount rate
#' @return Numeric value (population EVPI)
#' @export
calculate_population_evpi <- function(evpi_result, population_size, time_horizon = 10,
                                     discount_rate = 0.035) {

  per_patient_evpi <- evpi_result$evpi

  # Calculate present value of population EVPI
  years <- seq_len(time_horizon)
  discount_factors <- 1 / (1 + discount_rate)^years

  annual_evpi <- population_size * per_patient_evpi
  discounted_annual_evpi <- annual_evpi * discount_factors

  population_evpi <- sum(discounted_annual_evpi)

  message(sprintf("Per-patient EVPI: $%.2f", per_patient_evpi))
  message(sprintf("Population EVPI (%d patients/year over %d years): $%.0f",
                 population_size, time_horizon, population_evpi))

  population_evpi
}

#' Expected Value of Perfect Parameter Information (EVPPI)
#'
#' Calculate EVPI for specific parameters to prioritize further research.
#'
#' @param nma_result NMA result
#' @param cost_data Treatment costs
#' @param threshold Willingness-to-pay
#' @param parameters Which parameters to calculate EVPPI for
#' @param n_sim Number of simulations
#' @return evppi object
#' @export
calculate_evppi <- function(nma_result, cost_data, threshold, parameters,
                           n_sim = 5000) {

  message(sprintf("Calculating EVPPI for %d parameter(s)...", length(parameters)))

  # Placeholder for simplified EVPPI calculation
  # Full implementation would use methods from Strong et al. (2014)

  # For now, calculate single-parameter EVPPI
  treatments <- rownames(nma_result$TE.random)
  reference <- nma_result$reference.group

  evppi_results <- lapply(parameters, function(param) {

    # Identify which treatment this parameter corresponds to
    if (param %in% treatments) {

      message(sprintf("  Calculating EVPPI for: %s", param))

      # Simulate with and without uncertainty in this parameter
      # With uncertainty (full model)
      effects <- nma_result$TE.random[, reference]
      vcov_full <- nma_result$Cov.random

      sim_full <- MASS::mvrnorm(n = n_sim, mu = effects, Sigma = vcov_full)

      # Without uncertainty in this parameter (hold constant)
      vcov_partial <- vcov_full
      param_idx <- which(treatments == param)
      vcov_partial[param_idx, ] <- 0
      vcov_partial[, param_idx] <- 0

      sim_partial <- MASS::mvrnorm(n = n_sim, mu = effects, Sigma = vcov_partial)

      # Calculate EVPI for full model
      evpi_full <- calculate_evpi_from_simulations(sim_full, cost_data, threshold,
                                                   treatments, reference)

      # Calculate EVPI for partial model
      evpi_partial <- calculate_evpi_from_simulations(sim_partial, cost_data, threshold,
                                                      treatments, reference)

      # EVPPI = difference
      evppi <- evpi_full - evppi_partial

      data.frame(
        Parameter = param,
        EVPPI = evppi,
        Percent_of_EVPI = (evppi / evpi_full) * 100,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Parameter = param,
        EVPPI = NA,
        Percent_of_EVPI = NA,
        stringsAsFactors = FALSE
      )
    }
  })

  evppi_df <- do.call(rbind, evppi_results)
  evppi_df <- evppi_df[order(-evppi_df$EVPPI), ]

  result <- list(
    evppi_table = evppi_df,
    threshold = threshold,
    parameters = parameters
  )

  class(result) <- c("evppi", "list")
  result
}

#' Helper: Calculate EVPI from simulations
#' @keywords internal
calculate_evpi_from_simulations <- function(simulated_effects, cost_data, threshold,
                                           treatments, reference) {

  # Merge costs
  treatment_data <- data.frame(
    Treatment = treatments,
    stringsAsFactors = FALSE
  )

  treatment_data <- merge(treatment_data, cost_data,
                         by.x = "Treatment", by.y = "treatment",
                         all.x = TRUE)

  ref_idx <- which(treatments == reference)
  ref_cost <- treatment_data$cost[ref_idx]

  # Calculate net benefits for each simulation
  simulation_nb <- matrix(NA, nrow = nrow(simulated_effects), ncol = length(treatments))

  for (i in seq_len(nrow(simulated_effects))) {
    sim_effects <- simulated_effects[i, ]
    nb <- (sim_effects - sim_effects[ref_idx]) * threshold -
          (treatment_data$cost - ref_cost)
    simulation_nb[i, ] <- nb
  }

  max_nb_per_sim <- apply(simulation_nb, 1, max)
  expected_nb_perfect <- mean(max_nb_per_sim)

  mean_nb <- colMeans(simulation_nb)
  expected_nb_current <- max(mean_nb)

  evpi <- expected_nb_perfect - expected_nb_current
  evpi
}

#' Cost-Effectiveness Acceptability Curve (CEAC)
#'
#' Plot probability each treatment is cost-effective across thresholds.
#'
#' @param nma_result NMA result
#' @param cost_data Treatment costs
#' @param threshold_range Range of thresholds
#' @param n_sim Number of simulations
#' @return ggplot2 object
#' @export
plot_ceac <- function(nma_result, cost_data, threshold_range = seq(0, 50000, by = 2000),
                     n_sim = 10000) {

  message(sprintf("Generating CEAC with %d simulations...", n_sim))

  treatments <- rownames(nma_result$TE.random)
  reference <- nma_result$reference.group

  effects <- nma_result$TE.random[, reference]
  vcov_matrix <- nma_result$Cov.random

  # Merge costs
  treatment_data <- data.frame(
    Treatment = treatments,
    stringsAsFactors = FALSE
  )

  treatment_data <- merge(treatment_data, cost_data,
                         by.x = "Treatment", by.y = "treatment",
                         all.x = TRUE)

  ref_idx <- which(treatments == reference)
  ref_cost <- treatment_data$cost[ref_idx]

  # Simulate effects
  simulated_effects <- MASS::mvrnorm(n = n_sim, mu = effects, Sigma = vcov_matrix)

  # For each threshold, calculate probability each treatment is optimal
  ceac_results <- lapply(threshold_range, function(threshold) {

    # Calculate net benefits for all simulations
    simulation_nb <- matrix(NA, nrow = n_sim, ncol = length(treatments))

    for (i in seq_len(n_sim)) {
      sim_effects <- simulated_effects[i, ]
      nb <- (sim_effects - sim_effects[ref_idx]) * threshold -
            (treatment_data$cost - ref_cost)
      simulation_nb[i, ] <- nb
    }

    # Find optimal treatment in each simulation
    optimal_per_sim <- apply(simulation_nb, 1, which.max)

    # Probability each treatment is optimal
    prob_optimal <- table(factor(optimal_per_sim, levels = seq_along(treatments))) / n_sim

    data.frame(
      Threshold = threshold,
      Treatment = treatments,
      Probability = as.vector(prob_optimal),
      stringsAsFactors = FALSE
    )
  })

  ceac_df <- do.call(rbind, ceac_results)

  # Plot
  p <- ggplot2::ggplot(ceac_df, ggplot2::aes(x = Threshold, y = Probability,
                                             color = Treatment)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_x_continuous(labels = scales::dollar) +
    ggplot2::labs(
      title = "Cost-Effectiveness Acceptability Curve (CEAC)",
      subtitle = sprintf("Probability of cost-effectiveness by willingness-to-pay (n=%d simulations)", n_sim),
      x = "Willingness-to-Pay Threshold",
      y = "Probability Cost-Effective",
      color = "Treatment"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  print(p)
  invisible(p)
}

#' Print Threshold Analysis
#'
#' @param x threshold_analysis object
#' @param ... Additional arguments
#' @export
print.threshold_analysis <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Threshold Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Threshold range: $%.0f to $%.0f\n",
             min(x$threshold_range), max(x$threshold_range)))
  cat(sprintf("Reference treatment: %s\n\n", x$reference))

  if (nrow(x$switch_points) > 0) {
    cat("Treatment Selection Switch Points:\n\n")
    print(x$switch_points, row.names = FALSE)
    cat("\n")
  } else {
    cat("No treatment switches detected in threshold range.\n")
    cat(sprintf("Optimal treatment throughout: %s\n\n",
               x$threshold_results$Optimal_Treatment[1]))
  }

  cat("Treatment Data:\n\n")
  print(x$treatment_data[, c("Treatment", "Benefit", "cost", "Incremental_Cost", "Incremental_Benefit")],
        row.names = FALSE, digits = 2)

  cat("\n")
  cat("Use plot() to visualize threshold analysis\n\n")

  invisible(x)
}

#' Plot Threshold Analysis
#'
#' @param x threshold_analysis object
#' @param ... Additional arguments
#' @export
plot.threshold_analysis <- function(x, ...) {

  p <- ggplot2::ggplot(x$threshold_results,
                      ggplot2::aes(x = Threshold, y = Net_Benefit, color = Optimal_Treatment)) +
    ggplot2::geom_line(size = 1.5) +
    ggplot2::geom_vline(xintercept = x$switch_points$Threshold, linetype = "dashed",
                       alpha = 0.5) +
    ggplot2::scale_x_continuous(labels = scales::dollar) +
    ggplot2::labs(
      title = "Threshold Analysis: Optimal Treatment by Willingness-to-Pay",
      subtitle = sprintf("Reference: %s", x$reference),
      x = "Willingness-to-Pay Threshold",
      y = "Net Benefit",
      color = "Optimal Treatment"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  if (nrow(x$switch_points) > 0) {
    # Add annotations for switch points
    for (i in seq_len(nrow(x$switch_points))) {
      p <- p + ggplot2::annotate("text",
                                 x = x$switch_points$Threshold[i],
                                 y = max(x$threshold_results$Net_Benefit) * 0.9,
                                 label = sprintf("Switch at\n$%.0f", x$switch_points$Threshold[i]),
                                 size = 3, alpha = 0.7)
    }
  }

  print(p)
  invisible(p)
}

#' Print EVPI Results
#'
#' @param x evpi object
#' @param ... Additional arguments
#' @export
print.evpi <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Expected Value of Perfect Information (EVPI)\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Willingness-to-pay threshold: $%.0f\n", x$threshold))
  cat(sprintf("Number of simulations: %d\n\n", x$n_sim))

  cat(sprintf("EVPI per patient: $%.2f\n\n", x$evpi))

  cat("Current Decision:\n")
  cat(sprintf("  Optimal treatment: %s\n", x$optimal_treatment))
  cat(sprintf("  Expected net benefit: $%.2f\n\n", x$expected_nb_current_info))

  cat("Probability Each Treatment is Optimal:\n")
  for (i in seq_along(x$treatment_names)) {
    cat(sprintf("  %s: %.1f%%\n", x$treatment_names[i], x$prob_optimal[i] * 100))
  }

  cat("\n")
  cat("Interpretation:\n")
  cat(sprintf("  We would pay up to $%.2f per patient to eliminate all uncertainty\n", x$evpi))
  cat("  about treatment effects and make the perfect decision.\n\n")

  cat("Use calculate_population_evpi() to scale to population level\n\n")

  invisible(x)
}

#' Print EVPPI Results
#'
#' @param x evppi object
#' @param ... Additional arguments
#' @export
print.evppi <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Expected Value of Perfect Parameter Information (EVPPI)\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Willingness-to-pay threshold: $%.0f\n\n", x$threshold))

  cat("EVPPI by Parameter:\n\n")
  print(x$evppi_table, row.names = FALSE, digits = 2)

  cat("\n")
  cat("Interpretation: These values represent how much we would pay per patient\n")
  cat("to eliminate uncertainty about each specific parameter.\n")
  cat("Parameters with higher EVPPI should be prioritized for further research.\n\n")

  invisible(x)
}
