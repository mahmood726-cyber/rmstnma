#' Simulation and Power Analysis for Network Meta-Analysis
#'
#' Tools for simulating NMA data and conducting power analysis to plan
#' future network meta-analyses.
#'
#' @name simulation_power
#' @references
#' Thorlund K, Mills EJ (2012). Sample size and power considerations in
#' network meta-analysis. Systematic Reviews, 1:41.
#'
#' Nikolakopoulou A, et al. (2014). Characteristics of networks of interventions:
#' a description of a database of 186 published networks. PLOS ONE, 9(1):e86754.
NULL

#' Simulate Network Meta-Analysis Data
#'
#' Generate realistic NMA data for power analysis, method comparison, or
#' demonstration purposes.
#'
#' @param n_treatments Number of treatments
#' @param n_studies Number of studies
#' @param network_structure Type: "star", "line", "complete", "random"
#' @param effect_size True effect sizes (vector or single value)
#' @param heterogeneity Heterogeneity (tau)
#' @param sample_size_range Range of sample sizes per study
#' @param outcome_type "continuous" or "binary"
#' @param inconsistency Include inconsistency?
#' @param inconsistency_magnitude Magnitude of inconsistency
#' @return Simulated NMA data frame
#' @export
#' @examples
#' \dontrun{
#' # Star network
#' data_star <- simulate_nma_data(
#'   n_treatments = 5,
#'   n_studies = 20,
#'   network_structure = "star",
#'   effect_size = 0.5
#' )
#'
#' # Complete network with heterogeneity
#' data_complete <- simulate_nma_data(
#'   n_treatments = 4,
#'   n_studies = 30,
#'   network_structure = "complete",
#'   heterogeneity = 0.3
#' )
#' }
simulate_nma_data <- function(n_treatments = 5,
                              n_studies = 20,
                              network_structure = c("star", "line", "complete", "random"),
                              effect_size = 0.5,
                              heterogeneity = 0.2,
                              sample_size_range = c(50, 200),
                              outcome_type = c("continuous", "binary"),
                              inconsistency = FALSE,
                              inconsistency_magnitude = 0.3) {

  network_structure <- match.arg(network_structure)
  outcome_type <- match.arg(outcome_type)

  message(sprintf("Simulating %s network with %d treatments and %d studies...",
                 network_structure, n_treatments, n_studies))

  # Generate treatment names
  treatments <- c("Placebo", paste0("Treatment_", LETTERS[1:(n_treatments - 1)]))

  # Generate network structure
  comparisons <- generate_network_structure(treatments, network_structure, n_studies)

  # True effect sizes
  if (length(effect_size) == 1) {
    true_effects <- rep(effect_size, n_treatments - 1)
  } else {
    true_effects <- effect_size
  }

  # Generate study data
  study_list <- list()

  for (i in seq_len(nrow(comparisons))) {
    study_id <- sprintf("Study_%03d", i)

    treat1 <- comparisons$treat1[i]
    treat2 <- comparisons$treat2[i]

    # Sample size
    n1 <- sample(sample_size_range[1]:sample_size_range[2], 1)
    n2 <- sample(sample_size_range[1]:sample_size_range[2], 1)

    # True effect for this comparison
    idx1 <- which(treatments == treat1)
    idx2 <- which(treatments == treat2)

    if (idx1 == 1) {
      true_effect <- true_effects[idx2 - 1]
    } else if (idx2 == 1) {
      true_effect <- -true_effects[idx1 - 1]
    } else {
      true_effect <- true_effects[idx2 - 1] - true_effects[idx1 - 1]
    }

    # Add study-level random effect
    study_effect <- true_effect + rnorm(1, 0, heterogeneity)

    # Add inconsistency if requested
    if (inconsistency && runif(1) < 0.2) {  # 20% of comparisons have inconsistency
      study_effect <- study_effect + rnorm(1, 0, inconsistency_magnitude)
    }

    # Generate outcome data
    if (outcome_type == "continuous") {
      TE <- study_effect
      seTE <- sqrt(1 / n1 + 1 / n2)
    } else {
      # Binary outcome
      p1 <- plogis(0)  # Baseline risk
      p2 <- plogis(study_effect)
      events1 <- rbinom(1, n1, p1)
      events2 <- rbinom(1, n2, p2)

      # Calculate log OR
      TE <- log((events2 / (n2 - events2)) / (events1 / (n1 - events1)))
      if (!is.finite(TE)) TE <- 0  # Handle zero events

      seTE <- sqrt(1 / events1 + 1 / (n1 - events1) + 1 / events2 + 1 / (n2 - events2))
      if (!is.finite(seTE)) seTE <- 1  # Handle zero events
    }

    study_list[[i]] <- data.frame(
      studlab = study_id,
      treat1 = treat1,
      treat2 = treat2,
      TE = TE,
      seTE = seTE,
      n1 = n1,
      n2 = n2,
      stringsAsFactors = FALSE
    )
  }

  data <- do.call(rbind, study_list)

  message(sprintf("✓ Simulated %d comparisons", nrow(data)))

  attr(data, "true_effects") <- true_effects
  attr(data, "heterogeneity") <- heterogeneity
  attr(data, "network_structure") <- network_structure

  data
}

#' Generate Network Structure
#' @keywords internal
generate_network_structure <- function(treatments, structure, n_studies) {

  n_treatments <- length(treatments)

  if (structure == "star") {
    # All treatments compared to reference (placebo)
    reference <- treatments[1]
    other_treatments <- treatments[-1]

    comparisons <- data.frame(
      treat1 = rep(reference, n_studies),
      treat2 = sample(other_treatments, n_studies, replace = TRUE),
      stringsAsFactors = FALSE
    )

  } else if (structure == "line") {
    # Treatments in a line: A-B-C-D-E
    pairs <- data.frame(
      treat1 = treatments[-n_treatments],
      treat2 = treatments[-1],
      stringsAsFactors = FALSE
    )

    comparisons <- pairs[sample(nrow(pairs), n_studies, replace = TRUE), ]

  } else if (structure == "complete") {
    # All possible pairwise comparisons
    pairs <- expand.grid(
      treat1 = treatments,
      treat2 = treatments,
      stringsAsFactors = FALSE
    )
    pairs <- pairs[as.character(pairs$treat1) < as.character(pairs$treat2), ]

    comparisons <- pairs[sample(nrow(pairs), n_studies, replace = TRUE), ]

  } else {
    # Random structure
    comparisons <- data.frame(
      treat1 = sample(treatments, n_studies, replace = TRUE),
      treat2 = sample(treatments, n_studies, replace = TRUE),
      stringsAsFactors = FALSE
    )

    # Remove self-comparisons
    comparisons <- comparisons[comparisons$treat1 != comparisons$treat2, ]

    # Ensure at least one study per treatment
    for (trt in treatments) {
      if (!any(comparisons$treat1 == trt | comparisons$treat2 == trt)) {
        other_trt <- sample(treatments[treatments != trt], 1)
        comparisons <- rbind(comparisons, data.frame(
          treat1 = trt,
          treat2 = other_trt,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  comparisons
}

#' Power Analysis for Network Meta-Analysis
#'
#' Calculate statistical power for detecting treatment effects in NMA.
#'
#' @param n_treatments Number of treatments
#' @param n_studies_per_comparison Expected studies per comparison
#' @param effect_size Expected effect size to detect
#' @param heterogeneity Expected heterogeneity
#' @param alpha Significance level
#' @param network_structure Network structure
#' @param n_simulations Number of simulations
#' @return power_analysis object
#' @export
#' @examples
#' \dontrun{
#' power_result <- nma_power_analysis(
#'   n_treatments = 5,
#'   n_studies_per_comparison = 3,
#'   effect_size = 0.5,
#'   heterogeneity = 0.2
#' )
#' print(power_result)
#' }
nma_power_analysis <- function(n_treatments = 5,
                               n_studies_per_comparison = 3,
                               effect_size = 0.5,
                               heterogeneity = 0.2,
                               alpha = 0.05,
                               network_structure = "star",
                               n_simulations = 1000) {

  message(sprintf("Running power analysis with %d simulations...", n_simulations))

  n_studies <- n_studies_per_comparison * (n_treatments - 1)

  # Run simulations
  significant_results <- 0
  effect_estimates <- numeric(n_simulations)
  ci_widths <- numeric(n_simulations)

  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)

  for (i in seq_len(n_simulations)) {
    # Simulate data
    data <- simulate_nma_data(
      n_treatments = n_treatments,
      n_studies = n_studies,
      network_structure = network_structure,
      effect_size = effect_size,
      heterogeneity = heterogeneity
    )

    # Run NMA
    if (requireNamespace("netmeta", quietly = TRUE)) {
      nma <- netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = data
      )

      # Check if treatment effect is significant
      # Test Treatment_A vs Placebo
      effect_estimate <- nma$TE.random["Treatment_A", "Placebo"]
      se_estimate <- nma$seTE.random["Treatment_A", "Placebo"]

      z_score <- abs(effect_estimate / se_estimate)
      p_value <- 2 * (1 - pnorm(z_score))

      if (p_value < alpha) {
        significant_results <- significant_results + 1
      }

      effect_estimates[i] <- effect_estimate
      ci_widths[i] <- 2 * 1.96 * se_estimate
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  power <- significant_results / n_simulations

  # Calculate bias
  bias <- mean(effect_estimates) - effect_size

  result <- list(
    power = power,
    n_treatments = n_treatments,
    n_studies = n_studies,
    n_studies_per_comparison = n_studies_per_comparison,
    effect_size = effect_size,
    heterogeneity = heterogeneity,
    alpha = alpha,
    network_structure = network_structure,
    n_simulations = n_simulations,
    mean_effect_estimate = mean(effect_estimates),
    bias = bias,
    mean_ci_width = mean(ci_widths),
    coverage = sum(effect_estimates - ci_widths / 2 <= effect_size &
                  effect_estimates + ci_widths / 2 >= effect_size) / n_simulations
  )

  class(result) <- c("nma_power", "list")

  message(sprintf("\n✓ Power = %.1f%%", power * 100))

  result
}

#' Sample Size Calculation for NMA
#'
#' Calculate required sample size to achieve desired power.
#'
#' @param n_treatments Number of treatments
#' @param effect_size Effect size to detect
#' @param heterogeneity Expected heterogeneity
#' @param desired_power Desired statistical power
#' @param alpha Significance level
#' @param network_structure Network structure
#' @return Sample size recommendation
#' @export
calculate_nma_sample_size <- function(n_treatments = 5,
                                      effect_size = 0.5,
                                      heterogeneity = 0.2,
                                      desired_power = 0.80,
                                      alpha = 0.05,
                                      network_structure = "star") {

  message("Calculating required sample size...")

  # Try different sample sizes
  studies_per_comparison <- seq(2, 20, by = 2)
  powers <- numeric(length(studies_per_comparison))

  for (i in seq_along(studies_per_comparison)) {
    message(sprintf("  Testing %d studies per comparison...",
                   studies_per_comparison[i]))

    power_result <- nma_power_analysis(
      n_treatments = n_treatments,
      n_studies_per_comparison = studies_per_comparison[i],
      effect_size = effect_size,
      heterogeneity = heterogeneity,
      alpha = alpha,
      network_structure = network_structure,
      n_simulations = 500  # Reduced for speed
    )

    powers[i] <- power_result$power

    if (powers[i] >= desired_power) {
      break
    }
  }

  # Interpolate if needed
  if (max(powers) < desired_power) {
    required_studies <- NA
    message("\n⚠ Desired power not achieved with tested sample sizes")
  } else {
    idx <- which(powers >= desired_power)[1]
    required_studies <- studies_per_comparison[idx]
  }

  result <- list(
    required_studies_per_comparison = required_studies,
    total_studies = required_studies * (n_treatments - 1),
    achieved_power = powers[which.min(abs(powers - desired_power))],
    desired_power = desired_power,
    effect_size = effect_size,
    heterogeneity = heterogeneity,
    tested_samples = studies_per_comparison[seq_along(powers)],
    achieved_powers = powers
  )

  class(result) <- c("nma_sample_size", "list")

  if (!is.na(required_studies)) {
    message(sprintf("\n✓ Required: %d studies per comparison (%d total studies)",
                   required_studies, result$total_studies))
  }

  result
}

#' Optimal Allocation of Studies
#'
#' Determine optimal allocation of studies across comparisons to maximize power.
#'
#' @param n_treatments Number of treatments
#' @param total_studies Total number of available studies
#' @param effect_sizes Expected effect sizes for each treatment vs reference
#' @param reference Reference treatment
#' @return Optimal study allocation
#' @export
optimize_study_allocation <- function(n_treatments = 5,
                                     total_studies = 30,
                                     effect_sizes = NULL,
                                     reference = "Placebo") {

  message("Optimizing study allocation across comparisons...")

  treatments <- c(reference, paste0("Treatment_", LETTERS[1:(n_treatments - 1)]))

  if (is.null(effect_sizes)) {
    # Equal allocation
    allocation <- rep(floor(total_studies / (n_treatments - 1)), n_treatments - 1)
    remaining <- total_studies - sum(allocation)
    if (remaining > 0) {
      allocation[1:remaining] <- allocation[1:remaining] + 1
    }

  } else {
    # Allocate proportional to precision needed (inverse of effect size)
    # Smaller effects need more studies
    weights <- 1 / abs(effect_sizes)
    weights <- weights / sum(weights)

    allocation <- round(weights * total_studies)

    # Adjust to match exactly
    diff <- total_studies - sum(allocation)
    if (diff != 0) {
      idx <- which.max(weights)
      allocation[idx] <- allocation[idx] + diff
    }
  }

  comparison_names <- paste(treatments[-1], "vs", reference)

  result <- data.frame(
    Comparison = comparison_names,
    Allocated_Studies = allocation,
    Percentage = round(allocation / total_studies * 100, 1),
    stringsAsFactors = FALSE
  )

  message("\nOptimal allocation:")
  print(result, row.names = FALSE)

  result
}

#' Print Power Analysis Results
#' @param x nma_power object
#' @param ... Additional arguments
#' @export
print.nma_power <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Network Meta-Analysis Power Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Study Design:\n")
  cat(sprintf("  • Treatments: %d\n", x$n_treatments))
  cat(sprintf("  • Total studies: %d\n", x$n_studies))
  cat(sprintf("  • Studies per comparison: %d\n", x$n_studies_per_comparison))
  cat(sprintf("  • Network structure: %s\n\n", x$network_structure))

  cat("Parameters:\n")
  cat(sprintf("  • True effect size: %.3f\n", x$effect_size))
  cat(sprintf("  • Heterogeneity (τ): %.3f\n", x$heterogeneity))
  cat(sprintf("  • Significance level: %.3f\n\n", x$alpha))

  cat("Results:\n")
  cat(sprintf("  • Statistical power: %.1f%%\n", x$power * 100))
  cat(sprintf("  • Mean effect estimate: %.3f\n", x$mean_effect_estimate))
  cat(sprintf("  • Bias: %.3f\n", x$bias))
  cat(sprintf("  • Mean CI width: %.3f\n", x$mean_ci_width))
  cat(sprintf("  • Coverage: %.1f%%\n", x$coverage * 100))
  cat(sprintf("  • Simulations: %d\n\n", x$n_simulations))

  if (x$power < 0.80) {
    cat("⚠ Power is below the conventional threshold of 80%\n")
    cat("  Consider increasing sample size or number of studies\n\n")
  } else {
    cat("✓ Power is adequate (≥80%)\n\n")
  }

  invisible(x)
}

#' Print Sample Size Results
#' @param x nma_sample_size object
#' @param ... Additional arguments
#' @export
print.nma_sample_size <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  NMA Sample Size Calculation\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Target Parameters:\n")
  cat(sprintf("  • Desired power: %.1f%%\n", x$desired_power * 100))
  cat(sprintf("  • Effect size to detect: %.3f\n", x$effect_size))
  cat(sprintf("  • Expected heterogeneity: %.3f\n\n", x$heterogeneity))

  if (!is.na(x$required_studies_per_comparison)) {
    cat("Required Sample Size:\n")
    cat(sprintf("  • Studies per comparison: %d\n", x$required_studies_per_comparison))
    cat(sprintf("  • Total studies needed: %d\n", x$total_studies))
    cat(sprintf("  • Achieved power: %.1f%%\n\n", x$achieved_power * 100))
  } else {
    cat("⚠ Could not achieve desired power with tested sample sizes\n")
    cat("  Consider:\n")
    cat("    - Increasing effect size expectations\n")
    cat("    - Accepting lower power\n")
    cat("    - Reducing heterogeneity\n\n")
  }

  cat("Power Curve:\n")
  power_df <- data.frame(
    Studies_Per_Comparison = x$tested_samples,
    Power = sprintf("%.1f%%", x$achieved_powers * 100)
  )
  print(power_df, row.names = FALSE)

  cat("\n")

  invisible(x)
}

#' Plot Power Analysis
#' @param x nma_power or nma_sample_size object
#' @param ... Additional arguments
#' @export
plot.nma_power <- function(x, ...) {

  if (inherits(x, "nma_sample_size")) {
    # Power curve
    plot(x$tested_samples, x$achieved_powers * 100,
         type = "b",
         pch = 19,
         col = "#2c7fb8",
         lwd = 2,
         xlab = "Studies per Comparison",
         ylab = "Power (%)",
         main = "Power Curve for NMA Design",
         ylim = c(0, 100))

    abline(h = x$desired_power * 100, col = "red", lty = 2, lwd = 2)
    abline(h = 80, col = "gray50", lty = 3)

    if (!is.na(x$required_studies_per_comparison)) {
      abline(v = x$required_studies_per_comparison, col = "red", lty = 2, lwd = 2)
      text(x$required_studies_per_comparison, 10,
           sprintf("Required: %d", x$required_studies_per_comparison),
           pos = 4, col = "red")
    }

    legend("bottomright",
           legend = c(sprintf("Target: %.0f%%", x$desired_power * 100),
                     "Conventional: 80%"),
           col = c("red", "gray50"),
           lty = c(2, 3),
           lwd = 2)

  } else {
    # Power analysis visualization
    par(mfrow = c(1, 2))

    # Power gauge
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "",
         main = "Statistical Power")

    rect(0.2, 0.3, 0.8, 0.5, col = "gray90", border = NA)
    rect(0.2, 0.3, 0.2 + x$power * 0.6, 0.5,
         col = ifelse(x$power >= 0.80, "#1a9850", "#d73027"),
         border = NA)

    text(0.5, 0.6, sprintf("%.1f%%", x$power * 100),
         cex = 2, font = 2)

    text(0.5, 0.15, sprintf("%d studies | %.2f effect size",
                           x$n_studies, x$effect_size),
         cex = 0.8)

    # Coverage
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "",
         main = "95% CI Coverage")

    rect(0.2, 0.3, 0.8, 0.5, col = "gray90", border = NA)
    rect(0.2, 0.3, 0.2 + x$coverage * 0.6, 0.5,
         col = ifelse(x$coverage >= 0.90, "#1a9850", "#d73027"),
         border = NA)

    text(0.5, 0.6, sprintf("%.1f%%", x$coverage * 100),
         cex = 2, font = 2)

    text(0.5, 0.15, sprintf("Bias: %.3f", x$bias),
         cex = 0.8)

    par(mfrow = c(1, 1))
  }
}
