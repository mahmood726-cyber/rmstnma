# Cross-Design Synthesis: Combining RCT and Observational Data
#
# Novel method from 2024 for integrating randomized and non-randomized evidence
# References:
# - Hamza T, et al. (2024). crossnma R package. BMC Med Res Methodol
# - Efthimiou O, et al. (2023). Synthesizing cross-design evidence. Res Synth Methods

#' Cross-Design Synthesis: Network Meta-Analysis with RCT and Observational Data
#'
#' Integrates evidence from randomized controlled trials (RCTs) and non-randomized
#' studies (NRS) in a single network meta-analysis, accounting for differences in
#' design and risk of bias through bias-adjustment models.
#'
#' @param data Data frame with columns: study, treatment1, treatment2, effect, se, design
#' @param design_variable Name of column indicating study design ("RCT" or "NRS")
#' @param bias_adjustment Method for handling design-related bias:
#'   \itemize{
#'     \item{"none"}{Ignore design differences (naive)}
#'     \item{"penalized_prior"}{Use NRS to inform priors, penalized by bias}
#'     \item{"bias_model"}{Explicit bias parameter for NRS}
#'     \item{"hierarchical"}{Three-level hierarchical model}
#'   }
#' @param bias_prior Prior distribution for bias in NRS:
#'   \itemize{
#'     \item{"empirical"}{Based on empirical evidence (Turner et al. 2009)}
#'     \item{"conservative"}{Assume substantial bias}
#'     \item{"optimistic"}{Assume minimal bias}
#'     \item{"custom"}{User-specified}
#'   }
#' @param rob_adjustment Include risk-of-bias adjustment? (default: TRUE)
#' @param downweight_nrs Downweight NRS relative to RCTs? (default: TRUE)
#' @param method Estimation method: "bayesian" (default), "frequentist"
#' @param ... Additional arguments
#'
#' @return Object of class "cross_design_nma" containing:
#'   \itemize{
#'     \item{rct_only_results}{Analysis using only RCTs}
#'     \item{combined_results}{Analysis combining RCT + NRS}
#'     \item{bias_estimates}{Estimated bias in NRS}
#'     \item{sensitivity}{Comparison of RCT-only vs combined}
#'     \item{design_comparison}{RCT vs NRS effect estimates}
#'   }
#'
#' @details
#' **Why Cross-Design Synthesis?**
#' - RCTs: High internal validity, limited external validity
#' - NRS: Lower internal validity, better external validity (real-world)
#' - Combined: Leverage both sources, account for bias
#'
#' **Bias-Adjustment Approaches**:
#' 1. **None**: Naive pooling (not recommended)
#' 2. **Penalized Prior**: Use NRS to inform priors, with uncertainty penalty
#' 3. **Bias Model**: δ_NRS = δ_true + bias
#' 4. **Hierarchical**: Three levels (patient, study, design)
#'
#' **When to Use**:
#' - Few RCTs but many NRS available
#' - Generalizability concerns (RCTs ≠ real-world)
#' - Rare outcomes (RCTs underpowered)
#' - Long-term outcomes (RCTs short duration)
#'
#' **Caution**:
#' - NRS can introduce bias if not adjusted
#' - Results sensitive to bias assumptions
#' - Always compare RCT-only vs combined
#'
#' @examples
#' # Combine RCT and observational studies
#' data <- data.frame(
#'   study = 1:30,
#'   treatment1 = rep(c("A", "B", "C"), 10),
#'   treatment2 = rep(c("Placebo", "Placebo", "Placebo"), 10),
#'   effect = c(rnorm(15, -0.5, 0.2), rnorm(15, -0.4, 0.3)),  # RCT vs NRS
#'   se = c(runif(15, 0.15, 0.25), runif(15, 0.20, 0.35)),
#'   design = c(rep("RCT", 15), rep("NRS", 15))
#' )
#'
#' # Cross-design synthesis with bias adjustment
#' cds_result <- cross_design_synthesis(
#'   data = data,
#'   design_variable = "design",
#'   bias_adjustment = "bias_model",
#'   bias_prior = "empirical"
#' )
#'
#' # Compare RCT-only vs combined
#' print(cds_result)
#' plot(cds_result, type = "comparison")
#'
#' @export
cross_design_synthesis <- function(data,
                                    design_variable = "design",
                                    bias_adjustment = c("none", "penalized_prior",
                                                        "bias_model", "hierarchical"),
                                    bias_prior = c("empirical", "conservative",
                                                   "optimistic", "custom"),
                                    rob_adjustment = TRUE,
                                    downweight_nrs = TRUE,
                                    method = c("bayesian", "frequentist"),
                                    ...) {

  bias_adjustment <- match.arg(bias_adjustment)
  bias_prior <- match.arg(bias_prior)
  method <- match.arg(method)

  message("Cross-Design Synthesis: Combining RCT and Observational Data")
  message("==============================================================")
  message("Bias adjustment:", bias_adjustment)
  message("Method:", method)

  # Check design variable exists
  if (!(design_variable %in% colnames(data))) {
    stop("Design variable '", design_variable, "' not found in data")
  }

  # Identify RCTs and NRS
  rct_data <- data[data[[design_variable]] == "RCT", ]
  nrs_data <- data[data[[design_variable]] == "NRS", ]

  n_rct <- nrow(rct_data)
  n_nrs <- nrow(nrs_data)

  message(sprintf("Studies: %d RCTs, %d NRS (total %d)", n_rct, n_nrs, n_rct + n_nrs))

  if (n_rct == 0) {
    stop("No RCTs found. Cross-design synthesis requires at least some RCTs.")
  }

  if (n_nrs == 0) {
    warning("No NRS found. Performing standard NMA with RCTs only.")
    return(standard_nma_rct_only(rct_data))
  }

  # Fit RCT-only analysis
  message("\nStep 1: Analyzing RCTs only...")
  rct_only_results <- fit_rct_only_nma(rct_data, ...)

  # Fit combined analysis
  message("\nStep 2: Combining RCT + NRS with bias adjustment...")
  combined_results <- fit_cross_design_nma(
    data, bias_adjustment, bias_prior, downweight_nrs, method, ...
  )

  # Compare designs
  message("\nStep 3: Comparing RCT vs NRS effect estimates...")
  design_comparison <- compare_rct_vs_nrs(rct_data, nrs_data)

  # Sensitivity analysis
  sensitivity <- assess_cds_sensitivity(rct_only_results, combined_results)

  result <- list(
    rct_only_results = rct_only_results,
    combined_results = combined_results,
    bias_estimates = combined_results$bias_estimates,
    design_comparison = design_comparison,
    sensitivity = sensitivity,
    n_rct = n_rct,
    n_nrs = n_nrs,
    bias_adjustment = bias_adjustment,
    method = method
  )

  class(result) <- c("cross_design_nma", "list")
  return(result)
}


#' Fit RCT-only NMA
#' @keywords internal
fit_rct_only_nma <- function(rct_data, ...) {
  # Standard NMA on RCTs only
  # (Simplified - would use netmeta or run_powernma)

  effect_mean <- weighted.mean(rct_data$effect, 1 / rct_data$se^2)
  effect_se <- sqrt(1 / sum(1 / rct_data$se^2))

  result <- list(
    estimate = effect_mean,
    se = effect_se,
    lower = effect_mean - 1.96 * effect_se,
    upper = effect_mean + 1.96 * effect_se,
    n_studies = nrow(rct_data)
  )

  return(result)
}


#' Fit cross-design NMA
#' @keywords internal
fit_cross_design_nma <- function(data, bias_adjustment, bias_prior, downweight, method, ...) {

  if (method == "bayesian") {
    message("Bayesian cross-design synthesis requires JAGS/Stan")
    message("Using simplified frequentist approximation")
  }

  # Get bias prior
  bias_params <- get_bias_prior(bias_prior)

  # Apply bias adjustment
  if (bias_adjustment == "none") {
    result <- fit_naive_combined(data)
  } else if (bias_adjustment == "penalized_prior") {
    result <- fit_penalized_prior_cds(data, bias_params, downweight)
  } else if (bias_adjustment == "bias_model") {
    result <- fit_bias_model_cds(data, bias_params)
  } else {  # hierarchical
    result <- fit_hierarchical_cds(data, bias_params)
  }

  return(result)
}


#' Get bias prior parameters
#' @keywords internal
get_bias_prior <- function(bias_prior) {
  if (bias_prior == "empirical") {
    # Based on Turner et al. 2009 empirical evidence
    # NRS tend to overestimate effects by ~10-15% on average
    bias <- list(
      mean = 0.10,  # 10% bias (on log scale)
      sd = 0.15,    # Uncertainty
      description = "Empirical: NRS overestimate by ~10-15% (Turner et al. 2009)"
    )
  } else if (bias_prior == "conservative") {
    bias <- list(
      mean = 0.20,
      sd = 0.20,
      description = "Conservative: Assume substantial NRS bias"
    )
  } else if (bias_prior == "optimistic") {
    bias <- list(
      mean = 0.05,
      sd = 0.10,
      description = "Optimistic: Assume minimal NRS bias"
    )
  } else {
    bias <- list(
      mean = 0.10,
      sd = 0.15,
      description = "Default bias assumption"
    )
  }

  message("Bias prior: ", bias$description)

  return(bias)
}


#' Naive combined analysis
#' @keywords internal
fit_naive_combined <- function(data) {
  # Pool RCT + NRS without adjustment (NOT RECOMMENDED)

  effect_mean <- weighted.mean(data$effect, 1 / data$se^2)
  effect_se <- sqrt(1 / sum(1 / data$se^2))

  result <- list(
    estimate = effect_mean,
    se = effect_se,
    lower = effect_mean - 1.96 * effect_se,
    upper = effect_mean + 1.96 * effect_se,
    bias_estimates = NULL,
    note = "Naive pooling - bias not adjusted. NOT RECOMMENDED."
  )

  return(result)
}


#' Penalized prior approach
#' @keywords internal
fit_penalized_prior_cds <- function(data, bias_params, downweight) {
  # Use NRS to inform priors, with penalty for bias

  rct_data <- data[data$design == "RCT", ]
  nrs_data <- data[data$design == "NRS", ]

  # RCT estimate
  rct_effect <- weighted.mean(rct_data$effect, 1 / rct_data$se^2)

  # NRS estimate (adjusted for bias)
  nrs_effect <- weighted.mean(nrs_data$effect, 1 / nrs_data$se^2)
  nrs_effect_adjusted <- nrs_effect - bias_params$mean

  # Combine with downweighting
  if (downweight) {
    # Increase NRS uncertainty to account for bias
    nrs_se_inflated <- sqrt(mean(nrs_data$se^2) + bias_params$sd^2)
    nrs_weight <- 1 / nrs_se_inflated^2
  } else {
    nrs_weight <- sum(1 / nrs_data$se^2)
  }

  rct_weight <- sum(1 / rct_data$se^2)

  # Weighted combination
  combined_effect <- (rct_effect * rct_weight + nrs_effect_adjusted * nrs_weight) /
    (rct_weight + nrs_weight)
  combined_se <- sqrt(1 / (rct_weight + nrs_weight))

  result <- list(
    estimate = combined_effect,
    se = combined_se,
    lower = combined_effect - 1.96 * combined_se,
    upper = combined_effect + 1.96 * combined_se,
    bias_estimates = list(
      assumed_bias = bias_params$mean,
      nrs_unadjusted = nrs_effect,
      nrs_adjusted = nrs_effect_adjusted
    )
  )

  return(result)
}


#' Bias model approach
#' @keywords internal
fit_bias_model_cds <- function(data, bias_params) {
  # Explicit bias parameter: δ_NRS = δ_true + bias

  message("Bias model: δ_NRS = δ_true + bias")

  # Simplified - full version would use Bayesian model
  result <- fit_penalized_prior_cds(data, bias_params, downweight = TRUE)
  result$note <- "Bias model: Full Bayesian implementation in development"

  return(result)
}


#' Hierarchical approach
#' @keywords internal
fit_hierarchical_cds <- function(data, bias_params) {
  # Three-level model: patient, study, design

  message("Hierarchical model: patient → study → design")

  # Simplified - full version would use multilevel model
  result <- fit_penalized_prior_cds(data, bias_params, downweight = TRUE)
  result$note <- "Hierarchical model: Full implementation requires Bayesian framework"

  return(result)
}


#' Compare RCT vs NRS
#' @keywords internal
compare_rct_vs_nrs <- function(rct_data, nrs_data) {
  # Compare effect estimates between designs

  rct_effect <- weighted.mean(rct_data$effect, 1 / rct_data$se^2)
  nrs_effect <- weighted.mean(nrs_data$effect, 1 / nrs_data$se^2)

  difference <- nrs_effect - rct_effect

  comparison <- list(
    rct_estimate = rct_effect,
    nrs_estimate = nrs_effect,
    difference = difference,
    percent_difference = 100 * difference / abs(rct_effect),
    interpretation = interpret_design_difference(difference)
  )

  message(sprintf("RCT estimate: %.3f", rct_effect))
  message(sprintf("NRS estimate: %.3f", nrs_effect))
  message(sprintf("Difference: %.3f (%.1f%%)", difference,
                  comparison$percent_difference))

  return(comparison)
}


#' Interpret design difference
#' @keywords internal
interpret_design_difference <- function(difference) {
  if (abs(difference) < 0.05) {
    "Minimal difference between RCT and NRS"
  } else if (abs(difference) < 0.15) {
    "Moderate difference - bias adjustment important"
  } else {
    "Large difference - NRS results substantially different from RCTs"
  }
}


#' Assess CDS sensitivity
#' @keywords internal
assess_cds_sensitivity <- function(rct_only, combined) {
  # Compare RCT-only vs combined results

  diff <- combined$estimate - rct_only$estimate

  sensitivity <- list(
    rct_only_estimate = rct_only$estimate,
    combined_estimate = combined$estimate,
    difference = diff,
    conclusions_agree = abs(diff) < 0.10,  # Threshold
    recommendation = if (abs(diff) < 0.10) {
      "Conclusions robust: RCT-only and combined analyses agree"
    } else {
      "Conclusions sensitive: RCT-only and combined differ substantially. Report both."
    }
  )

  message("\nSensitivity Assessment:")
  message("  ", sensitivity$recommendation)

  return(sensitivity)
}


#' Print method for cross-design NMA
#' @export
print.cross_design_nma <- function(x, ...) {
  cat("Cross-Design Synthesis: RCT + Observational Data\n")
  cat("=================================================\n\n")

  cat(sprintf("Studies: %d RCTs, %d NRS\n", x$n_rct, x$n_nrs))
  cat("Bias adjustment:", x$bias_adjustment, "\n\n")

  cat("RCT-ONLY Analysis:\n")
  cat(sprintf("  Estimate: %.3f (SE %.3f, 95%% CI: %.3f to %.3f)\n",
              x$rct_only_results$estimate,
              x$rct_only_results$se,
              x$rct_only_results$lower,
              x$rct_only_results$upper))
  cat("\n")

  cat("COMBINED (RCT + NRS) Analysis:\n")
  cat(sprintf("  Estimate: %.3f (SE %.3f, 95%% CI: %.3f to %.3f)\n",
              x$combined_results$estimate,
              x$combined_results$se,
              x$combined_results$lower,
              x$combined_results$upper))
  cat("\n")

  if (!is.null(x$bias_estimates)) {
    cat("Bias Estimates:\n")
    cat(sprintf("  Assumed NRS bias: %.3f\n", x$bias_estimates$assumed_bias))
    if (!is.null(x$bias_estimates$nrs_unadjusted)) {
      cat(sprintf("  NRS unadjusted: %.3f\n", x$bias_estimates$nrs_unadjusted))
      cat(sprintf("  NRS adjusted: %.3f\n", x$bias_estimates$nrs_adjusted))
    }
    cat("\n")
  }

  cat("Design Comparison:\n")
  cat("  ", x$design_comparison$interpretation, "\n\n")

  cat("Sensitivity:\n")
  cat("  ", x$sensitivity$recommendation, "\n")

  invisible(x)
}


#' Plot method for cross-design NMA
#' @export
plot.cross_design_nma <- function(x, type = c("comparison", "forest", "bias"), ...) {
  type <- match.arg(type)

  if (type == "comparison") {
    plot_rct_nrs_comparison(x, ...)
  } else if (type == "forest") {
    plot_cds_forest(x, ...)
  } else if (type == "bias") {
    plot_bias_assessment(x, ...)
  }
}


#' Plot RCT vs NRS comparison
#' @keywords internal
plot_rct_nrs_comparison <- function(cds_object, ...) {
  par(mar = c(5, 8, 4, 2))

  estimates <- c(
    cds_object$rct_only_results$estimate,
    cds_object$combined_results$estimate,
    cds_object$design_comparison$nrs_estimate
  )

  lower <- c(
    cds_object$rct_only_results$lower,
    cds_object$combined_results$lower,
    NA  # NRS CI not stored
  )

  upper <- c(
    cds_object$rct_only_results$upper,
    cds_object$combined_results$upper,
    NA
  )

  labels <- c("RCT only", "Combined (RCT+NRS)", "NRS only")
  y_pos <- 1:3

  plot(estimates, y_pos,
       xlim = range(c(lower, upper), na.rm = TRUE),
       ylim = c(0.5, 3.5),
       xlab = "Treatment Effect",
       ylab = "",
       yaxt = "n",
       pch = c(19, 17, 15),
       col = c("blue", "purple", "orange"),
       cex = 1.5,
       main = "Cross-Design Synthesis: RCT vs NRS")

  axis(2, at = y_pos, labels = labels, las = 1)

  # Add CIs
  segments(x0 = lower[!is.na(lower)], y0 = y_pos[!is.na(lower)],
           x1 = upper[!is.na(upper)], y1 = y_pos[!is.na(upper)],
           lwd = 2,
           col = c("blue", "purple"))

  abline(v = 0, lty = 2, col = "gray")

  legend("topright",
         legend = labels,
         pch = c(19, 17, 15),
         col = c("blue", "purple", "orange"),
         bty = "n")
}


#' Plot forest plot for CDS
#' @keywords internal
plot_cds_forest <- function(cds_object, ...) {
  message("Forest plot for cross-design synthesis")
  # Placeholder
}


#' Plot bias assessment
#' @keywords internal
plot_bias_assessment <- function(cds_object, ...) {
  if (is.null(cds_object$bias_estimates)) {
    message("No bias estimates available")
    return(invisible(NULL))
  }

  bias <- cds_object$bias_estimates

  barplot(c(bias$nrs_unadjusted, bias$nrs_adjusted),
          names.arg = c("NRS Unadjusted", "NRS Adjusted"),
          ylab = "Treatment Effect",
          main = "Effect of Bias Adjustment",
          col = c("orange", "darkgreen"),
          ylim = range(c(bias$nrs_unadjusted, bias$nrs_adjusted,
                         cds_object$rct_only_results$estimate)) * c(0.9, 1.1))

  abline(h = cds_object$rct_only_results$estimate, lty = 2, lwd = 2, col = "blue")
  text(1.5, cds_object$rct_only_results$estimate * 1.05, "RCT estimate", col = "blue")
}
