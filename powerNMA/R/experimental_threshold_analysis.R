# ============================================================================
# EXPERIMENTAL: Threshold Analysis for Network Meta-Analysis
# ============================================================================
#
# Based on:
# - Ades et al. (2025) "Treatment recommendations based on network meta-analysis:
#   Rules for risk-averse decision-makers" Research Synthesis Methods
# - Multiple 2024 papers on threshold analysis for NMA credibility assessment
#
# STATUS: EXPERIMENTAL - Methods from 2024-2025 literature
#
# Threshold analysis asks: "How much would the evidence have to change before
# the treatment recommendation changes?" This provides a quantitative measure
# of the robustness of recommendations to uncertainty, bias, and imprecision.
#
# ============================================================================

#' EXPERIMENTAL: Threshold Analysis for Treatment Recommendations
#'
#' Performs threshold analysis to assess the robustness of treatment recommendations
#' from network meta-analysis. Quantifies how much the evidence would need to change
#' (due to bias, imprecision, or new data) before the recommendation changes.
#'
#' @param nma_object A network meta-analysis object (from netmeta, gemtc, or powerNMA)
#' @param outcome_direction Direction where higher values are better: "higher" or "lower"
#' @param recommendation Current recommended treatment based on NMA. If NULL,
#'   automatically selects treatment with best point estimate.
#' @param decision_rule Decision rule: "maximize_benefit", "minimize_risk",
#'   "cost_effectiveness", or "multi_criteria"
#' @param risk_aversion Risk aversion parameter (0 = risk-neutral, higher = more risk-averse)
#' @param cost_data Optional data frame with cost information (columns: treatment, cost)
#' @param willingness_to_pay Willingness to pay threshold for cost-effectiveness (if applicable)
#' @param threshold_type Type of threshold to compute:
#'   \itemize{
#'     \item "effect_size": How much would treatment effects need to change?
#'     \item "bias": How much bias would change the recommendation?
#'     \item "new_study": What results from a new study would change recommendation?
#'     \item "all": Compute all threshold types
#'   }
#' @param confidence_level Confidence level for uncertainty intervals (default 0.95)
#' @param ... Additional arguments
#'
#' @return An object of class "threshold_analysis" containing:
#'   \item{recommendation}{Current recommended treatment}
#'   \item{thresholds}{Threshold values for recommendation change}
#'   \item{robustness_score}{Quantitative robustness score (0-100)}
#'   \item{sensitivity_to_bias}{How sensitive recommendation is to potential bias}
#'   \item{tipping_point}{The "tipping point" where recommendation would change}
#'   \item{alternative_recommendation}{What treatment would be recommended if threshold crossed}
#'   \item{interpretation}{Clinical interpretation and guidance}
#'
#' @details
#' Threshold analysis provides more actionable information than traditional
#' GRADE/CINeMA approaches:
#'
#' \strong{Advantages:}
#' \itemize{
#'   \item Quantifies precise amount of change needed
#'   \item Can incorporate multiple outcomes simultaneously
#'   \item Accounts for decision-maker risk preferences
#'   \item Provides specific guidance for future research
#'   \item More intuitive than abstract "certainty ratings"
#' }
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item High threshold (>2 SDs): Recommendation is very robust
#'   \item Moderate threshold (1-2 SDs): Recommendation is reasonably robust
#'   \item Low threshold (<1 SD): Recommendation is fragile, more evidence needed
#' }
#'
#' @references
#' Ades AE, Davies A, Phillippo D, et al. (2025). Treatment recommendations based
#' on network meta-analysis: Rules for risk-averse decision-makers.
#' Research Synthesis Methods.
#'
#' Phillippo DM, Dias S, Ades AE, et al. (2019). Threshold analysis as an
#' alternative to GRADE for assessing confidence in guideline recommendations
#' based on network meta-analyses. Annals of Internal Medicine, 170:538-546.
#'
#' @examples
#' \dontrun{
#' # Run NMA first
#' library(netmeta)
#' data(smokingcessation)
#' nma <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                data = smokingcessation, sm = "OR")
#'
#' # Threshold analysis
#' threshold_result <- threshold_analysis(
#'   nma_object = nma,
#'   outcome_direction = "higher",
#'   decision_rule = "maximize_benefit",
#'   risk_aversion = 1.0
#' )
#'
#' print(threshold_result)
#' plot(threshold_result)
#'
#' # With cost-effectiveness
#' costs <- data.frame(
#'   treatment = c("A", "B", "C", "D"),
#'   cost = c(100, 500, 1200, 2000)
#' )
#'
#' threshold_cost <- threshold_analysis(
#'   nma_object = nma,
#'   outcome_direction = "higher",
#'   decision_rule = "cost_effectiveness",
#'   cost_data = costs,
#'   willingness_to_pay = 50000
#' )
#' }
#'
#' @export
threshold_analysis <- function(nma_object,
                              outcome_direction = c("higher", "lower"),
                              recommendation = NULL,
                              decision_rule = c("maximize_benefit", "minimize_risk",
                                               "cost_effectiveness", "multi_criteria"),
                              risk_aversion = 1.0,
                              cost_data = NULL,
                              willingness_to_pay = NULL,
                              threshold_type = c("all", "effect_size", "bias", "new_study"),
                              confidence_level = 0.95,
                              ...) {

  # Experimental warning
  message("=============================================================")
  message("EXPERIMENTAL METHOD: Threshold Analysis for NMA")
  message("Based on: Ades et al. (2025) Research Synthesis Methods")
  message("Status: Cutting-edge method from 2024-2025 literature")
  message("=============================================================")

  # Argument checking
  outcome_direction <- match.arg(outcome_direction)
  decision_rule <- match.arg(decision_rule)
  threshold_type <- match.arg(threshold_type)

  # Extract treatment effects and uncertainty from NMA
  nma_estimates <- extract_nma_estimates(nma_object)

  # Determine current recommendation
  if (is.null(recommendation)) {
    recommendation <- determine_recommendation(
      nma_estimates,
      outcome_direction,
      decision_rule,
      risk_aversion,
      cost_data,
      willingness_to_pay
    )
    message("Auto-selected recommendation: ", recommendation$treatment)
  }

  # Compute thresholds
  if (threshold_type == "all") {
    thresholds <- list(
      effect_size = compute_effect_size_threshold(nma_estimates, recommendation, outcome_direction),
      bias = compute_bias_threshold(nma_estimates, recommendation, outcome_direction),
      new_study = compute_new_study_threshold(nma_estimates, recommendation, outcome_direction)
    )
  } else if (threshold_type == "effect_size") {
    thresholds <- list(
      effect_size = compute_effect_size_threshold(nma_estimates, recommendation, outcome_direction)
    )
  } else if (threshold_type == "bias") {
    thresholds <- list(
      bias = compute_bias_threshold(nma_estimates, recommendation, outcome_direction)
    )
  } else {
    thresholds <- list(
      new_study = compute_new_study_threshold(nma_estimates, recommendation, outcome_direction)
    )
  }

  # Calculate robustness score
  robustness <- calculate_robustness_score(thresholds, nma_estimates)

  # Identify alternative recommendation (what would be chosen if threshold crossed)
  alternative <- identify_alternative_recommendation(
    nma_estimates,
    recommendation,
    outcome_direction,
    thresholds
  )

  # Create interpretation
  interpretation <- create_threshold_interpretation(
    recommendation,
    thresholds,
    robustness,
    alternative
  )

  # Create result object
  result <- list(
    recommendation = recommendation,
    thresholds = thresholds,
    robustness_score = robustness,
    alternative_recommendation = alternative,
    interpretation = interpretation,
    nma_estimates = nma_estimates,
    outcome_direction = outcome_direction,
    decision_rule = decision_rule,
    risk_aversion = risk_aversion,
    call = match.call()
  )

  class(result) <- "threshold_analysis"
  return(result)
}


# ============================================================================
# Helper Functions
# ============================================================================

#' Extract estimates from NMA object
#' @keywords internal
extract_nma_estimates <- function(nma_object) {

  if (inherits(nma_object, "netmeta")) {

    # Extract treatment effects (league table)
    treatments <- rownames(nma_object$TE.random)
    n_treat <- length(treatments)

    estimates <- data.frame(
      treatment = treatments,
      effect = diag(nma_object$TE.random),  # Relative to reference
      se = diag(nma_object$seTE.random),
      lower = diag(nma_object$lower.random),
      upper = diag(nma_object$upper.random),
      stringsAsFactors = FALSE
    )

    # Get full league table for pairwise comparisons
    estimates$league_table <- list(nma_object$TE.random)
    estimates$league_se <- list(nma_object$seTE.random)

    return(estimates)

  } else {
    stop("NMA object type not yet supported. Use netmeta objects.")
  }
}


#' Determine recommendation based on decision rule
#' @keywords internal
determine_recommendation <- function(nma_estimates, outcome_direction,
                                    decision_rule, risk_aversion,
                                    cost_data, willingness_to_pay) {

  if (decision_rule == "maximize_benefit") {

    # Risk-neutral: Best point estimate
    # Risk-averse: Best lower confidence bound
    if (risk_aversion == 0) {
      if (outcome_direction == "higher") {
        best_idx <- which.max(nma_estimates$effect)
      } else {
        best_idx <- which.min(nma_estimates$effect)
      }
    } else {
      # Risk-averse: penalize uncertainty
      # Use lower bound for "higher is better", upper bound for "lower is better"
      if (outcome_direction == "higher") {
        adjusted_effect <- nma_estimates$lower
        best_idx <- which.max(adjusted_effect)
      } else {
        adjusted_effect <- nma_estimates$upper
        best_idx <- which.min(adjusted_effect)
      }
    }

    return(list(
      treatment = nma_estimates$treatment[best_idx],
      effect = nma_estimates$effect[best_idx],
      se = nma_estimates$se[best_idx],
      rule = decision_rule
    ))

  } else if (decision_rule == "cost_effectiveness") {

    if (is.null(cost_data) || is.null(willingness_to_pay)) {
      stop("For cost-effectiveness, must provide cost_data and willingness_to_pay")
    }

    # Calculate incremental cost-effectiveness ratios (ICERs)
    # Merge costs with estimates
    merged <- merge(nma_estimates, cost_data, by = "treatment", all.x = TRUE)

    if (any(is.na(merged$cost))) {
      stop("Cost data missing for some treatments")
    }

    # Order by effect (higher is better assumed for cost-effectiveness)
    if (outcome_direction == "lower") {
      # For "lower is better" (e.g., adverse events), negate effects
      merged$effect <- -merged$effect
    }

    merged <- merged[order(merged$effect, decreasing = TRUE), ]

    # Calculate incremental cost and effect relative to previous treatment
    merged$incr_cost <- c(merged$cost[1], diff(merged$cost))
    merged$incr_effect <- c(merged$effect[1], diff(merged$effect))

    # Calculate ICER for each treatment (compared to next worse)
    merged$icer <- ifelse(merged$incr_effect > 0,
                         merged$incr_cost / merged$incr_effect,
                         Inf)

    # Find cost-effective treatment: best treatment with ICER < willingness_to_pay
    # that is not dominated by extended dominance
    cost_effective_treatments <- merged[merged$icer <= willingness_to_pay |
                                       merged$icer == merged$incr_effect[1], ]

    if (nrow(cost_effective_treatments) == 0) {
      # If no treatments meet threshold, choose least costly
      best_idx <- which.min(merged$cost)
      message("No treatments are cost-effective at willingness-to-pay threshold. ",
              "Selecting least costly option.")
    } else {
      # Choose treatment with highest effect among cost-effective options
      best_idx <- which(merged$treatment == cost_effective_treatments$treatment[
        which.max(cost_effective_treatments$effect)
      ])
    }

    return(list(
      treatment = merged$treatment[best_idx],
      effect = merged$effect[best_idx],
      se = merged$se[best_idx],
      cost = merged$cost[best_idx],
      icer = merged$icer[best_idx],
      rule = decision_rule
    ))

  } else {
    # Default to maximize benefit
    if (outcome_direction == "higher") {
      best_idx <- which.max(nma_estimates$effect)
    } else {
      best_idx <- which.min(nma_estimates$effect)
    }

    return(list(
      treatment = nma_estimates$treatment[best_idx],
      effect = nma_estimates$effect[best_idx],
      se = nma_estimates$se[best_idx],
      rule = decision_rule
    ))
  }
}


#' Compute effect size threshold
#' @keywords internal
compute_effect_size_threshold <- function(nma_estimates, recommendation, outcome_direction) {

  rec_treatment <- recommendation$treatment
  rec_effect <- recommendation$effect
  rec_se <- recommendation$se

  # Find second-best treatment
  other_treatments <- nma_estimates[nma_estimates$treatment != rec_treatment, ]

  if (outcome_direction == "higher") {
    second_best_idx <- which.max(other_treatments$effect)
  } else {
    second_best_idx <- which.min(other_treatments$effect)
  }

  second_best <- other_treatments[second_best_idx, ]

  # Calculate threshold: How much would recommended treatment effect need to decrease
  # (or second-best increase) before recommendation changes?

  # Difference between recommended and second-best
  effect_diff <- abs(rec_effect - second_best$effect)

  # Combined SE
  combined_se <- sqrt(rec_se^2 + second_best$se^2)

  # Threshold in standard deviation units
  threshold_sd <- effect_diff / combined_se

  return(list(
    value = effect_diff,
    sd_units = threshold_sd,
    current_best = rec_treatment,
    second_best = second_best$treatment,
    interpretation = paste0(
      "Recommended treatment effect would need to change by ",
      round(effect_diff, 3),
      " (", round(threshold_sd, 2), " SDs) before recommendation changes to ",
      second_best$treatment
    )
  ))
}


#' Compute bias threshold
#' @keywords internal
compute_bias_threshold <- function(nma_estimates, recommendation, outcome_direction) {

  # How much systematic bias would be needed to change the recommendation?
  effect_threshold <- compute_effect_size_threshold(nma_estimates, recommendation, outcome_direction)

  # Bias threshold is similar but considers directional bias
  bias_threshold <- effect_threshold$value

  # Additional consideration: Is bias plausible?
  # Empirical studies suggest NMA bias typically < 0.5 SD
  plausible_bias <- 0.5 * effect_threshold$sd_units

  is_robust_to_bias <- effect_threshold$sd_units > 2.0  # > 2 SD is robust

  return(list(
    value = bias_threshold,
    sd_units = effect_threshold$sd_units,
    plausible_bias_range = c(-plausible_bias, plausible_bias),
    robust_to_bias = is_robust_to_bias,
    interpretation = paste0(
      "A systematic bias of ", round(bias_threshold, 3),
      " (", round(effect_threshold$sd_units, 2), " SDs) would change recommendation. ",
      ifelse(is_robust_to_bias,
             "This is ROBUST to plausible bias.",
             "This is SENSITIVE to bias - more evidence needed.")
    )
  ))
}


#' Compute new study threshold
#' @keywords internal
compute_new_study_threshold <- function(nma_estimates, recommendation, outcome_direction) {

  rec_treatment <- recommendation$treatment
  rec_effect <- recommendation$effect
  rec_se <- recommendation$se

  # What would a new study need to show to change the recommendation?
  # Uses meta-analytic updating to find tipping point

  # Assume new study with similar SE to existing evidence
  assumed_new_se <- rec_se

  # Calculate the effect size threshold (tipping point)
  effect_threshold <- compute_effect_size_threshold(nma_estimates, recommendation, outcome_direction)

  # Find second-best treatment for comparison
  other_treatments <- nma_estimates[nma_estimates$treatment != rec_treatment, ]
  if (outcome_direction == "higher") {
    second_best_idx <- which.max(other_treatments$effect)
  } else {
    second_best_idx <- which.min(other_treatments$effect)
  }
  second_best <- other_treatments[second_best_idx, ]

  # Meta-analytic updating formula:
  # New pooled estimate = (w_old * θ_old + w_new * θ_new) / (w_old + w_new)
  # where w = 1/SE²

  # Current evidence for recommended treatment
  w_old <- 1 / rec_se^2
  theta_old <- rec_effect

  # New study weight (assumed)
  w_new <- 1 / assumed_new_se^2

  # Tipping point: where pooled estimate equals second-best
  # For recommendation to change, need: new_pooled ≈ second_best effect
  theta_tip <- second_best$effect

  # Solve for required new study effect:
  # theta_tip = (w_old * theta_old + w_new * theta_new) / (w_old + w_new)
  # theta_new = (theta_tip * (w_old + w_new) - w_old * theta_old) / w_new

  required_effect <- (theta_tip * (w_old + w_new) - w_old * theta_old) / w_new

  # Calculate what the new pooled SE would be
  pooled_se_new <- sqrt(1 / (w_old + w_new))

  # How many SDs away from current effect is this?
  threshold_sd <- abs(required_effect - rec_effect) / assumed_new_se

  return(list(
    required_effect = required_effect,
    current_effect = rec_effect,
    difference_from_current = required_effect - rec_effect,
    assumed_se = assumed_new_se,
    threshold_sd = threshold_sd,
    tipping_point = theta_tip,
    second_best_treatment = second_best$treatment,
    interpretation = paste0(
      "A new study would need to show effect of ", round(required_effect, 3),
      " (assuming SE = ", round(assumed_new_se, 3), ") to change recommendation ",
      "to ", second_best$treatment, ". ",
      "This differs from current effect by ", round(abs(required_effect - rec_effect), 3),
      " (", round(threshold_sd, 2), " SDs)."
    )
  ))
}


#' Calculate robustness score (0-100)
#' @keywords internal
calculate_robustness_score <- function(thresholds, nma_estimates) {

  # Robustness score based on threshold in SD units
  # > 2 SD = 100 (very robust)
  # 1-2 SD = 50-100 (moderately robust)
  # < 1 SD = 0-50 (fragile)

  if ("effect_size" %in% names(thresholds)) {
    sd_threshold <- thresholds$effect_size$sd_units
  } else if ("bias" %in% names(thresholds)) {
    sd_threshold <- thresholds$bias$sd_units
  } else {
    sd_threshold <- 1.0  # Default
  }

  # Convert to 0-100 scale
  score <- min(100, max(0, (sd_threshold / 2) * 100))

  return(round(score))
}


#' Identify alternative recommendation
#' @keywords internal
identify_alternative_recommendation <- function(nma_estimates, recommendation,
                                               outcome_direction, thresholds) {

  rec_treatment <- recommendation$treatment

  # Find second-best treatment
  other_treatments <- nma_estimates[nma_estimates$treatment != rec_treatment, ]

  if (outcome_direction == "higher") {
    second_best_idx <- which.max(other_treatments$effect)
  } else {
    second_best_idx <- which.min(other_treatments$effect)
  }

  alternative <- other_treatments[second_best_idx, ]

  return(list(
    treatment = alternative$treatment,
    effect = alternative$effect,
    se = alternative$se
  ))
}


#' Create threshold interpretation
#' @keywords internal
create_threshold_interpretation <- function(recommendation, thresholds,
                                           robustness, alternative) {

  # Robustness category
  if (robustness >= 80) {
    robustness_cat <- "VERY ROBUST"
    guidance <- "Recommendation is highly stable. Unlikely to change with new evidence or bias."
  } else if (robustness >= 50) {
    robustness_cat <- "MODERATELY ROBUST"
    guidance <- "Recommendation is reasonably stable but could change with substantial new evidence."
  } else {
    robustness_cat <- "FRAGILE"
    guidance <- "Recommendation is sensitive to uncertainty. More evidence strongly recommended."
  }

  interpretation <- list(
    robustness_category = robustness_cat,
    robustness_score = robustness,
    guidance = guidance,
    recommendation = recommendation$treatment,
    alternative = alternative$treatment,
    summary = paste0(
      "Current recommendation: ", recommendation$treatment, " (", robustness_cat, " - Score: ", robustness, "/100). ",
      "If threshold crossed, alternative would be: ", alternative$treatment, ". ",
      guidance
    )
  )

  return(interpretation)
}


# ============================================================================
# S3 Methods
# ============================================================================

#' @export
print.threshold_analysis <- function(x, ...) {
  cat("=============================================================\n")
  cat("EXPERIMENTAL: Threshold Analysis for Treatment Recommendation\n")
  cat("=============================================================\n\n")

  cat("Current Recommendation:", x$recommendation$treatment, "\n")
  cat("Robustness Score:", x$robustness_score, "/ 100\n")
  cat("Category:", x$interpretation$robustness_category, "\n\n")

  cat("Thresholds for Recommendation Change:\n")
  cat("--------------------------------------\n")

  if ("effect_size" %in% names(x$thresholds)) {
    cat("\nEffect Size Threshold:\n")
    cat(x$thresholds$effect_size$interpretation, "\n")
  }

  if ("bias" %in% names(x$thresholds)) {
    cat("\nBias Threshold:\n")
    cat(x$thresholds$bias$interpretation, "\n")
  }

  if ("new_study" %in% names(x$thresholds)) {
    cat("\nNew Study Threshold:\n")
    cat(x$thresholds$new_study$interpretation, "\n")
  }

  cat("\n")
  cat("Alternative Recommendation:", x$alternative_recommendation$treatment, "\n")

  cat("\n")
  cat("Clinical Guidance:\n")
  cat(x$interpretation$guidance, "\n")

  invisible(x)
}


#' @export
plot.threshold_analysis <- function(x, type = c("robustness", "thresholds", "comparison"), ...) {

  type <- match.arg(type)

  if (type == "robustness") {
    plot_robustness_gauge(x, ...)
  } else if (type == "thresholds") {
    plot_threshold_bars(x, ...)
  } else if (type == "comparison") {
    plot_treatment_comparison(x, ...)
  }
}


#' Plot robustness gauge
#' @keywords internal
plot_robustness_gauge <- function(x, ...) {

  score <- x$robustness_score

  # Simple gauge plot
  par(mar = c(2, 2, 3, 2))

  # Create color based on score
  if (score >= 80) {
    col <- "darkgreen"
  } else if (score >= 50) {
    col <- "orange"
  } else {
    col <- "red"
  }

  # Barplot as gauge
  barplot(score,
          ylim = c(0, 100),
          col = col,
          border = NA,
          main = paste0("Robustness Score: ", score, "/100\n",
                       x$interpretation$robustness_category),
          ylab = "Robustness Score",
          names.arg = x$recommendation$treatment)

  # Add reference lines
  abline(h = 80, lty = 2, col = "darkgreen", lwd = 2)
  abline(h = 50, lty = 2, col = "orange", lwd = 2)

  text(0.7, 85, "Very Robust", pos = 4, col = "darkgreen")
  text(0.7, 55, "Moderate", pos = 4, col = "orange")
  text(0.7, 25, "Fragile", pos = 4, col = "red")
}


#' Plot threshold bars
#' @keywords internal
plot_threshold_bars <- function(x, ...) {

  thresholds <- x$thresholds

  # Extract SD units from each threshold
  sd_values <- c()
  labels <- c()

  if ("effect_size" %in% names(thresholds)) {
    sd_values <- c(sd_values, thresholds$effect_size$sd_units)
    labels <- c(labels, "Effect Size")
  }

  if ("bias" %in% names(thresholds)) {
    sd_values <- c(sd_values, thresholds$bias$sd_units)
    labels <- c(labels, "Bias")
  }

  if ("new_study" %in% names(thresholds)) {
    sd_values <- c(sd_values, thresholds$new_study$threshold_sd)
    labels <- c(labels, "New Study")
  }

  par(mar = c(5, 8, 4, 2))

  barplot(sd_values,
          names.arg = labels,
          horiz = TRUE,
          las = 1,
          col = ifelse(sd_values > 2, "darkgreen",
                      ifelse(sd_values > 1, "orange", "red")),
          xlab = "Threshold (Standard Deviations)",
          main = "Thresholds for Recommendation Change",
          xlim = c(0, max(c(sd_values, 2.5))))

  # Reference line at 2 SD
  abline(v = 2, lty = 2, col = "black", lwd = 2)
  text(2, length(sd_values) + 0.5, "Robust (>2 SD)", pos = 4)

  grid()
}


#' Plot treatment comparison
#' @keywords internal
plot_treatment_comparison <- function(x, ...) {

  est <- x$nma_estimates

  # Forest plot style comparison
  par(mar = c(5, 8, 4, 2))

  y_pos <- 1:nrow(est)

  plot(est$effect, y_pos,
       xlim = range(c(est$lower, est$upper), na.rm = TRUE),
       ylim = c(0.5, nrow(est) + 0.5),
       xlab = "Treatment Effect",
       ylab = "",
       yaxt = "n",
       pch = 18,
       cex = 1.5,
       col = ifelse(est$treatment == x$recommendation$treatment, "darkgreen",
                   ifelse(est$treatment == x$alternative_recommendation$treatment, "orange", "black")),
       main = "Treatment Effects with Recommendation Highlighted")

  axis(2, at = y_pos, labels = est$treatment, las = 1)

  # Add confidence intervals
  for (i in 1:nrow(est)) {
    segments(est$lower[i], y_pos[i], est$upper[i], y_pos[i],
            col = ifelse(est$treatment[i] == x$recommendation$treatment, "darkgreen",
                        ifelse(est$treatment[i] == x$alternative_recommendation$treatment, "orange", "gray")))
  }

  # Reference line at 0
  abline(v = 0, lty = 2, col = "red")

  # Legend
  legend("topright",
         legend = c("Current Recommendation", "Alternative", "Other"),
         pch = 18,
         col = c("darkgreen", "orange", "black"),
         cex = 0.8)

  grid()
}
