#' Heterogeneity Assessment for Network Meta-Analysis
#'
#' Comprehensive tools for assessing and reporting heterogeneity including
#' I², τ², prediction intervals, and variance decomposition.
#'
#' @name heterogeneity
NULL

#' Comprehensive Heterogeneity Report
#'
#' Generate a comprehensive report of heterogeneity statistics for NMA.
#' Includes I², τ², between-study variance, prediction intervals, and
#' interpretations.
#'
#' @param nma_result A netmeta object
#' @param reference Reference treatment for reporting
#' @return heterogeneity_report object
#' @export
#' @references
#' Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis.
#' Statistics in Medicine. 2002;21(11):1539-1558.
#'
#' Riley RD, et al. Interpretation of random effects meta-analyses.
#' BMJ. 2011;342:d549.
#'
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' het_report <- heterogeneity_report(nma)
#' print(het_report)
#' }
heterogeneity_report <- function(nma_result, reference = NULL) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  # Extract heterogeneity statistics
  tau <- nma_result$tau
  tau2 <- tau^2
  I2 <- nma_result$I2

  # Q statistic
  Q <- nma_result$Q
  df <- nma_result$df.Q
  p_Q <- pchisq(Q, df, lower.tail = FALSE)

  # Get treatment effects for prediction intervals
  treatments <- rownames(nma_result$TE.random)
  treatments_nonref <- setdiff(treatments, reference)

  # Calculate prediction intervals for each treatment vs reference
  pred_intervals <- calculate_prediction_intervals_all(nma_result, reference)

  # Heterogeneity interpretation
  i2_interpretation <- interpret_i2(I2)
  tau_interpretation <- interpret_tau(tau, nma_result$sm)

  # Summary
  report <- list(
    tau = tau,
    tau2 = tau2,
    I2 = I2,
    I2_percent = I2 * 100,
    Q_statistic = Q,
    Q_df = df,
    Q_pvalue = p_Q,
    prediction_intervals = pred_intervals,
    i2_interpretation = i2_interpretation,
    tau_interpretation = tau_interpretation,
    sm = nma_result$sm,
    reference = reference,
    n_treatments = length(treatments),
    n_studies = length(unique(nma_result$studlab))
  )

  class(report) <- c("heterogeneity_report", "list")
  report
}

#' Calculate Prediction Intervals for All Treatments
#'
#' @param nma_result netmeta object
#' @param reference Reference treatment
#' @return Data frame with prediction intervals
#' @keywords internal
calculate_prediction_intervals_all <- function(nma_result, reference) {

  treatments <- rownames(nma_result$TE.random)
  treatments_nonref <- setdiff(treatments, reference)

  tau <- nma_result$tau

  pred_int_list <- lapply(treatments_nonref, function(trt) {
    effect <- nma_result$TE.random[trt, reference]
    se <- nma_result$seTE.random[trt, reference]
    ci_lower <- nma_result$lower.random[trt, reference]
    ci_upper <- nma_result$upper.random[trt, reference]

    # Prediction interval
    pred_se <- sqrt(se^2 + tau^2)
    pi_lower <- effect - 1.96 * pred_se
    pi_upper <- effect + 1.96 * pred_se

    # Width comparison
    ci_width <- ci_upper - ci_lower
    pi_width <- pi_upper - pi_lower
    width_ratio <- pi_width / ci_width

    data.frame(
      Treatment = trt,
      Effect = effect,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      PI_Lower = pi_lower,
      PI_Upper = pi_upper,
      CI_Width = ci_width,
      PI_Width = pi_width,
      Width_Ratio = width_ratio,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, pred_int_list)
}

#' Interpret I² Statistic
#'
#' @param I2 I-squared value (0 to 1)
#' @return Character interpretation
#' @keywords internal
interpret_i2 <- function(I2) {
  I2_percent <- I2 * 100

  if (I2_percent < 25) {
    "Low heterogeneity (I² < 25%)"
  } else if (I2_percent < 50) {
    "Moderate heterogeneity (25% ≤ I² < 50%)"
  } else if (I2_percent < 75) {
    "Substantial heterogeneity (50% ≤ I² < 75%)"
  } else {
    "Considerable heterogeneity (I² ≥ 75%)"
  }
}

#' Interpret Tau Statistic
#'
#' @param tau Tau value (between-study SD)
#' @param sm Summary measure type
#' @return Character interpretation
#' @keywords internal
interpret_tau <- function(tau, sm) {

  # Context-specific interpretation
  if (sm %in% c("OR", "RR", "HR")) {
    # Log scale
    if (tau < 0.1) {
      "Minimal between-study heterogeneity (τ < 0.1)"
    } else if (tau < 0.3) {
      "Low to moderate heterogeneity (0.1 ≤ τ < 0.3)"
    } else if (tau < 0.5) {
      "Moderate to substantial heterogeneity (0.3 ≤ τ < 0.5)"
    } else {
      "Substantial heterogeneity (τ ≥ 0.5)"
    }
  } else {
    # Mean difference scale
    sprintf("Between-study standard deviation: τ = %.3f", tau)
  }
}

#' Calculate I² for Specific Comparison
#'
#' Calculate I² for a direct pairwise comparison within the network
#'
#' @param data Pairwise data
#' @param treat1 First treatment
#' @param treat2 Second treatment
#' @return List with I², Q, and other statistics
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' i2_result <- calculate_comparison_i2(data, "A", "B")
#' print(i2_result)
#' }
calculate_comparison_i2 <- function(data, treat1, treat2) {

  # Extract data for this comparison
  comp_data <- data[
    (data$treat1 == treat1 & data$treat2 == treat2) |
    (data$treat1 == treat2 & data$treat2 == treat1),
  ]

  if (nrow(comp_data) < 2) {
    return(list(
      comparison = paste(treat1, "vs", treat2),
      n_studies = nrow(comp_data),
      I2 = NA,
      message = "Insufficient studies for I² calculation (need ≥2)"
    ))
  }

  # Run meta-analysis
  ma <- meta::metagen(
    TE = comp_data$TE,
    seTE = comp_data$seTE,
    studlab = comp_data$studlab
  )

  list(
    comparison = paste(treat1, "vs", treat2),
    n_studies = nrow(comp_data),
    I2 = ma$I2,
    I2_percent = ma$I2 * 100,
    tau = ma$tau,
    tau2 = ma$tau^2,
    Q = ma$Q,
    Q_df = ma$df.Q,
    Q_pvalue = ma$pval.Q,
    interpretation = interpret_i2(ma$I2)
  )
}

#' Calculate I² for All Direct Comparisons
#'
#' Calculate heterogeneity statistics for all direct comparisons in the network
#'
#' @param data Pairwise data
#' @return Data frame with I² for each comparison
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' all_i2 <- calculate_all_comparison_i2(data)
#' print(all_i2)
#' }
calculate_all_comparison_i2 <- function(data) {

  # Identify unique comparisons
  comparisons <- unique(data.frame(
    t1 = pmin(as.character(data$treat1), as.character(data$treat2)),
    t2 = pmax(as.character(data$treat1), as.character(data$treat2)),
    stringsAsFactors = FALSE
  ))

  results <- lapply(1:nrow(comparisons), function(i) {
    result <- calculate_comparison_i2(data, comparisons$t1[i], comparisons$t2[i])
    as.data.frame(result, stringsAsFactors = FALSE)
  })

  result_df <- do.call(rbind, results)
  result_df <- result_df[order(-result_df$I2_percent), ]

  class(result_df) <- c("comparison_i2", "data.frame")
  result_df
}

#' Variance Decomposition
#'
#' Decompose total variance into within-study and between-study components
#'
#' @param nma_result A netmeta object
#' @return List with variance components
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' variance_decomp <- variance_decomposition(nma)
#' print(variance_decomp)
#' }
variance_decomposition <- function(nma_result) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Between-study variance
  tau2 <- nma_result$tau^2

  # Within-study variance (average of sampling variances)
  within_var <- mean(nma_result$seTE^2, na.rm = TRUE)

  # Total variance
  total_var <- tau2 + within_var

  # Proportions
  prop_between <- tau2 / total_var
  prop_within <- within_var / total_var

  list(
    between_study_variance = tau2,
    within_study_variance = within_var,
    total_variance = total_var,
    proportion_between = prop_between,
    proportion_within = prop_within,
    I2 = nma_result$I2,
    interpretation = sprintf(
      "%.1f%% of total variance is due to between-study heterogeneity",
      prop_between * 100
    )
  )
}

#' Prediction Interval Calculation
#'
#' Calculate prediction interval for a specific treatment comparison.
#' Prediction intervals account for heterogeneity and predict where future
#' study effects might fall.
#'
#' @param nma_result A netmeta object
#' @param treatment Treatment to compare to reference
#' @param reference Reference treatment
#' @param level Confidence level (default: 0.95)
#' @return List with effect, CI, and PI
#' @export
#' @references
#' IntHout J, et al. Plea for routinely presenting prediction intervals in
#' meta-analysis. BMJ Open. 2016;6(7):e010247.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' pi <- prediction_interval(nma, "DrugA", "Placebo")
#' print(pi)
#' }
prediction_interval <- function(nma_result, treatment, reference = NULL, level = 0.95) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  # Get effect estimate and SE
  effect <- nma_result$TE.random[treatment, reference]
  se <- nma_result$seTE.random[treatment, reference]
  tau <- nma_result$tau

  # Confidence interval (from nma_result)
  ci_lower <- nma_result$lower.random[treatment, reference]
  ci_upper <- nma_result$upper.random[treatment, reference]

  # Prediction interval
  z <- qnorm((1 + level) / 2)
  pred_se <- sqrt(se^2 + tau^2)
  pi_lower <- effect - z * pred_se
  pi_upper <- effect + z * pred_se

  result <- list(
    treatment = treatment,
    reference = reference,
    effect = effect,
    se = se,
    tau = tau,
    level = level,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    pi_lower = pi_lower,
    pi_upper = pi_upper,
    ci_width = ci_upper - ci_lower,
    pi_width = pi_upper - pi_lower,
    interpretation = sprintf(
      "We are %.0f%% confident the true effect is between %.3f and %.3f. ",
      level * 100, ci_lower, ci_upper
    ) + sprintf(
      "We predict that %.0f%% of future similar studies will have effects between %.3f and %.3f.",
      level * 100, pi_lower, pi_upper
    )
  )

  class(result) <- c("prediction_interval", "list")
  result
}

#' Plot Heterogeneity Statistics
#'
#' Create visualization of heterogeneity across comparisons
#'
#' @param comparison_i2_result Result from calculate_all_comparison_i2()
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' all_i2 <- calculate_all_comparison_i2(data)
#' plot_heterogeneity(all_i2)
#' }
plot_heterogeneity <- function(comparison_i2_result) {

  if (!inherits(comparison_i2_result, "comparison_i2")) {
    stop("Input must be from calculate_all_comparison_i2()")
  }

  # Remove NA values
  plot_data <- comparison_i2_result[!is.na(comparison_i2_result$I2_percent), ]

  if (nrow(plot_data) == 0) {
    stop("No valid I² values to plot")
  }

  # Order by I²
  plot_data <- plot_data[order(plot_data$I2_percent), ]
  plot_data$comparison <- factor(plot_data$comparison, levels = plot_data$comparison)

  # Color by interpretation
  plot_data$color_group <- cut(
    plot_data$I2_percent,
    breaks = c(0, 25, 50, 75, 100),
    labels = c("Low", "Moderate", "Substantial", "Considerable"),
    include.lowest = TRUE
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = comparison, y = I2_percent, fill = color_group)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = c(25, 50, 75), linetype = "dashed", color = "gray40") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("Low" = "#1a9850", "Moderate" = "#91cf60",
                "Substantial" = "#fee08b", "Considerable" = "#d73027"),
      name = "Heterogeneity Level"
    ) +
    ggplot2::labs(
      title = "Heterogeneity (I²) Across Direct Comparisons",
      subtitle = "Higher I² indicates greater between-study heterogeneity",
      x = "Comparison",
      y = "I² (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank()
    )
}

#' Print Heterogeneity Report
#'
#' @param x heterogeneity_report object
#' @param ... Additional arguments
#' @export
print.heterogeneity_report <- function(x, ...) {
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("  Comprehensive Heterogeneity Assessment\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Network: %d treatments, %d studies\n", x$n_treatments, x$n_studies))
  cat(sprintf("Reference: %s\n", x$reference))
  cat(sprintf("Effect measure: %s\n\n", x$sm))

  cat("══ Global Heterogeneity Statistics ══\n\n")

  cat(sprintf("τ (between-study SD):    %.3f\n", x$tau))
  cat(sprintf("τ² (between-study var):  %.3f\n", x$tau2))
  cat(sprintf("I² (heterogeneity):      %.1f%%\n", x$I2_percent))
  cat(sprintf("  Interpretation: %s\n\n", x$i2_interpretation))

  cat(sprintf("Q statistic:             %.2f (df = %d, p = %.4f)\n",
             x$Q_statistic, x$Q_df, x$Q_pvalue))
  if (x$Q_pvalue < 0.05) {
    cat("  ⚠ Significant heterogeneity detected (p < 0.05)\n")
  } else {
    cat("  ✓ No significant heterogeneity detected (p ≥ 0.05)\n")
  }

  cat("\n")
  cat(x$tau_interpretation)
  cat("\n\n")

  cat("══ Prediction Intervals ══\n\n")
  cat("Prediction intervals account for heterogeneity and predict where\n")
  cat("future study effects might fall (wider than confidence intervals).\n\n")

  if (!is.null(x$prediction_intervals) && nrow(x$prediction_intervals) > 0) {
    pi_summary <- x$prediction_intervals[, c("Treatment", "Effect", "CI_Lower", "CI_Upper",
                                             "PI_Lower", "PI_Upper", "Width_Ratio")]
    names(pi_summary) <- c("Treatment", "Effect", "CI_Low", "CI_High",
                          "PI_Low", "PI_High", "PI/CI Ratio")

    print(pi_summary, row.names = FALSE, digits = 3)

    cat("\nNote: PI/CI Ratio shows how much wider prediction intervals are\n")
    cat("      due to heterogeneity (higher = more heterogeneity impact)\n")
  }

  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")

  invisible(x)
}

#' Print Prediction Interval
#'
#' @param x prediction_interval object
#' @param ... Additional arguments
#' @export
print.prediction_interval <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Prediction Interval\n")
  cat("══════════════════════════════════════════════════\n\n")

  cat(sprintf("Comparison: %s vs %s\n", x$treatment, x$reference))
  cat(sprintf("Effect estimate: %.3f (SE = %.3f)\n", x$effect, x$se))
  cat(sprintf("Between-study SD (τ): %.3f\n\n", x$tau))

  cat(sprintf("%.0f%% Confidence Interval:  [%.3f, %.3f]  (width: %.3f)\n",
             x$level * 100, x$ci_lower, x$ci_upper, x$ci_width))
  cat(sprintf("%.0f%% Prediction Interval:  [%.3f, %.3f]  (width: %.3f)\n\n",
             x$level * 100, x$pi_lower, x$pi_upper, x$pi_width))

  cat("Interpretation:\n")
  cat(strwrap(x$interpretation, width = 60, prefix = "  "), sep = "\n")
  cat("\n")

  invisible(x)
}
