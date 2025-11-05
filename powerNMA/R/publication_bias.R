#' Publication Bias Assessment for Network Meta-Analysis
#'
#' Tools for detecting and assessing publication bias including Egger's test,
#' comparison-adjusted funnel plots, and small-study effects testing.
#'
#' @name publication_bias
NULL

#' Comprehensive Publication Bias Assessment
#'
#' Perform multiple tests for publication bias and small-study effects in NMA.
#' Includes comparison-adjusted Egger's test, funnel plot asymmetry test,
#' and visual diagnostics.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @param method Test method: "Egger", "Begg", "both"
#' @return publication_bias object with test results
#' @export
#' @references
#' Chaimani A, Salanti G. Using network meta-analysis to evaluate the existence
#' of small-study effects in a network of interventions.
#' Research Synthesis Methods. 2012;3(2):161-176.
#'
#' Egger M, et al. Bias in meta-analysis detected by a simple, graphical test.
#' BMJ. 1997;315(7109):629-634.
#'
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' pub_bias <- assess_publication_bias(nma, data)
#' print(pub_bias)
#' plot(pub_bias)
#' }
assess_publication_bias <- function(nma_result,
                                   data,
                                   method = c("Egger", "Begg", "both")) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  method <- match.arg(method)

  results <- list()

  # Comparison-adjusted Egger's test
  if (method %in% c("Egger", "both")) {
    message("Running comparison-adjusted Egger's test...")
    results$egger <- comparison_adjusted_egger(nma_result, data)
  }

  # Comparison-adjusted Begg's test (rank correlation)
  if (method %in% c("Begg", "both")) {
    message("Running comparison-adjusted Begg's test...")
    results$begg <- comparison_adjusted_begg(nma_result, data)
  }

  # Design-by-comparison interaction test
  message("Testing for design-by-comparison interaction...")
  results$design_comparison <- test_design_comparison_interaction(nma_result, data)

  # Summary
  results$summary <- summarize_bias_tests(results)

  class(results) <- c("publication_bias", "list")
  results
}

#' Comparison-Adjusted Egger's Test
#'
#' Performs Egger's regression test adjusted for different comparisons.
#' Tests whether small studies show different effects than large studies.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @return List with Egger test results
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' egger_result <- comparison_adjusted_egger(nma, data)
#' print(egger_result)
#' }
comparison_adjusted_egger <- function(nma_result, data) {

  # Calculate comparison-adjusted effects
  studlab <- nma_result$studlab
  treat1 <- nma_result$treat1
  treat2 <- nma_result$treat2
  TE <- nma_result$TE
  seTE <- nma_result$seTE

  # Get comparison-specific effects from NMA (random effects)
  comparison_effects <- nma_result$TE.random[cbind(treat2, treat1)]

  # Adjusted effects (observed - expected)
  adjusted_effects <- TE - comparison_effects

  # Precision (1/SE) and inverse variance weights
  precision <- 1 / seTE
  weights <- 1 / seTE^2

  # Egger's regression: adjusted_effect ~ precision
  # Using weighted least squares
  reg_data <- data.frame(
    adj_effect = adjusted_effects,
    precision = precision,
    weights = weights
  )

  # Remove any NA values
  reg_data <- reg_data[complete.cases(reg_data), ]

  if (nrow(reg_data) < 5) {
    warning("Insufficient data for Egger's test (need ≥5 studies)")
    return(list(
      test = "Egger",
      n_studies = nrow(reg_data),
      insufficient_data = TRUE
    ))
  }

  # Fit weighted regression
  fit <- lm(adj_effect ~ precision, data = reg_data, weights = weights)

  # Extract results
  coef_summary <- summary(fit)$coefficients
  intercept <- coef_summary["(Intercept)", "Estimate"]
  intercept_se <- coef_summary["(Intercept)", "Std. Error"]
  intercept_t <- coef_summary["(Intercept)", "t value"]
  intercept_p <- coef_summary["(Intercept)", "Pr(>|t|)"]

  slope <- coef_summary["precision", "Estimate"]
  slope_p <- coef_summary["precision", "Pr(>|t|)"]

  list(
    test = "Comparison-adjusted Egger",
    n_studies = nrow(reg_data),
    intercept = intercept,
    intercept_se = intercept_se,
    intercept_t = intercept_t,
    intercept_p = intercept_p,
    slope = slope,
    slope_p = slope_p,
    interpretation = interpret_egger(intercept_p),
    model = fit,
    data = reg_data
  )
}

#' Comparison-Adjusted Begg's Test
#'
#' Rank correlation test for funnel plot asymmetry, adjusted for comparisons.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @return List with Begg test results
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' begg_result <- comparison_adjusted_begg(nma, data)
#' }
comparison_adjusted_begg <- function(nma_result, data) {

  # Calculate comparison-adjusted effects
  treat1 <- nma_result$treat1
  treat2 <- nma_result$treat2
  TE <- nma_result$TE
  seTE <- nma_result$seTE

  comparison_effects <- nma_result$TE.random[cbind(treat2, treat1)]
  adjusted_effects <- TE - comparison_effects

  # Remove NA
  valid_idx <- !is.na(adjusted_effects) & !is.na(seTE)
  adj_eff <- adjusted_effects[valid_idx]
  se <- seTE[valid_idx]

  if (length(adj_eff) < 5) {
    warning("Insufficient data for Begg's test (need ≥5 studies)")
    return(list(
      test = "Begg",
      n_studies = length(adj_eff),
      insufficient_data = TRUE
    ))
  }

  # Kendall's tau correlation between adjusted effects and their variances
  cor_test <- cor.test(adj_eff, se^2, method = "kendall")

  list(
    test = "Comparison-adjusted Begg",
    n_studies = length(adj_eff),
    tau = cor_test$estimate,
    z_score = cor_test$statistic,
    p_value = cor_test$p.value,
    interpretation = interpret_begg(cor_test$p.value)
  )
}

#' Test Design-by-Comparison Interaction
#'
#' Test whether small-study effects differ across different types of comparisons.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @return List with test results
#' @keywords internal
test_design_comparison_interaction <- function(nma_result, data) {

  # Create comparison identifiers
  comparison_id <- paste(
    pmin(as.character(nma_result$treat1), as.character(nma_result$treat2)),
    pmax(as.character(nma_result$treat1), as.character(nma_result$treat2)),
    sep = " vs "
  )

  # Only test if we have multiple comparisons with multiple studies each
  comp_counts <- table(comparison_id)
  if (length(comp_counts) < 2 || max(comp_counts) < 2) {
    return(list(
      test = "Design-by-comparison interaction",
      insufficient_data = TRUE
    ))
  }

  # Get adjusted effects and precision
  treat1 <- nma_result$treat1
  treat2 <- nma_result$treat2
  TE <- nma_result$TE
  seTE <- nma_result$seTE

  comparison_effects <- nma_result$TE.random[cbind(treat2, treat1)]
  adjusted_effects <- TE - comparison_effects
  precision <- 1 / seTE

  # Test for interaction
  test_data <- data.frame(
    adj_effect = adjusted_effects,
    precision = precision,
    comparison = comparison_id,
    weights = 1 / seTE^2
  )

  test_data <- test_data[complete.cases(test_data), ]

  if (nrow(test_data) < 10) {
    return(list(
      test = "Design-by-comparison interaction",
      insufficient_data = TRUE
    ))
  }

  # Fit model with interaction
  fit_interaction <- lm(adj_effect ~ precision * comparison,
                       data = test_data, weights = weights)

  # ANOVA to test interaction
  anova_result <- anova(fit_interaction)

  # Extract p-value for interaction term
  interaction_p <- anova_result["precision:comparison", "Pr(>F)"]

  list(
    test = "Design-by-comparison interaction",
    n_studies = nrow(test_data),
    n_comparisons = length(unique(test_data$comparison)),
    interaction_p = interaction_p,
    interpretation = if (!is.na(interaction_p) && interaction_p < 0.10) {
      "Significant interaction: Small-study effects differ across comparisons"
    } else {
      "No significant interaction detected"
    }
  )
}

#' Interpret Egger Test Result
#'
#' @param p_value P-value from Egger test
#' @return Character interpretation
#' @keywords internal
interpret_egger <- function(p_value) {
  if (p_value < 0.05) {
    "Significant asymmetry detected (p < 0.05) - potential small-study effects or publication bias"
  } else if (p_value < 0.10) {
    "Weak evidence of asymmetry (0.05 ≤ p < 0.10) - borderline concern"
  } else {
    "No significant asymmetry detected (p ≥ 0.10)"
  }
}

#' Interpret Begg Test Result
#'
#' @param p_value P-value from Begg test
#' @return Character interpretation
#' @keywords internal
interpret_begg <- function(p_value) {
  if (p_value < 0.05) {
    "Significant rank correlation detected (p < 0.05) - potential small-study effects"
  } else if (p_value < 0.10) {
    "Weak evidence of rank correlation (0.05 ≤ p < 0.10)"
  } else {
    "No significant rank correlation detected (p ≥ 0.10)"
  }
}

#' Summarize Bias Tests
#'
#' @param results List of test results
#' @return Summary data frame
#' @keywords internal
summarize_bias_tests <- function(results) {

  summary_list <- list()

  if (!is.null(results$egger)) {
    if (!isTRUE(results$egger$insufficient_data)) {
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Test = "Egger's test",
        Statistic = sprintf("t = %.3f", results$egger$intercept_t),
        P_value = results$egger$intercept_p,
        Significant = results$egger$intercept_p < 0.10,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.null(results$begg)) {
    if (!isTRUE(results$begg$insufficient_data)) {
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Test = "Begg's test",
        Statistic = sprintf("τ = %.3f", results$begg$tau),
        P_value = results$begg$p_value,
        Significant = results$begg$p_value < 0.10,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.null(results$design_comparison)) {
    if (!isTRUE(results$design_comparison$insufficient_data)) {
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Test = "Design×Comparison",
        Statistic = "F-test",
        P_value = results$design_comparison$interaction_p,
        Significant = results$design_comparison$interaction_p < 0.10,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(summary_list) > 0) {
    do.call(rbind, summary_list)
  } else {
    data.frame(Test = character(0), Statistic = character(0),
              P_value = numeric(0), Significant = logical(0))
  }
}

#' Contour-Enhanced Funnel Plot
#'
#' Create a funnel plot with significance contours to distinguish asymmetry
#' due to publication bias from asymmetry due to heterogeneity.
#'
#' @param nma_result A netmeta object
#' @param contours Significance levels for contours (default: c(0.01, 0.05, 0.10))
#' @return ggplot2 object
#' @export
#' @references
#' Peters JL, et al. Contour-enhanced meta-analysis funnel plots help
#' distinguish publication bias from other causes of asymmetry.
#' Journal of Clinical Epidemiology. 2008;61(10):991-996.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' contour_funnel_plot(nma)
#' }
contour_funnel_plot <- function(nma_result, contours = c(0.01, 0.05, 0.10)) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Get comparison-adjusted effects
  treat1 <- nma_result$treat1
  treat2 <- nma_result$treat2
  TE <- nma_result$TE
  seTE <- nma_result$seTE

  comparison_effects <- nma_result$TE.random[cbind(treat2, treat1)]
  adjusted_effects <- TE - comparison_effects

  plot_data <- data.frame(
    Effect = adjusted_effects,
    SE = seTE,
    Precision = 1/seTE
  )

  plot_data <- plot_data[complete.cases(plot_data), ]

  # Create contour lines
  se_range <- seq(0, max(plot_data$SE, na.rm = TRUE), length.out = 100)

  contour_data <- lapply(contours, function(alpha) {
    z <- qnorm(1 - alpha/2)
    data.frame(
      SE = se_range,
      Lower = -z * se_range,
      Upper = z * se_range,
      Alpha = sprintf("p = %.2f", alpha)
    )
  })

  contour_df <- do.call(rbind, contour_data)

  # Create plot
  p <- ggplot2::ggplot() +
    # Contour regions
    ggplot2::geom_ribbon(
      data = contour_df,
      ggplot2::aes(x = SE, ymin = Lower, ymax = Upper, fill = Alpha),
      alpha = 0.2
    ) +
    # Zero line
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    # Study points
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = SE, y = Effect),
      size = 2, alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = c("#d73027", "#fee08b", "#91cf60"),
      name = "Significance"
    ) +
    ggplot2::scale_x_reverse() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Contour-Enhanced Funnel Plot",
      subtitle = "Comparison-adjusted effects with significance contours",
      x = "Standard Error (SE)",
      y = "Comparison-Adjusted Effect"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Generate Publication Bias Report
#'
#' Create a comprehensive report combining multiple bias assessments
#'
#' @param pub_bias A publication_bias object
#' @param nma_result Original netmeta object
#' @param output_file Optional file path to save report
#' @return Character vector with report text
#' @export
#' @examples
#' \dontrun{
#' pub_bias <- assess_publication_bias(nma, data)
#' report <- generate_bias_report(pub_bias, nma)
#' cat(report, sep = "\n")
#' }
generate_bias_report <- function(pub_bias, nma_result, output_file = NULL) {

  if (!inherits(pub_bias, "publication_bias")) {
    stop("pub_bias must be a publication_bias object")
  }

  report <- c(
    "",
    "═══════════════════════════════════════════════════════════",
    "  Publication Bias Assessment Report",
    "═══════════════════════════════════════════════════════════",
    "",
    sprintf("Network: %d treatments, %d studies",
           length(unique(c(nma_result$treat1, nma_result$treat2))),
           length(unique(nma_result$studlab))),
    ""
  )

  # Overall summary
  if (!is.null(pub_bias$summary) && nrow(pub_bias$summary) > 0) {
    report <- c(report, "Test Summary:", "")
    summary_text <- capture.output(print(pub_bias$summary, row.names = FALSE))
    report <- c(report, summary_text, "")

    n_sig <- sum(pub_bias$summary$Significant, na.rm = TRUE)
    if (n_sig > 0) {
      report <- c(report,
                 sprintf("⚠ WARNING: %d of %d tests suggest potential bias", n_sig, nrow(pub_bias$summary)),
                 "")
    } else {
      report <- c(report,
                 "✓ No tests indicate significant publication bias",
                 "")
    }
  }

  # Egger details
  if (!is.null(pub_bias$egger) && !isTRUE(pub_bias$egger$insufficient_data)) {
    report <- c(report,
               "─────────────────────────────────────────────────────────",
               "Egger's Test Details:",
               "",
               sprintf("  Intercept: %.4f (SE = %.4f)", pub_bias$egger$intercept, pub_bias$egger$intercept_se),
               sprintf("  t-statistic: %.3f", pub_bias$egger$intercept_t),
               sprintf("  P-value: %.4f", pub_bias$egger$intercept_p),
               sprintf("  Interpretation: %s", pub_bias$egger$interpretation),
               "")
  }

  # Begg details
  if (!is.null(pub_bias$begg) && !isTRUE(pub_bias$begg$insufficient_data)) {
    report <- c(report,
               "─────────────────────────────────────────────────────────",
               "Begg's Test Details:",
               "",
               sprintf("  Kendall's τ: %.4f", pub_bias$begg$tau),
               sprintf("  P-value: %.4f", pub_bias$begg$p_value),
               sprintf("  Interpretation: %s", pub_bias$begg$interpretation),
               "")
  }

  # Recommendations
  report <- c(report,
             "─────────────────────────────────────────────────────────",
             "Recommendations:",
             "")

  if (!is.null(pub_bias$summary) && any(pub_bias$summary$Significant, na.rm = TRUE)) {
    report <- c(report,
               "• Evidence of potential publication bias detected",
               "• Consider sensitivity analyses excluding small studies",
               "• Assess risk of bias in individual studies",
               "• Search for unpublished studies and grey literature",
               "• Interpret network meta-analysis results with caution",
               "")
  } else {
    report <- c(report,
               "• No strong evidence of publication bias",
               "• Results appear robust to small-study effects",
               "• Continue with standard NMA interpretation",
               "")
  }

  report <- c(report,
             "═══════════════════════════════════════════════════════════",
             "")

  # Save to file if requested
  if (!is.null(output_file)) {
    writeLines(report, output_file)
    message(sprintf("Report saved to: %s", output_file))
  }

  report
}

#' Print Publication Bias Results
#'
#' @param x publication_bias object
#' @param ... Additional arguments
#' @export
print.publication_bias <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Publication Bias Assessment\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  if (!is.null(x$summary) && nrow(x$summary) > 0) {
    cat("Test Summary:\n\n")
    print(x$summary, row.names = FALSE)
    cat("\n")

    n_sig <- sum(x$summary$Significant, na.rm = TRUE)
    if (n_sig > 0) {
      cat(sprintf("⚠ %d of %d tests suggest potential publication bias\n\n", n_sig, nrow(x$summary)))
    } else {
      cat("✓ No significant publication bias detected\n\n")
    }
  }

  # Egger
  if (!is.null(x$egger) && !isTRUE(x$egger$insufficient_data)) {
    cat("Egger's Test:\n")
    cat(sprintf("  P-value: %.4f\n", x$egger$intercept_p))
    cat(sprintf("  %s\n\n", x$egger$interpretation))
  }

  # Begg
  if (!is.null(x$begg) && !isTRUE(x$begg$insufficient_data)) {
    cat("Begg's Test:\n")
    cat(sprintf("  P-value: %.4f\n", x$begg$p_value))
    cat(sprintf("  %s\n\n", x$begg$interpretation))
  }

  cat("Use plot() to visualize funnel plot asymmetry\n")
  cat("Use generate_bias_report() for detailed report\n\n")

  invisible(x)
}

#' Plot Publication Bias Object
#'
#' @param x publication_bias object
#' @param type Plot type: "funnel" or "contour"
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.publication_bias <- function(x, type = c("funnel", "contour"), ...) {
  type <- match.arg(type)

  # Extract nma_result from Egger model data if available
  if (type == "funnel" && !is.null(x$egger) && !is.null(x$egger$data)) {
    plot_data <- x$egger$data

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = adj_effect, y = precision)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = "Comparison-Adjusted Funnel Plot",
        x = "Comparison-Adjusted Effect",
        y = "Precision (1/SE)"
      ) +
      ggplot2::theme_minimal()

    return(p)
  }

  message("Use contour_funnel_plot(nma_result) for contour-enhanced plot")
  invisible(NULL)
}
