#' CINeMA: Confidence in Network Meta-Analysis
#'
#' Implementation of the CINeMA framework for assessing confidence in NMA results.
#' Evaluates 6 domains: within-study bias, reporting bias, indirectness,
#' imprecision, heterogeneity, and incoherence.
#'
#' @name cinema
#' @references
#' Nikolakopoulou A, et al. (2020). CINeMA: An approach for assessing confidence
#' in the results of a network meta-analysis. PLOS Medicine, 17(4):e1003082.
NULL

#' Assess Confidence using CINeMA Framework
#'
#' Comprehensive confidence assessment for NMA results following CINeMA methodology.
#' Evaluates 6 domains and provides an overall confidence rating.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @param risk_of_bias Data frame with risk of bias assessments (optional)
#' @param inconsistency_result Results from inconsistency assessment (optional)
#' @param heterogeneity_result Results from heterogeneity assessment (optional)
#' @param comparison Specific comparison to assess (NULL for all)
#' @return cinema_assessment object
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' cinema <- assess_cinema_confidence(nma, data)
#' print(cinema)
#' }
assess_cinema_confidence <- function(nma_result,
                                    data,
                                    risk_of_bias = NULL,
                                    inconsistency_result = NULL,
                                    heterogeneity_result = NULL,
                                    comparison = NULL) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Get all comparisons or specific one
  treatments <- rownames(nma_result$TE.random)
  reference <- nma_result$reference.group

  if (is.null(comparison)) {
    comparisons <- setdiff(treatments, reference)
  } else {
    comparisons <- comparison
  }

  # Assess each comparison
  assessments <- lapply(comparisons, function(comp) {
    assess_comparison_cinema(
      nma_result = nma_result,
      data = data,
      comparison = comp,
      reference = reference,
      risk_of_bias = risk_of_bias,
      inconsistency_result = inconsistency_result,
      heterogeneity_result = heterogeneity_result
    )
  })

  names(assessments) <- comparisons

  # Create summary table
  summary_table <- do.call(rbind, lapply(names(assessments), function(comp) {
    a <- assessments[[comp]]
    data.frame(
      Comparison = paste(comp, "vs", reference),
      WithinStudyBias = a$domain_ratings$within_study_bias,
      ReportingBias = a$domain_ratings$reporting_bias,
      Indirectness = a$domain_ratings$indirectness,
      Imprecision = a$domain_ratings$imprecision,
      Heterogeneity = a$domain_ratings$heterogeneity,
      Incoherence = a$domain_ratings$incoherence,
      OverallConfidence = a$overall_confidence,
      stringsAsFactors = FALSE
    )
  }))

  result <- list(
    assessments = assessments,
    summary = summary_table,
    n_comparisons = length(comparisons),
    reference = reference
  )

  class(result) <- c("cinema_assessment", "list")
  result
}

#' Assess Single Comparison using CINeMA
#'
#' @param nma_result netmeta object
#' @param data Pairwise data
#' @param comparison Treatment to assess
#' @param reference Reference treatment
#' @param risk_of_bias ROB data
#' @param inconsistency_result Inconsistency results
#' @param heterogeneity_result Heterogeneity results
#' @return List with domain assessments
#' @keywords internal
assess_comparison_cinema <- function(nma_result, data, comparison, reference,
                                     risk_of_bias, inconsistency_result,
                                     heterogeneity_result) {

  # Domain 1: Within-study bias
  within_bias <- assess_within_study_bias(data, comparison, reference, risk_of_bias)

  # Domain 2: Reporting bias (publication bias)
  reporting_bias <- assess_reporting_bias(nma_result, data)

  # Domain 3: Indirectness
  indirectness <- assess_indirectness(nma_result, data, comparison, reference)

  # Domain 4: Imprecision
  imprecision <- assess_imprecision(nma_result, comparison, reference)

  # Domain 5: Heterogeneity
  heterogeneity <- assess_heterogeneity_cinema(nma_result, heterogeneity_result)

  # Domain 6: Incoherence (inconsistency)
  incoherence <- assess_incoherence(nma_result, inconsistency_result, comparison, reference)

  # Overall confidence rating
  overall <- calculate_overall_confidence(
    within_bias, reporting_bias, indirectness,
    imprecision, heterogeneity, incoherence
  )

  list(
    comparison = paste(comparison, "vs", reference),
    domain_ratings = list(
      within_study_bias = within_bias$rating,
      reporting_bias = reporting_bias$rating,
      indirectness = indirectness$rating,
      imprecision = imprecision$rating,
      heterogeneity = heterogeneity$rating,
      incoherence = incoherence$rating
    ),
    domain_details = list(
      within_study_bias = within_bias,
      reporting_bias = reporting_bias,
      indirectness = indirectness,
      imprecision = imprecision,
      heterogeneity = heterogeneity,
      incoherence = incoherence
    ),
    overall_confidence = overall$rating,
    overall_justification = overall$justification
  )
}

#' Domain 1: Within-Study Bias Assessment
#'
#' @keywords internal
assess_within_study_bias <- function(data, comparison, reference, risk_of_bias) {

  if (is.null(risk_of_bias)) {
    return(list(
      rating = "No concerns (not assessed)",
      justification = "Risk of bias data not provided. Assume studies are well-conducted.",
      studies_high_rob = NA
    ))
  }

  # Match studies contributing to this comparison
  # This is simplified - in practice, need contribution matrix
  contributing_studies <- unique(data$studlab)

  if (!"studlab" %in% names(risk_of_bias) || !"overall_rob" %in% names(risk_of_bias)) {
    return(list(
      rating = "No concerns (format issue)",
      justification = "Risk of bias data format not recognized.",
      studies_high_rob = NA
    ))
  }

  # Count studies with high ROB
  rob_studies <- risk_of_bias[risk_of_bias$studlab %in% contributing_studies, ]
  if (nrow(rob_studies) == 0) {
    return(list(
      rating = "No concerns (no data)",
      justification = "No risk of bias data for contributing studies.",
      studies_high_rob = NA
    ))
  }

  high_rob <- sum(rob_studies$overall_rob == "High", na.rm = TRUE)
  total_studies <- nrow(rob_studies)
  pct_high <- high_rob / total_studies

  # Rating based on percentage with high ROB
  if (pct_high >= 0.50) {
    rating <- "Major concerns"
    justification <- sprintf("%d/%d (%.0f%%) contributing studies at high risk of bias",
                            high_rob, total_studies, pct_high * 100)
  } else if (pct_high >= 0.25) {
    rating <- "Some concerns"
    justification <- sprintf("%d/%d (%.0f%%) contributing studies at high risk of bias",
                            high_rob, total_studies, pct_high * 100)
  } else {
    rating <- "No concerns"
    justification <- sprintf("Only %d/%d (%.0f%%) studies at high risk of bias",
                            high_rob, total_studies, pct_high * 100)
  }

  list(
    rating = rating,
    justification = justification,
    studies_high_rob = high_rob,
    total_studies = total_studies,
    percent_high_rob = pct_high
  )
}

#' Domain 2: Reporting Bias Assessment
#'
#' @keywords internal
assess_reporting_bias <- function(nma_result, data) {

  n_studies <- length(unique(data$studlab))

  # Simple heuristic - needs actual publication bias tests
  if (n_studies < 10) {
    rating <- "Some concerns"
    justification <- sprintf("Small number of studies (n=%d); publication bias difficult to assess", n_studies)
  } else {
    # In practice, use results from publication bias tests
    rating <- "No concerns (not formally tested)"
    justification <- sprintf("Sufficient studies (n=%d); formal publication bias testing recommended", n_studies)
  }

  list(
    rating = rating,
    justification = justification,
    n_studies = n_studies
  )
}

#' Domain 3: Indirectness Assessment
#'
#' @keywords internal
assess_indirectness <- function(nma_result, data, comparison, reference) {

  # Check if direct evidence exists
  direct_evidence <- data[
    (data$treat1 == comparison & data$treat2 == reference) |
    (data$treat1 == reference & data$treat2 == comparison),
  ]

  has_direct <- nrow(direct_evidence) > 0

  if (has_direct) {
    rating <- "No concerns"
    justification <- sprintf("Direct evidence available (%d studies)", nrow(direct_evidence))
  } else {
    rating <- "Some concerns"
    justification <- "No direct evidence; effect estimated indirectly through network"
  }

  list(
    rating = rating,
    justification = justification,
    has_direct_evidence = has_direct,
    n_direct_studies = nrow(direct_evidence)
  )
}

#' Domain 4: Imprecision Assessment
#'
#' @keywords internal
assess_imprecision <- function(nma_result, comparison, reference) {

  # Get confidence interval
  effect <- nma_result$TE.random[comparison, reference]
  lower <- nma_result$lower.random[comparison, reference]
  upper <- nma_result$upper.random[comparison, reference]
  se <- nma_result$seTE.random[comparison, reference]

  ci_width <- upper - lower

  # Check if CI crosses null or includes clinically unimportant effects
  # For effect measures on log scale (OR, RR, HR)
  sm <- nma_result$sm

  if (sm %in% c("OR", "RR", "HR")) {
    # On log scale, null is 0
    crosses_null <- (lower < 0 & upper > 0)

    # Wide CI if > 2 on original scale (log(2) ≈ 0.693)
    is_wide <- ci_width > 1.386  # 2 * log(2)
  } else {
    # Mean difference scale
    crosses_null <- (lower < 0 & upper > 0)
    is_wide <- ci_width > (2 * se)  # Arbitrary threshold
  }

  if (crosses_null && is_wide) {
    rating <- "Major concerns"
    justification <- "Wide confidence interval crossing null effect"
  } else if (crosses_null || is_wide) {
    rating <- "Some concerns"
    justification <- if (crosses_null) {
      "Confidence interval includes null effect"
    } else {
      "Confidence interval is wide"
    }
  } else {
    rating <- "No concerns"
    justification <- "Confidence interval is precise and does not include null"
  }

  list(
    rating = rating,
    justification = justification,
    ci_width = ci_width,
    crosses_null = crosses_null,
    effect = effect,
    lower = lower,
    upper = upper
  )
}

#' Domain 5: Heterogeneity Assessment for CINeMA
#'
#' @keywords internal
assess_heterogeneity_cinema <- function(nma_result, heterogeneity_result) {

  if (!is.null(heterogeneity_result)) {
    I2 <- heterogeneity_result$I2
    tau <- heterogeneity_result$tau
  } else {
    I2 <- nma_result$I2
    tau <- nma_result$tau
  }

  I2_percent <- I2 * 100

  if (I2_percent >= 75) {
    rating <- "Major concerns"
    justification <- sprintf("Considerable heterogeneity (I² = %.0f%%, τ = %.3f)", I2_percent, tau)
  } else if (I2_percent >= 50) {
    rating <- "Some concerns"
    justification <- sprintf("Substantial heterogeneity (I² = %.0f%%, τ = %.3f)", I2_percent, tau)
  } else {
    rating <- "No concerns"
    justification <- sprintf("Low to moderate heterogeneity (I² = %.0f%%, τ = %.3f)", I2_percent, tau)
  }

  list(
    rating = rating,
    justification = justification,
    I2 = I2,
    tau = tau
  )
}

#' Domain 6: Incoherence Assessment
#'
#' @keywords internal
assess_incoherence <- function(nma_result, inconsistency_result, comparison, reference) {

  if (is.null(inconsistency_result)) {
    return(list(
      rating = "No concerns (not assessed)",
      justification = "Inconsistency not formally tested"
    ))
  }

  # Check if this comparison shows inconsistency
  comp_name <- paste(
    pmin(comparison, reference),
    pmax(comparison, reference),
    sep = ":"
  )

  if (inherits(inconsistency_result, "node_split") &&
      !is.null(inconsistency_result$summary)) {

    comp_result <- inconsistency_result$summary[
      inconsistency_result$summary$comparison == comp_name,
    ]

    if (nrow(comp_result) > 0 && !is.na(comp_result$p_value[1])) {
      p_val <- comp_result$p_value[1]

      if (p_val < 0.05) {
        rating <- "Major concerns"
        justification <- sprintf("Significant inconsistency detected (p = %.3f)", p_val)
      } else if (p_val < 0.10) {
        rating <- "Some concerns"
        justification <- sprintf("Borderline inconsistency (p = %.3f)", p_val)
      } else {
        rating <- "No concerns"
        justification <- sprintf("No inconsistency detected (p = %.3f)", p_val)
      }

      return(list(
        rating = rating,
        justification = justification,
        p_value = p_val
      ))
    }
  }

  # Default if no specific test available
  list(
    rating = "No concerns (insufficient data)",
    justification = "Inconsistency test not available for this comparison"
  )
}

#' Calculate Overall Confidence Rating
#'
#' @keywords internal
calculate_overall_confidence <- function(within_bias, reporting_bias, indirectness,
                                        imprecision, heterogeneity, incoherence) {

  ratings <- c(
    within_bias$rating,
    reporting_bias$rating,
    indirectness$rating,
    imprecision$rating,
    heterogeneity$rating,
    incoherence$rating
  )

  # Count concerns
  n_major <- sum(grepl("Major concerns", ratings))
  n_some <- sum(grepl("Some concerns", ratings))
  n_no <- sum(grepl("No concerns", ratings))

  # CINeMA algorithm for overall rating
  if (n_major >= 2) {
    overall <- "Very low"
    justification <- sprintf("%d domains with major concerns", n_major)
  } else if (n_major == 1 && n_some >= 2) {
    overall <- "Low"
    justification <- sprintf("1 domain with major concerns, %d with some concerns", n_some)
  } else if (n_major == 1 || n_some >= 3) {
    overall <- "Moderate"
    justification <- if (n_major == 1) {
      "1 domain with major concerns"
    } else {
      sprintf("%d domains with some concerns", n_some)
    }
  } else {
    overall <- "High"
    justification <- "No major concerns; few domains with some concerns"
  }

  list(
    rating = overall,
    justification = justification,
    n_major_concerns = n_major,
    n_some_concerns = n_some,
    n_no_concerns = n_no
  )
}

#' Create CINeMA Contribution Matrix
#'
#' Calculate the contribution of each direct comparison to each network estimate.
#' This is a simplified version - full implementation requires netmeta internals.
#'
#' @param nma_result A netmeta object
#' @return Matrix of contributions (percentage)
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' contrib <- cinema_contribution_matrix(nma)
#' }
cinema_contribution_matrix <- function(nma_result) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Try to use netmeta's contribution matrix if available
  if (!is.null(nma_result$A.matrix)) {
    return(nma_result$A.matrix)
  }

  # Simplified calculation
  message("Full contribution matrix calculation requires netmeta internals")
  message("Returning simplified approximation")

  NULL
}

#' Visualize CINeMA Assessment
#'
#' Create a traffic light plot showing confidence ratings across comparisons.
#'
#' @param cinema_result A cinema_assessment object
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' cinema <- assess_cinema_confidence(nma, data)
#' plot_cinema_assessment(cinema)
#' }
plot_cinema_assessment <- function(cinema_result) {

  if (!inherits(cinema_result, "cinema_assessment")) {
    stop("cinema_result must be a cinema_assessment object")
  }

  # Prepare data for plotting
  summary <- cinema_result$summary

  # Convert to long format
  plot_data <- tidyr::pivot_longer(
    summary,
    cols = c(WithinStudyBias, ReportingBias, Indirectness,
            Imprecision, Heterogeneity, Incoherence),
    names_to = "Domain",
    values_to = "Rating"
  )

  # Clean domain names
  plot_data$Domain <- gsub("([a-z])([A-Z])", "\\1 \\2", plot_data$Domain)

  # Color mapping
  color_map <- c(
    "No concerns" = "#1a9850",
    "No concerns (not assessed)" = "#91cf60",
    "No concerns (format issue)" = "#91cf60",
    "No concerns (no data)" = "#91cf60",
    "No concerns (insufficient data)" = "#91cf60",
    "Some concerns" = "#fee08b",
    "Major concerns" = "#d73027"
  )

  # Simplify ratings for color mapping
  plot_data$Rating_Simple <- ifelse(
    grepl("No concerns", plot_data$Rating), "No concerns",
    ifelse(grepl("Some concerns", plot_data$Rating), "Some concerns",
          "Major concerns")
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Domain, y = Comparison, fill = Rating_Simple)) +
    ggplot2::geom_tile(color = "white", size = 1) +
    ggplot2::scale_fill_manual(
      values = c(
        "No concerns" = "#1a9850",
        "Some concerns" = "#fee08b",
        "Major concerns" = "#d73027"
      ),
      name = "Confidence"
    ) +
    ggplot2::labs(
      title = "CINeMA Confidence Assessment",
      subtitle = "Traffic light plot showing confidence across domains",
      x = "Domain",
      y = "Comparison"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Print CINeMA Assessment
#'
#' @param x cinema_assessment object
#' @param ... Additional arguments
#' @export
print.cinema_assessment <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════════════════\n")
  cat("  CINeMA: Confidence in Network Meta-Analysis\n")
  cat("══════════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Comparisons assessed: %d\n", x$n_comparisons))
  cat(sprintf("Reference treatment: %s\n\n", x$reference))

  cat("Summary of Confidence Ratings:\n\n")
  print(x$summary, row.names = FALSE)

  cat("\n")
  cat("Overall Confidence Levels:\n")
  conf_counts <- table(x$summary$OverallConfidence)
  for (level in names(conf_counts)) {
    cat(sprintf("  • %s: %d comparison(s)\n", level, conf_counts[level]))
  }

  cat("\n")
  cat("CINeMA Domains:\n")
  cat("  1. Within-study bias - Risk of bias in contributing studies\n")
  cat("  2. Reporting bias - Publication and selective reporting bias\n")
  cat("  3. Indirectness - Availability of direct vs indirect evidence\n")
  cat("  4. Imprecision - Width of confidence intervals\n")
  cat("  5. Heterogeneity - Between-study variability\n")
  cat("  6. Incoherence - Statistical inconsistency in the network\n")

  cat("\n")
  cat("Use plot() to visualize the traffic light plot\n")
  cat("Access detailed assessments via result$assessments\n\n")

  invisible(x)
}

#' Plot CINeMA Assessment
#'
#' @param x cinema_assessment object
#' @param ... Additional arguments
#' @export
plot.cinema_assessment <- function(x, ...) {
  plot_cinema_assessment(x)
}
