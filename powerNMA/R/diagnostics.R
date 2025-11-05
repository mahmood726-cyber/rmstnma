# ============================================================================
# Diagnostic Utilities for powerNMA
# ============================================================================
#
# Provides comprehensive diagnostic functions for assessing data quality,
# network properties, and analysis results.
#

#' Comprehensive Data Quality Diagnostics
#'
#' Performs a thorough check of NMA data quality and provides actionable
#' recommendations for improvement.
#'
#' @param data Data frame with NMA data (pairwise format)
#' @param return_details Logical; if TRUE, returns detailed diagnostic object
#' @return If return_details=FALSE, prints summary. If TRUE, returns list with diagnostics.
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' diagnose_nma_data(data)
#'
#' # Get detailed diagnostics
#' diag <- diagnose_nma_data(data, return_details = TRUE)
#' names(diag)
#' }
diagnose_nma_data <- function(data, return_details = FALSE) {
  # Validate basic structure first
  if (!is.data.frame(data)) {
    stop("[diagnose_nma_data] Input must be a data frame", call. = FALSE)
  }

  diagnostics <- list()
  issues <- character()
  warnings_list <- character()
  recommendations <- character()

  # 1. Check required columns
  required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    issues <- c(issues, sprintf("Missing required columns: %s",
                               paste(missing_cols, collapse = ", ")))
  }

  diagnostics$has_required_columns <- length(missing_cols) == 0

  if (!diagnostics$has_required_columns) {
    # Can't continue without required columns
    diagnostics$issues <- issues
    diagnostics$status <- "FAILED"
    if (!return_details) {
      cat("\n❌ DIAGNOSTIC FAILED\n")
      cat("Issues:\n")
      for (issue in issues) cat("  •", issue, "\n")
    }
    return(if (return_details) diagnostics else invisible(NULL))
  }

  # 2. Sample size
  diagnostics$n_studies <- length(unique(data$studlab))
  diagnostics$n_comparisons <- nrow(data)
  diagnostics$n_treatments <- length(unique(c(data$treat1, data$treat2)))

  if (diagnostics$n_studies < 3) {
    warnings_list <- c(warnings_list, sprintf("Only %d studies - results may be unreliable",
                                              diagnostics$n_studies))
    recommendations <- c(recommendations, "Consider adding more studies or using caution in interpretation")
  }

  if (diagnostics$n_comparisons < .MIN_COMPARISONS_FOR_NMA) {
    issues <- c(issues, sprintf("Only %d comparison(s) - need at least %d",
                               diagnostics$n_comparisons, .MIN_COMPARISONS_FOR_NMA))
  }

  # 3. Data quality checks
  diagnostics$has_na_TE <- any(is.na(data$TE))
  diagnostics$has_na_seTE <- any(is.na(data$seTE))
  diagnostics$has_inf_TE <- any(!is.finite(data$TE))
  diagnostics$has_inf_seTE <- any(!is.finite(data$seTE))
  diagnostics$has_negative_seTE <- any(data$seTE <= 0, na.rm = TRUE)

  if (diagnostics$has_na_TE) {
    issues <- c(issues, sprintf("Found %d NA values in TE column", sum(is.na(data$TE))))
  }
  if (diagnostics$has_na_seTE) {
    issues <- c(issues, sprintf("Found %d NA values in seTE column", sum(is.na(data$seTE))))
  }
  if (diagnostics$has_inf_TE) {
    issues <- c(issues, sprintf("Found %d non-finite values in TE column",
                               sum(!is.finite(data$TE))))
  }
  if (diagnostics$has_inf_seTE) {
    issues <- c(issues, sprintf("Found %d non-finite values in seTE column",
                               sum(!is.finite(data$seTE))))
  }
  if (diagnostics$has_negative_seTE) {
    issues <- c(issues, sprintf("Found %d non-positive seTE values", sum(data$seTE <= 0, na.rm = TRUE)))
  }

  # 4. Effect size distribution
  diagnostics$TE_range <- range(data$TE[is.finite(data$TE)])
  diagnostics$TE_mean <- mean(data$TE[is.finite(data$TE)])
  diagnostics$TE_sd <- sd(data$TE[is.finite(data$TE)])

  # Check for extreme values
  if (diagnostics$TE_range[2] - diagnostics$TE_range[1] > 10) {
    warnings_list <- c(warnings_list, sprintf("Very wide effect size range: [%.2f, %.2f]",
                                              diagnostics$TE_range[1], diagnostics$TE_range[2]))
    recommendations <- c(recommendations, "Check for outliers or data entry errors")
  }

  # 5. Standard error distribution
  diagnostics$seTE_range <- range(data$seTE[is.finite(data$seTE)])
  diagnostics$seTE_mean <- mean(data$seTE[is.finite(data$seTE)])
  diagnostics$seTE_median <- median(data$seTE[is.finite(data$seTE)])

  # Check for very large or very small SE
  if (diagnostics$seTE_range[2] > 5) {
    warnings_list <- c(warnings_list, sprintf("Some very large standard errors (max: %.2f)",
                                              diagnostics$seTE_range[2]))
    recommendations <- c(recommendations, "Check for small sample sizes or data quality issues")
  }

  # 6. Network connectivity
  treatment_counts <- table(c(data$treat1, data$treat2))
  diagnostics$treatment_counts <- as.data.frame(treatment_counts)
  names(diagnostics$treatment_counts) <- c("Treatment", "N_Comparisons")

  # Find treatments with few comparisons
  sparse_treatments <- names(treatment_counts)[treatment_counts < 2]
  if (length(sparse_treatments) > 0) {
    warnings_list <- c(warnings_list, sprintf("%d treatment(s) with < 2 comparisons: %s",
                                              length(sparse_treatments),
                                              paste(sparse_treatments, collapse = ", ")))
    recommendations <- c(recommendations, "Consider removing poorly connected treatments")
  }

  # 7. Multi-arm trials detection
  comparisons_per_study <- table(data$studlab)
  multi_arm_studies <- names(comparisons_per_study)[comparisons_per_study >= 3]
  diagnostics$n_multi_arm_studies <- length(multi_arm_studies)
  diagnostics$has_multi_arm_trials <- diagnostics$n_multi_arm_studies > 0

  if (diagnostics$has_multi_arm_trials) {
    diagnostics$multi_arm_studies <- multi_arm_studies
  }

  # 8. Reference treatment detection
  most_common_treatment <- names(sort(treatment_counts, decreasing = TRUE))[1]
  diagnostics$suggested_reference <- most_common_treatment
  diagnostics$reference_n_comparisons <- as.numeric(treatment_counts[most_common_treatment])

  # Overall status
  diagnostics$n_issues <- length(issues)
  diagnostics$n_warnings <- length(warnings_list)
  diagnostics$issues <- issues
  diagnostics$warnings <- warnings_list
  diagnostics$recommendations <- recommendations

  if (length(issues) == 0) {
    diagnostics$status <- "PASS"
  } else if (length(issues) < 3) {
    diagnostics$status <- "PASS_WITH_WARNINGS"
  } else {
    diagnostics$status <- "FAILED"
  }

  # Print summary if not returning details
  if (!return_details) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("  NMA Data Quality Diagnostics\n")
    cat("══════════════════════════════════════════════════════════\n\n")

    cat("📊 Data Summary:\n")
    cat(sprintf("  • Studies: %d\n", diagnostics$n_studies))
    cat(sprintf("  • Comparisons: %d\n", diagnostics$n_comparisons))
    cat(sprintf("  • Treatments: %d\n", diagnostics$n_treatments))
    if (diagnostics$has_multi_arm_trials) {
      cat(sprintf("  • Multi-arm trials: %d\n", diagnostics$n_multi_arm_studies))
    }
    cat("\n")

    cat("📈 Effect Sizes (TE):\n")
    cat(sprintf("  • Range: [%.3f, %.3f]\n", diagnostics$TE_range[1], diagnostics$TE_range[2]))
    cat(sprintf("  • Mean (SD): %.3f (%.3f)\n", diagnostics$TE_mean, diagnostics$TE_sd))
    cat("\n")

    cat("🎯 Standard Errors (seTE):\n")
    cat(sprintf("  • Range: [%.3f, %.3f]\n", diagnostics$seTE_range[1], diagnostics$seTE_range[2]))
    cat(sprintf("  • Median: %.3f\n", diagnostics$seTE_median))
    cat("\n")

    cat("🔗 Network Properties:\n")
    cat(sprintf("  • Suggested reference: %s (%d comparisons)\n",
                diagnostics$suggested_reference, diagnostics$reference_n_comparisons))
    cat("\n")

    # Status
    status_symbol <- switch(diagnostics$status,
                           "PASS" = "✅",
                           "PASS_WITH_WARNINGS" = "⚠️ ",
                           "FAILED" = "❌",
                           "❓")
    cat(sprintf("%s Overall Status: %s\n\n", status_symbol, diagnostics$status))

    # Issues
    if (length(issues) > 0) {
      cat("❌ Issues Found:\n")
      for (issue in issues) {
        cat("  •", issue, "\n")
      }
      cat("\n")
    }

    # Warnings
    if (length(warnings_list) > 0) {
      cat("⚠️  Warnings:\n")
      for (warn in warnings_list) {
        cat("  •", warn, "\n")
      }
      cat("\n")
    }

    # Recommendations
    if (length(recommendations) > 0) {
      cat("💡 Recommendations:\n")
      for (rec in recommendations) {
        cat("  •", rec, "\n")
      }
      cat("\n")
    }

    if (diagnostics$status == "PASS") {
      cat("✅ Data is ready for network meta-analysis!\n\n")
    }

    cat("══════════════════════════════════════════════════════════\n\n")
  }

  if (return_details) {
    return(diagnostics)
  } else {
    invisible(diagnostics)
  }
}

#' Quick Data Summary
#'
#' Provides a one-line summary of NMA data.
#'
#' @param data Data frame with NMA data
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 20)
#' quick_summary(data)
quick_summary <- function(data) {
  n_studies <- length(unique(data$studlab))
  n_comparisons <- nrow(data)
  n_treatments <- length(unique(c(data$treat1, data$treat2)))

  cat(sprintf("📊 %d studies, %d comparisons, %d treatments\n",
              n_studies, n_comparisons, n_treatments))
  invisible(list(n_studies = n_studies,
                n_comparisons = n_comparisons,
                n_treatments = n_treatments))
}

#' Diagnose Network Connectivity
#'
#' Checks if the network is connected and identifies disconnected components.
#'
#' @param data Data frame with NMA data
#' @return List with connectivity diagnostics
#' @export
diagnose_network_connectivity <- function(data) {
  if (!has_pkg("igraph")) {
    stop("[diagnose_network_connectivity] Package 'igraph' required", call. = FALSE)
  }

  # Build network graph
  edges <- data.frame(
    from = as.character(data$treat1),
    to = as.character(data$treat2),
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  # Check connectivity
  is_connected <- igraph::is_connected(g)
  n_components <- igraph::count_components(g)
  components <- igraph::components(g)

  result <- list(
    is_connected = is_connected,
    n_components = n_components,
    component_sizes = components$csize,
    treatments_by_component = split(names(components$membership), components$membership)
  )

  # Print summary
  if (is_connected) {
    cat("✅ Network is fully connected\n")
  } else {
    cat(sprintf("⚠️  Network has %d disconnected components\n", n_components))
    for (i in seq_along(result$treatments_by_component)) {
      comp_treatments <- result$treatments_by_component[[i]]
      cat(sprintf("   Component %d (%d treatments): %s\n",
                  i, length(comp_treatments),
                  paste(comp_treatments, collapse = ", ")))
    }
  }

  invisible(result)
}

#' Print Data Quality Report Card
#'
#' Creates a visual report card for data quality.
#'
#' @param diagnostics Output from diagnose_nma_data(return_details = TRUE)
#' @export
print_quality_report_card <- function(diagnostics) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║          NMA Data Quality Report Card                   ║\n")
  cat("╠══════════════════════════════════════════════════════════╣\n")

  # Grading criteria
  grade_required_cols <- if(diagnostics$has_required_columns) "A" else "F"
  grade_sample_size <- if(diagnostics$n_studies >= 10) "A" else if(diagnostics$n_studies >= 5) "B" else "C"
  grade_data_quality <- if(diagnostics$n_issues == 0) "A" else if(diagnostics$n_issues <= 2) "B" else "F"
  grade_connectivity <- if(diagnostics$n_treatments >= 3) "A" else "C"

  cat(sprintf("║  Required Columns          [%s]                          ║\n", grade_required_cols))
  cat(sprintf("║  Sample Size               [%s]  (%d studies)            ║\n",
              grade_sample_size, diagnostics$n_studies))
  cat(sprintf("║  Data Quality              [%s]  (%d issues)             ║\n",
              grade_data_quality, diagnostics$n_issues))
  cat(sprintf("║  Network Connectivity      [%s]  (%d treatments)         ║\n",
              grade_connectivity, diagnostics$n_treatments))
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # Overall grade
  grades <- c(grade_required_cols, grade_sample_size, grade_data_quality, grade_connectivity)
  if (all(grades == "A")) {
    cat("🏆 Overall Grade: A - Excellent!\n")
  } else if ("F" %in% grades) {
    cat("❌ Overall Grade: F - Needs significant improvement\n")
  } else if (any(c("B", "C") %in% grades)) {
    cat("⚠️  Overall Grade: B - Good, with room for improvement\n")
  }
  cat("\n")
}
