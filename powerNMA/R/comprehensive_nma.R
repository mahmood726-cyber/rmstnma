#' Comprehensive Network Meta-Analysis with All Advanced Features
#'
#' Perform a complete NMA analysis including all Phase 5 enhancements:
#' rankings, inconsistency assessment, heterogeneity, publication bias,
#' visualization, and comprehensive reporting.
#'
#' @name comprehensive_nma
NULL

#' Run Comprehensive Network Meta-Analysis
#'
#' Execute a complete NMA pipeline with all advanced statistical assessments
#' and publication-ready outputs. This is the main function for conducting
#' state-of-the-art network meta-analyses.
#'
#' @param data Pairwise NMA data frame with columns: treat1, treat2, TE, seTE, studlab
#' @param reference Reference treatment (auto-detected if NULL)
#' @param sm Summary measure: "OR", "RR", "HR", "MD", "SMD"
#' @param small_values_good For rankings: are smaller values better? (e.g., TRUE for adverse events)
#' @param assess_inconsistency Run inconsistency tests (node-splitting, etc.)
#' @param assess_publication_bias Run publication bias tests
#' @param generate_plots Create all visualization plots
#' @param output_dir Directory for saving results (NULL for no saving)
#' @return comprehensive_nma_result object with all analyses
#' @export
#' @examples
#' \dontrun{
#' # Complete analysis with all features
#' data <- simulate_nma_data(n_studies = 30)
#' results <- run_comprehensive_nma(
#'   data = data,
#'   sm = "OR",
#'   small_values_good = TRUE,
#'   output_dir = "nma_results"
#' )
#'
#' # View results
#' print(results)
#' summary(results)
#'
#' # Access individual components
#' results$rankings        # SUCRA scores
#' results$inconsistency   # Node-splitting results
#' results$heterogeneity   # I², τ² statistics
#' results$pub_bias        # Publication bias tests
#' results$league_table    # League table
#' }
run_comprehensive_nma <- function(data,
                                 reference = NULL,
                                 sm = c("OR", "RR", "HR", "MD", "SMD"),
                                 small_values_good = TRUE,
                                 assess_inconsistency = TRUE,
                                 assess_publication_bias = TRUE,
                                 generate_plots = TRUE,
                                 output_dir = NULL) {

  sm <- match.arg(sm)

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Comprehensive Network Meta-Analysis\n")
  cat("  Phase 5: Advanced Statistical Methods & Publication Tools\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Validate data
  required_cols <- c("treat1", "treat2", "TE", "seTE", "studlab")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  cat("✓ Data validated\n")
  cat(sprintf("  • %d studies\n", length(unique(data$studlab))))
  cat(sprintf("  • %d treatments\n", length(unique(c(data$treat1, data$treat2)))))
  cat(sprintf("  • %d comparisons\n\n", nrow(data)))

  # Auto-detect reference
  if (is.null(reference)) {
    counts <- table(c(as.character(data$treat1), as.character(data$treat2)))
    reference <- names(sort(counts, decreasing = TRUE))[1]
    cat(sprintf("✓ Reference treatment auto-detected: '%s'\n\n", reference))
  }

  results <- list()

  # ═══════════════════════════════════════════════════════════════
  # 1. Core Network Meta-Analysis
  # ═══════════════════════════════════════════════════════════════

  cat("STEP 1: Running core network meta-analysis...\n")

  nma <- netmeta::netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = treat1,
    treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference,
    comb.fixed = TRUE,
    comb.random = TRUE
  )

  results$nma <- nma
  cat("  ✓ Network meta-analysis complete\n\n")

  # ═══════════════════════════════════════════════════════════════
  # 2. Treatment Rankings & SUCRA
  # ═══════════════════════════════════════════════════════════════

  cat("STEP 2: Calculating treatment rankings and SUCRA scores...\n")

  sucra <- tryCatch({
    calculate_sucra(nma, small_values = ifelse(small_values_good, "good", "bad"))
  }, error = function(e) {
    warning(sprintf("SUCRA calculation failed: %s", e$message))
    NULL
  })

  results$sucra <- sucra
  results$rankings <- sucra

  if (!is.null(sucra)) {
    cat("  ✓ SUCRA scores calculated\n")
    cat(sprintf("    • Best treatment: %s (SUCRA: %.1f%%)\n",
               sucra$summary$Treatment[1],
               sucra$summary$SUCRA[1]))
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 3. Heterogeneity Assessment
  # ═══════════════════════════════════════════════════════════════

  cat("STEP 3: Assessing heterogeneity...\n")

  het <- tryCatch({
    heterogeneity_report(nma, reference)
  }, error = function(e) {
    warning(sprintf("Heterogeneity assessment failed: %s", e$message))
    NULL
  })

  results$heterogeneity <- het

  if (!is.null(het)) {
    cat(sprintf("  ✓ Heterogeneity assessed: I² = %.1f%%, τ = %.3f\n",
               het$I2_percent, het$tau))
    cat(sprintf("    • %s\n", het$i2_interpretation))
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 4. Inconsistency Assessment
  # ═══════════════════════════════════════════════════════════════

  if (assess_inconsistency) {
    cat("STEP 4: Assessing network inconsistency...\n")

    # Node-splitting
    node_split <- tryCatch({
      node_splitting(nma, data, p_threshold = 0.10)
    }, error = function(e) {
      warning(sprintf("Node-splitting failed: %s", e$message))
      NULL
    })

    results$node_splitting <- node_split

    if (!is.null(node_split) && node_split$n_tests > 0) {
      cat(sprintf("  ✓ Node-splitting: %d comparisons tested\n", node_split$n_tests))
      if (node_split$n_inconsistent > 0) {
        cat(sprintf("    ⚠ %d comparisons show inconsistency\n", node_split$n_inconsistent))
      } else {
        cat("    ✓ No significant inconsistency detected\n")
      }
    }

    # Design inconsistency
    design_inc <- tryCatch({
      design_inconsistency(nma, data)
    }, error = function(e) {
      warning(sprintf("Design inconsistency test failed: %s", e$message))
      NULL
    })

    results$design_inconsistency <- design_inc

    if (!is.null(design_inc)) {
      cat(sprintf("  ✓ Design inconsistency: %s\n", design_inc$interpretation))
    }

    # Loop inconsistency
    loop_inc <- tryCatch({
      loop_inconsistency(nma, data)
    }, error = function(e) {
      warning(sprintf("Loop inconsistency test failed: %s", e$message))
      NULL
    })

    results$loop_inconsistency <- loop_inc

    cat("\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # 5. Publication Bias Assessment
  # ═══════════════════════════════════════════════════════════════

  if (assess_publication_bias) {
    cat("STEP 5: Assessing publication bias...\n")

    pub_bias <- tryCatch({
      assess_publication_bias(nma, data, method = "both")
    }, error = function(e) {
      warning(sprintf("Publication bias assessment failed: %s", e$message))
      NULL
    })

    results$publication_bias <- pub_bias
    results$pub_bias <- pub_bias

    if (!is.null(pub_bias) && !is.null(pub_bias$summary)) {
      n_sig <- sum(pub_bias$summary$Significant, na.rm = TRUE)
      if (n_sig > 0) {
        cat(sprintf("  ⚠ %d of %d tests suggest potential bias\n", n_sig, nrow(pub_bias$summary)))
      } else {
        cat("  ✓ No significant publication bias detected\n")
      }
    }
    cat("\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # 6. League Table
  # ═══════════════════════════════════════════════════════════════

  cat("STEP 6: Generating league table...\n")

  league <- tryCatch({
    create_league_table(nma, sucra, digits = 2, format = "effect_ci", sort_by = "sucra")
  }, error = function(e) {
    warning(sprintf("League table generation failed: %s", e$message))
    NULL
  })

  results$league_table <- league

  if (!is.null(league)) {
    cat("  ✓ League table generated\n")
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 7. Network Geometry
  # ═══════════════════════════════════════════════════════════════

  cat("STEP 7: Calculating network geometry...\n")

  geometry <- tryCatch({
    network_geometry(data)
  }, error = function(e) {
    warning(sprintf("Network geometry calculation failed: %s", e$message))
    NULL
  })

  results$network_geometry <- geometry

  if (!is.null(geometry)) {
    cat(sprintf("  ✓ Network density: %.2f\n", geometry$density))
    cat(sprintf("    • Mean degree: %.1f\n", geometry$mean_degree))
  }
  cat("\n")

  # ═══════════════════════════════════════════════════════════════
  # 8. Visualization
  # ═══════════════════════════════════════════════════════════════

  if (generate_plots) {
    cat("STEP 8: Generating visualization plots...\n")

    plots <- tryCatch({
      generate_nma_plots(nma, data, sucra, output_dir = if (!is.null(output_dir)) {
        file.path(output_dir, "plots")
      } else NULL)
    }, error = function(e) {
      warning(sprintf("Plot generation failed: %s", e$message))
      list()
    })

    results$plots <- plots
    cat(sprintf("  ✓ Generated %d plots\n", length(plots)))
    cat("\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # 9. Save Results
  # ═══════════════════════════════════════════════════════════════

  if (!is.null(output_dir)) {
    cat("STEP 9: Saving results...\n")

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save league table
    if (!is.null(league)) {
      tryCatch({
        write_league_table(league, file.path(output_dir, "league_table.csv"))
        cat("  ✓ League table saved\n")
      }, error = function(e) {
        warning(sprintf("Failed to save league table: %s", e$message))
      })
    }

    # Save SUCRA rankings
    if (!is.null(sucra)) {
      tryCatch({
        utils::write.csv(sucra$summary,
                        file.path(output_dir, "sucra_rankings.csv"),
                        row.names = FALSE)
        cat("  ✓ SUCRA rankings saved\n")
      }, error = function(e) {
        warning(sprintf("Failed to save SUCRA rankings: %s", e$message))
      })
    }

    # Save heterogeneity report
    if (!is.null(het)) {
      tryCatch({
        sink(file.path(output_dir, "heterogeneity_report.txt"))
        print(het)
        sink()
        cat("  ✓ Heterogeneity report saved\n")
      }, error = function(e) {
        sink()
        warning(sprintf("Failed to save heterogeneity report: %s", e$message))
      })
    }

    # Save publication bias report
    if (!is.null(pub_bias)) {
      tryCatch({
        report <- generate_bias_report(pub_bias, nma,
                                      output_file = file.path(output_dir, "publication_bias_report.txt"))
        cat("  ✓ Publication bias report saved\n")
      }, error = function(e) {
        warning(sprintf("Failed to save publication bias report: %s", e$message))
      })
    }

    # Save comprehensive summary
    tryCatch({
      sink(file.path(output_dir, "comprehensive_summary.txt"))
      print_comprehensive_summary(results)
      sink()
      cat("  ✓ Comprehensive summary saved\n")
    }, error = function(e) {
      sink()
      warning(sprintf("Failed to save comprehensive summary: %s", e$message))
    })

    cat(sprintf("\n  All results saved to: %s\n", output_dir))
    cat("\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Finalize
  # ═══════════════════════════════════════════════════════════════

  results$data <- data
  results$reference <- reference
  results$sm <- sm

  class(results) <- c("comprehensive_nma_result", "list")

  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  ✓ Comprehensive NMA analysis complete!\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  cat("Access results:\n")
  cat("  • results$nma              - Core NMA results\n")
  cat("  • results$rankings         - SUCRA scores\n")
  cat("  • results$heterogeneity    - Heterogeneity statistics\n")
  cat("  • results$inconsistency    - Inconsistency tests\n")
  cat("  • results$pub_bias         - Publication bias tests\n")
  cat("  • results$league_table     - League table\n")
  cat("  • results$plots            - All visualization plots\n\n")

  cat("Use print(results) or summary(results) for overview\n\n")

  results
}

#' Print Comprehensive Summary
#'
#' @param results comprehensive_nma_result object
#' @keywords internal
print_comprehensive_summary <- function(results) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Comprehensive Network Meta-Analysis Summary\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Network characteristics
  cat("NETWORK CHARACTERISTICS\n")
  cat("───────────────────────────────────────────────────────────────\n")
  if (!is.null(results$network_geometry)) {
    cat(sprintf("Studies: %d\n", results$network_geometry$n_studies))
    cat(sprintf("Treatments: %d\n", results$network_geometry$n_treatments))
    cat(sprintf("Comparisons: %d\n", results$network_geometry$n_comparisons))
    cat(sprintf("Network density: %.2f\n", results$network_geometry$density))
  }
  cat("\n")

  # Treatment rankings
  if (!is.null(results$rankings)) {
    cat("TOP 5 TREATMENTS (by SUCRA)\n")
    cat("───────────────────────────────────────────────────────────────\n")
    top5 <- head(results$rankings$summary, 5)
    for (i in 1:nrow(top5)) {
      cat(sprintf("%d. %s - SUCRA: %.1f%%, P(best): %.1f%%\n",
                 i, top5$Treatment[i], top5$SUCRA[i], top5$ProbBest[i]))
    }
    cat("\n")
  }

  # Heterogeneity
  if (!is.null(results$heterogeneity)) {
    cat("HETEROGENEITY\n")
    cat("───────────────────────────────────────────────────────────────\n")
    cat(sprintf("I² = %.1f%% - %s\n",
               results$heterogeneity$I2_percent,
               results$heterogeneity$i2_interpretation))
    cat(sprintf("τ = %.3f - %s\n",
               results$heterogeneity$tau,
               results$heterogeneity$tau_interpretation))
    cat("\n")
  }

  # Inconsistency
  if (!is.null(results$node_splitting)) {
    cat("INCONSISTENCY\n")
    cat("───────────────────────────────────────────────────────────────\n")
    cat(sprintf("Node-splitting: %d of %d comparisons inconsistent\n",
               results$node_splitting$n_inconsistent,
               results$node_splitting$n_tests))

    if (results$node_splitting$n_inconsistent > 0) {
      cat("⚠ Review inconsistent comparisons\n")
    } else {
      cat("✓ Network appears consistent\n")
    }
    cat("\n")
  }

  # Publication bias
  if (!is.null(results$pub_bias) && !is.null(results$pub_bias$summary)) {
    cat("PUBLICATION BIAS\n")
    cat("───────────────────────────────────────────────────────────────\n")
    n_sig <- sum(results$pub_bias$summary$Significant, na.rm = TRUE)
    cat(sprintf("%d of %d tests significant\n",
               n_sig, nrow(results$pub_bias$summary)))

    if (n_sig > 0) {
      cat("⚠ Potential publication bias detected\n")
    } else {
      cat("✓ No evidence of publication bias\n")
    }
    cat("\n")
  }

  cat("═══════════════════════════════════════════════════════════════\n")
}

#' Print Comprehensive NMA Result
#'
#' @param x comprehensive_nma_result object
#' @param ... Additional arguments
#' @export
print.comprehensive_nma_result <- function(x, ...) {
  print_comprehensive_summary(x)
  invisible(x)
}

#' Summary of Comprehensive NMA Result
#'
#' @param object comprehensive_nma_result object
#' @param ... Additional arguments
#' @export
summary.comprehensive_nma_result <- function(object, ...) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Comprehensive NMA: Detailed Summary\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Print each component
  if (!is.null(object$nma)) {
    cat("CORE NMA RESULTS:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(object$nma)
    cat("\n\n")
  }

  if (!is.null(object$rankings)) {
    cat("TREATMENT RANKINGS:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(object$rankings)
    cat("\n\n")
  }

  if (!is.null(object$heterogeneity)) {
    cat("HETEROGENEITY ASSESSMENT:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(object$heterogeneity)
    cat("\n\n")
  }

  if (!is.null(object$node_splitting)) {
    cat("INCONSISTENCY ASSESSMENT (Node-Splitting):\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(object$node_splitting)
    cat("\n\n")
  }

  if (!is.null(object$pub_bias)) {
    cat("PUBLICATION BIAS ASSESSMENT:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(object$pub_bias)
    cat("\n\n")
  }

  invisible(object)
}
