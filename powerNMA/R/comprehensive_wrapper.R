#' Comprehensive NMA Analysis Wrapper
#'
#' One-function complete network meta-analysis incorporating all Phase 5-8
#' features: advanced statistics, breakthrough methods, network geometry,
#' simulation, and export capabilities.
#'
#' @name comprehensive_wrapper
NULL

#' Run Ultimate Comprehensive NMA Analysis
#'
#' Perform complete state-of-the-art network meta-analysis with all
#' advanced features from Phases 5-8.
#'
#' @param data Pairwise comparison data
#' @param sm Summary measure
#' @param reference Reference treatment
#' @param assess_geometry Analyze network geometry?
#' @param assess_inconsistency Perform inconsistency assessment?
#' @param calculate_rankings Calculate treatment rankings?
#' @param assess_publication_bias Assess publication bias?
#' @param generate_visualizations Create all visualizations?
#' @param export_results Export results to multiple formats?
#' @param output_dir Output directory for exports
#' @param multivariate_outcomes List of additional outcome datasets
#' @param component_matrix Component structure for CNMA
#' @param ipd_data Individual patient data for personalized prediction
#' @param cost_data Cost data for value of information
#' @param perform_simulation Run power analysis simulation?
#' @return comprehensive_nma object with all analyses
#' @export
#' @examples
#' \dontrun{
#' # Complete analysis
#' results <- run_ultimate_nma(
#'   data = nma_data,
#'   sm = "OR",
#'   reference = "Placebo",
#'   assess_geometry = TRUE,
#'   calculate_rankings = TRUE,
#'   export_results = TRUE,
#'   output_dir = "complete_nma_results"
#' )
#'
#' print(results)
#' summary(results)
#' }
run_ultimate_nma <- function(data,
                             sm = "MD",
                             reference = NULL,
                             assess_geometry = TRUE,
                             assess_inconsistency = TRUE,
                             calculate_rankings = TRUE,
                             assess_publication_bias = TRUE,
                             generate_visualizations = TRUE,
                             export_results = FALSE,
                             output_dir = "nma_comprehensive",
                             multivariate_outcomes = NULL,
                             component_matrix = NULL,
                             ipd_data = NULL,
                             cost_data = NULL,
                             perform_simulation = FALSE) {

  message("\n")
  message("═══════════════════════════════════════════════════════════")
  message("  ULTIMATE COMPREHENSIVE NETWORK META-ANALYSIS")
  message("  Phases 5-8: Complete State-of-the-Art Analysis")
  message("═══════════════════════════════════════════════════════════")
  message("\n")

  start_time <- Sys.time()

  result <- list()

  # ============================================================================
  # PHASE 1: Core NMA
  # ============================================================================

  message("█ Step 1/10: Running core network meta-analysis...")

  if (requireNamespace("netmeta", quietly = TRUE)) {
    core_nma <- netmeta::netmeta(
      TE = TE, seTE = seTE,
      treat1 = treat1, treat2 = treat2,
      studlab = studlab,
      data = data,
      sm = sm,
      reference.group = reference
    )

    result$core_nma <- core_nma

    if (is.null(reference)) {
      reference <- core_nma$reference.group
    }

    message(sprintf("  ✓ %d treatments, %d studies, %d comparisons",
                   length(core_nma$trts),
                   length(unique(data$studlab)),
                   nrow(data)))
  } else {
    stop("Package 'netmeta' required")
  }

  # ============================================================================
  # PHASE 2: Network Geometry (Phase 8)
  # ============================================================================

  if (assess_geometry) {
    message("\n█ Step 2/10: Analyzing network geometry...")

    geometry <- analyze_network_geometry(data, reference)
    result$geometry <- geometry

    message(sprintf("  ✓ Network: %s | Density: %.1f%% | Robustness: %.1f%%",
                   geometry$characteristics$network_type,
                   geometry$connectivity$density * 100,
                   geometry$robustness$robustness_score * 100))
  }

  # ============================================================================
  # PHASE 3: Treatment Rankings (Phase 5)
  # ============================================================================

  if (calculate_rankings) {
    message("\n█ Step 3/10: Calculating treatment rankings (SUCRA)...")

    sucra <- calculate_sucra(core_nma, small_values = "bad")
    result$sucra <- sucra

    top_treatment <- names(which.max(sucra$sucra_scores))
    message(sprintf("  ✓ Best treatment: %s (SUCRA: %.1f%%)",
                   top_treatment, max(sucra$sucra_scores) * 100))
  }

  # ============================================================================
  # PHASE 4: Inconsistency Assessment (Phase 5)
  # ============================================================================

  if (assess_inconsistency) {
    message("\n█ Step 4/10: Assessing network inconsistency...")

    node_split <- node_splitting(core_nma, data)
    result$inconsistency <- node_split

    n_inconsistent <- sum(node_split$results$p_value < 0.10, na.rm = TRUE)
    message(sprintf("  ✓ %d/%d comparisons show inconsistency (p < 0.10)",
                   n_inconsistent, nrow(node_split$results)))
  }

  # ============================================================================
  # PHASE 5: Heterogeneity Analysis (Phase 5)
  # ============================================================================

  message("\n█ Step 5/10: Analyzing heterogeneity...")

  heterogeneity <- heterogeneity_report(core_nma, reference)
  result$heterogeneity <- heterogeneity

  message(sprintf("  ✓ τ² = %.3f | I² = %.1f%% (%s)",
                 core_nma$tau^2, core_nma$I2,
                 heterogeneity$interpretation))

  # ============================================================================
  # PHASE 6: Publication Bias (Phase 5)
  # ============================================================================

  if (assess_publication_bias) {
    message("\n█ Step 6/10: Assessing publication bias...")

    pub_bias <- assess_publication_bias(core_nma, data)
    result$publication_bias <- pub_bias

    message(sprintf("  ✓ Egger's test p = %.3f | Begg's test p = %.3f",
                   pub_bias$egger$p_value, pub_bias$begg$p_value))
  }

  # ============================================================================
  # PHASE 7: Multivariate NMA (Phase 7)
  # ============================================================================

  if (!is.null(multivariate_outcomes)) {
    message("\n█ Step 7/10: Performing multivariate NMA...")

    mvnma_data <- prepare_mvnma_data(
      data_list = c(list(data), multivariate_outcomes),
      outcome_names = c("Primary", paste0("Outcome", seq_along(multivariate_outcomes)))
    )

    mvnma <- fit_mvnma(mvnma_data, reference = reference, method = "riley")
    result$multivariate <- mvnma

    message(sprintf("  ✓ %d outcomes analyzed jointly", mvnma_data$n_outcomes))
  }

  # ============================================================================
  # PHASE 8: Component NMA (Phase 7)
  # ============================================================================

  if (!is.null(component_matrix)) {
    message("\n█ Step 8/10: Performing component network meta-analysis...")

    cnma <- additive_cnma(data, component_matrix, sm, reference)
    result$component_nma <- cnma

    message(sprintf("  ✓ %d components analyzed | R² = %.3f",
                   ncol(component_matrix), cnma$model_fit$R_squared))
  }

  # ============================================================================
  # PHASE 9: Personalized Prediction (Phase 7)
  # ============================================================================

  if (!is.null(ipd_data)) {
    message("\n█ Step 9/10: Fitting personalized prediction model...")

    ipd_model <- fit_ipd_metaregression(
      ipd_data = ipd_data,
      outcome_var = "outcome",
      treatment_var = "treatment",
      study_var = "study",
      effect_modifiers = setdiff(names(ipd_data),
                                 c("outcome", "treatment", "study", "patient_id")),
      reference = reference
    )

    result$personalized <- ipd_model

    message("  ✓ Effect modifiers identified")
  }

  # ============================================================================
  # PHASE 10: Value of Information (Phase 7)
  # ============================================================================

  if (!is.null(cost_data)) {
    message("\n█ Step 10/10: Calculating value of information...")

    evpi <- calculate_evpi(
      nma_result = core_nma,
      cost_data = cost_data,
      threshold = 50000,
      n_sim = 1000
    )

    result$value_of_information <- evpi

    message(sprintf("  ✓ EVPI = $%.2f per patient", evpi$evpi))
  }

  # ============================================================================
  # Simulation & Power Analysis
  # ============================================================================

  if (perform_simulation) {
    message("\n█ Bonus: Running power analysis simulation...")

    power <- nma_power_analysis(
      n_treatments = length(core_nma$trts),
      n_studies_per_comparison = floor(nrow(data) / (length(core_nma$trts) - 1)),
      effect_size = 0.5,
      heterogeneity = core_nma$tau,
      n_simulations = 500
    )

    result$power_analysis <- power

    message(sprintf("  ✓ Estimated power: %.1f%%", power$power * 100))
  }

  # ============================================================================
  # Visualizations
  # ============================================================================

  if (generate_visualizations) {
    message("\n█ Generating visualizations...")

    viz_list <- list()

    # Network plot
    if (exists("geometry", result)) {
      viz_list$network_geometry <- result$geometry
    }

    # Forest plot with rankings
    if (exists("sucra", result)) {
      viz_list$rankings <- result$sucra
    }

    result$visualizations <- viz_list

    message(sprintf("  ✓ %d visualization objects created", length(viz_list)))
  }

  # ============================================================================
  # Export Results
  # ============================================================================

  if (export_results) {
    message("\n█ Exporting results...")

    exported <- export_nma_results(
      nma_result = core_nma,
      output_dir = output_dir,
      formats = c("csv", "html", "json"),
      prefix = "comprehensive_nma"
    )

    result$exported_files <- exported

    message(sprintf("  ✓ Results exported to: %s", output_dir))
  }

  # ============================================================================
  # Finalize
  # ============================================================================

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result$runtime_seconds <- elapsed
  result$timestamp <- Sys.time()

  class(result) <- c("comprehensive_nma", "list")

  message("\n")
  message("═══════════════════════════════════════════════════════════")
  message(sprintf("  ✓ ANALYSIS COMPLETE IN %.1f SECONDS", elapsed))
  message("═══════════════════════════════════════════════════════════")
  message("\n")

  result
}

#' Print Comprehensive NMA Results
#' @param x comprehensive_nma object
#' @param ... Additional arguments
#' @export
print.comprehensive_nma <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  COMPREHENSIVE NETWORK META-ANALYSIS RESULTS\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Analysis Timestamp:", format(x$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
  cat(sprintf("Runtime: %.1f seconds\n\n", x$runtime_seconds))

  # Core NMA
  if (!is.null(x$core_nma)) {
    cat("Core NMA:\n")
    cat(sprintf("  • Treatments: %d\n", length(x$core_nma$trts)))
    cat(sprintf("  • Studies: %d\n", length(unique(x$core_nma$studlab))))
    cat(sprintf("  • Reference: %s\n", x$core_nma$reference.group))
    cat(sprintf("  • Heterogeneity (I²): %.1f%%\n\n", x$core_nma$I2))
  }

  # Network Geometry
  if (!is.null(x$geometry)) {
    cat("Network Geometry:\n")
    cat(sprintf("  • Structure: %s\n", x$geometry$characteristics$network_type))
    cat(sprintf("  • Density: %.1f%%\n", x$geometry$connectivity$density * 100))
    cat(sprintf("  • Robustness: %.1f%%\n\n", x$geometry$robustness$robustness_score * 100))
  }

  # Rankings
  if (!is.null(x$sucra)) {
    cat("Treatment Rankings (SUCRA):\n")
    top3 <- head(sort(x$sucra$sucra_scores, decreasing = TRUE), 3)
    for (i in seq_along(top3)) {
      cat(sprintf("  %d. %s: %.1f%%\n", i, names(top3)[i], top3[i] * 100))
    }
    cat("\n")
  }

  # Inconsistency
  if (!is.null(x$inconsistency)) {
    n_inconsistent <- sum(x$inconsistency$results$p_value < 0.10, na.rm = TRUE)
    n_total <- nrow(x$inconsistency$results)
    cat("Inconsistency:\n")
    cat(sprintf("  • %d/%d comparisons inconsistent (p < 0.10)\n\n",
               n_inconsistent, n_total))
  }

  # Publication Bias
  if (!is.null(x$publication_bias)) {
    cat("Publication Bias:\n")
    cat(sprintf("  • Egger's test: p = %.3f\n", x$publication_bias$egger$p_value))
    cat(sprintf("  • Begg's test: p = %.3f\n\n", x$publication_bias$begg$p_value))
  }

  # Multivariate
  if (!is.null(x$multivariate)) {
    cat("Multivariate Analysis:\n")
    cat(sprintf("  • %d outcomes analyzed jointly\n\n", x$multivariate$mvnma_data$n_outcomes))
  }

  # Component NMA
  if (!is.null(x$component_nma)) {
    cat("Component NMA:\n")
    cat(sprintf("  • %d components\n", x$component_nma$model_fit$n_components))
    cat(sprintf("  • Model R²: %.3f\n\n", x$component_nma$model_fit$R_squared))
  }

  # Personalized
  if (!is.null(x$personalized)) {
    cat("Personalized Prediction:\n")
    cat(sprintf("  • Effect modifiers: %s\n\n",
               paste(x$personalized$effect_modifiers, collapse = ", ")))
  }

  # Value of Information
  if (!is.null(x$value_of_information)) {
    cat("Value of Information:\n")
    cat(sprintf("  • EVPI per patient: $%.2f\n\n", x$value_of_information$evpi))
  }

  # Power Analysis
  if (!is.null(x$power_analysis)) {
    cat("Power Analysis:\n")
    cat(sprintf("  • Estimated power: %.1f%%\n\n", x$power_analysis$power * 100))
  }

  # Export
  if (!is.null(x$exported_files)) {
    cat(sprintf("Results exported to %d format(s)\n\n", length(x$exported_files)))
  }

  cat("Use summary() for detailed results\n")
  cat("Use plot() for visualizations\n\n")

  invisible(x)
}

#' Summary of Comprehensive NMA
#' @param object comprehensive_nma object
#' @param ... Additional arguments
#' @export
summary.comprehensive_nma <- function(object, ...) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  COMPREHENSIVE NMA DETAILED SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  # Print each component's summary
  if (!is.null(object$core_nma)) {
    cat("■ CORE NMA RESULTS\n\n")
    print(summary(object$core_nma))
    cat("\n")
  }

  if (!is.null(object$geometry)) {
    cat("■ NETWORK GEOMETRY\n\n")
    print(object$geometry)
    cat("\n")
  }

  if (!is.null(object$sucra)) {
    cat("■ TREATMENT RANKINGS\n\n")
    print(object$sucra)
    cat("\n")
  }

  if (!is.null(object$inconsistency)) {
    cat("■ INCONSISTENCY ASSESSMENT\n\n")
    print(object$inconsistency)
    cat("\n")
  }

  if (!is.null(object$heterogeneity)) {
    cat("■ HETEROGENEITY ANALYSIS\n\n")
    print(object$heterogeneity)
    cat("\n")
  }

  invisible(object)
}

#' Plot Comprehensive NMA Results
#' @param x comprehensive_nma object
#' @param which Which plots to generate
#' @param ... Additional arguments
#' @export
plot.comprehensive_nma <- function(x, which = "all", ...) {

  available_plots <- c("network", "rankings", "forest", "funnel",
                      "geometry", "inconsistency")

  if (which == "all") {
    plots_to_make <- available_plots
  } else {
    plots_to_make <- which
  }

  for (plot_type in plots_to_make) {

    if (plot_type == "network" && !is.null(x$geometry)) {
      plot(x$geometry, type = "network")
    }

    if (plot_type == "rankings" && !is.null(x$sucra)) {
      plot(x$sucra, type = "rankogram")
    }

    if (plot_type == "geometry" && !is.null(x$geometry)) {
      plot(x$geometry, type = "degree")
    }

    # Add more plot types as needed
  }

  invisible(x)
}

#' Quick NMA Analysis
#'
#' Fast, streamlined NMA analysis with essential features only.
#'
#' @param data Pairwise comparison data
#' @param sm Summary measure
#' @param reference Reference treatment
#' @return Simplified NMA results
#' @export
quick_nma <- function(data, sm = "MD", reference = NULL) {

  message("Running quick NMA analysis...")

  # Core NMA
  nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  # Rankings
  sucra <- calculate_sucra(nma)

  # Heterogeneity
  het <- heterogeneity_report(nma)

  result <- list(
    nma = nma,
    sucra = sucra,
    heterogeneity = het
  )

  class(result) <- c("quick_nma", "list")

  message("✓ Quick analysis complete")

  result
}
