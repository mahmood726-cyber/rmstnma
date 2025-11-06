#' Check if Analysis Mode is Valid
#'
#' @param mode Character: "standard" or "experimental"
#' @return Logical
#' @keywords internal
is_valid_mode <- function(mode) {
  mode %in% c("standard", "experimental")
}

#' Get Mode Description
#'
#' @param mode Character: "standard" or "experimental"
#' @return Character description
#' @keywords internal
get_mode_description <- function(mode) {
  switch(mode,
    "standard" = "Validated methods suitable for systematic reviews and clinical guidelines",
    "experimental" = "Novel methods requiring validation - research use only",
    "Unknown mode"
  )
}

#' Display Experimental Mode Warning
#'
#' @keywords internal
warn_experimental_mode <- function() {
  warning("
\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
              EXPERIMENTAL MODE ENABLED
\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550

These methods are NOT validated for clinical decision-making.

SUITABLE FOR:
  \u2022 Academic methods research
  \u2022 Exploratory analysis
  \u2022 Hypothesis generation
  \u2022 Sensitivity analyses (complementing standard NMA)

NOT SUITABLE FOR:
  \u2717 Cochrane systematic reviews
  \u2717 Clinical practice guidelines (WHO, NICE, etc.)
  \u2717 Regulatory submissions
  \u2717 Clinical decision-making

VALIDATION STATUS:
  See VALIDATION_PLAN.md for detailed status of each method.
  Methods paper in preparation.

PROCEED WITH CAUTION - Research use only.
\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
  ", immediate. = TRUE, call. = FALSE)
}

#' Run Standard NMA (Validated Methods Only)
#'
#' Executes network meta-analysis using only validated methods based on
#' netmeta and gemtc packages. Suitable for systematic reviews and clinical
#' guidelines.
#'
#' @param data Pairwise comparison data
#' @param config Analysis configuration from setup_powernma()
#' @param sm Summary measure (HR, OR, RR, MD, SMD)
#' @param reference Reference treatment name
#' @param use_bayesian Logical: run Bayesian NMA with gemtc
#' @param use_sensitivity Logical: run sensitivity analyses
#' @param ... Additional arguments
#'
#' @return List of class powernma_result with standard NMA results
#' @keywords internal
run_standard_nma <- function(data,
                             config = NULL,
                             sm = "HR",
                             reference = NULL,
                             use_bayesian = FALSE,
                             use_sensitivity = TRUE,
                             ...) {

  msg("Running STANDARD MODE analysis")
  msg("Using validated methods from netmeta and gemtc packages")

  # Initialize results list
  results <- list()

  # 1. Core network meta-analysis (netmeta wrapper)
  msg("Step 1/5: Running frequentist NMA (netmeta)")
  results$network <- robust_netmeta(
    data = data,
    sm = sm,
    reference = reference,
    ...
  )

  # 2. Network geometry
  msg("Step 2/5: Calculating network geometry")
  results$geometry <- network_geometry(data)

  # 3. Inconsistency assessment
  msg("Step 3/5: Assessing inconsistency")
  results$inconsistency <- tryCatch({
    check_inconsistency(results$network)
  }, error = function(e) {
    warning("Inconsistency assessment failed: ", e$message)
    NULL
  })

  # 4. Sensitivity analyses
  if (use_sensitivity) {
    msg("Step 4/5: Running sensitivity analyses")
    results$sensitivity <- loo_sensitivity_simple(data, sm = sm, reference = reference)
  } else {
    results$sensitivity <- NULL
  }

  # 5. Bayesian NMA (optional, requires gemtc)
  if (use_bayesian) {
    msg("Step 5/5: Running Bayesian NMA (gemtc)")
    results$bayesian <- tryCatch({
      run_bayesian_nma_simple(data, sm = sm, reference = reference)
    }, error = function(e) {
      warning("Bayesian NMA failed: ", e$message, "\nInstall 'gemtc' package for Bayesian analysis.")
      NULL
    })
  } else {
    msg("Step 5/5: Skipping Bayesian NMA (not requested)")
    results$bayesian <- NULL
  }

  # Metadata
  results$mode <- "standard"
  results$validated <- TRUE
  results$methods <- "Standard frequentist network meta-analysis"
  results$validation_status <- "VALIDATED - suitable for publication"
  results$based_on <- sprintf("netmeta %s", utils::packageVersion("netmeta"))

  if (use_bayesian && !is.null(results$bayesian)) {
    results$based_on <- paste0(results$based_on,
                               sprintf(", gemtc %s", utils::packageVersion("gemtc")))
  }

  results$timestamp <- Sys.time()
  results$config <- config

  class(results) <- c("powernma_result", "powernma_standard", "list")
  results
}

#' Run Experimental NMA (Novel Methods - Research Use)
#'
#' Executes network meta-analysis using experimental methods including
#' time-varying approaches (RMST, milestone survival), transportability
#' weighting, and advanced sensitivity analyses.
#'
#' WARNING: These methods require validation. Not suitable for clinical
#' decision-making or systematic reviews.
#'
#' @param data Data (pairwise or IPD)
#' @param data_type Character: "pairwise" or "ipd"
#' @param config Analysis configuration
#' @param sm Summary measure
#' @param reference Reference treatment
#' @param use_timevarying Logical: run RMST and milestone analysis (IPD only)
#' @param use_transportability Logical: run transportability weighting
#' @param use_advanced Logical: run advanced sensitivity analyses
#' @param target_population Data frame with target population characteristics (for transportability)
#' @param rmst_tau Time horizon for RMST (days)
#' @param milestone_times Vector of milestone times (days)
#' @param ... Additional arguments
#'
#' @return List of class powernma_result with experimental results
#' @keywords internal
run_experimental_nma <- function(data,
                                 data_type = c("pairwise", "ipd"),
                                 config = NULL,
                                 sm = "HR",
                                 reference = NULL,
                                 use_timevarying = TRUE,
                                 use_transportability = FALSE,
                                 use_advanced = TRUE,
                                 target_population = NULL,
                                 rmst_tau = 365,
                                 milestone_times = c(90, 180, 365),
                                 ...) {

  data_type <- match.arg(data_type)

  msg("Running EXPERIMENTAL MODE analysis")
  msg("WARNING: Methods not validated for clinical use")

  # Initialize results
  results <- list()
  warnings_list <- character()

  # Branch by data type
  if (data_type == "pairwise") {

    # Standard NMA first
    msg("Step 1/4: Running standard NMA")
    results$network <- robust_netmeta(data, sm = sm, reference = reference, ...)

    # Experimental: Transportability weighting
    if (use_transportability && !is.null(target_population)) {
      msg("Step 2/4: Computing transportability weights (EXPERIMENTAL)")
      results$transportability <- tryCatch({
        weights <- compute_transport_weights(data, target_population)
        diagnostics <- transportability_diagnostics(data, target_population, weights)
        list(weights = weights, diagnostics = diagnostics)
      }, error = function(e) {
        warning("Transportability weighting failed: ", e$message)
        NULL
      })
      warnings_list <- c(warnings_list,
        "Transportability requires positivity assumption - check diagnostics")
    } else {
      results$transportability <- NULL
    }

    # Experimental: Advanced sensitivity analyses
    if (use_advanced) {
      msg("Step 3/4: Running advanced sensitivity analyses (EXPERIMENTAL)")

      # PET-PEESE publication bias
      results$publication_bias <- tryCatch({
        pet_peese_analysis(data)
      }, error = function(e) {
        warning("PET-PEESE failed: ", e$message)
        NULL
      })

      # Leave-one-treatment-out
      results$loto <- tryCatch({
        loto_sensitivity(data, reference = reference, sm = sm)
      }, error = function(e) {
        warning("LOTO sensitivity failed: ", e$message)
        NULL
      })

      # Multiverse robustness
      results$multiverse <- tryCatch({
        run_multiverse(data, reference = reference)
      }, error = function(e) {
        warning("Multiverse analysis failed: ", e$message)
        NULL
      })

      warnings_list <- c(warnings_list,
        "PET-PEESE for NMA requires adaptation from pairwise MA",
        "Multiverse results should be interpreted cautiously")
    }

    msg("Step 4/4: Finalizing results")

  } else {
    # IPD data: Time-varying methods

    msg("Step 1/3: Validating IPD structure")
    validate_ipd(data)

    if (use_timevarying) {
      msg("Step 2/3: Running time-varying NMA (EXPERIMENTAL)")

      # RMST NMA
      results$rmst <- tryCatch({
        rmst_nma(data, tau = rmst_tau, reference = reference)
      }, error = function(e) {
        warning("RMST NMA failed: ", e$message)
        NULL
      })

      # Milestone NMA
      results$milestone <- tryCatch({
        milestone_nma(data, times = milestone_times, reference = reference)
      }, error = function(e) {
        warning("Milestone NMA failed: ", e$message)
        NULL
      })

      warnings_list <- c(warnings_list,
        "RMST sign convention requires mathematical validation (see VALIDATION_PLAN.md)",
        "Milestone analysis uses extend=FALSE by default - check follow-up adequacy",
        "Multi-arm trial handling may be incomplete - verify results",
        "Continuity correction uses standard method - check sparse data handling")
    }

    msg("Step 3/3: Finalizing results")
  }

  # Metadata
  results$mode <- "experimental"
  results$validated <- FALSE
  results$data_type <- data_type
  results$methods <- if (data_type == "ipd") {
    "Experimental time-varying network meta-analysis (RMST, Milestone)"
  } else {
    "Experimental advanced NMA methods (transportability, publication bias)"
  }
  results$validation_status <- "EXPERIMENTAL - research use only"
  results$warnings <- warnings_list
  results$timestamp <- Sys.time()
  results$config <- config

  # Critical warnings per method
  if (!is.null(results$rmst)) {
    results$rmst$warning <- "RMST sign convention not yet mathematically validated"
  }
  if (!is.null(results$milestone)) {
    results$milestone$warning <- "Milestone extrapolation behavior requires careful checking"
  }

  class(results) <- c("powernma_result", "powernma_experimental", "list")
  results
}

#' Check Inconsistency in Network
#'
#' @param network_result Result from netmeta
#' @return List with inconsistency assessment
#' @keywords internal
check_inconsistency <- function(network_result) {
  if (!inherits(network_result, "netmeta")) {
    warning("Input is not a netmeta object")
    return(NULL)
  }

  # Design-by-treatment interaction test
  tryCatch({
    list(
      test = "Design-by-treatment interaction",
      Q = network_result$Q.decomp$Q,
      df = network_result$Q.decomp$df,
      pval = network_result$Q.decomp$pval,
      significant = network_result$Q.decomp$pval < 0.05,
      interpretation = if (network_result$Q.decomp$pval < 0.05) {
        "Significant inconsistency detected (p < 0.05)"
      } else {
        "No significant inconsistency detected"
      }
    )
  }, error = function(e) {
    warning("Could not extract inconsistency statistics: ", e$message)
    NULL
  })
}

#' Validate Mode-Data Type Combination
#'
#' @param mode Character: "standard" or "experimental"
#' @param data_type Character: "pairwise" or "ipd"
#' @return Logical (invisibly), stops if invalid
#' @keywords internal
validate_mode_datatype <- function(mode, data_type) {
  if (mode == "standard" && data_type == "ipd") {
    stop("
IPD-based time-varying methods (RMST, milestone) are EXPERIMENTAL.

These methods have NOT been validated for clinical decision-making.

To use these methods, set mode = 'experimental':
  run_powernma(data, data_type = 'ipd', mode = 'experimental')

For STANDARD validated NMA, provide pairwise data:
  run_powernma(data, data_type = 'pairwise', mode = 'standard')

See VALIDATION_PLAN.md for details on validation status.
    ", call. = FALSE)
  }
  invisible(TRUE)
}
