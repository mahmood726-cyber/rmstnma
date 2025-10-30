#' Run Comprehensive powerNMA Analysis
#'
#' Execute a complete network meta-analysis pipeline combining standard NMA,
#' time-varying methods, Bayesian inference, and advanced diagnostics.
#'
#' @param data Either pairwise NMA data or individual patient data (IPD)
#' @param data_type Type of input data: "pairwise" or "ipd"
#' @param ref_treatment Reference treatment (auto-detected if NULL)
#' @param target_population Named list of covariate values for transportability
#' @param config Configuration object from setup_powernma()
#' @return powernma_result object with comprehensive results
#' @export
#' @examples
#' \dontrun{
#' # Standard NMA with simulated data
#' data <- simulate_nma_data(n_studies = 30)
#' results <- run_powernma(data, data_type = "pairwise")
#'
#' # Time-varying NMA with IPD
#' ipd <- generate_example_ipd(n_trials = 10)
#' results <- run_powernma(ipd, data_type = "ipd",
#'                         config = setup_powernma(use_timevarying = TRUE))
#'
#' # View results
#' print(results)
#' summary(results)
#' plot(results)
#' }
run_powernma <- function(data,
                        data_type = c("pairwise", "ipd"),
                        ref_treatment = NULL,
                        target_population = NULL,
                        config = setup_powernma()) {

  data_type <- match.arg(data_type)

  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  powerNMA: Comprehensive Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════\n\n")

  set.seed(config$seed)

  # Initialize results container
  results <- list()

  # ========== DATA PREPARATION ==========

  if (data_type == "pairwise") {
    data <- clean_nma_data(data)
    validate_nma_input(data)
    msg("Validated pairwise data: %d comparisons from %d studies",
        nrow(data), length(unique(data$studlab)))
  } else {
    validate_ipd(data)
    msg("Validated IPD: %d patients from %d trials",
        nrow(data), length(unique(data$trial)))
  }

  # Auto-detect reference treatment
  if (is.null(ref_treatment)) {
    if (data_type == "pairwise") {
      counts <- table(c(data$treat1, data$treat2))
    } else {
      counts <- table(data$treatment)
    }
    ref_treatment <- names(sort(counts, decreasing = TRUE))[1]
    msg("Reference treatment auto-detected: '%s'", ref_treatment)
  }

  # ========== STANDARD NMA (if pairwise data) ==========

  if (data_type == "pairwise") {
    msg("Running standard network meta-analysis...")

    results$main_nma <- robust_netmeta(
      data = data,
      ref_treatment = ref_treatment,
      sm = config$sm,
      random = TRUE
    )

    results$main_nma$sm <- config$sm

    cat(sprintf("\nHeterogeneity τ: %.4f\n", results$main_nma$tau))
    cat(sprintf("I²: %.1f%%\n", results$main_nma$I2.random * 100))

    # Network geometry
    results$geometry <- network_geometry(data)
    msg("Network: %d treatments, %d comparisons (density: %.2f)",
        results$geometry$n_treatments,
        results$geometry$n_comparisons,
        results$geometry$density)

    # Inconsistency checks
    msg("Checking for inconsistency...")
    results$global_inconsistency <- .safe_try(
      netmeta::decomp.design(results$main_nma),
      "Global inconsistency test"
    )

    if (results$main_nma$k >= 10) {
      results$local_inconsistency <- .safe_try(
        netmeta::netsplit(results$main_nma),
        "Node-splitting"
      )
    } else {
      msg("Skipping node-splitting (k < 10)")
    }
  }

  # ========== TIME-VARYING ANALYSES (if IPD) ==========

  if (data_type == "ipd" && config$use_timevarying) {
    msg("Running time-varying analyses...")

    # RMST NMA
    if (length(config$tau_list) > 0) {
      msg("Computing RMST NMA at tau = %s",
          paste(config$tau_list, collapse = ", "))

      results$rmst_nma <- rmst_nma(
        ipd = data,
        tau_list = config$tau_list,
        reference = ref_treatment
      )
    }

    # Milestone NMA
    if (length(config$milestone_times) > 0) {
      msg("Computing Milestone NMA at times = %s",
          paste(config$milestone_times, collapse = ", "))

      results$milestone_nma <- milestone_nma(
        ipd = data,
        times = config$milestone_times,
        reference = ref_treatment
      )
    }
  }

  # ========== BAYESIAN NMA ==========

  if (config$use_bayesian && data_type == "pairwise") {
    if (has_pkg("gemtc") && has_pkg("coda")) {
      msg("Running Bayesian NMA...")
      results$bayesian <- .safe_try(
        run_bayesian_nma_simple(data, ref_treatment, config$sm),
        "Bayesian NMA"
      )
    } else {
      msg("Skipping Bayesian NMA (gemtc/coda not available)")
    }
  }

  # ========== SENSITIVITY ANALYSES ==========

  if (config$run_sensitivity && data_type == "pairwise") {
    msg("Running sensitivity analyses...")

    # Leave-one-out
    msg("Leave-one-study-out...")
    results$loo <- loo_sensitivity_simple(data, ref_treatment, config$sm)
  }

  # ========== RANKING ==========

  if (!is.null(results$main_nma)) {
    results$ranking <- .safe_try(
      netmeta::netrank(results$main_nma, small.values = "good"),
      "Treatment ranking"
    )
  }

  # ========== EXPORT RESULTS ==========

  if (config$export_results) {
    if (!dir.exists(config$output_dir)) {
      dir.create(config$output_dir, recursive = TRUE)
    }

    output_file <- file.path(config$output_dir, "powernma_results.rds")
    saveRDS(results, output_file)
    msg("Results saved to: %s", output_file)

    # Export key tables
    export_powernma_tables(results, config$output_dir)
  }

  # ========== CREATE SUMMARY ==========

  results$summary <- create_powernma_summary(results, ref_treatment, config$sm)

  # ========== FINALIZE ==========

  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  Analysis Complete!\n")
  cat("═══════════════════════════════════════════════════\n\n")

  output <- list(
    results = results,
    data = data,
    ref_treatment = ref_treatment,
    config = config
  )

  class(output) <- c("powernma_result", "list")
  output
}

#' Simple Bayesian NMA wrapper
#' @keywords internal
run_bayesian_nma_simple <- function(data, ref_treatment, sm) {
  if (!has_pkg("gemtc")) return(NULL)

  # Prepare data
  net_data <- data %>%
    dplyr::select(studlab, treat1, treat2, TE, seTE) %>%
    dplyr::mutate(
      studlab = as.character(studlab),
      treat1 = gsub("[^A-Za-z0-9_]", "_", as.character(treat1)),
      treat2 = gsub("[^A-Za-z0-9_]", "_", as.character(treat2))
    )

  net <- gemtc::mtc.network(data.re = net_data)
  model <- gemtc::mtc.model(net, type = "consistency",
                            likelihood = "normal", link = "identity")
  results <- gemtc::mtc.run(model, n.adapt = 5000, n.iter = 10000,
                           n.chains = 3)

  list(network = net, model = model, results = results)
}

#' Simple LOO sensitivity
#' @keywords internal
loo_sensitivity_simple <- function(data, ref_treatment, sm) {
  studs <- unique(data$studlab)

  if (length(studs) < 5) {
    msg("Too few studies for LOO")
    return(NULL)
  }

  tau_loo <- vapply(studs, function(s) {
    dd <- data %>% dplyr::filter(studlab != s)

    fit <- tryCatch(
      robust_netmeta(dd, ref_treatment, sm = sm, random = TRUE),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NA_real_)
    fit$tau
  }, numeric(1))

  tibble::tibble(
    studlab = studs,
    tau = tau_loo
  )
}

#' Create summary statistics
#' @keywords internal
create_powernma_summary <- function(results, ref_treatment, sm) {
  summary <- list(
    ref_treatment = ref_treatment,
    sm = sm
  )

  if (!is.null(results$main_nma)) {
    summary$tau <- results$main_nma$tau
    summary$I2 <- results$main_nma$I2.random

    # Extract P-scores
    if (!is.null(results$ranking)) {
      ps <- if (!is.null(results$ranking$Pscore)) {
        tibble::tibble(
          treatment = names(results$ranking$Pscore),
          Pscore = as.numeric(results$ranking$Pscore)
        )
      } else if (!is.null(results$ranking$random$Pscore)) {
        tibble::tibble(
          treatment = names(results$ranking$random$Pscore),
          Pscore = as.numeric(results$ranking$random$Pscore)
        )
      } else NULL

      if (!is.null(ps)) {
        summary$top5_pscore <- ps %>%
          dplyr::arrange(dplyr::desc(Pscore)) %>%
          head(5)
      }
    }
  }

  summary
}

#' Export key tables
#' @keywords internal
export_powernma_tables <- function(results, dir) {
  # P-scores
  if (!is.null(results$summary$top5_pscore)) {
    readr::write_csv(
      results$summary$top5_pscore,
      file.path(dir, "pscores.csv")
    )
  }

  # LOO results
  if (!is.null(results$loo)) {
    readr::write_csv(
      results$loo,
      file.path(dir, "loo_sensitivity.csv")
    )
  }

  msg("Tables exported to: %s", dir)
}
