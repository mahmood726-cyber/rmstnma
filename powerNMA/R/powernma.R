#' Run Comprehensive powerNMA Analysis
#'
#' Execute a complete network meta-analysis pipeline. Supports two modes:
#'
#' **STANDARD MODE** (default): Validated methods suitable for systematic reviews
#' and clinical guidelines. Uses netmeta and gemtc packages. Cochrane-approved.
#'
#' **EXPERIMENTAL MODE**: Novel methods requiring validation (time-varying NMA,
#' transportability). Research use only. NOT suitable for clinical decision-making.
#'
#' @param data Either pairwise NMA data or individual patient data (IPD)
#' @param data_type Type of input data: "pairwise" or "ipd"
#' @param mode Analysis mode: "standard" (default, validated) or "experimental" (research)
#' @param ref_treatment Reference treatment (auto-detected if NULL)
#' @param target_population Named list of covariate values for transportability (experimental mode)
#' @param config Configuration object from setup_powernma()
#' @return powernma_result object with comprehensive results
#' @export
#' @examples
#' \dontrun{
#' # STANDARD MODE: Validated NMA for systematic reviews
#' data <- simulate_nma_data(n_studies = 30)
#' results <- run_powernma(data, data_type = "pairwise", mode = "standard")
#'
#' # EXPERIMENTAL MODE: Time-varying NMA (research use)
#' ipd <- generate_example_ipd(n_trials = 10)
#' results <- run_powernma(ipd, data_type = "ipd", mode = "experimental")
#'
#' # View results
#' print(results)
#' summary(results)
#' }
run_powernma <- function(data,
                        data_type = c("pairwise", "ipd"),
                        mode = c("standard", "experimental"),
                        ref_treatment = NULL,
                        target_population = NULL,
                        config = setup_powernma()) {

  data_type <- match.arg(data_type)
  mode <- match.arg(mode)

  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  powerNMA: Comprehensive Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("Mode:", toupper(mode), "\n")
  cat("═══════════════════════════════════════════════════\n\n")

  # Validate mode-datatype combination
  validate_mode_datatype(mode, data_type)

  # Show experimental warning if needed
  if (mode == "experimental") {
    warn_experimental_mode()
  }

  set.seed(config$seed)

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

  # ========== ROUTE TO APPROPRIATE MODE ==========

  results <- if (mode == "standard") {
    # STANDARD MODE: Validated methods only
    run_standard_nma(
      data = data,
      config = config,
      sm = config$sm,
      reference = ref_treatment,
      use_bayesian = config$use_bayesian,
      use_sensitivity = config$run_sensitivity
    )
  } else {
    # EXPERIMENTAL MODE: Novel methods with warnings
    run_experimental_nma(
      data = data,
      data_type = data_type,
      config = config,
      sm = config$sm,
      reference = ref_treatment,
      use_timevarying = config$use_timevarying,
      use_transportability = !is.null(target_population),
      use_advanced = TRUE,
      target_population = target_population,
      rmst_tau = if (length(config$tau_list) > 0) config$tau_list[1] else 365,
      milestone_times = if (length(config$milestone_times) > 0) config$milestone_times else c(90, 180, 365)
    )
  }

  # ========== EXPORT RESULTS (if requested) ==========

  if (config$export_results) {
    if (!dir.exists(config$output_dir)) {
      dir.create(config$output_dir, recursive = TRUE)
    }

    output_file <- file.path(config$output_dir, "powernma_results.rds")
    saveRDS(results, output_file)
    msg("Results saved to: %s", output_file)

    # Export tables using unified function
    export_powernma_tables_unified(results, config$output_dir, verbose = TRUE)
  }

  # ========== FINALIZE ==========

  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  Analysis Complete!\n")
  cat("═══════════════════════════════════════════════════\n\n")

  # Add metadata
  results$data <- data
  results$ref_treatment <- ref_treatment
  results$config <- config

  results
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
  results <- gemtc::mtc.run(model,
                           n.adapt = .DEFAULT_N_ADAPT,
                           n.iter = .DEFAULT_N_ITER,
                           n.chains = .DEFAULT_N_CHAINS)

  list(network = net, model = model, results = results)
}

#' Simple LOO sensitivity
#' @keywords internal
loo_sensitivity_simple <- function(data, ref_treatment, sm) {
  studs <- unique(data$studlab)

  if (length(studs) < .MIN_STUDIES_FOR_LOO) {
    msg("Too few studies for LOO (need at least %d)", .MIN_STUDIES_FOR_LOO)
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

#' Export Results Tables (Unified)
#'
#' Consolidated export function handling both standard and experimental modes.
#' Replaces the old export_powernma_tables() and export_powernma_tables_v2()
#' functions with a single, more maintainable implementation.
#'
#' @param results powerNMA results object
#' @param dir Output directory path
#' @param verbose Print export messages
#' @keywords internal
export_powernma_tables_unified <- function(results, dir, verbose = TRUE) {
  # Validate inputs
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  mode <- results$mode %||% "standard"
  exported_files <- character()

  # Helper to safely export a table
  safe_export <- function(data, filename, description) {
    if (!is.null(data)) {
      filepath <- file.path(dir, filename)
      tryCatch({
        readr::write_csv(data, filepath)
        exported_files <<- c(exported_files, description)
        if (verbose) msg("Exported: %s", description)
      }, error = function(e) {
        warning(sprintf("Failed to export %s: %s", description, e$message))
      })
    }
  }

  # Common exports (both modes)
  if (!is.null(results$summary$top5_pscore)) {
    safe_export(results$summary$top5_pscore, "top_treatments.csv", "Top 5 treatments (P-scores)")
  }

  # Mode-specific exports
  if (mode == "standard") {
    # Standard NMA exports
    if (!is.null(results$network)) {
      estimates <- tryCatch({
        data.frame(
          treatment = results$network$trts,
          TE_fixed = results$network$TE.fixed,
          seTE_fixed = results$network$seTE.fixed,
          TE_random = results$network$TE.random,
          seTE_random = results$network$seTE.random,
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)

      if (!is.null(estimates)) {
        safe_export(estimates, "network_estimates.csv", "Network effect estimates")
      }
    }

    if (!is.null(results$sensitivity)) {
      safe_export(results$sensitivity, "loo_sensitivity.csv", "Leave-one-out sensitivity")
    }

    if (!is.null(results$geometry)) {
      safe_export(results$geometry, "network_geometry.csv", "Network geometry metrics")
    }

  } else {
    # Experimental mode exports
    if (!is.null(results$rmst)) {
      safe_export(results$rmst$data, "rmst_results.csv", "RMST NMA results")
    }

    if (!is.null(results$milestone)) {
      safe_export(results$milestone$data, "milestone_results.csv", "Milestone survival results")
    }

    if (!is.null(results$transportability)) {
      weights_df <- data.frame(
        study = names(results$transportability$weights),
        weight = as.numeric(results$transportability$weights),
        stringsAsFactors = FALSE
      )
      safe_export(weights_df, "transportability_weights.csv", "Transportability weights")
    }

    if (!is.null(results$publication_bias)) {
      safe_export(results$publication_bias, "publication_bias.csv", "Publication bias assessment")
    }

    if (!is.null(results$loto)) {
      safe_export(results$loto, "loto_sensitivity.csv", "Leave-one-treatment-out sensitivity")
    }
  }

  # Summary message
  if (verbose && length(exported_files) > 0) {
    msg("Successfully exported %d table(s) to: %s", length(exported_files), dir)
  }

  invisible(exported_files)
}

#' @keywords internal
#' @deprecated Use export_powernma_tables_unified() instead
export_powernma_tables <- function(results, dir) {
  .Deprecated("export_powernma_tables_unified")
  export_powernma_tables_unified(results, dir, verbose = TRUE)
}

#' @keywords internal
#' @deprecated Use export_powernma_tables_unified() instead
export_powernma_tables_v2 <- function(results, dir) {
  .Deprecated("export_powernma_tables_unified")
  export_powernma_tables_unified(results, dir, verbose = TRUE)
}
