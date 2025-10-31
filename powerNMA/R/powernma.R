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

    # Export mode-specific tables
    export_powernma_tables_v2(results, config$output_dir)
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

#' Export mode-specific tables (v2)
#' @keywords internal
export_powernma_tables_v2 <- function(results, dir) {
  mode <- results$mode %||% "standard"

  if (mode == "standard") {
    # Export standard NMA results
    if (!is.null(results$network)) {
      # Network estimates
      tryCatch({
        estimates <- data.frame(
          treatment = results$network$trts,
          TE_fixed = results$network$TE.fixed,
          seTE_fixed = results$network$seTE.fixed,
          TE_random = results$network$TE.random,
          seTE_random = results$network$seTE.random
        )
        readr::write_csv(estimates, file.path(dir, "network_estimates.csv"))
      }, error = function(e) NULL)
    }

    # Sensitivity results
    if (!is.null(results$sensitivity)) {
      readr::write_csv(results$sensitivity, file.path(dir, "loo_sensitivity.csv"))
    }

  } else {
    # Export experimental results
    if (!is.null(results$rmst)) {
      readr::write_csv(results$rmst$data, file.path(dir, "rmst_results.csv"))
    }

    if (!is.null(results$milestone)) {
      readr::write_csv(results$milestone$data, file.path(dir, "milestone_results.csv"))
    }

    if (!is.null(results$transportability)) {
      readr::write_csv(
        data.frame(weight = results$transportability$weights),
        file.path(dir, "transportability_weights.csv")
      )
    }

    if (!is.null(results$publication_bias)) {
      readr::write_csv(results$publication_bias,
                      file.path(dir, "publication_bias.csv"))
    }
  }

  msg("Tables exported to: %s", dir)
}
