#' Bayesian Model Selection for Network Meta-Analysis
#'
#' Tools for comparing Bayesian NMA models using information criteria:
#' DIC (Deviance Information Criterion), WAIC (Watanabe-Akaike Information
#' Criterion), and LOOIC (Leave-One-Out Information Criterion).
#'
#' @name bayesian_model_selection
#' @references
#' Spiegelhalter DJ, et al. (2002). Bayesian measures of model complexity and fit.
#' Journal of the Royal Statistical Society B, 64(4):583-639.
#'
#' Vehtari A, et al. (2017). Practical Bayesian model evaluation using leave-one-out
#' cross-validation and WAIC. Statistics and Computing, 27(5):1413-1432.
NULL

#' Calculate DIC for Bayesian NMA Model
#'
#' Deviance Information Criterion penalizes model complexity with effective
#' number of parameters. Lower DIC indicates better model fit.
#'
#' @param bayesian_result Bayesian NMA result object (from gemtc or rjags)
#' @param type Type of Bayesian result: "gemtc", "rjags", "stan"
#' @return List with DIC components
#' @export
#' @references
#' Spiegelhalter DJ, et al. (2002). Bayesian measures of model complexity and fit.
#'
#' @examples
#' \dontrun{
#' # From gemtc
#' network <- gemtc::mtc.network(data.re(smoking))
#' model <- gemtc::mtc.model(network)
#' result <- gemtc::mtc.run(model)
#' dic <- calculate_dic(result, type = "gemtc")
#' print(dic)
#' }
calculate_dic <- function(bayesian_result, type = c("gemtc", "rjags", "stan")) {

  type <- match.arg(type)

  if (type == "gemtc") {
    dic_result <- extract_dic_gemtc(bayesian_result)
  } else if (type == "rjags") {
    dic_result <- extract_dic_rjags(bayesian_result)
  } else if (type == "stan") {
    dic_result <- extract_dic_stan(bayesian_result)
  } else {
    stop("Unsupported Bayesian model type")
  }

  class(dic_result) <- c("dic_result", "list")
  dic_result
}

#' Extract DIC from gemtc Object
#'
#' @param result gemtc result (mtc.result object)
#' @return List with DIC components
#' @keywords internal
extract_dic_gemtc <- function(result) {

  if (!requireNamespace("gemtc", quietly = TRUE)) {
    stop("Package 'gemtc' required for this function")
  }

  if (!inherits(result, "mtc.result")) {
    stop("result must be an mtc.result object from gemtc")
  }

  # Extract DIC from gemtc
  dic_obj <- try(gemtc::dic.samples(result$model, n.iter = 10000), silent = TRUE)

  if (inherits(dic_obj, "try-error")) {
    warning("Failed to calculate DIC: ", dic_obj)
    return(list(
      DIC = NA,
      pD = NA,
      Dbar = NA,
      available = FALSE,
      message = "DIC calculation failed"
    ))
  }

  # Extract components
  deviance <- dic_obj$deviance
  penalty <- dic_obj$penalty
  dic_value <- sum(deviance) + sum(penalty)

  list(
    DIC = dic_value,
    pD = sum(penalty),  # Effective number of parameters
    Dbar = sum(deviance),  # Posterior mean deviance
    deviance_components = deviance,
    penalty_components = penalty,
    available = TRUE,
    interpretation = interpret_dic(dic_value, sum(penalty))
  )
}

#' Extract DIC from rjags Object
#'
#' @param result rjags.samples object
#' @return List with DIC components
#' @keywords internal
extract_dic_rjags <- function(result) {

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' required")
  }

  # Try to extract DIC monitor
  if ("dic" %in% names(result)) {
    dic_summary <- summary(result$dic)

    list(
      DIC = dic_summary$statistics["deviance", "Mean"] + dic_summary$statistics["penalty", "Mean"],
      pD = dic_summary$statistics["penalty", "Mean"],
      Dbar = dic_summary$statistics["deviance", "Mean"],
      available = TRUE,
      interpretation = interpret_dic(
        dic_summary$statistics["deviance", "Mean"] + dic_summary$statistics["penalty", "Mean"],
        dic_summary$statistics["penalty", "Mean"]
      )
    )
  } else {
    list(
      DIC = NA,
      pD = NA,
      Dbar = NA,
      available = FALSE,
      message = "DIC not available in rjags result. Ensure dic=TRUE in jags.samples()"
    )
  }
}

#' Extract DIC from Stan Object
#'
#' @param result Stan fit object
#' @return List with DIC approximation
#' @keywords internal
extract_dic_stan <- function(result) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' required")
  }

  # Stan doesn't directly compute DIC, but we can approximate
  # Extract log-likelihood
  log_lik <- rstan::extract(result, "log_lik")$log_lik

  if (is.null(log_lik)) {
    return(list(
      DIC = NA,
      available = FALSE,
      message = "log_lik not found in Stan model. Add generated quantities block."
    ))
  }

  # Calculate components
  dev <- -2 * log_lik
  Dbar <- mean(rowSums(dev))
  D_theta_bar <- -2 * sum(colMeans(log_lik))
  pD <- Dbar - D_theta_bar
  DIC <- Dbar + pD

  list(
    DIC = DIC,
    pD = pD,
    Dbar = Dbar,
    available = TRUE,
    note = "DIC approximated from Stan log-likelihood",
    interpretation = interpret_dic(DIC, pD)
  )
}

#' Calculate WAIC for Bayesian NMA
#'
#' Watanabe-Akaike Information Criterion is fully Bayesian and asymptotically
#' equivalent to Bayesian cross-validation. Preferred over DIC for hierarchical models.
#'
#' @param bayesian_result Bayesian NMA result
#' @param type Type of Bayesian result
#' @return List with WAIC components
#' @export
#' @references
#' Watanabe S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' Journal of Machine Learning Research, 11:3571-3594.
#'
#' @examples
#' \dontrun{
#' waic <- calculate_waic(bayesian_result, type = "gemtc")
#' print(waic)
#' }
calculate_waic <- function(bayesian_result, type = c("gemtc", "rjags", "stan")) {

  type <- match.arg(type)

  # Extract log-likelihood samples
  log_lik <- extract_log_likelihood(bayesian_result, type)

  if (is.null(log_lik)) {
    return(list(
      WAIC = NA,
      pWAIC = NA,
      lppd = NA,
      available = FALSE,
      message = "Could not extract log-likelihood for WAIC calculation"
    ))
  }

  # Calculate WAIC components
  # log_lik should be: [n_iterations x n_observations]

  # Computed log pointwise predictive density
  lppd <- sum(log(colMeans(exp(log_lik))))

  # Effective number of parameters (method 1: variance)
  pWAIC1 <- sum(apply(log_lik, 2, var))

  # Effective number of parameters (method 2: more stable)
  pWAIC2 <- sum(log(colMeans(exp(log_lik))) - colMeans(log_lik))

  # WAIC
  WAIC1 <- -2 * (lppd - pWAIC1)
  WAIC2 <- -2 * (lppd - pWAIC2)

  # Standard error of WAIC
  waic_vec <- -2 * (log(colMeans(exp(log_lik))) - apply(log_lik, 2, var))
  se_WAIC <- sqrt(ncol(log_lik) * var(waic_vec))

  list(
    WAIC1 = WAIC1,
    WAIC2 = WAIC2,  # Preferred version
    WAIC = WAIC2,   # Default to method 2
    pWAIC1 = pWAIC1,
    pWAIC2 = pWAIC2,
    lppd = lppd,
    se_WAIC = se_WAIC,
    available = TRUE,
    interpretation = interpret_waic(WAIC2, pWAIC2)
  )
}

#' Calculate LOOIC using Leave-One-Out Cross-Validation
#'
#' LOOIC (LOO-IC) using Pareto Smoothed Importance Sampling (PSIS-LOO).
#' Most reliable for model comparison in Bayesian hierarchical models.
#'
#' @param bayesian_result Bayesian NMA result
#' @param type Type of Bayesian result
#' @return List with LOOIC components
#' @export
#' @references
#' Vehtari A, et al. (2017). Practical Bayesian model evaluation using
#' leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5):1413-1432.
#'
#' @examples
#' \dontrun{
#' looic <- calculate_looic(bayesian_result, type = "stan")
#' print(looic)
#' }
calculate_looic <- function(bayesian_result, type = c("gemtc", "rjags", "stan")) {

  type <- match.arg(type)

  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' required for LOOIC calculation. Install with: install.packages('loo')")
  }

  # Extract log-likelihood
  log_lik <- extract_log_likelihood(bayesian_result, type)

  if (is.null(log_lik)) {
    return(list(
      LOOIC = NA,
      pLOO = NA,
      available = FALSE,
      message = "Could not extract log-likelihood for LOOIC calculation"
    ))
  }

  # Calculate PSIS-LOO
  loo_result <- tryCatch({
    loo::loo(log_lik)
  }, error = function(e) {
    warning("LOOIC calculation failed: ", e$message)
    NULL
  })

  if (is.null(loo_result)) {
    return(list(
      LOOIC = NA,
      pLOO = NA,
      available = FALSE,
      message = "LOOIC calculation failed"
    ))
  }

  # Extract components
  list(
    LOOIC = loo_result$estimates["looic", "Estimate"],
    pLOO = loo_result$estimates["p_loo", "Estimate"],
    se_LOOIC = loo_result$estimates["looic", "SE"],
    elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
    n_eff = sum(loo_result$diagnostics$n_eff),
    pareto_k = loo_result$diagnostics$pareto_k,
    available = TRUE,
    full_result = loo_result,
    interpretation = interpret_looic(loo_result)
  )
}

#' Compare Multiple Bayesian Models
#'
#' Compare multiple Bayesian NMA models using DIC, WAIC, and LOOIC.
#'
#' @param models Named list of Bayesian NMA results
#' @param type Type of Bayesian results
#' @param criteria Which criteria to calculate: "all", "DIC", "WAIC", "LOOIC"
#' @return Data frame with model comparison
#' @export
#' @examples
#' \dontrun{
#' models <- list(
#'   "Fixed Effects" = fe_model,
#'   "Random Effects" = re_model,
#'   "Inconsistency" = inconsistency_model
#' )
#' comparison <- compare_bayesian_models(models, type = "gemtc")
#' print(comparison)
#' }
compare_bayesian_models <- function(models,
                                   type = c("gemtc", "rjags", "stan"),
                                   criteria = c("all", "DIC", "WAIC", "LOOIC")) {

  type <- match.arg(type)
  criteria <- match.arg(criteria)

  if (!is.list(models) || length(models) < 2) {
    stop("models must be a list with at least 2 models")
  }

  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- paste0("Model", seq_along(models))
  }

  results <- list()

  for (i in seq_along(models)) {
    model <- models[[i]]
    name <- model_names[i]

    result_row <- list(Model = name)

    # DIC
    if (criteria %in% c("all", "DIC")) {
      dic <- tryCatch({
        calculate_dic(model, type)
      }, error = function(e) {
        list(DIC = NA, pD = NA, available = FALSE)
      })

      result_row$DIC <- dic$DIC
      result_row$pD <- dic$pD
    }

    # WAIC
    if (criteria %in% c("all", "WAIC")) {
      waic <- tryCatch({
        calculate_waic(model, type)
      }, error = function(e) {
        list(WAIC = NA, pWAIC = NA, available = FALSE)
      })

      result_row$WAIC <- waic$WAIC
      result_row$pWAIC <- waic$pWAIC2
    }

    # LOOIC
    if (criteria %in% c("all", "LOOIC")) {
      looic <- tryCatch({
        calculate_looic(model, type)
      }, error = function(e) {
        list(LOOIC = NA, pLOO = NA, available = FALSE)
      })

      result_row$LOOIC <- looic$LOOIC
      result_row$pLOO <- looic$pLOO
    }

    results[[i]] <- as.data.frame(result_row, stringsAsFactors = FALSE)
  }

  comparison_df <- do.call(rbind, results)

  # Add rankings and differences
  if ("DIC" %in% names(comparison_df)) {
    comparison_df$DIC_rank <- rank(comparison_df$DIC, na.last = "keep")
    comparison_df$DIC_diff <- comparison_df$DIC - min(comparison_df$DIC, na.rm = TRUE)
  }

  if ("WAIC" %in% names(comparison_df)) {
    comparison_df$WAIC_rank <- rank(comparison_df$WAIC, na.last = "keep")
    comparison_df$WAIC_diff <- comparison_df$WAIC - min(comparison_df$WAIC, na.rm = TRUE)
  }

  if ("LOOIC" %in% names(comparison_df)) {
    comparison_df$LOOIC_rank <- rank(comparison_df$LOOIC, na.last = "keep")
    comparison_df$LOOIC_diff <- comparison_df$LOOIC - min(comparison_df$LOOIC, na.rm = TRUE)
  }

  class(comparison_df) <- c("model_comparison", "data.frame")
  comparison_df
}

#' Extract Log-Likelihood from Bayesian Result
#'
#' @param result Bayesian result object
#' @param type Type of result
#' @return Matrix of log-likelihood values
#' @keywords internal
extract_log_likelihood <- function(result, type) {

  if (type == "gemtc") {
    # gemtc doesn't directly expose log-likelihood
    # Would need to extract from JAGS model
    warning("Log-likelihood extraction from gemtc not directly supported")
    return(NULL)

  } else if (type == "rjags") {
    # Extract from rjags monitor
    if ("log_lik" %in% names(result)) {
      return(result$log_lik)
    } else {
      warning("log_lik not found in rjags result")
      return(NULL)
    }

  } else if (type == "stan") {
    if (requireNamespace("rstan", quietly = TRUE)) {
      log_lik <- rstan::extract(result, "log_lik")$log_lik
      return(log_lik)
    }
  }

  NULL
}

#' Interpret DIC Value
#'
#' @param dic DIC value
#' @param pD Effective number of parameters
#' @return Character interpretation
#' @keywords internal
interpret_dic <- function(dic, pD) {
  sprintf(
    "DIC = %.1f (pD = %.1f effective parameters). Lower values indicate better fit. Compare models by DIC differences > 5.",
    dic, pD
  )
}

#' Interpret WAIC Value
#'
#' @param waic WAIC value
#' @param pWAIC Effective number of parameters
#' @return Character interpretation
#' @keywords internal
interpret_waic <- function(waic, pWAIC) {
  sprintf(
    "WAIC = %.1f (pWAIC = %.1f). Fully Bayesian criterion. Differences > 2 SE are meaningful.",
    waic, pWAIC
  )
}

#' Interpret LOOIC Result
#'
#' @param loo_result loo object from loo package
#' @return Character interpretation
#' @keywords internal
interpret_looic <- function(loo_result) {

  # Check Pareto k diagnostics
  k_threshold <- 0.7
  bad_k <- sum(loo_result$diagnostics$pareto_k > k_threshold, na.rm = TRUE)

  interpretation <- sprintf(
    "LOOIC = %.1f. Most reliable criterion for hierarchical models.",
    loo_result$estimates["looic", "Estimate"]
  )

  if (bad_k > 0) {
    interpretation <- paste0(
      interpretation,
      sprintf("\nWarning: %d observations with Pareto k > %.1f (unreliable estimates).",
             bad_k, k_threshold)
    )
  }

  interpretation
}

#' Print Model Comparison
#'
#' @param x model_comparison object
#' @param ... Additional arguments
#' @export
print.model_comparison <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Bayesian Model Comparison\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Information Criteria (lower is better):\n\n")
  print(as.data.frame(x), row.names = FALSE)

  cat("\n")
  cat("Interpretation Guidelines:\n")
  cat("  • DIC/WAIC/LOOIC differences < 2: Models essentially equivalent\n")
  cat("  • Differences 2-5: Weak preference for lower value\n")
  cat("  • Differences 5-10: Strong preference for lower value\n")
  cat("  • Differences > 10: Decisive preference for lower value\n\n")

  cat("Preferred Criteria:\n")
  cat("  • LOOIC: Most reliable for hierarchical models (gold standard)\n")
  cat("  • WAIC: Fully Bayesian, good for most models\n")
  cat("  • DIC: Traditional, may be unstable for complex models\n\n")

  # Identify best model by each criterion
  if ("DIC_rank" %in% names(x)) {
    best_dic <- x$Model[which.min(x$DIC)]
    cat(sprintf("Best model by DIC: %s\n", best_dic))
  }

  if ("WAIC_rank" %in% names(x)) {
    best_waic <- x$Model[which.min(x$WAIC)]
    cat(sprintf("Best model by WAIC: %s\n", best_waic))
  }

  if ("LOOIC_rank" %in% names(x)) {
    best_looic <- x$Model[which.min(x$LOOIC)]
    cat(sprintf("Best model by LOOIC: %s\n", best_looic))
  }

  cat("\n")

  invisible(x)
}

#' Print DIC Result
#'
#' @param x dic_result object
#' @param ... Additional arguments
#' @export
print.dic_result <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Deviance Information Criterion (DIC)\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  if (!x$available) {
    cat("DIC not available\n")
    if (!is.null(x$message)) {
      cat(sprintf("Reason: %s\n", x$message))
    }
    return(invisible(x))
  }

  cat(sprintf("DIC:   %.2f\n", x$DIC))
  cat(sprintf("pD:    %.2f  (effective number of parameters)\n", x$pD))
  cat(sprintf("D̄:     %.2f  (posterior mean deviance)\n\n", x$Dbar))

  if (!is.null(x$interpretation)) {
    cat("Interpretation:\n")
    cat(strwrap(x$interpretation, width = 60, prefix = "  "), sep = "\n")
  }

  cat("\n")

  invisible(x)
}
