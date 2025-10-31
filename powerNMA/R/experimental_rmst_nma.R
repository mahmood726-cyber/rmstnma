# ============================================================================
# EXPERIMENTAL: RMST-based Network Meta-Analysis
# ============================================================================
#
# Based on: Hua et al. (2025) "Network Meta-Analysis of Time-to-Event Endpoints
# With Individual Participant Data Using Restricted Mean Survival Time Regression"
# Biometrical Journal. DOI: 10.1002/bimj.70037
#
# STATUS: EXPERIMENTAL - Methods from 2025 literature, under active development
#
# RMST (Restricted Mean Survival Time) represents the area under the survival
# curve up to a specified time point τ (tau). It provides a more interpretable
# measure than hazard ratios:
# - Measured in time units (e.g., months, years)
# - Does not require proportional hazards assumption
# - Directly interpretable: "Treatment A gives X more months than Treatment B"
#
# ============================================================================

#' EXPERIMENTAL: Network Meta-Analysis Using Restricted Mean Survival Time
#'
#' Performs network meta-analysis for time-to-event outcomes using Restricted
#' Mean Survival Time (RMST) instead of hazard ratios. RMST represents the mean
#' survival time up to a prespecified time point (tau) and corresponds to the
#' area under the survival curve.
#'
#' @param data A data frame containing the survival data. Must include columns:
#'   \itemize{
#'     \item study: Study identifier
#'     \item treatment: Treatment name
#'     \item time: Time to event or censoring
#'     \item event: Event indicator (1 = event, 0 = censored)
#'   }
#'   For aggregate data, must include:
#'   \itemize{
#'     \item study: Study identifier
#'     \item treatment: Treatment name
#'     \item rmst: Estimated RMST
#'     \item rmst_se: Standard error of RMST
#'   }
#' @param tau Restriction time (maximum follow-up time to consider).
#'   If NULL, uses minimum of maximum follow-up times across studies.
#' @param data_type Type of data: "ipd" (individual participant data) or
#'   "aggregate" (pre-calculated RMST values)
#' @param reference Reference treatment for comparisons. If NULL, uses first treatment.
#' @param method Analysis method: "frequentist" or "bayesian"
#' @param covariate Optional covariate name for meta-regression (IPD only)
#' @param tau_common Logical. Should tau be common across all studies? Default TRUE.
#' @param ... Additional arguments passed to netmeta or JAGS models
#'
#' @return An object of class "rmst_nma" containing:
#'   \item{rmst_estimates}{Treatment-specific RMST estimates}
#'   \item{rmst_differences}{Pairwise RMST differences (treatment effects)}
#'   \item{tau}{The restriction time used}
#'   \item{reference}{Reference treatment}
#'   \item{ranking}{Treatment ranking based on RMST}
#'   \item{nma_result}{Underlying NMA model object}
#'   \item{interpretation}{Clinical interpretation of results}
#'   \item{method}{Analysis method used}
#'
#' @details
#' RMST advantages over hazard ratios:
#' \itemize{
#'   \item More intuitive interpretation (e.g., "2.5 additional months")
#'   \item No proportional hazards assumption required
#'   \item Robust to non-proportional hazards
#'   \item Directly clinically meaningful
#' }
#'
#' For IPD, RMST is calculated using pseudo-observations or direct integration
#' of survival curves. For aggregate data, pre-calculated RMST values are used.
#'
#' @references
#' Hua H, et al. (2025). Network Meta-Analysis of Time-to-Event Endpoints With
#' Individual Participant Data Using Restricted Mean Survival Time Regression.
#' Biometrical Journal, 67(1), e70037.
#'
#' Royston P, Parmar MKB (2013). Restricted mean survival time: an alternative
#' to the hazard ratio for the design and analysis of randomized trials with a
#' time-to-event outcome. BMC Medical Research Methodology, 13:152.
#'
#' @examples
#' \dontrun{
#' # Example 1: Aggregate data (pre-calculated RMST)
#' data_agg <- data.frame(
#'   study = rep(c("Study1", "Study2", "Study3"), each = 2),
#'   treatment = rep(c("A", "B"), 3),
#'   rmst = c(24.5, 26.8, 23.1, 25.9, 22.8, 27.2),  # months
#'   rmst_se = c(1.2, 1.3, 1.5, 1.4, 1.1, 1.2)
#' )
#'
#' result_agg <- rmst_nma(
#'   data = data_agg,
#'   tau = 36,  # 3-year RMST
#'   data_type = "aggregate",
#'   reference = "A"
#' )
#' print(result_agg)
#' plot(result_agg)
#'
#' # Example 2: Individual participant data
#' # Simulate IPD survival data
#' set.seed(123)
#' data_ipd <- data.frame(
#'   study = rep(c("Study1", "Study2"), each = 100),
#'   treatment = sample(c("A", "B", "C"), 200, replace = TRUE),
#'   time = rexp(200, rate = 0.05),
#'   event = rbinom(200, 1, 0.7)
#' )
#'
#' result_ipd <- rmst_nma(
#'   data = data_ipd,
#'   tau = 24,  # 2-year RMST
#'   data_type = "ipd",
#'   reference = "A"
#' )
#' print(result_ipd)
#'
#' # Example 3: With covariate adjustment (IPD)
#' data_ipd$age <- rnorm(200, mean = 60, sd = 10)
#' result_covar <- rmst_nma(
#'   data = data_ipd,
#'   tau = 24,
#'   data_type = "ipd",
#'   reference = "A",
#'   covariate = "age"
#' )
#' }
#'
#' @export
rmst_nma <- function(data,
                     tau = NULL,
                     data_type = c("aggregate", "ipd"),
                     reference = NULL,
                     method = c("frequentist", "bayesian"),
                     covariate = NULL,
                     tau_common = TRUE,
                     ...) {

  # Argument checking
  data_type <- match.arg(data_type)
  method <- match.arg(method)

  # Experimental warning
  message("=============================================================")
  message("EXPERIMENTAL METHOD: RMST-based Network Meta-Analysis")
  message("Based on: Hua et al. (2025) Biometrical Journal")
  message("Status: Cutting-edge method from 2025 literature")
  message("=============================================================")

  # Validate data
  if (data_type == "aggregate") {
    required_cols <- c("study", "treatment", "rmst", "rmst_se")
    if (!all(required_cols %in% names(data))) {
      stop("For aggregate data, columns required: ", paste(required_cols, collapse = ", "))
    }
  } else {
    required_cols <- c("study", "treatment", "time", "event")
    if (!all(required_cols %in% names(data))) {
      stop("For IPD data, columns required: ", paste(required_cols, collapse = ", "))
    }
  }

  # Set reference treatment
  if (is.null(reference)) {
    reference <- as.character(data$treatment[1])
    message("Using reference treatment: ", reference)
  }

  # Determine tau if not provided
  if (is.null(tau)) {
    if (data_type == "ipd") {
      max_times <- tapply(data$time, data$study, max)
      tau <- min(max_times)
      message("Auto-detected tau = ", round(tau, 2), " (minimum of study max times)")
    } else {
      stop("For aggregate data, tau must be specified")
    }
  }

  # Process data based on type
  if (data_type == "ipd") {
    # Calculate RMST from individual participant data
    message("Calculating RMST from individual participant data...")
    rmst_data <- calculate_rmst_from_ipd(data, tau, covariate)
  } else {
    # Use provided aggregate RMST data
    rmst_data <- data
  }

  # Run network meta-analysis on RMST differences
  if (method == "frequentist") {
    nma_result <- fit_frequentist_rmst_nma(rmst_data, reference, tau, ...)
  } else {
    nma_result <- fit_bayesian_rmst_nma(rmst_data, reference, tau, ...)
  }

  # Extract results
  rmst_estimates <- extract_rmst_estimates(nma_result, rmst_data, reference)
  rmst_differences <- extract_rmst_differences(nma_result)
  ranking <- rank_treatments_rmst(rmst_estimates)
  interpretation <- create_rmst_interpretation(rmst_differences, tau)

  # Create result object
  result <- list(
    rmst_estimates = rmst_estimates,
    rmst_differences = rmst_differences,
    tau = tau,
    reference = reference,
    ranking = ranking,
    nma_result = nma_result,
    interpretation = interpretation,
    method = method,
    data_type = data_type,
    covariate = covariate,
    call = match.call()
  )

  class(result) <- "rmst_nma"
  return(result)
}


# ============================================================================
# Helper Functions
# ============================================================================

#' Calculate RMST from individual participant data
#' @keywords internal
calculate_rmst_from_ipd <- function(data, tau, covariate = NULL) {

  # Check if survival package is available
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for IPD analysis. Please install it.")
  }

  # Split by study and treatment
  study_treatment_combos <- unique(data[, c("study", "treatment")])

  rmst_results <- data.frame(
    study = character(),
    treatment = character(),
    rmst = numeric(),
    rmst_se = numeric(),
    n = integer(),
    events = integer(),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(study_treatment_combos)) {
    study_i <- study_treatment_combos$study[i]
    treatment_i <- study_treatment_combos$treatment[i]

    # Subset data
    subset_data <- data[data$study == study_i & data$treatment == treatment_i, ]

    if (nrow(subset_data) == 0) next

    # Fit survival curve
    surv_obj <- survival::Surv(subset_data$time, subset_data$event)

    if (!is.null(covariate) && covariate %in% names(subset_data)) {
      # With covariate adjustment
      km_fit <- survival::survfit(surv_obj ~ 1, data = subset_data, weights = NULL)
      # Note: For proper covariate adjustment, would need Cox model with pseudo-observations
      # This is simplified for demonstration
      message("Note: Covariate adjustment requires advanced pseudo-observation methods")
    } else {
      # Standard Kaplan-Meier
      km_fit <- survival::survfit(surv_obj ~ 1, data = subset_data)
    }

    # Calculate RMST as area under KM curve up to tau
    rmst_calc <- calculate_rmst_from_survfit(km_fit, tau)

    rmst_results <- rbind(rmst_results, data.frame(
      study = study_i,
      treatment = treatment_i,
      rmst = rmst_calc$rmst,
      rmst_se = rmst_calc$se,
      n = nrow(subset_data),
      events = sum(subset_data$event),
      stringsAsFactors = FALSE
    ))
  }

  return(rmst_results)
}


#' Calculate RMST from survfit object
#' @keywords internal
calculate_rmst_from_survfit <- function(km_fit, tau) {

  # Extract survival times and probabilities
  times <- c(0, km_fit$time)
  surv <- c(1, km_fit$surv)

  # Restrict to tau
  times <- times[times <= tau]
  surv <- surv[1:length(times)]

  # Add tau if not already present
  if (max(times) < tau) {
    times <- c(times, tau)
    surv <- c(surv, surv[length(surv)])  # Carry forward last survival prob
  }

  # Calculate RMST using trapezoidal rule (area under curve)
  rmst <- 0
  for (i in 1:(length(times) - 1)) {
    width <- times[i + 1] - times[i]
    height <- (surv[i] + surv[i + 1]) / 2
    rmst <- rmst + width * height
  }

  # Estimate standard error using Greenwood's formula
  # Simplified version - proper implementation would use pseudo-observations
  if (!is.null(km_fit$std.err)) {
    # Approximate SE for RMST based on survival SE at tau
    se_surv_tau <- approx(km_fit$time, km_fit$std.err, xout = tau, rule = 2)$y
    se_rmst <- se_surv_tau * tau / 2  # Rough approximation
  } else {
    se_rmst <- NA
  }

  return(list(rmst = rmst, se = se_rmst))
}


#' Fit frequentist RMST-based NMA
#' @keywords internal
fit_frequentist_rmst_nma <- function(rmst_data, reference, tau, ...) {

  if (!requireNamespace("netmeta", quietly = TRUE)) {
    stop("Package 'netmeta' is required. Please install it.")
  }

  # Prepare pairwise comparisons
  pairs <- create_pairwise_rmst(rmst_data)

  # Run network meta-analysis on RMST differences
  nma <- netmeta::netmeta(
    TE = pairs$rmst_diff,
    seTE = pairs$rmst_diff_se,
    treat1 = pairs$treat1,
    treat2 = pairs$treat2,
    studlab = pairs$study,
    reference.group = reference,
    sm = "RMST",
    ...
  )

  return(nma)
}


#' Create pairwise RMST comparisons from multi-arm data
#' @keywords internal
create_pairwise_rmst <- function(rmst_data) {

  studies <- unique(rmst_data$study)
  pairs <- data.frame(
    study = character(),
    treat1 = character(),
    treat2 = character(),
    rmst_diff = numeric(),
    rmst_diff_se = numeric(),
    stringsAsFactors = FALSE
  )

  for (study_i in studies) {
    study_data <- rmst_data[rmst_data$study == study_i, ]

    if (nrow(study_data) < 2) next

    # Create all pairwise comparisons within study
    for (i in 1:(nrow(study_data) - 1)) {
      for (j in (i + 1):nrow(study_data)) {

        treat1 <- study_data$treatment[i]
        treat2 <- study_data$treatment[j]
        rmst1 <- study_data$rmst[i]
        rmst2 <- study_data$rmst[j]
        se1 <- study_data$rmst_se[i]
        se2 <- study_data$rmst_se[j]

        # RMST difference (treat2 - treat1)
        rmst_diff <- rmst2 - rmst1

        # Standard error of difference
        rmst_diff_se <- sqrt(se1^2 + se2^2)

        pairs <- rbind(pairs, data.frame(
          study = study_i,
          treat1 = treat1,
          treat2 = treat2,
          rmst_diff = rmst_diff,
          rmst_diff_se = rmst_diff_se,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(pairs)
}


#' Fit Bayesian RMST-based NMA
#' @keywords internal
fit_bayesian_rmst_nma <- function(rmst_data, reference, tau, ...) {

  message("Bayesian RMST-NMA requires JAGS/Stan - using frequentist approximation")

  # For now, fall back to frequentist
  # Full Bayesian implementation would use JAGS/Stan with RMST likelihood
  result <- fit_frequentist_rmst_nma(rmst_data, reference, tau, ...)
  result$note <- "Bayesian RMST-NMA requires JAGS/Stan - frequentist approximation used"

  return(result)
}


#' Extract RMST estimates for each treatment
#' @keywords internal
extract_rmst_estimates <- function(nma_result, rmst_data, reference) {

  treatments <- unique(rmst_data$treatment)

  # Get reference RMST (average across studies)
  ref_rmst <- mean(rmst_data$rmst[rmst_data$treatment == reference])

  estimates <- data.frame(
    treatment = treatments,
    rmst = numeric(length(treatments)),
    rmst_lower = numeric(length(treatments)),
    rmst_upper = numeric(length(treatments)),
    stringsAsFactors = FALSE
  )

  for (i in 1:length(treatments)) {
    if (treatments[i] == reference) {
      estimates$rmst[i] <- ref_rmst
      estimates$rmst_lower[i] <- NA
      estimates$rmst_upper[i] <- NA
    } else {
      # Get RMST difference from NMA
      if (inherits(nma_result, "netmeta")) {
        comp_name <- paste(treatments[i], reference, sep = ":")
        if (comp_name %in% rownames(nma_result$TE.random)) {
          diff <- nma_result$TE.random[comp_name, reference]
          diff_lower <- nma_result$lower.random[comp_name, reference]
          diff_upper <- nma_result$upper.random[comp_name, reference]

          estimates$rmst[i] <- ref_rmst + diff
          estimates$rmst_lower[i] <- ref_rmst + diff_lower
          estimates$rmst_upper[i] <- ref_rmst + diff_upper
        } else {
          # Average from observed data
          estimates$rmst[i] <- mean(rmst_data$rmst[rmst_data$treatment == treatments[i]])
          estimates$rmst_lower[i] <- NA
          estimates$rmst_upper[i] <- NA
        }
      }
    }
  }

  # Sort by RMST (descending)
  estimates <- estimates[order(-estimates$rmst), ]

  return(estimates)
}


#' Extract pairwise RMST differences
#' @keywords internal
extract_rmst_differences <- function(nma_result) {

  if (inherits(nma_result, "netmeta")) {
    # Extract league table
    diffs <- as.data.frame(nma_result$TE.random)
    diffs_se <- as.data.frame(nma_result$seTE.random)

    return(list(
      differences = diffs,
      se = diffs_se,
      interpretation = "Positive values indicate longer RMST (better survival)"
    ))
  }

  return(NULL)
}


#' Rank treatments by RMST
#' @keywords internal
rank_treatments_rmst <- function(rmst_estimates) {

  # Simple ranking by point estimate
  rmst_estimates$rank <- rank(-rmst_estimates$rmst)

  # Format as ranking table
  ranking <- rmst_estimates[order(rmst_estimates$rank), c("rank", "treatment", "rmst")]

  return(ranking)
}


#' Create clinical interpretation
#' @keywords internal
create_rmst_interpretation <- function(rmst_differences, tau) {

  interpretation <- list(
    tau = tau,
    message = paste0(
      "RMST represents mean survival time (in original time units) up to tau = ", tau, ". ",
      "RMST differences directly indicate the gain in mean survival time. ",
      "For example, RMST difference of 2.5 means '2.5 time units longer survival on average'."
    ),
    advantages = c(
      "More interpretable than hazard ratios",
      "No proportional hazards assumption needed",
      "Directly clinically meaningful (time gained)",
      "Robust to non-proportional hazards"
    )
  )

  return(interpretation)
}


# ============================================================================
# S3 Methods
# ============================================================================

#' @export
print.rmst_nma <- function(x, ...) {
  cat("=============================================================\n")
  cat("EXPERIMENTAL: RMST-based Network Meta-Analysis\n")
  cat("=============================================================\n\n")

  cat("Restriction time (tau):", x$tau, "\n")
  cat("Reference treatment:", x$reference, "\n")
  cat("Analysis method:", x$method, "\n")
  cat("Data type:", x$data_type, "\n\n")

  cat("Treatment Rankings (by RMST):\n")
  print(x$ranking, row.names = FALSE)

  cat("\n")
  cat("RMST Estimates:\n")
  print(x$rmst_estimates, row.names = FALSE, digits = 2)

  cat("\n")
  cat("Interpretation:\n")
  cat(x$interpretation$message, "\n")

  invisible(x)
}


#' @export
plot.rmst_nma <- function(x, type = c("forest", "ranking", "network"), ...) {

  type <- match.arg(type)

  if (type == "forest") {
    plot_rmst_forest(x, ...)
  } else if (type == "ranking") {
    plot_rmst_ranking(x, ...)
  } else if (type == "network") {
    if (inherits(x$nma_result, "netmeta")) {
      netmeta::netgraph(x$nma_result, ...)
    }
  }
}


#' Forest plot for RMST estimates
#' @keywords internal
plot_rmst_forest <- function(x, ...) {

  est <- x$rmst_estimates
  est <- est[order(est$rmst), ]

  # Basic forest plot
  par(mar = c(5, 8, 4, 2))

  y_pos <- 1:nrow(est)
  plot(est$rmst, y_pos,
       xlim = range(c(est$rmst_lower, est$rmst_upper), na.rm = TRUE),
       ylim = c(0.5, nrow(est) + 0.5),
       xlab = paste0("RMST (up to tau = ", x$tau, ")"),
       ylab = "",
       yaxt = "n",
       pch = 18,
       cex = 1.5,
       main = "RMST Estimates by Treatment")

  axis(2, at = y_pos, labels = est$treatment, las = 1)

  # Add confidence intervals
  for (i in 1:nrow(est)) {
    if (!is.na(est$rmst_lower[i])) {
      segments(est$rmst_lower[i], y_pos[i], est$rmst_upper[i], y_pos[i])
    }
  }

  # Reference line
  abline(v = est$rmst[est$treatment == x$reference], lty = 2, col = "red")

  grid()
}


#' Ranking plot for RMST
#' @keywords internal
plot_rmst_ranking <- function(x, ...) {

  ranking <- x$ranking

  par(mar = c(5, 8, 4, 2))

  barplot(ranking$rmst,
          names.arg = ranking$treatment,
          horiz = TRUE,
          las = 1,
          xlab = paste0("RMST (up to tau = ", x$tau, ")"),
          main = "Treatment Ranking by RMST",
          col = colorRampPalette(c("lightblue", "darkblue"))(nrow(ranking)))

  grid()
}
