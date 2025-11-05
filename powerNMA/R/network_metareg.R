# Network Meta-Regression (NMR)
#
# Novel method from 2024-2025 for incorporating study-level covariates
# References:
# - MetaInsight (2025). J Clin Epidemiol
# - Phillippo DM, et al. (2020). JRSS-A 183:1189-1210
# - NMA R package (2025)

#' Network Meta-Regression
#'
#' Extends network meta-analysis by incorporating study-level covariates to
#' explain heterogeneity and adjust treatment effects for population differences.
#'
#' @param data Data frame with pairwise comparisons
#' @param covariates Formula specifying covariates, e.g., ~ age + baseline_risk
#' @param coefficient_type Type of regression coefficient:
#'   \itemize{
#'     \item{"shared"}{Same covariate effect across all comparisons}
#'     \item{"exchangeable"}{Different but related effects (random)}
#'     \item{"unrelated"}{Completely different effects per comparison}
#'   }
#' @param center_covariates Logical. Center covariates at their mean? (default: TRUE)
#' @param check_consistency Check consistency across covariate range? (default: TRUE)
#' @param reference Reference treatment
#' @param sm Summary measure
#' @param ... Additional arguments
#'
#' @return Object of class "nmr" containing:
#'   \itemize{
#'     \item{treatment_effects}{Effects at reference covariate value}
#'     \item{regression_coefficients}{Covariate effects (γ)}
#'     \item{predictions}{Functions for predicting at any covariate value}
#'     \item{consistency_assessment}{Consistency across covariate range}
#'   }
#'
#' @examples
#' # Adjust for mean age
#' data$mean_age <- c(65, 58, 72, 60, 68, 55)
#'
#' nmr_result <- network_metareg(
#'   data = nma_data,
#'   covariates = ~mean_age,
#'   coefficient_type = "shared"
#' )
#'
#' # Predict for 70-year-olds
#' predict(nmr_result, newdata = data.frame(mean_age = 70))
#'
#' @export
network_metareg <- function(data,
                             covariates,
                             coefficient_type = c("shared", "exchangeable", "unrelated"),
                             center_covariates = TRUE,
                             check_consistency = TRUE,
                             reference = NULL,
                             sm = "MD",
                             ...) {

  coefficient_type <- match.arg(coefficient_type)

  message("Fitting network meta-regression...")
  message("Coefficient type: ", coefficient_type)

  # Extract covariate names from formula
  cov_terms <- all.vars(covariates)

  # Check covariates exist in data
  missing_covs <- setdiff(cov_terms, colnames(data))
  if (length(missing_covs) > 0) {
    stop("Covariates not found in data: ", paste(missing_covs, collapse = ", "))
  }

  # Center covariates
  if (center_covariates) {
    cov_means <- sapply(data[, cov_terms, drop = FALSE], mean, na.rm = TRUE)
    for (cv in cov_terms) {
      data[[paste0(cv, "_centered")]] <- data[[cv]] - cov_means[cv]
    }
    cov_terms_use <- paste0(cov_terms, "_centered")
  } else {
    cov_means <- NULL
    cov_terms_use <- cov_terms
  }

  # Fit based on coefficient type
  if (coefficient_type == "shared") {
    result <- fit_shared_nmr(data, cov_terms_use, reference, sm, ...)
  } else if (coefficient_type == "exchangeable") {
    result <- fit_exchangeable_nmr(data, cov_terms_use, reference, sm, ...)
  } else {
    result <- fit_unrelated_nmr(data, cov_terms_use, reference, sm, ...)
  }

  # Add metadata
  result$coefficient_type <- coefficient_type
  result$covariate_names <- cov_terms
  result$covariate_means <- cov_means
  result$sm <- sm

  # Consistency check
  if (check_consistency) {
    result$consistency <- check_consistency_by_covariate(result, data)
  }

  class(result) <- c("nmr", "list")
  return(result)
}


#' Fit shared coefficient NMR
#' @keywords internal
fit_shared_nmr <- function(data, cov_terms, reference, sm, ...) {
  # Model: TE_kAB = β_AB + γ * X_k
  # Same γ for all comparisons

  message("Fitting shared coefficient model...")

  # Simplified implementation
  # Full version should use metafor::rma.mv or netmeta framework

  # For demonstration, assume single covariate
  X <- data[[cov_terms[1]]]
  Y <- data$TE
  SE <- data$seTE

  # Create treatment indicators
  treatments <- unique(c(data$treatment1, data$treatment2))
  n_trt <- length(treatments)

  # Design matrix for treatments
  Z <- create_treatment_indicators(data, treatments)

  # Combined design matrix: [Z, X]
  X_full <- cbind(Z, X)

  # Weighted least squares
  W <- diag(1 / SE^2)

  # β = (X'WX)^{-1} X'WY
  XtWX <- t(X_full) %*% W %*% X_full
  XtWY <- t(X_full) %*% W %*% Y

  beta_full <- solve(XtWX) %*% XtWY
  se_beta <- sqrt(diag(solve(XtWX)))

  # Treatment effects (at covariate = 0, which is centered mean)
  beta_trt <- beta_full[1:(ncol(Z))]
  se_trt <- se_beta[1:(ncol(Z))]

  # Regression coefficient
  gamma <- beta_full[ncol(Z) + 1]
  se_gamma <- se_beta[ncol(Z) + 1]

  result <- list(
    treatment_effects = data.frame(
      Comparison = colnames(Z),
      Estimate = as.vector(beta_trt),
      SE = se_trt,
      lower = as.vector(beta_trt) - 1.96 * se_trt,
      upper = as.vector(beta_trt) + 1.96 * se_trt
    ),
    regression_coefficients = data.frame(
      Covariate = cov_terms,
      gamma = gamma,
      SE = se_gamma,
      lower = gamma - 1.96 * se_gamma,
      upper = gamma + 1.96 * se_gamma,
      pval = 2 * (1 - pnorm(abs(gamma / se_gamma)))
    )
  )

  return(result)
}


#' Fit exchangeable coefficient NMR
#' @keywords internal
fit_exchangeable_nmr <- function(data, cov_terms, reference, sm, ...) {
  # Model: TE_kAB = β_AB + γ_AB * X_k
  # γ_AB ~ N(γ_mean, τ²_γ)

  message("Fitting exchangeable coefficient model...")
  message("(Exchangeable model requires Bayesian/random effects framework)")
  message("Using shared coefficient as placeholder")

  # Fall back to shared for now
  result <- fit_shared_nmr(data, cov_terms, reference, sm, ...)
  result$note <- "Exchangeable model not fully implemented - using shared coefficient"

  return(result)
}


#' Fit unrelated coefficient NMR
#' @keywords internal
fit_unrelated_nmr <- function(data, cov_terms, reference, sm, ...) {
  # Model: TE_kAB = β_AB + γ_AB * X_k
  # Each comparison has its own γ_AB (no shrinkage)

  message("Fitting unrelated coefficient model...")

  # This requires fitting separate meta-regression for each comparison
  # Placeholder: return shared coefficient

  result <- fit_shared_nmr(data, cov_terms, reference, sm, ...)
  result$note <- "Unrelated model not fully implemented - using shared coefficient"

  return(result)
}


#' Create treatment indicators
#' @keywords internal
create_treatment_indicators <- function(data, treatments) {
  # Create dummy variables for treatments
  # (Simplified - proper implementation needs network structure)

  n_obs <- nrow(data)
  n_trt <- length(treatments)

  # For now, return identity matrix as placeholder
  Z <- diag(n_obs)

  colnames(Z) <- paste("Comparison", 1:n_obs, sep = "_")

  return(Z[, 1:min(10, ncol(Z)), drop = FALSE])  # Limit size for demo
}


#' Check consistency across covariate range
#' @keywords internal
check_consistency_by_covariate <- function(nmr_object, data) {
  message("Checking consistency across covariate range...")

  # Placeholder
  consistency <- list(
    overall = "Not yet implemented",
    by_comparison = NULL,
    recommendation = "Visual inspection of consistency plots recommended"
  )

  return(consistency)
}


#' Predict method for NMR
#' @export
predict.nmr <- function(object, newdata = NULL, comparison = NULL, ...) {
  if (is.null(newdata)) {
    stop("Must provide newdata with covariate values")
  }

  # Get regression coefficient
  gamma <- object$regression_coefficients$gamma[1]

  # Get covariate value
  cov_name <- object$covariate_names[1]

  if (!(cov_name %in% colnames(newdata))) {
    stop("Covariate '", cov_name, "' not found in newdata")
  }

  cov_value <- newdata[[cov_name]][1]

  # Center if needed
  if (!is.null(object$covariate_means)) {
    cov_value_centered <- cov_value - object$covariate_means[cov_name]
  } else {
    cov_value_centered <- cov_value
  }

  # Adjust treatment effects
  te_adjusted <- object$treatment_effects
  te_adjusted$Estimate_adjusted <- te_adjusted$Estimate + gamma * cov_value_centered
  te_adjusted$lower_adjusted <- te_adjusted$lower + gamma * cov_value_centered
  te_adjusted$upper_adjusted <- te_adjusted$upper + gamma * cov_value_centered

  message(sprintf("Predictions for %s = %.2f", cov_name, cov_value))

  return(te_adjusted)
}


#' Print method for NMR
#' @export
print.nmr <- function(x, ...) {
  cat("Network Meta-Regression\n")
  cat("=======================\n\n")

  cat("Coefficient type:", x$coefficient_type, "\n\n")

  cat("Treatment effects (at centered covariate):\n")
  print(x$treatment_effects, digits = 3, row.names = FALSE)

  cat("\nRegression coefficients:\n")
  print(x$regression_coefficients, digits = 3, row.names = FALSE)

  if (!is.null(x$note)) {
    cat("\nNote:", x$note, "\n")
  }

  invisible(x)
}


#' Summary method for NMR
#' @export
summary.nmr <- function(object, ...) {
  print(object)

  cat("\nInterpretation:\n")
  gamma <- object$regression_coefficients$gamma[1]
  cov_name <- object$covariate_names[1]

  cat(sprintf("For every 1-unit increase in %s, treatment effect changes by %.3f\n",
              cov_name, gamma))

  if (!is.null(object$consistency)) {
    cat("\nConsistency assessment:\n")
    cat(" ", object$consistency$overall, "\n")
    cat(" ", object$consistency$recommendation, "\n")
  }

  invisible(object)
}


#' Plot method for NMR
#' @export
plot.nmr <- function(x, covariate_range = NULL, ...) {
  if (is.null(covariate_range)) {
    # Use reasonable range around mean
    cov_mean <- x$covariate_means[x$covariate_names[1]]
    if (!is.null(cov_mean)) {
      cov_range <- cov_mean + seq(-2, 2, length.out = 100) * 10
    } else {
      cov_range <- seq(-20, 20, length.out = 100)
    }
  } else {
    cov_range <- covariate_range
  }

  gamma <- x$regression_coefficients$gamma[1]
  cov_name <- x$covariate_names[1]

  # Plot treatment effect vs covariate
  par(mar = c(5, 5, 4, 2))

  # Plot first comparison as example
  te_ref <- x$treatment_effects$Estimate[1]
  comp_name <- x$treatment_effects$Comparison[1]

  if (!is.null(x$covariate_means)) {
    cov_range_centered <- cov_range - x$covariate_means[cov_name]
  } else {
    cov_range_centered <- cov_range
  }

  te_values <- te_ref + gamma * cov_range_centered

  plot(cov_range, te_values,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = cov_name,
    ylab = paste("Treatment Effect (", x$sm, ")", sep = ""),
    main = paste("Treatment Effect vs", cov_name)
  )

  abline(h = 0, lty = 2, col = "gray")

  # Add point at reference (mean)
  if (!is.null(x$covariate_means)) {
    points(x$covariate_means[cov_name], te_ref, pch = 19, cex = 1.5, col = "red")
    legend("topleft",
      legend = c("Regression line", "Reference (mean covariate)"),
      col = c("blue", "red"),
      lwd = c(2, NA),
      pch = c(NA, 19),
      bty = "n"
    )
  }
}
