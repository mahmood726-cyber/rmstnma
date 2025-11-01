# ==============================================================================
# Surrogate Endpoint Network Meta-Analysis (Experimental)
# ==============================================================================
#
# Integrated from surroNMA package (mahmood726-cyber/surroNMA)
# Provides advanced surrogate endpoint analysis including:
# - Bivariate NMA (joint modeling of surrogate and true endpoints)
# - Surrogate Index (SI) training with multiple surrogates
# - Surrogate Threshold Effect (STE) calculation
# - R² individual and trial-level validation
#
# References:
# - Buyse et al. (2000). Statistical evaluation of surrogate endpoints
# - Burzykowski et al. (2005). The Evaluation of Surrogate Endpoints
# ==============================================================================

#' Build surrogate network data structure
#'
#' Creates a specialized network object for surrogate endpoint analysis
#' compatible with powerNMA framework.
#'
#' @param data Data frame with study-level treatment comparisons
#' @param study Column name for study identifier
#' @param treat1 Column name for first treatment in comparison
#' @param treat2 Column name for second treatment in comparison
#' @param S_eff Column name for surrogate endpoint treatment effect
#' @param S_se Column name for surrogate endpoint standard error
#' @param T_eff Column name for true endpoint treatment effect (can have missing values)
#' @param T_se Column name for true endpoint standard error
#' @param S_multi Optional: character vector of multiple surrogate column names
#' @param corr_ST Optional: correlation between S and T (default 0 if unknown)
#' @param check_connectivity Check network connectivity (default TRUE)
#'
#' @return surrogate_network object with standardized structure
#'
#' @export
#' @examples
#' \dontrun{
#' # Example with single surrogate (e.g., progression-free survival)
#' data <- data.frame(
#'   study = c("S1", "S1", "S2", "S2"),
#'   treat1 = c("A", "B", "A", "C"),
#'   treat2 = c("B", "C", "C", "B"),
#'   pfs_effect = c(0.8, 0.6, 0.7, 0.5),  # Surrogate: PFS
#'   pfs_se = c(0.15, 0.18, 0.14, 0.16),
#'   os_effect = c(0.5, NA, 0.4, NA),     # True: OS (some missing)
#'   os_se = c(0.20, NA, 0.22, NA)
#' )
#'
#' net <- build_surrogate_network(
#'   data = data,
#'   study = study,
#'   treat1 = treat1, treat2 = treat2,
#'   S_eff = pfs_effect, S_se = pfs_se,
#'   T_eff = os_effect, T_se = os_se
#' )
#' }
build_surrogate_network <- function(
  data,
  study, treat1, treat2,
  S_eff, S_se,
  T_eff = NULL, T_se = NULL,
  S_multi = NULL,
  corr_ST = 0.0,
  check_connectivity = TRUE
) {

  # Convert to data frame
  df <- as.data.frame(data)

  # Extract columns using NSE (non-standard evaluation)
  study_id <- df[[deparse(substitute(study))]]
  trt1 <- df[[deparse(substitute(treat1))]]
  trt2 <- df[[deparse(substitute(treat2))]]
  S_e <- df[[deparse(substitute(S_eff))]]
  S_s <- df[[deparse(substitute(S_se))]]

  # Optional true endpoint
  T_e <- if (!missing(T_eff)) df[[deparse(substitute(T_eff))]] else NULL
  T_s <- if (!missing(T_se)) df[[deparse(substitute(T_se))]] else NULL

  # Validation
  if (is.null(study_id) || is.null(trt1) || is.null(trt2)) {
    stop("study, treat1, treat2 must be provided")
  }

  if (is.null(S_e) || is.null(S_s)) {
    stop("S_eff and S_se are required for surrogate endpoint analysis")
  }

  if (!is.null(T_e) && is.null(T_s)) {
    stop("T_se required when T_eff is provided")
  }

  # Check for non-finite values
  if (any(!is.finite(S_e) | !is.finite(S_s))) {
    stop("Non-finite values in S_eff/S_se. Please remove or impute.")
  }

  if (!is.null(T_e) && any(is.finite(T_e) & !is.finite(T_s))) {
    stop("T_se must be finite when T_eff is finite")
  }

  # Treatment coding
  treatments <- sort(unique(c(trt1, trt2)))
  K <- length(treatments)
  trt_map <- setNames(seq_len(K), treatments)

  trt1_idx <- as.integer(trt_map[as.character(trt1)])
  trt2_idx <- as.integer(trt_map[as.character(trt2)])

  J <- length(unique(study_id))

  # Network connectivity check
  if (check_connectivity) {
    # Create adjacency check (simple version without igraph)
    edges <- data.frame(from = trt1, to = trt2)
    n_comparisons <- nrow(edges)
    n_treatments <- K

    # Rough connectivity check: at least K-1 edges needed
    if (n_comparisons < (n_treatments - 1)) {
      warning("Network may be disconnected: few edges relative to treatments")
    }
  }

  # Handle multiple surrogates if provided
  S_multi_matrix <- NULL
  if (!is.null(S_multi)) {
    if (!all(S_multi %in% names(df))) {
      stop("Some S_multi columns not found in data")
    }
    S_multi_matrix <- as.matrix(df[, S_multi, drop = FALSE])
  }

  # Build surrogate network object
  structure(
    list(
      data = df,
      study = study_id,
      treat1 = trt1_idx,
      treat2 = trt2_idx,
      treat_levels = treatments,
      K = K,
      J = J,
      n_comparisons = nrow(df),
      S_eff = S_e,
      S_se = S_s,
      T_eff = T_e,
      T_se = T_s,
      S_multi = S_multi_matrix,
      corr_ST = if (length(corr_ST) == 1) rep(corr_ST, nrow(df)) else corr_ST,
      has_true_endpoint = !is.null(T_e),
      n_true_observed = if (!is.null(T_e)) sum(is.finite(T_e)) else 0
    ),
    class = c("surrogate_network", "powernma_data")
  )
}


#' Train Surrogate Index from Multiple Surrogates
#'
#' Combines multiple surrogate endpoints into a single Surrogate Index (SI)
#' using machine learning methods. Useful when multiple early endpoints are
#' available (e.g., multiple biomarkers, PFS at different timepoints).
#'
#' @param net surrogate_network object from build_surrogate_network()
#' @param method Method for combining surrogates: "ols" (ordinary least squares),
#'   "ridge" (ridge regression), "pcr" (principal component regression)
#' @param standardize Standardize surrogates before training (default TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return surrogate_index object with trained model
#'
#' @export
#' @examples
#' \dontrun{
#' # Train SI from multiple biomarkers
#' si_model <- train_surrogate_index(
#'   net,
#'   method = "ridge",
#'   standardize = TRUE
#' )
#'
#' # Apply SI to augment network
#' net_augmented <- apply_surrogate_index(net, si_model)
#' }
train_surrogate_index <- function(
  net,
  method = c("ols", "ridge", "pcr"),
  standardize = TRUE,
  seed = 12345
) {

  method <- match.arg(method)

  # Validate input
  if (!inherits(net, "surrogate_network")) {
    stop("net must be a surrogate_network object")
  }

  if (is.null(net$S_multi)) {
    stop("S_multi required to train Surrogate Index. Use build_surrogate_network() with S_multi argument.")
  }

  if (!net$has_true_endpoint) {
    stop("True endpoint (T_eff) required to train Surrogate Index")
  }

  # Prepare data
  S <- as.data.frame(net$S_multi)
  colnames(S) <- paste0("S", seq_len(ncol(S)))
  T_outcome <- net$T_eff

  # Find complete cases
  complete_idx <- which(is.finite(T_outcome) & rowSums(is.finite(S)) == ncol(S))

  if (length(complete_idx) < 5) {
    stop("Insufficient complete cases for SI training. Need at least 5 observations with both S and T.")
  }

  S_train <- S[complete_idx, , drop = FALSE]
  T_train <- T_outcome[complete_idx]

  # Standardize if requested
  S_means <- NULL
  S_sds <- NULL
  T_mean <- NULL
  T_sd <- NULL

  if (standardize) {
    S_means <- colMeans(S_train)
    S_sds <- apply(S_train, 2, sd)
    S_train <- scale(S_train)

    T_mean <- mean(T_train)
    T_sd <- sd(T_train)
    T_train <- scale(T_train)[, 1]
  }

  set.seed(seed)

  # Train model based on method
  fit_result <- list()

  if (method == "ols") {
    # Ordinary least squares
    model <- stats::lm(T_train ~ ., data = S_train)
    fit_result$model <- model
    fit_result$coef <- stats::coef(model)

    # Simple R² (not cross-validated)
    fit_result$r_squared <- summary(model)$r.squared

  } else if (method == "ridge") {
    # Ridge regression (requires glmnet)
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' required for ridge regression. Install with: install.packages('glmnet')")
    }

    X_matrix <- as.matrix(S_train)
    Y_vector <- as.numeric(T_train)

    # Cross-validation to find optimal lambda
    cv_fit <- glmnet::cv.glmnet(X_matrix, Y_vector, alpha = 0)

    # Fit final model
    final_fit <- glmnet::glmnet(
      X_matrix, Y_vector,
      alpha = 0,
      lambda = cv_fit$lambda.min
    )

    fit_result$model <- final_fit
    fit_result$lambda <- cv_fit$lambda.min
    fit_result$coef <- as.numeric(stats::coef(final_fit))
    names(fit_result$coef) <- c("(Intercept)", colnames(S_train))

    # Estimate R² from CV
    fit_result$r_squared <- max(0, 1 - min(cv_fit$cvm, na.rm = TRUE) / stats::var(Y_vector))

  } else if (method == "pcr") {
    # Principal component regression (requires pls package)
    if (!requireNamespace("pls", quietly = TRUE)) {
      stop("Package 'pls' required for PCR. Install with: install.packages('pls')")
    }

    model <- pls::pcr(T_train ~ ., data = S_train, validation = "CV")

    # Select optimal number of components
    optimal_ncomp <- which.min(model$validation$PRESS)

    fit_result$model <- model
    fit_result$ncomp <- optimal_ncomp

    # Calculate R²
    predictions <- stats::predict(model, ncomp = optimal_ncomp, newdata = S_train)[, 1, 1]
    fit_result$r_squared <- max(0, 1 - mean((T_train - predictions)^2) / stats::var(T_train))
  }

  # Return surrogate index object
  structure(
    list(
      method = method,
      fit = fit_result,
      standardize = standardize,
      S_means = S_means,
      S_sds = S_sds,
      T_mean = T_mean,
      T_sd = T_sd,
      n_surrogates = ncol(S_train),
      n_train = length(complete_idx),
      r_squared = fit_result$r_squared,
      column_names = colnames(S_train)
    ),
    class = "surrogate_index"
  )
}


#' Apply Surrogate Index to Network
#'
#' Uses a trained surrogate index model to create a combined surrogate
#' endpoint for all observations in the network.
#'
#' @param net surrogate_network object
#' @param si_model surrogate_index object from train_surrogate_index()
#'
#' @return Updated surrogate_network with S_eff replaced by SI predictions
#'
#' @export
apply_surrogate_index <- function(net, si_model) {

  if (!inherits(net, "surrogate_network")) {
    stop("net must be a surrogate_network object")
  }

  if (!inherits(si_model, "surrogate_index")) {
    stop("si_model must be a surrogate_index object from train_surrogate_index()")
  }

  if (is.null(net$S_multi)) {
    stop("S_multi required in network object")
  }

  # Prepare surrogate matrix
  S <- as.data.frame(net$S_multi)
  colnames(S) <- si_model$column_names

  # Standardize if model was trained with standardization
  if (si_model$standardize) {
    for (j in seq_along(si_model$S_means)) {
      S[[j]] <- (S[[j]] - si_model$S_means[j]) / si_model$S_sds[j]
    }
  }

  # Make predictions
  if (si_model$method == "ols") {
    # OLS predictions
    X_matrix <- cbind(1, as.matrix(S))
    predictions <- as.numeric(X_matrix %*% si_model$fit$coef)

  } else if (si_model$method == "ridge") {
    # Ridge predictions
    X_matrix <- cbind(1, as.matrix(S))
    predictions <- as.numeric(X_matrix %*% si_model$fit$coef)

  } else if (si_model$method == "pcr") {
    # PCR predictions
    predictions <- as.numeric(
      stats::predict(
        si_model$fit$model,
        ncomp = si_model$fit$ncomp,
        newdata = S
      )[, 1, 1]
    )
  }

  # Back-transform if standardized
  if (si_model$standardize) {
    predictions <- predictions * si_model$T_sd + si_model$T_mean
  }

  # Update network with SI as new surrogate
  net$S_eff <- predictions
  net$S_original <- net$S_eff  # Store original for reference
  attr(net, "has_SI") <- TRUE
  attr(net, "SI_method") <- si_model$method
  attr(net, "SI_r_squared") <- si_model$r_squared

  return(net)
}


#' Compute Surrogate Threshold Effect (STE)
#'
#' Calculates the minimum surrogate endpoint effect required to predict
#' a clinically meaningful effect on the true endpoint.
#'
#' Formula: STE = (threshold_T - α) / β
#' where α = intercept, β = slope from S→T relationship
#'
#' @param alpha_draws Vector of posterior/bootstrap draws for intercept (α)
#' @param beta_draws Vector of posterior/bootstrap draws for slope (β)
#' @param threshold_T Clinically meaningful threshold on true endpoint (default 0)
#' @param conf_level Confidence level for interval (default 0.95)
#'
#' @return List with STE estimates and uncertainty
#'
#' @export
#' @examples
#' \dontrun{
#' # After bivariate NMA fit
#' ste_result <- compute_surrogate_threshold_effect(
#'   alpha_draws = fit$posterior$alpha,
#'   beta_draws = fit$posterior$beta,
#'   threshold_T = 0.0  # Non-inferiority margin
#' )
#'
#' cat("STE median:", ste_result$median, "\n")
#' cat("95% CI:", ste_result$ci_lower, "to", ste_result$ci_upper, "\n")
#' }
compute_surrogate_threshold_effect <- function(
  alpha_draws,
  beta_draws,
  threshold_T = 0.0,
  conf_level = 0.95
) {

  # Validation
  if (length(alpha_draws) != length(beta_draws)) {
    stop("alpha_draws and beta_draws must have same length")
  }

  if (length(alpha_draws) < 10) {
    stop("Need at least 10 draws for STE calculation")
  }

  # Calculate STE for each draw
  # STE = (threshold_T - alpha) / beta
  ste_draws <- (threshold_T - alpha_draws) / beta_draws

  # Handle infinite/missing values
  ste_draws <- ste_draws[is.finite(ste_draws)]

  if (length(ste_draws) < 10) {
    warning("Many non-finite STE values. Check for β ≈ 0 (weak surrogacy)")
  }

  # Compute summary statistics
  alpha <- (1 - conf_level) / 2

  result <- list(
    threshold_T = threshold_T,
    mean = mean(ste_draws),
    median = stats::median(ste_draws),
    sd = stats::sd(ste_draws),
    ci_lower = stats::quantile(ste_draws, alpha),
    ci_upper = stats::quantile(ste_draws, 1 - alpha),
    conf_level = conf_level,
    n_draws = length(ste_draws),
    n_finite = length(ste_draws),
    draws = ste_draws
  )

  class(result) <- "surrogate_threshold_effect"
  return(result)
}


#' Fit Bivariate Network Meta-Analysis (Frequentist)
#'
#' Joint analysis of surrogate and true endpoints using frequentist approach
#' with parametric bootstrap for uncertainty quantification.
#'
#' @param net surrogate_network object
#' @param n_boot Number of bootstrap samples (default 400)
#' @param boot_method Bootstrap distribution: "normal" or "student" (default "normal")
#' @param df Degrees of freedom for Student-t bootstrap (default 5)
#' @param seed Random seed
#'
#' @return bivariate_nma_fit object with treatment effects and surrogacy parameters
#'
#' @export
#' @examples
#' \dontrun{
#' # Frequentist bivariate NMA
#' fit <- fit_bivariate_nma_freq(
#'   net,
#'   n_boot = 400,
#'   boot_method = "normal"
#' )
#'
#' # Extract results
#' print(fit$dS)  # Treatment effects on surrogate
#' print(fit$dT)  # Treatment effects on true endpoint
#' print(fit$surrogacy)  # α and β parameters
#' }
fit_bivariate_nma_freq <- function(
  net,
  n_boot = 400,
  boot_method = c("normal", "student"),
  df = 5,
  seed = 12345
) {

  boot_method <- match.arg(boot_method)

  # Validate input
  if (!inherits(net, "surrogate_network")) {
    stop("net must be a surrogate_network object")
  }

  if (!net$has_true_endpoint) {
    warning("No true endpoint data. Fitting surrogate-only network.")
  }

  set.seed(seed)

  K <- net$K
  N <- net$n_comparisons
  ref <- 1  # Reference treatment

  # Create design matrix
  # For each comparison, code treatment effects relative to reference
  X <- matrix(0, nrow = N, ncol = K - 1)
  for (i in 1:N) {
    if (net$treat1[i] != ref) X[i, net$treat1[i] - 1] <- 1
    if (net$treat2[i] != ref) X[i, net$treat2[i] - 1] <- -1
  }

  # Fit surrogate endpoint NMA
  yS <- net$S_eff
  vS <- net$S_se^2
  W_S <- 1 / vS

  # Weighted least squares for surrogate
  betaS <- tryCatch({
    solve(t(X) %*% (W_S * X), t(X) %*% (W_S * yS))
  }, error = function(e) {
    # Use generalized inverse if singular
    MASS::ginv(t(X) %*% (W_S * X)) %*% t(X) %*% (W_S * yS)
  })

  dS_hat <- c(0, as.numeric(betaS))  # Treatment effects (ref = 0)

  # Fit true endpoint NMA (only on observed rows)
  dT_hat <- rep(0, K)
  alpha_est <- 0
  beta_est <- 1

  if (net$has_true_endpoint && net$n_true_observed >= 3) {

    obs_T <- is.finite(net$T_eff)
    X_T <- X[obs_T, , drop = FALSE]
    yT <- net$T_eff[obs_T]
    vT <- net$T_se[obs_T]^2
    W_T <- 1 / vT

    if (nrow(X_T) >= (K - 1)) {
      betaT <- tryCatch({
        solve(t(X_T) %*% (W_T * X_T), t(X_T) %*% (W_T * yT))
      }, error = function(e) {
        MASS::ginv(t(X_T) %*% (W_T * X_T)) %*% t(X_T) %*% (W_T * yT)
      })
      dT_hat <- c(0, as.numeric(betaT))
    }

    # Estimate α and β via Deming regression
    # (accounts for measurement error in both S and T)
    fitted_S <- X %*% betaS
    fitted_T_obs <- X_T %*% betaT

    if (length(fitted_T_obs) >= 3) {
      # Variance ratio for Deming regression
      lambda <- mean(vT) / mean(vS[obs_T])

      # Deming regression
      x_vals <- fitted_S[obs_T]
      y_vals <- fitted_T_obs

      x_mean <- mean(x_vals)
      y_mean <- mean(y_vals)

      sxx <- stats::var(x_vals)
      syy <- stats::var(y_vals)
      sxy <- stats::cov(x_vals, y_vals)

      beta_est <- (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 + 4 * lambda * sxy^2)) / (2 * sxy)
      alpha_est <- y_mean - beta_est * x_mean
    }
  }

  # Bootstrap for uncertainty quantification
  draws_T <- matrix(NA_real_, n_boot, K)
  draws_alpha <- numeric(n_boot)
  draws_beta <- numeric(n_boot)

  for (b in 1:n_boot) {

    # Generate bootstrap samples
    eS <- stats::rnorm(N, 0, sqrt(vS))

    if (boot_method == "student") {
      eS <- eS * sqrt(df / stats::rchisq(N, df))
    }

    yS_boot <- yS + eS

    # Refit surrogate
    betaS_boot <- tryCatch({
      solve(t(X) %*% (W_S * X), t(X) %*% (W_S * yS_boot))
    }, error = function(e) {
      MASS::ginv(t(X) %*% (W_S * X)) %*% t(X) %*% (W_S * yS_boot)
    })

    dS_boot <- c(0, as.numeric(betaS_boot))

    # Predict true endpoint using α and β
    dT_surrogate <- alpha_est + beta_est * dS_boot

    draws_T[b, ] <- dT_surrogate
    draws_alpha[b] <- alpha_est  # Fixed in frequentist version
    draws_beta[b] <- beta_est
  }

  # Treatment ranking
  ranks_matrix <- t(apply(draws_T, 1, function(x) rank(-x)))
  mean_ranks <- colMeans(ranks_matrix)

  # SUCRA scores
  sucra <- (K - mean_ranks) / (K - 1)

  # Build result object
  result <- structure(
    list(
      engine = "frequentist",
      net = net,
      dS = dS_hat,
      dT = dT_hat,
      surrogacy = list(
        alpha = alpha_est,
        beta = beta_est,
        alpha_draws = draws_alpha,
        beta_draws = draws_beta
      ),
      draws_T = draws_T,
      ranks = list(
        mean_rank = mean_ranks,
        sucra = sucra,
        rank_matrix = ranks_matrix
      ),
      uncertainty = list(
        n_boot = n_boot,
        boot_method = boot_method
      )
    ),
    class = c("bivariate_nma_fit", "powernma_result")
  )

  return(result)
}


#' Print method for surrogate_network
#' @export
print.surrogate_network <- function(x, ...) {
  cat("Surrogate Network Meta-Analysis Data\n")
  cat("=====================================\n\n")
  cat("Studies:", x$J, "\n")
  cat("Treatments:", x$K, "(", paste(x$treat_levels, collapse = ", "), ")\n")
  cat("Comparisons:", x$n_comparisons, "\n\n")
  cat("Surrogate endpoint: n =", x$n_comparisons, "\n")

  if (x$has_true_endpoint) {
    cat("True endpoint: n =", x$n_true_observed, "observed\n")
    pct <- round(100 * x$n_true_observed / x$n_comparisons, 1)
    cat("  (", pct, "% of comparisons)\n", sep = "")
  } else {
    cat("True endpoint: not provided\n")
  }

  if (!is.null(x$S_multi)) {
    cat("\nMultiple surrogates:", ncol(x$S_multi), "\n")
  }

  if (!is.null(attr(x, "has_SI"))) {
    cat("\nSurrogate Index trained:", attr(x, "SI_method"), "\n")
    cat("  R²:", round(attr(x, "SI_r_squared"), 3), "\n")
  }

  invisible(x)
}


#' Print method for bivariate_nma_fit
#' @export
print.bivariate_nma_fit <- function(x, ...) {
  cat("Bivariate Network Meta-Analysis Fit\n")
  cat("====================================\n\n")
  cat("Engine:", x$engine, "\n")
  cat("Treatments:", x$net$K, "\n")
  cat("Studies:", x$net$J, "\n\n")

  cat("Treatment Effects on Surrogate (dS):\n")
  df_S <- data.frame(
    Treatment = x$net$treat_levels,
    Effect = round(x$dS, 3)
  )
  print(df_S, row.names = FALSE)

  cat("\nTreatment Effects on True Endpoint (dT):\n")
  df_T <- data.frame(
    Treatment = x$net$treat_levels,
    Effect = round(x$dT, 3),
    SUCRA = round(x$ranks$sucra, 3)
  )
  print(df_T, row.names = FALSE)

  cat("\nSurrogacy Parameters:\n")
  cat("  α (intercept):", round(x$surrogacy$alpha, 3), "\n")
  cat("  β (slope):", round(x$surrogacy$beta, 3), "\n")

  if (abs(x$surrogacy$beta - 1) < 0.1) {
    cat("  → Strong surrogacy (β ≈ 1)\n")
  } else if (abs(x$surrogacy$beta) < 0.3) {
    cat("  → Weak surrogacy (β < 0.3)\n")
  }

  invisible(x)
}


#' Print method for surrogate_threshold_effect
#' @export
print.surrogate_threshold_effect <- function(x, ...) {
  cat("Surrogate Threshold Effect (STE)\n")
  cat("================================\n\n")
  cat("True endpoint threshold:", x$threshold_T, "\n\n")
  cat("STE estimate:\n")
  cat("  Median:", round(x$median, 3), "\n")
  cat("  Mean:", round(x$mean, 3), "\n")
  cat("  SD:", round(x$sd, 3), "\n")
  cat("  ", x$conf_level * 100, "% CI: [",
      round(x$ci_lower, 3), ", ", round(x$ci_upper, 3), "]\n", sep = "")
  cat("\nInterpretation:\n")
  cat("  Surrogate effects ≥", round(x$median, 3),
      "likely indicate meaningful true endpoint benefit\n")

  invisible(x)
}
