#' Advanced Weighting Schemes
#'
#' Compute transportability weights based on covariate distances
#'
#' @param data Study data with covariates
#' @param target_population Named list of target covariate values
#' @param metric Distance metric ("mahalanobis" or "euclidean")
#' @param kernel Kernel function ("gaussian" or "tricube")
#' @param truncation Truncation proportion for extreme weights
#' @param min_weight Minimum weight floor
#' @return Vector of study weights
#' @export
compute_transport_weights <- function(data, target_population,
                                     metric = c("mahalanobis", "euclidean"),
                                     kernel = c("gaussian", "tricube"),
                                     truncation = 0.02,
                                     min_weight = 1e-6) {
  metric <- match.arg(metric)
  kernel <- match.arg(kernel)

  study_df <- data %>% dplyr::distinct(studlab, .keep_all = TRUE)
  covars <- intersect(names(target_population), names(study_df))

  if (!length(covars)) {
    msg("Transportability: no overlapping covariates; using equal weights")
    return(rep(1, nrow(study_df)))
  }

  X <- as.matrix(study_df[, covars, drop = FALSE])
  target_vec <- unlist(target_population[covars])

  # Compute distances
  if (metric == "mahalanobis") {
    S <- stats::cov(stats::na.omit(X))
    if (any(is.na(S)) || det(as.matrix(S)) <= .Machine$double.eps) {
      S <- diag(ncol(X))
    }
    diff_mat <- X - matrix(target_vec, nrow(X), length(target_vec), byrow = TRUE)
    d2 <- diag(diff_mat %*% solve(S) %*% t(diff_mat))
    d <- sqrt(pmax(0, d2))
  } else {
    d <- sqrt(rowSums((X - matrix(target_vec, nrow(X), length(target_vec),
                                  byrow = TRUE))^2))
  }

  # Apply kernel
  s <- stats::median(d, na.rm = TRUE) + 1e-8
  if (kernel == "gaussian") {
    w <- exp(-0.5 * (d / s)^2)
  } else {
    u <- safe_clip(d / s, 0, 1)
    w <- (1 - u^3)^3
  }

  w[!is.finite(w)] <- 0

  # Truncate extreme weights
  if (truncation > 0) {
    lo <- stats::quantile(w, probs = truncation, na.rm = TRUE)
    hi <- stats::quantile(w, probs = 1 - truncation, na.rm = TRUE)
    w <- safe_clip(w, lo, hi)
  }

  w <- pmax(min_weight, w)
  w <- w / mean(w, na.rm = TRUE)

  tibble::tibble(studlab = study_df$studlab, weight = w)
}

#' PET-PEESE Publication Bias Analysis
#'
#' Precision-Effect Test and Precision-Effect Estimate with Standard Error
#'
#' @param data Pairwise NMA data
#' @param min_k Minimum number of studies per comparison
#' @return Data frame with PET-PEESE results
#' @export
pet_peese_analysis <- function(data, min_k = 10) {
  df <- data %>%
    dplyr::mutate(
      cmp = paste(pmin(treat1, treat2), pmax(treat1, treat2), sep = " vs ")
    )

  comps <- df %>%
    dplyr::count(cmp) %>%
    dplyr::filter(n >= min_k) %>%
    dplyr::pull(cmp)

  if (!length(comps)) {
    msg("No comparison with >= %d studies", min_k)
    return(NULL)
  }

  results_list <- lapply(comps, function(cmp) {
    d <- df %>% dplyr::filter(cmp == !!cmp)
    w <- 1 / (d$seTE^2)

    # PET
    pet <- .safe_try(
      stats::lm(TE ~ seTE, data = d, weights = w),
      "PET regression"
    )

    pet_int <- NA_real_
    pet_p <- NA_real_
    use_peese <- FALSE

    if (!inherits(pet, "try-error")) {
      pet_summ <- summary(pet)
      pet_int <- stats::coef(pet_summ)["(Intercept)", "Estimate"]
      pet_p <- stats::coef(pet_summ)["seTE", "Pr(>|t|)"]
      use_peese <- is.finite(pet_p) && pet_p < 0.10
    }

    peese_int <- NA_real_
    if (use_peese) {
      peese <- .safe_try(
        stats::lm(TE ~ I(seTE^2), data = d, weights = w),
        "PEESE regression"
      )
      if (!inherits(peese, "try-error")) {
        peese_int <- stats::coef(summary(peese))["(Intercept)", "Estimate"]
      }
    }

    tibble::tibble(
      comparison = cmp,
      k = nrow(d),
      PET_intercept = pet_int,
      PET_p_se = pet_p,
      PEESE_intercept = peese_int,
      used_PEESE = use_peese
    )
  })

  dplyr::bind_rows(results_list)
}

#' Leave-One-Treatment-Out Sensitivity
#'
#' @param data Pairwise data
#' @param ref_treatment Reference treatment
#' @param sm Summary measure
#' @return LOTO results
#' @export
loto_sensitivity <- function(data, ref_treatment, sm = "HR") {
  trts <- sort(unique(c(data$treat1, data$treat2)))

  results_list <- lapply(trts, function(t) {
    dd <- data %>% dplyr::filter(treat1 != t & treat2 != t)

    if (nrow(dd) < 5) {
      return(tibble::tibble(
        treatment = t,
        tau = NA_real_,
        I2 = NA_real_
      ))
    }

    new_ref <- ifelse(ref_treatment == t,
                     setdiff(trts, t)[1],
                     ref_treatment)

    fit <- .safe_try(
      robust_netmeta(dd, new_ref, sm = sm, random = TRUE),
      paste0("LOTO: ", t)
    )

    if (inherits(fit, "try-error")) {
      tibble::tibble(treatment = t, tau = NA_real_, I2 = NA_real_)
    } else {
      tibble::tibble(
        treatment = t,
        tau = fit$tau,
        I2 = fit$I2.random
      )
    }
  })

  dplyr::bind_rows(results_list)
}

#' Multiverse Analysis
#'
#' Test robustness across modeling choices
#'
#' @param data Pairwise data
#' @param ref_treatment Reference treatment
#' @param sm Summary measure
#' @return Multiverse results
#' @export
run_multiverse <- function(data, ref_treatment, sm = "HR") {
  msg("Running multiverse sensitivity analysis...")

  grid <- tibble::tribble(
    ~scenario, ~random, ~fixed,
    "Random effects only", TRUE, FALSE,
    "Fixed effects only", FALSE, TRUE,
    "Both models", TRUE, TRUE
  )

  results_list <- lapply(seq_len(nrow(grid)), function(i) {
    row <- grid[i, ]

    fit <- .safe_try(
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = data,
        sm = sm,
        fixed = row$fixed,
        random = row$random,
        reference.group = ref_treatment
      ),
      paste0("Multiverse: ", row$scenario)
    )

    if (inherits(fit, "try-error")) {
      tibble::tibble(
        scenario = row$scenario,
        tau = NA_real_,
        I2 = NA_real_
      )
    } else {
      tibble::tibble(
        scenario = row$scenario,
        tau = ifelse(row$random, fit$tau, NA_real_),
        I2 = ifelse(row$random, fit$I2.random, NA_real_)
      )
    }
  })

  dplyr::bind_rows(results_list)
}

#' Transportability Diagnostics
#'
#' Assess effective sample size and weight distribution
#'
#' @param weights Named vector of study weights
#' @return Diagnostic statistics
#' @export
transportability_diagnostics <- function(data, target_population, weights,
                                         covariates = NULL) {
  # Extract weights if it's a tibble from compute_transport_weights()
  if (is.data.frame(weights) && "weight" %in% names(weights)) {
    weight_vec <- weights$weight
    study_labels <- weights$studlab
  } else {
    weight_vec <- as.numeric(weights)
    study_labels <- if (!is.null(names(weights))) names(weights) else seq_along(weights)
  }

  # Identify covariates for balance assessment
  if (is.null(covariates)) {
    covariates <- intersect(names(target_population), names(data))
  }

  # 1. EFFECTIVE SAMPLE SIZE
  ess <- (sum(weight_vec))^2 / sum(weight_vec^2)
  ess_ratio <- ess / length(weight_vec)

  # 2. COVARIATE BALANCE (SMD before and after weighting)
  balance <- NULL
  if (length(covariates) > 0) {
    study_df <- data %>% dplyr::distinct(studlab, .keep_all = TRUE)

    # Match weights to studies
    study_df$weight <- weight_vec[match(study_df$studlab, study_labels)]

    balance <- purrr::map_dfr(covariates, function(cov) {
      study_vals <- study_df[[cov]]
      target_val <- target_population[[cov]]

      if (is.numeric(study_vals) && is.numeric(target_val)) {
        # Standardized Mean Difference (SMD)
        # SMD = (mean_study - mean_target) / pooled_sd

        # Unweighted
        mean_unweighted <- mean(study_vals, na.rm = TRUE)
        sd_unweighted <- stats::sd(study_vals, na.rm = TRUE)
        smd_unweighted <- (mean_unweighted - target_val) / (sd_unweighted + 1e-8)

        # Weighted
        mean_weighted <- stats::weighted.mean(study_vals, w = study_df$weight, na.rm = TRUE)
        # Weighted SD
        var_weighted <- stats::weighted.mean((study_vals - mean_weighted)^2,
                                            w = study_df$weight, na.rm = TRUE)
        sd_weighted <- sqrt(var_weighted)
        smd_weighted <- (mean_weighted - target_val) / (sd_weighted + 1e-8)

        improvement <- abs(smd_unweighted) - abs(smd_weighted)

        tibble::tibble(
          covariate = cov,
          target_value = target_val,
          mean_unweighted = mean_unweighted,
          mean_weighted = mean_weighted,
          smd_unweighted = smd_unweighted,
          smd_weighted = smd_weighted,
          improvement = improvement,
          balanced = abs(smd_weighted) < 0.1
        )
      } else {
        # Categorical variable - skip for now
        tibble::tibble(
          covariate = cov,
          target_value = as.character(target_val),
          mean_unweighted = NA,
          mean_weighted = NA,
          smd_unweighted = NA,
          smd_weighted = NA,
          improvement = NA,
          balanced = NA
        )
      }
    })
  }

  # 3. POSITIVITY CHECK (convex hull for continuous covariates)
  positivity <- list(all_within_hull = NA, message = "")

  if (length(covariates) >= 2) {
    numeric_covs <- covariates[sapply(covariates, function(c) is.numeric(data[[c]]))]

    if (length(numeric_covs) >= 2) {
      study_df <- data %>% dplyr::distinct(studlab, .keep_all = TRUE)
      study_matrix <- as.matrix(study_df[, numeric_covs[1:min(2, length(numeric_covs))]])
      target_vec <- unlist(target_population[numeric_covs[1:min(2, length(numeric_covs))]])

      # Simple positivity check: are target values within min/max range?
      within_range <- all(sapply(seq_along(target_vec), function(i) {
        val <- target_vec[i]
        col_vals <- study_matrix[, i]
        val >= min(col_vals, na.rm = TRUE) && val <= max(col_vals, na.rm = TRUE)
      }))

      positivity$all_within_hull <- within_range
      positivity$message <- if (within_range) {
        "Target population characteristics within study data range (positivity likely satisfied)"
      } else {
        "WARNING: Target population outside study data range (positivity violated)"
      }
    }
  }

  # 4. WEIGHT DISTRIBUTION SUMMARY
  weight_summary <- list(
    min = min(weight_vec),
    q25 = stats::quantile(weight_vec, 0.25),
    median = stats::median(weight_vec),
    q75 = stats::quantile(weight_vec, 0.75),
    max = max(weight_vec),
    mean = mean(weight_vec),
    sd = stats::sd(weight_vec),
    cv = stats::sd(weight_vec) / mean(weight_vec)
  )

  # Count extreme weights
  extreme_low <- sum(weight_vec < stats::quantile(weight_vec, 0.05))
  extreme_high <- sum(weight_vec > stats::quantile(weight_vec, 0.95))

  # 5. GENERATE RECOMMENDATION
  recommendation <- generate_transportability_recommendation(
    ess_ratio = ess_ratio,
    balance = balance,
    positivity = positivity,
    weight_summary = weight_summary
  )

  # Return comprehensive diagnostics
  structure(
    list(
      effective_sample_size = ess,
      ess_ratio = ess_ratio,
      n_studies = length(weight_vec),
      balance_table = balance,
      positivity = positivity,
      weight_distribution = weight_summary,
      extreme_weights = list(low = extreme_low, high = extreme_high),
      recommendation = recommendation
    ),
    class = c("transportability_diagnostics", "list")
  )
}

#' Generate Transportability Recommendation
#' @keywords internal
generate_transportability_recommendation <- function(ess_ratio, balance,
                                                    positivity, weight_summary) {
  issues <- character()
  severity <- "OK"

  # Check ESS
  if (ess_ratio < 0.5) {
    issues <- c(issues,
      sprintf("Low effective sample size (%.1f%% of original)", ess_ratio * 100))
    severity <- "CAUTION"
  }

  if (ess_ratio < 0.3) {
    severity <- "WARNING"
  }

  # Check balance
  if (!is.null(balance) && nrow(balance) > 0) {
    unbalanced <- balance %>%
      dplyr::filter(!is.na(balanced), !balanced)

    if (nrow(unbalanced) > 0) {
      issues <- c(issues,
        sprintf("%d covariate(s) still imbalanced after weighting (SMD > 0.1)",
                nrow(unbalanced)))
      severity <- ifelse(severity == "WARNING", "WARNING", "CAUTION")
    }
  }

  # Check positivity
  if (!is.na(positivity$all_within_hull) && !positivity$all_within_hull) {
    issues <- c(issues,
      "Target population outside study covariate range (positivity violation)")
    severity <- "ERROR"
  }

  # Check weight variability
  if (weight_summary$cv > 1.5) {
    issues <- c(issues,
      sprintf("High weight variability (CV = %.2f)", weight_summary$cv))
    if (severity == "OK") severity <- "CAUTION"
  }

  # Generate recommendation text
  if (length(issues) == 0) {
    rec_text <- "Transportability weights appear reasonable. Proceed with caution and report sensitivity analyses."
  } else {
    rec_text <- paste0(
      severity, ": ",
      paste(issues, collapse = "; "),
      ". ",
      switch(severity,
        "CAUTION" = "Consider sensitivity analyses.",
        "WARNING" = "Use with extreme caution. Consider additional covariates or alternative methods.",
        "ERROR" = "Transportability not valid. Do not use these weights.",
        ""
      )
    )
  }

  list(
    severity = severity,
    issues = issues,
    text = rec_text
  )
}

#' Print method for transportability diagnostics
#' @export
print.transportability_diagnostics <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  Transportability Diagnostics\n")
  cat("═══════════════════════════════════════════════════\n\n")

  cat(sprintf("Number of studies: %d\n", x$n_studies))
  cat(sprintf("Effective sample size: %.1f (%.1f%% of original)\n",
              x$effective_sample_size, x$ess_ratio * 100))
  cat("\n")

  cat("Weight distribution:\n")
  cat(sprintf("  Min:    %.3f\n", x$weight_distribution$min))
  cat(sprintf("  Q1:     %.3f\n", x$weight_distribution$q25))
  cat(sprintf("  Median: %.3f\n", x$weight_distribution$median))
  cat(sprintf("  Q3:     %.3f\n", x$weight_distribution$q75))
  cat(sprintf("  Max:    %.3f\n", x$weight_distribution$max))
  cat(sprintf("  CV:     %.2f\n", x$weight_distribution$cv))
  cat("\n")

  if (!is.null(x$balance_table) && nrow(x$balance_table) > 0) {
    cat("Covariate balance (SMD):\n")
    print(x$balance_table %>%
            dplyr::select(covariate, smd_unweighted, smd_weighted, balanced), n = Inf)
    cat("\n")
  }

  cat("Positivity:\n")
  cat(sprintf("  %s\n", x$positivity$message))
  cat("\n")

  cat("RECOMMENDATION:\n")
  cat(sprintf("  %s\n", x$recommendation$text))
  cat("\n")

  invisible(x)
}

#' Component Network Meta-Analysis
#'
#' Analyze additive components of combination treatments
#'
#' @param nma_fit netmeta object
#' @return Component NMA results or NULL
#' @export
run_component_nma <- function(nma_fit) {
  trts <- rownames(nma_fit$TE.random %||% nma_fit$TE.fixed)

  if (!any(grepl("\\+", trts))) {
    msg("No combination treatments detected (e.g., 'A+B')")
    return(NULL)
  }

  if (!has_pkg("netmeta")) {
    msg("netmeta required for component NMA")
    return(NULL)
  }

  out <- .safe_try(
    netmeta::netcomb(nma_fit),
    "Component NMA"
  )

  if (inherits(out, "try-error")) {
    msg("Component NMA failed")
    return(NULL)
  }

  out
}
