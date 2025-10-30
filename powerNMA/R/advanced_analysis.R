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
transportability_diagnostics <- function(weights) {
  w <- as.numeric(weights)

  tibble::tibble(
    effective_sample_size = sum(w)^2 / sum(w^2),
    min_weight = min(w),
    max_weight = max(w),
    mean_weight = mean(w),
    sd_weight = stats::sd(w),
    cv_weight = stats::sd(w) / mean(w)
  )
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
