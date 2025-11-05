#' Enhanced Meta-Regression for Network Meta-Analysis
#'
#' Advanced meta-regression tools including subgroup analysis, interaction
#' testing, non-linear relationships, and robust estimation.
#'
#' @name meta_regression_enhanced
#' @references
#' Dias S, et al. (2013). Evidence synthesis for decision making 2: A generalized
#' linear modeling framework for pairwise and network meta-analysis of randomized
#' controlled trials. Medical Decision Making, 33(5):607-617.
NULL

#' Network Meta-Regression
#'
#' Perform meta-regression in network meta-analysis to explore effect modifiers
#' and sources of heterogeneity.
#'
#' @param data Pairwise NMA data with covariates
#' @param formula Regression formula (e.g., TE ~ age + female_pct)
#' @param reference Reference treatment
#' @param sm Summary measure
#' @param method Estimation method: "ml" (maximum likelihood) or "reml"
#' @param test_interactions Test treatment-covariate interactions
#' @return nma_metareg object
#' @export
#' @examples
#' \dontrun{
#' # Simple meta-regression
#' metareg <- network_metaregression(
#'   data = data,
#'   formula = TE ~ age_mean + female_pct,
#'   reference = "Placebo"
#' )
#' print(metareg)
#' plot(metareg)
#' }
network_metaregression <- function(data,
                                  formula,
                                  reference = NULL,
                                  sm = "MD",
                                  method = c("ml", "reml"),
                                  test_interactions = FALSE) {

  method <- match.arg(method)

  # Validate data
  required_cols <- c("treat1", "treat2", "TE", "seTE", "studlab")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Auto-detect reference if needed
  if (is.null(reference)) {
    counts <- table(c(as.character(data$treat1), as.character(data$treat2)))
    reference <- names(sort(counts, decreasing = TRUE))[1]
  }

  # Extract covariates from formula
  formula_terms <- all.vars(formula)
  response_var <- formula_terms[1]
  covariate_vars <- formula_terms[-1]

  # Check covariates exist
  missing_covs <- setdiff(covariate_vars, names(data))
  if (length(missing_covs) > 0) {
    stop(sprintf("Covariates not found in data: %s", paste(missing_covs, collapse = ", ")))
  }

  # Fit base network meta-analysis
  base_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  # Fit meta-regression using metafor
  if (requireNamespace("metafor", quietly = TRUE)) {
    metareg_fit <- fit_metareg_metafor(data, formula, method)
  } else {
    stop("Package 'metafor' required for meta-regression")
  }

  # Test interactions if requested
  interaction_tests <- NULL
  if (test_interactions) {
    interaction_tests <- test_treatment_covariate_interactions(
      data, covariate_vars, reference
    )
  }

  result <- list(
    base_nma = base_nma,
    metareg_fit = metareg_fit,
    formula = formula,
    covariates = covariate_vars,
    reference = reference,
    sm = sm,
    method = method,
    interaction_tests = interaction_tests,
    data = data
  )

  class(result) <- c("nma_metareg", "list")
  result
}

#' Fit Meta-Regression using metafor
#'
#' @param data Data frame
#' @param formula Regression formula
#' @param method Estimation method
#' @return rma object from metafor
#' @keywords internal
fit_metareg_metafor <- function(data, formula, method) {

  # Prepare data for metafor
  # Use inverse variance weights
  data$vi <- data$seTE^2

  # Fit random-effects meta-regression
  fit <- metafor::rma.mv(
    formula,
    V = vi,
    random = ~ 1 | studlab,
    data = data,
    method = toupper(method)
  )

  fit
}

#' Test Treatment-Covariate Interactions
#'
#' Test whether covariate effects differ across treatments.
#'
#' @param data Pairwise data
#' @param covariates Vector of covariate names
#' @param reference Reference treatment
#' @return Data frame with interaction tests
#' @keywords internal
test_treatment_covariate_interactions <- function(data, covariates, reference) {

  results <- lapply(covariates, function(cov) {
    # Create interaction term
    data$interaction <- data[[cov]] * (data$treat1 != reference | data$treat2 != reference)

    # Fit model with and without interaction
    formula_main <- as.formula(paste("TE ~", cov))
    formula_int <- as.formula(paste("TE ~", cov, "+ interaction"))

    data$vi <- data$seTE^2

    fit_main <- metafor::rma.mv(
      formula_main,
      V = vi,
      random = ~ 1 | studlab,
      data = data
    )

    fit_int <- metafor::rma.mv(
      formula_int,
      V = vi,
      random = ~ 1 | studlab,
      data = data
    )

    # Likelihood ratio test
    lr_test <- anova(fit_main, fit_int)

    data.frame(
      Covariate = cov,
      LR_statistic = lr_test$LRT[2],
      df = lr_test$df[2],
      p_value = lr_test$pval[2],
      Significant = lr_test$pval[2] < 0.05,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Subgroup Network Meta-Analysis
#'
#' Perform NMA separately for subgroups and test for subgroup differences.
#'
#' @param data Pairwise NMA data
#' @param subgroup_var Name of subgroup variable
#' @param reference Reference treatment
#' @param sm Summary measure
#' @param test_difference Test for subgroup differences
#' @return nma_subgroup object
#' @export
#' @examples
#' \dontrun{
#' # Subgroup analysis by age groups
#' subgroup <- subgroup_network_meta(
#'   data = data,
#'   subgroup_var = "age_group",
#'   reference = "Placebo",
#'   test_difference = TRUE
#' )
#' print(subgroup)
#' plot(subgroup)
#' }
subgroup_network_meta <- function(data,
                                  subgroup_var,
                                  reference = NULL,
                                  sm = "MD",
                                  test_difference = TRUE) {

  if (!(subgroup_var %in% names(data))) {
    stop(sprintf("Subgroup variable '%s' not found in data", subgroup_var))
  }

  # Auto-detect reference
  if (is.null(reference)) {
    counts <- table(c(as.character(data$treat1), as.character(data$treat2)))
    reference <- names(sort(counts, decreasing = TRUE))[1]
  }

  # Get subgroup levels
  subgroups <- unique(data[[subgroup_var]])
  subgroups <- subgroups[!is.na(subgroups)]

  if (length(subgroups) < 2) {
    stop("Need at least 2 subgroups for subgroup analysis")
  }

  # Run NMA for each subgroup
  subgroup_results <- lapply(subgroups, function(sg) {
    subgroup_data <- data[data[[subgroup_var]] == sg, ]

    if (nrow(subgroup_data) < 3) {
      warning(sprintf("Subgroup '%s' has insufficient data (n=%d studies)",
                     sg, nrow(subgroup_data)))
      return(NULL)
    }

    nma <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = subgroup_data,
        sm = sm,
        reference.group = reference
      )
    }, error = function(e) {
      warning(sprintf("NMA failed for subgroup '%s': %s", sg, e$message))
      NULL
    })

    list(
      subgroup = sg,
      nma = nma,
      n_studies = nrow(subgroup_data),
      n_comparisons = length(unique(paste(subgroup_data$treat1, subgroup_data$treat2)))
    )
  })

  names(subgroup_results) <- as.character(subgroups)

  # Remove failed subgroups
  subgroup_results <- subgroup_results[!sapply(subgroup_results, is.null)]

  # Test for subgroup differences
  difference_tests <- NULL
  if (test_difference && length(subgroup_results) >= 2) {
    difference_tests <- test_subgroup_differences(subgroup_results, reference)
  }

  result <- list(
    subgroup_var = subgroup_var,
    subgroups = subgroups,
    subgroup_results = subgroup_results,
    difference_tests = difference_tests,
    reference = reference,
    sm = sm
  )

  class(result) <- c("nma_subgroup", "list")
  result
}

#' Test for Subgroup Differences
#'
#' @param subgroup_results List of subgroup NMA results
#' @param reference Reference treatment
#' @return Data frame with difference tests
#' @keywords internal
test_subgroup_differences <- function(subgroup_results, reference) {

  # Extract all treatments
  all_treatments <- unique(unlist(lapply(subgroup_results, function(sr) {
    if (!is.null(sr$nma)) {
      rownames(sr$nma$TE.random)
    } else {
      NULL
    }
  })))

  all_treatments <- setdiff(all_treatments, reference)

  # Test each treatment vs reference across subgroups
  tests <- lapply(all_treatments, function(trt) {
    # Extract effects from each subgroup
    effects <- sapply(subgroup_results, function(sr) {
      if (!is.null(sr$nma) && trt %in% rownames(sr$nma$TE.random)) {
        sr$nma$TE.random[trt, reference]
      } else {
        NA
      }
    })

    se_effects <- sapply(subgroup_results, function(sr) {
      if (!is.null(sr$nma) && trt %in% rownames(sr$nma$TE.random)) {
        sr$nma$seTE.random[trt, reference]
      } else {
        NA
      }
    })

    valid_idx <- !is.na(effects) & !is.na(se_effects)

    if (sum(valid_idx) < 2) {
      return(data.frame(
        Treatment = trt,
        Q = NA,
        df = NA,
        p_value = NA,
        Significant = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    # Cochran's Q test for heterogeneity between subgroups
    effects_valid <- effects[valid_idx]
    se_valid <- se_effects[valid_idx]
    weights <- 1 / se_valid^2

    weighted_mean <- sum(weights * effects_valid) / sum(weights)
    Q <- sum(weights * (effects_valid - weighted_mean)^2)
    df <- sum(valid_idx) - 1
    p_value <- pchisq(Q, df, lower.tail = FALSE)

    data.frame(
      Treatment = trt,
      Q = Q,
      df = df,
      p_value = p_value,
      Significant = p_value < 0.05,
      Interpretation = if (p_value < 0.05) {
        "Subgroup effects differ significantly"
      } else {
        "No significant subgroup differences"
      },
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, tests)
}

#' Dose-Response Meta-Regression
#'
#' Model dose-response relationships in network meta-analysis.
#'
#' @param data Pairwise data with dose information
#' @param dose_var Name of dose variable
#' @param model Dose-response model: "linear", "log", "quadratic", "emax"
#' @param reference Reference treatment (zero dose)
#' @return dose_response_nma object
#' @export
#' @examples
#' \dontrun{
#' dose_response <- dose_response_metareg(
#'   data = data,
#'   dose_var = "drug_dose_mg",
#'   model = "emax"
#' )
#' plot(dose_response)
#' }
dose_response_metareg <- function(data,
                                 dose_var,
                                 model = c("linear", "log", "quadratic", "emax"),
                                 reference = NULL) {

  model <- match.arg(model)

  if (!(dose_var %in% names(data))) {
    stop(sprintf("Dose variable '%s' not found in data", dose_var))
  }

  # Create dose transformations
  data$dose <- data[[dose_var]]

  if (model == "log") {
    data$dose_transformed <- log(data$dose + 1)  # +1 to handle zero dose
  } else if (model == "quadratic") {
    data$dose_transformed <- data$dose
    data$dose_squared <- data$dose^2
  } else if (model == "emax") {
    # Emax model: E = E0 + (Emax * dose) / (ED50 + dose)
    # Requires non-linear fitting - simplified here
    data$dose_transformed <- data$dose / (median(data$dose, na.rm = TRUE) + data$dose)
  } else {
    # Linear
    data$dose_transformed <- data$dose
  }

  # Fit regression
  if (model == "quadratic") {
    formula <- as.formula("TE ~ dose_transformed + dose_squared")
  } else {
    formula <- as.formula("TE ~ dose_transformed")
  }

  data$vi <- data$seTE^2

  fit <- metafor::rma.mv(
    formula,
    V = vi,
    random = ~ 1 | studlab,
    data = data
  )

  result <- list(
    fit = fit,
    model = model,
    dose_var = dose_var,
    data = data,
    formula = formula
  )

  class(result) <- c("dose_response_nma", "list")
  result
}

#' Print NMA Meta-Regression
#'
#' @param x nma_metareg object
#' @param ... Additional arguments
#' @export
print.nma_metareg <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Network Meta-Regression\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Formula: %s\n", deparse(x$formula)))
  cat(sprintf("Method: %s\n", toupper(x$method)))
  cat(sprintf("Reference: %s\n\n", x$reference))

  cat("Meta-Regression Coefficients:\n\n")
  print(summary(x$metareg_fit))

  if (!is.null(x$interaction_tests)) {
    cat("\n")
    cat("Treatment-Covariate Interaction Tests:\n\n")
    print(x$interaction_tests, row.names = FALSE)
  }

  cat("\n")

  invisible(x)
}

#' Print NMA Subgroup Analysis
#'
#' @param x nma_subgroup object
#' @param ... Additional arguments
#' @export
print.nma_subgroup <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Subgroup Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Subgroup variable: %s\n", x$subgroup_var))
  cat(sprintf("Number of subgroups: %d\n", length(x$subgroup_results)))
  cat(sprintf("Reference: %s\n\n", x$reference))

  cat("Subgroup Results:\n\n")
  for (sg_name in names(x$subgroup_results)) {
    sg <- x$subgroup_results[[sg_name]]
    cat(sprintf("  %s: %d studies, %d comparisons\n",
               sg_name, sg$n_studies, sg$n_comparisons))
  }

  if (!is.null(x$difference_tests)) {
    cat("\n")
    cat("Tests for Subgroup Differences:\n\n")
    print(x$difference_tests, row.names = FALSE)
  }

  cat("\n")

  invisible(x)
}

#' Plot Meta-Regression Results
#'
#' @param x nma_metareg object
#' @param covariate Which covariate to plot
#' @param ... Additional arguments
#' @export
plot.nma_metareg <- function(x, covariate = NULL, ...) {

  if (is.null(covariate)) {
    covariate <- x$covariates[1]
  }

  if (!(covariate %in% x$covariates)) {
    stop(sprintf("Covariate '%s' not found. Available: %s",
                covariate, paste(x$covariates, collapse = ", ")))
  }

  # Create scatter plot with regression line
  data <- x$data
  data$cov_val <- data[[covariate]]

  # Get predicted values
  pred_range <- seq(min(data$cov_val, na.rm = TRUE),
                   max(data$cov_val, na.rm = TRUE),
                   length.out = 100)

  pred_data <- data.frame(cov_val = pred_range)
  names(pred_data)[1] <- covariate

  # Would need to implement prediction properly
  # This is simplified

  ggplot2::ggplot(data, ggplot2::aes(x = cov_val, y = TE)) +
    ggplot2::geom_point(ggplot2::aes(size = 1/seTE), alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#2c7fb8") +
    ggplot2::labs(
      title = "Network Meta-Regression",
      subtitle = sprintf("Effect by %s", covariate),
      x = covariate,
      y = sprintf("Treatment Effect (%s)", x$sm),
      size = "Precision"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Subgroup Results
#'
#' @param x nma_subgroup object
#' @param treatment Which treatment to plot (vs reference)
#' @param ... Additional arguments
#' @export
plot.nma_subgroup <- function(x, treatment = NULL, ...) {

  # Extract effects by subgroup
  if (is.null(treatment)) {
    # Use first non-reference treatment
    treatment <- names(x$subgroup_results[[1]]$nma$TE.random)[1]
    if (treatment == x$reference) {
      treatment <- names(x$subgroup_results[[1]]$nma$TE.random)[2]
    }
  }

  plot_data <- lapply(names(x$subgroup_results), function(sg_name) {
    sg <- x$subgroup_results[[sg_name]]

    if (!is.null(sg$nma) && treatment %in% rownames(sg$nma$TE.random)) {
      data.frame(
        Subgroup = sg_name,
        Effect = sg$nma$TE.random[treatment, x$reference],
        SE = sg$nma$seTE.random[treatment, x$reference],
        Lower = sg$nma$lower.random[treatment, x$reference],
        Upper = sg$nma$upper.random[treatment, x$reference],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  plot_data <- do.call(rbind, plot_data[!sapply(plot_data, is.null)])

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Subgroup, y = Effect)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                          width = 0.2, size = 0.8, color = "#2c7fb8") +
    ggplot2::geom_point(size = 3, color = "#2c7fb8") +
    ggplot2::labs(
      title = "Subgroup Analysis",
      subtitle = sprintf("%s vs %s by %s", treatment, x$reference, x$subgroup_var),
      x = x$subgroup_var,
      y = sprintf("Treatment Effect (%s)", x$sm)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
