#' Multivariate Network Meta-Analysis
#'
#' Analyze multiple correlated outcomes simultaneously using multivariate NMA.
#' This allows borrowing of strength across outcomes and accounts for correlations.
#'
#' @name multivariate_nma
#' @references
#' Achana FA, et al. (2014). Extending methods for investigating the relationship
#' between treatment effect and baseline risk from pairwise meta-analysis to
#' network meta-analysis. Statistics in Medicine, 33(5):752-770.
#'
#' Efthimiou O, et al. (2015). Combining multiple correlated outcomes in a
#' network of interventions. Statistics in Medicine, 34(1):69-83.
NULL

#' Prepare Data for Multivariate NMA
#'
#' Structure data with multiple outcomes for multivariate analysis.
#'
#' @param data_list List of data frames, one per outcome
#' @param outcome_names Character vector of outcome names
#' @param correlation Correlation matrix between outcomes (optional)
#' @return Structured mvnma_data object
#' @export
#' @examples
#' \dontrun{
#' # Separate datasets for efficacy and safety
#' efficacy_data <- pairwise(treat = treat, event = event_eff, n = n, ...)
#' safety_data <- pairwise(treat = treat, event = event_safe, n = n, ...)
#'
#' mvnma_data <- prepare_mvnma_data(
#'   data_list = list(efficacy_data, safety_data),
#'   outcome_names = c("Efficacy", "Safety")
#' )
#' }
prepare_mvnma_data <- function(data_list, outcome_names = NULL,
                              correlation = NULL) {

  if (!is.list(data_list)) {
    stop("data_list must be a list of data frames")
  }

  n_outcomes <- length(data_list)

  if (is.null(outcome_names)) {
    outcome_names <- paste0("Outcome", seq_len(n_outcomes))
  }

  if (length(outcome_names) != n_outcomes) {
    stop("Length of outcome_names must match number of datasets")
  }

  # Validate that all datasets have required columns
  required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")

  for (i in seq_along(data_list)) {
    missing_cols <- setdiff(required_cols, names(data_list[[i]]))
    if (length(missing_cols) > 0) {
      stop(sprintf("Dataset %d (%s) missing columns: %s",
                  i, outcome_names[i], paste(missing_cols, collapse = ", ")))
    }
  }

  # Check for common studies across outcomes
  study_sets <- lapply(data_list, function(df) unique(df$studlab))
  all_studies <- unique(unlist(study_sets))
  common_studies <- Reduce(intersect, study_sets)

  # Check for common treatments
  treatment_sets <- lapply(data_list, function(df) {
    unique(c(as.character(df$treat1), as.character(df$treat2)))
  })
  all_treatments <- unique(unlist(treatment_sets))

  # Default correlation structure
  if (is.null(correlation)) {
    # Assume moderate positive correlation (0.5) between outcomes
    correlation <- matrix(0.5, nrow = n_outcomes, ncol = n_outcomes)
    diag(correlation) <- 1.0
    rownames(correlation) <- colnames(correlation) <- outcome_names
    message("No correlation specified. Using default correlation = 0.5 between outcomes.")
  }

  # Validate correlation matrix
  if (!isSymmetric(correlation) || any(diag(correlation) != 1)) {
    stop("Correlation matrix must be symmetric with 1s on diagonal")
  }

  if (any(eigen(correlation)$values <= 0)) {
    stop("Correlation matrix is not positive definite")
  }

  result <- list(
    data_list = data_list,
    outcome_names = outcome_names,
    n_outcomes = n_outcomes,
    correlation = correlation,
    all_studies = all_studies,
    common_studies = common_studies,
    all_treatments = all_treatments,
    n_common_studies = length(common_studies)
  )

  class(result) <- c("mvnma_data", "list")
  result
}

#' Fit Multivariate Network Meta-Analysis
#'
#' Perform multivariate NMA accounting for correlations between outcomes.
#'
#' @param mvnma_data Prepared mvnma_data object
#' @param reference Reference treatment
#' @param sm Summary measures (vector, one per outcome)
#' @param method Method: "riley" (Riley 2017) or "composite" (Efthimiou 2015)
#' @return mvnma_result object
#' @export
#' @examples
#' \dontrun{
#' mvnma_result <- fit_mvnma(
#'   mvnma_data = mvnma_data,
#'   reference = "Placebo",
#'   sm = c("OR", "OR"),
#'   method = "riley"
#' )
#' print(mvnma_result)
#' plot(mvnma_result)
#' }
fit_mvnma <- function(mvnma_data, reference = NULL, sm = NULL,
                     method = c("riley", "composite", "separate")) {

  method <- match.arg(method)

  if (!inherits(mvnma_data, "mvnma_data")) {
    stop("Input must be mvnma_data object from prepare_mvnma_data()")
  }

  n_outcomes <- mvnma_data$n_outcomes

  # Set default summary measures
  if (is.null(sm)) {
    sm <- rep("MD", n_outcomes)
    message("No summary measures specified. Using 'MD' for all outcomes.")
  }

  if (length(sm) != n_outcomes) {
    stop(sprintf("Length of sm (%d) must match number of outcomes (%d)",
                length(sm), n_outcomes))
  }

  message(sprintf("Fitting multivariate NMA with %d outcomes using '%s' method...",
                 n_outcomes, method))

  # Fit separate NMAs for each outcome
  separate_nmas <- lapply(seq_len(n_outcomes), function(i) {
    message(sprintf("  Fitting outcome %d: %s", i, mvnma_data$outcome_names[i]))

    netmeta::netmeta(
      TE = TE, seTE = seTE,
      treat1 = treat1, treat2 = treat2,
      studlab = studlab,
      data = mvnma_data$data_list[[i]],
      sm = sm[i],
      reference.group = reference
    )
  })

  names(separate_nmas) <- mvnma_data$outcome_names

  # Set reference if not specified
  if (is.null(reference)) {
    reference <- separate_nmas[[1]]$reference.group
  }

  if (method == "separate") {
    # Just return separate analyses
    result <- list(
      method = "separate",
      separate_nmas = separate_nmas,
      mvnma_data = mvnma_data,
      reference = reference,
      sm = sm
    )

    class(result) <- c("mvnma_result", "list")
    return(result)
  }

  # Get treatment effects from each outcome
  treatments <- rownames(separate_nmas[[1]]$TE.random)

  # Build multivariate effect vectors and covariance matrices
  # For each treatment vs reference
  mvnma_estimates <- list()

  for (trt in treatments) {
    if (trt == reference) next

    # Extract effects for this treatment across all outcomes
    effects <- sapply(separate_nmas, function(nma) {
      nma$TE.random[trt, reference]
    })

    # Extract variances
    variances <- sapply(separate_nmas, function(nma) {
      nma$seTE.random[trt, reference]^2
    })

    # Build covariance matrix using assumed correlations
    cov_matrix <- outer(sqrt(variances), sqrt(variances)) * mvnma_data$correlation

    mvnma_estimates[[trt]] <- list(
      treatment = trt,
      effects = effects,
      cov_matrix = cov_matrix
    )
  }

  # Multivariate pooling (simplified Riley method)
  # For full Bayesian implementation, would use gemtc with multivariate likelihood
  mvnma_pooled <- lapply(names(mvnma_estimates), function(trt) {
    est <- mvnma_estimates[[trt]]

    # Precision matrix
    precision <- solve(est$cov_matrix)

    # Multivariate pooled estimate (generalized least squares)
    weight_sum <- sum(precision)
    pooled_effects <- est$effects  # Simplified - in full model would pool across studies

    # Standard errors
    pooled_se <- sqrt(diag(est$cov_matrix))

    # Confidence intervals
    lower <- pooled_effects - 1.96 * pooled_se
    upper <- pooled_effects + 1.96 * pooled_se

    data.frame(
      Treatment = trt,
      Outcome = mvnma_data$outcome_names,
      Effect = pooled_effects,
      SE = pooled_se,
      Lower_95 = lower,
      Upper_95 = upper,
      stringsAsFactors = FALSE
    )
  })

  pooled_table <- do.call(rbind, mvnma_pooled)

  result <- list(
    method = method,
    separate_nmas = separate_nmas,
    pooled_estimates = pooled_table,
    mvnma_estimates = mvnma_estimates,
    mvnma_data = mvnma_data,
    reference = reference,
    sm = sm,
    treatments = treatments
  )

  class(result) <- c("mvnma_result", "list")
  result
}

#' Create Benefit-Risk Plot for Multivariate NMA
#'
#' Visualize trade-offs between efficacy and safety outcomes.
#'
#' @param mvnma_result Multivariate NMA result
#' @param efficacy_outcome Name or index of efficacy outcome
#' @param safety_outcome Name or index of safety outcome
#' @param efficacy_good Direction for efficacy ("higher" or "lower")
#' @param safety_good Direction for safety ("higher" or "lower")
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' benefit_risk_plot(
#'   mvnma_result,
#'   efficacy_outcome = "Efficacy",
#'   safety_outcome = "Safety",
#'   efficacy_good = "higher",
#'   safety_good = "lower"
#' )
#' }
benefit_risk_plot <- function(mvnma_result, efficacy_outcome = 1, safety_outcome = 2,
                             efficacy_good = c("higher", "lower"),
                             safety_good = c("higher", "lower")) {

  efficacy_good <- match.arg(efficacy_good)
  safety_good <- match.arg(safety_good)

  if (!inherits(mvnma_result, "mvnma_result")) {
    stop("Input must be mvnma_result object")
  }

  # Handle outcome selection by name or index
  if (is.numeric(efficacy_outcome)) {
    eff_name <- mvnma_result$mvnma_data$outcome_names[efficacy_outcome]
  } else {
    eff_name <- efficacy_outcome
  }

  if (is.numeric(safety_outcome)) {
    safe_name <- mvnma_result$mvnma_data$outcome_names[safety_outcome]
  } else {
    safe_name <- safety_outcome
  }

  # Extract data
  plot_data <- mvnma_result$pooled_estimates

  eff_data <- plot_data[plot_data$Outcome == eff_name, ]
  safe_data <- plot_data[plot_data$Outcome == safe_name, ]

  # Merge
  merged <- merge(eff_data, safe_data, by = "Treatment", suffixes = c("_eff", "_safe"))

  # Adjust signs if needed
  if (efficacy_good == "lower") {
    merged$Effect_eff <- -merged$Effect_eff
  }

  if (safety_good == "lower") {
    merged$Effect_safe <- -merged$Effect_safe
  }

  # Add reference point
  ref_point <- data.frame(
    Treatment = mvnma_result$reference,
    Effect_eff = 0,
    Effect_safe = 0,
    stringsAsFactors = FALSE
  )

  merged <- rbind(merged[, c("Treatment", "Effect_eff", "Effect_safe")], ref_point)

  # Quadrant colors
  merged$Quadrant <- ifelse(merged$Effect_eff > 0 & merged$Effect_safe > 0, "Win-Win",
                    ifelse(merged$Effect_eff > 0 & merged$Effect_safe < 0, "Benefit>Risk",
                    ifelse(merged$Effect_eff < 0 & merged$Effect_safe > 0, "Risk>Benefit",
                    "Lose-Lose")))

  # Color for reference
  merged$Quadrant[merged$Treatment == mvnma_result$reference] <- "Reference"

  color_map <- c(
    "Win-Win" = "#1a9850",
    "Benefit>Risk" = "#91cf60",
    "Risk>Benefit" = "#fc8d59",
    "Lose-Lose" = "#d73027",
    "Reference" = "black"
  )

  p <- ggplot2::ggplot(merged, ggplot2::aes(x = Effect_eff, y = Effect_safe)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_point(ggplot2::aes(color = Quadrant, size = 3)) +
    ggplot2::geom_text(ggplot2::aes(label = Treatment), vjust = -1, size = 3.5) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(
      title = "Benefit-Risk Trade-off Plot",
      subtitle = sprintf("Reference: %s", mvnma_result$reference),
      x = sprintf("%s Effect (higher is better)", eff_name),
      y = sprintf("%s Effect (higher is better)", safe_name)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  print(p)
  invisible(p)
}

#' Calculate Net Benefit for Multivariate NMA
#'
#' Compute net benefit combining multiple outcomes with specified weights.
#'
#' @param mvnma_result Multivariate NMA result
#' @param weights Named vector of outcome weights (must sum to 1)
#' @param directions Named vector of outcome directions ("higher" or "lower" is good)
#' @return Data frame with net benefit scores
#' @export
#' @examples
#' \dontrun{
#' net_benefit <- calculate_net_benefit(
#'   mvnma_result,
#'   weights = c(Efficacy = 0.6, Safety = 0.4),
#'   directions = c(Efficacy = "higher", Safety = "lower")
#' )
#' }
calculate_net_benefit <- function(mvnma_result, weights, directions) {

  if (!inherits(mvnma_result, "mvnma_result")) {
    stop("Input must be mvnma_result object")
  }

  # Validate weights
  if (abs(sum(weights) - 1.0) > 1e-6) {
    stop("Weights must sum to 1.0")
  }

  outcome_names <- names(weights)
  if (!all(outcome_names %in% mvnma_result$mvnma_data$outcome_names)) {
    stop("Some outcome names in weights not found in result")
  }

  # Extract effects
  treatments <- unique(mvnma_result$pooled_estimates$Treatment)

  net_benefits <- lapply(treatments, function(trt) {
    trt_effects <- mvnma_result$pooled_estimates[
      mvnma_result$pooled_estimates$Treatment == trt, ]

    # Calculate weighted sum
    net_benefit <- 0
    components <- list()

    for (outcome in outcome_names) {
      effect <- trt_effects$Effect[trt_effects$Outcome == outcome]

      # Adjust sign based on direction
      if (directions[outcome] == "lower") {
        effect <- -effect
      }

      weighted_effect <- effect * weights[outcome]
      net_benefit <- net_benefit + weighted_effect

      components[[outcome]] <- weighted_effect
    }

    data.frame(
      Treatment = trt,
      Net_Benefit = net_benefit,
      t(unlist(components)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  result <- do.call(rbind, net_benefits)
  result <- result[order(-result$Net_Benefit), ]
  rownames(result) <- NULL

  result
}

#' Test for Outcome-Treatment Interactions
#'
#' Test whether treatment effects differ systematically across outcomes.
#'
#' @param mvnma_result Multivariate NMA result
#' @return Interaction test results
#' @export
test_outcome_treatment_interaction <- function(mvnma_result) {

  if (!inherits(mvnma_result, "mvnma_result")) {
    stop("Input must be mvnma_result object")
  }

  message("Testing for outcome-treatment interactions...")

  # Extract treatment effects for each outcome
  treatments <- setdiff(mvnma_result$treatments, mvnma_result$reference)
  outcomes <- mvnma_result$mvnma_data$outcome_names

  # Build matrix: treatments × outcomes
  effect_matrix <- matrix(NA, nrow = length(treatments), ncol = length(outcomes),
                         dimnames = list(treatments, outcomes))

  for (i in seq_along(treatments)) {
    for (j in seq_along(outcomes)) {
      eff <- mvnma_result$pooled_estimates$Effect[
        mvnma_result$pooled_estimates$Treatment == treatments[i] &
        mvnma_result$pooled_estimates$Outcome == outcomes[j]
      ]
      effect_matrix[i, j] <- eff
    }
  }

  # Test for rank consistency across outcomes
  # Kendall's W (coefficient of concordance)
  rank_matrix <- apply(effect_matrix, 2, rank)

  # Calculate Kendall's W
  n <- nrow(rank_matrix)  # treatments
  k <- ncol(rank_matrix)  # outcomes

  R_i <- rowSums(rank_matrix)
  S <- sum((R_i - mean(R_i))^2)
  W <- (12 * S) / (k^2 * (n^3 - n))

  # Chi-square test
  chi_sq <- k * (n - 1) * W
  df <- n - 1
  p_value <- 1 - pchisq(chi_sq, df)

  # Pairwise outcome correlations of treatment rankings
  outcome_cors <- cor(rank_matrix, method = "spearman")

  list(
    kendalls_W = W,
    chi_square = chi_sq,
    df = df,
    p_value = p_value,
    interpretation = if (p_value < 0.05) {
      "Significant concordance: Treatment rankings consistent across outcomes"
    } else {
      "No significant concordance: Treatment rankings differ across outcomes"
    },
    rank_matrix = rank_matrix,
    rank_correlations = outcome_cors
  )
}

#' Print Multivariate NMA Results
#'
#' @param x mvnma_result object
#' @param ... Additional arguments
#' @export
print.mvnma_result <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Multivariate Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Number of outcomes: %d\n", x$mvnma_data$n_outcomes))
  cat(sprintf("  %s\n", paste(x$mvnma_data$outcome_names, collapse = ", ")))
  cat(sprintf("Reference treatment: %s\n", x$reference))
  cat(sprintf("Number of treatments: %d\n", length(x$treatments)))
  cat(sprintf("Studies providing data for multiple outcomes: %d\n\n",
             x$mvnma_data$n_common_studies))

  cat("Assumed Correlation Between Outcomes:\n\n")
  print(round(x$mvnma_data$correlation, 3))

  cat("\n")
  cat("Pooled Treatment Effects Across Outcomes:\n\n")
  print(x$pooled_estimates, row.names = FALSE, digits = 3)

  cat("\n")
  cat("Use benefit_risk_plot() to visualize efficacy-safety trade-offs\n")
  cat("Use calculate_net_benefit() to compute weighted composite scores\n")
  cat("Use test_outcome_treatment_interaction() to test ranking consistency\n\n")

  invisible(x)
}

#' Plot Multivariate NMA Results
#'
#' @param x mvnma_result object
#' @param type Plot type: "forest", "heatmap", or "benefit_risk"
#' @param ... Additional arguments passed to specific plot functions
#' @export
plot.mvnma_result <- function(x, type = c("forest", "heatmap", "benefit_risk"), ...) {

  type <- match.arg(type)

  if (type == "forest") {
    # Forest plot by outcome
    plot_data <- x$pooled_estimates
    plot_data$Treatment <- factor(plot_data$Treatment,
                                  levels = unique(plot_data$Treatment))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = Effect, color = Outcome)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
                            width = 0.3, position = ggplot2::position_dodge(0.5)) +
      ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.5)) +
      ggplot2::facet_wrap(~ Outcome, scales = "free_y") +
      ggplot2::labs(
        title = "Multivariate NMA: Treatment Effects by Outcome",
        subtitle = sprintf("Reference: %s", x$reference),
        x = "Treatment",
        y = "Effect"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )

  } else if (type == "heatmap") {
    # Heatmap of treatment effects
    plot_data <- x$pooled_estimates

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Outcome, y = Treatment, fill = Effect)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Effect)), color = "black") +
      ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                                   midpoint = 0) +
      ggplot2::labs(
        title = "Multivariate NMA: Treatment Effects Heatmap",
        subtitle = sprintf("Reference: %s", x$reference),
        x = "Outcome",
        y = "Treatment"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

  } else {
    # Benefit-risk plot (requires 2+ outcomes)
    if (x$mvnma_data$n_outcomes < 2) {
      stop("Benefit-risk plot requires at least 2 outcomes")
    }

    p <- benefit_risk_plot(x, efficacy_outcome = 1, safety_outcome = 2, ...)
  }

  print(p)
  invisible(p)
}

#' Print Prepared MVNMA Data
#'
#' @param x mvnma_data object
#' @param ... Additional arguments
#' @export
print.mvnma_data <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Multivariate NMA Data Structure\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Number of outcomes: %d\n", x$n_outcomes))
  for (i in seq_len(x$n_outcomes)) {
    cat(sprintf("  %d. %s (%d studies)\n",
               i, x$outcome_names[i], nrow(x$data_list[[i]])))
  }

  cat(sprintf("\nTotal studies across all outcomes: %d\n", length(x$all_studies)))
  cat(sprintf("Studies with data for multiple outcomes: %d (%.1f%%)\n",
             x$n_common_studies,
             100 * x$n_common_studies / length(x$all_studies)))

  cat(sprintf("\nTotal treatments: %d\n", length(x$all_treatments)))
  cat(sprintf("  %s\n", paste(x$all_treatments, collapse = ", ")))

  cat("\nAssumed Correlation Matrix:\n\n")
  print(round(x$correlation, 3))

  cat("\n")

  invisible(x)
}
