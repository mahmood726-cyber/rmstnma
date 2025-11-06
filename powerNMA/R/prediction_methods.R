# Enhanced Prediction and Ranking Methods
#
# Novel methods from 2021-2024 for predictive treatment ranking
# References:
# - Rosenberger KJ, et al. (2021). BMC Med Res Methodol 21:213
# - Nikolakopoulou A, et al. (2022). PMC 9071581

#' Predictive P-score and SUCRA
#'
#' Calculates predictive treatment rankings that account for heterogeneity when
#' predicting to future studies/populations. Goes beyond standard SUCRA/P-score
#' by adding heterogeneity (τ²) to uncertainty.
#'
#' @param nma_object Network meta-analysis object (from run_powernma or netmeta)
#' @param n_simulations Number of simulations for predictive distribution (default: 10000)
#' @param target_heterogeneity Heterogeneity assumption for predictions:
#'   \itemize{
#'     \item{"observed"}{Use observed τ² from NMA}
#'     \item{"conservative"}{Use τ² + buffer (1.5×observed)}
#'     \item{"optimistic"}{Use 0.5×τ²}
#'   }
#' @param return_samples Return simulated samples? (default: FALSE)
#'
#' @return Object of class "predictive_ranking" with:
#'   \itemize{
#'     \item{standard_pscore}{Standard P-scores}
#'     \item{predictive_pscore}{Predictive P-scores (accounting for τ²)}
#'     \item{standard_sucra}{Standard SUCRA}
#'     \item{predictive_sucra}{Predictive SUCRA}
#'     \item{prob_best_standard}{P(treatment is best) - standard}
#'     \item{prob_best_predictive}{P(treatment is best) - predictive}
#'     \item{rankograms}{Rank distribution plots data}
#'   }
#'
#' @details
#' Standard metrics answer: "Which treatment was best in existing studies?"
#' Predictive metrics answer: "Which treatment is likely best in a future study?"
#'
#' The key difference is heterogeneity:
#' - Standard: TE ~ N(μ, σ²)
#' - Predictive: TE_future ~ N(μ, σ² + τ²)
#'
#' @examples
#' # After running NMA
#' nma_result <- run_powernma(data, mode = "standard")
#'
#' # Predictive rankings
#' pred_rank <- predictive_ranking(nma_result)
#'
#' # Compare standard vs predictive
#' print(pred_rank)
#' plot(pred_rank)
#'
#' @export
predictive_ranking <- function(nma_object,
                                n_simulations = 10000,
                                target_heterogeneity = c("observed", "conservative", "optimistic"),
                                return_samples = FALSE) {

  target_heterogeneity <- match.arg(target_heterogeneity)

  message("Calculating predictive treatment rankings...")
  message("Simulations: ", n_simulations)
  message("Heterogeneity: ", target_heterogeneity)

  # Extract from NMA object
  if ("network" %in% names(nma_object)) {
    # powerNMA object
    nm <- nma_object$network
    treatments <- nm$trts
    te_random <- nm$TE.random
    se_random <- nm$seTE.random
    tau <- nm$tau
  } else {
    # Direct netmeta object
    treatments <- nma_object$trts
    te_random <- nma_object$TE.random
    se_random <- nma_object$seTE.random
    tau <- nma_object$tau
  }

  n_trt <- length(treatments)

  # Adjust tau based on target
  if (target_heterogeneity == "conservative") {
    tau_pred <- tau * 1.5
  } else if (target_heterogeneity == "optimistic") {
    tau_pred <- tau * 0.5
  } else {
    tau_pred <- tau
  }

  message(sprintf("Using tau = %.3f for predictions", tau_pred))

  # Standard rankings (no additional heterogeneity)
  standard_ranks <- simulate_rankings(te_random, se_random, tau = 0, n_sim = n_simulations)

  # Predictive rankings (with heterogeneity)
  predictive_ranks <- simulate_rankings(te_random, se_random, tau = tau_pred, n_sim = n_simulations)

  # Calculate P-scores
  standard_pscore <- calculate_pscore_from_ranks(standard_ranks, treatments)
  predictive_pscore <- calculate_pscore_from_ranks(predictive_ranks, treatments)

  # Calculate SUCRA
  standard_sucra <- calculate_sucra_from_ranks(standard_ranks, treatments)
  predictive_sucra <- calculate_sucra_from_ranks(predictive_ranks, treatments)

  # Probability of being best
  prob_best_standard <- colMeans(standard_ranks == 1)
  names(prob_best_standard) <- treatments

  prob_best_predictive <- colMeans(predictive_ranks == 1)
  names(prob_best_predictive) <- treatments

  # Rankograms data
  rankograms_standard <- create_rankogram_data(standard_ranks, treatments)
  rankograms_predictive <- create_rankogram_data(predictive_ranks, treatments)

  result <- list(
    standard_pscore = standard_pscore,
    predictive_pscore = predictive_pscore,
    standard_sucra = standard_sucra,
    predictive_sucra = predictive_sucra,
    prob_best_standard = prob_best_standard,
    prob_best_predictive = prob_best_predictive,
    rankograms_standard = rankograms_standard,
    rankograms_predictive = rankograms_predictive,
    treatments = treatments,
    tau = tau,
    tau_predictive = tau_pred
  )

  if (return_samples) {
    result$standard_ranks <- standard_ranks
    result$predictive_ranks = predictive_ranks
  }

  class(result) <- "predictive_ranking"
  return(result)
}


#' Simulate treatment rankings
#' @keywords internal
simulate_rankings <- function(te, se, tau, n_sim) {
  n_trt <- length(te)

  # Simulate treatment effects
  # TE_i ~ N(μ_i, σ_i² + τ²)
  sims <- matrix(NA, nrow = n_sim, ncol = n_trt)

  for (i in 1:n_trt) {
    total_var <- se[i]^2 + tau^2
    sims[, i] <- rnorm(n_sim, mean = te[i], sd = sqrt(total_var))
  }

  # Rank in each simulation (higher = better)
  # Assuming higher treatment effect is better; reverse if needed
  ranks <- t(apply(-sims, 1, rank, ties.method = "random"))  # Negative for descending

  return(ranks)
}


#' Calculate P-score from simulated ranks
#' @keywords internal
calculate_pscore_from_ranks <- function(ranks, treatments) {
  n_trt <- ncol(ranks)

  # P-score = (1/(n-1)) * Σ P(TE_i > TE_k)
  # Approximated from rank probabilities

  p_scores <- rep(NA, n_trt)

  for (i in 1:n_trt) {
    # Mean rank for treatment i
    mean_rank <- mean(ranks[, i])

    # Convert to P-score: (n - mean_rank + 0.5) / n
    # Higher rank number = worse, so invert
    p_scores[i] <- (n_trt - mean_rank + 0.5) / n_trt
  }

  names(p_scores) <- treatments
  return(p_scores)
}


#' Calculate SUCRA from simulated ranks
#' @keywords internal
calculate_sucra_from_ranks <- function(ranks, treatments) {
  n_trt <- ncol(ranks)

  sucra <- rep(NA, n_trt)

  for (i in 1:n_trt) {
    # SUCRA = (1/(n-1)) * Σ cumulative rank probabilities
    rank_probs <- table(factor(ranks[, i], levels = 1:n_trt)) / nrow(ranks)

    cumrank <- cumsum(rank_probs)

    sucra[i] <- sum(cumrank[1:(n_trt - 1)]) / (n_trt - 1)
  }

  names(sucra) <- treatments
  return(sucra)
}


#' Create rankogram data
#' @keywords internal
create_rankogram_data <- function(ranks, treatments) {
  n_trt <- ncol(ranks)

  rankogram_list <- list()

  for (i in 1:n_trt) {
    trt <- treatments[i]
    rank_probs <- table(factor(ranks[, i], levels = 1:n_trt)) / nrow(ranks)

    rankogram_list[[trt]] <- data.frame(
      treatment = trt,
      rank = 1:n_trt,
      probability = as.vector(rank_probs)
    )
  }

  rankogram_data <- do.call(rbind, rankogram_list)
  rownames(rankogram_data) <- NULL

  return(rankogram_data)
}


#' Treatment Hierarchy for Multiple Outcomes
#'
#' Creates treatment hierarchy considering multiple outcomes with optional weights.
#' Useful for benefit-risk analysis.
#'
#' @param nma_objects List of NMA objects for different outcomes
#' @param outcome_names Names of outcomes
#' @param weights Outcome weights (default: equal weights)
#' @param higher_better Logical vector: is higher value better for each outcome?
#' @param ... Additional arguments
#'
#' @return Object of class "treatment_hierarchy" with multi-outcome rankings
#'
#' @examples
#' # Efficacy and safety NMAs
#' nma_efficacy <- run_powernma(efficacy_data, mode = "standard")
#' nma_safety <- run_powernma(safety_data, mode = "standard")
#'
#' hierarchy <- treatment_hierarchy(
#'   nma_objects = list(nma_efficacy, nma_safety),
#'   outcome_names = c("Efficacy", "Safety"),
#'   weights = c(0.6, 0.4),
#'   higher_better = c(TRUE, FALSE)  # Higher efficacy good, higher adverse events bad
#' )
#'
#' @export
treatment_hierarchy <- function(nma_objects,
                                 outcome_names,
                                 weights = NULL,
                                 higher_better = NULL,
                                 ...) {

  n_outcomes <- length(nma_objects)

  if (is.null(weights)) {
    weights <- rep(1 / n_outcomes, n_outcomes)
  }

  if (is.null(higher_better)) {
    higher_better <- rep(TRUE, n_outcomes)
  }

  # Normalize weights
  weights <- weights / sum(weights)

  message("Creating treatment hierarchy with ", n_outcomes, " outcomes")
  message("Weights: ", paste(round(weights, 2), collapse = ", "))

  # Extract treatment effects for each outcome
  te_list <- list()
  se_list <- list()

  for (i in 1:n_outcomes) {
    if ("network" %in% names(nma_objects[[i]])) {
      te_list[[i]] <- nma_objects[[i]]$network$TE.random
      se_list[[i]] <- nma_objects[[i]]$network$seTE.random
    } else {
      te_list[[i]] <- nma_objects[[i]]$TE.random
      se_list[[i]] <- nma_objects[[i]]$seTE.random
    }

    # Reverse sign if lower is better
    if (!higher_better[i]) {
      te_list[[i]] <- -te_list[[i]]
    }
  }

  # Calculate weighted composite score
  treatments <- rownames(te_list[[1]])
  n_trt <- length(treatments)

  composite_score <- rep(0, n_trt)

  for (i in 1:n_outcomes) {
    # Standardize to mean=0, sd=1
    te_std <- (te_list[[i]] - mean(te_list[[i]], na.rm = TRUE)) / sd(te_list[[i]], na.rm = TRUE)

    composite_score <- composite_score + weights[i] * te_std
  }

  names(composite_score) <- treatments

  # Rank by composite score
  composite_rank <- rank(-composite_score, ties.method = "average")

  result <- list(
    composite_score = composite_score,
    composite_rank = composite_rank,
    outcome_weights = weights,
    outcome_names = outcome_names,
    treatments = treatments
  )

  class(result) <- "treatment_hierarchy"
  return(result)
}


#' Print method for predictive ranking
#' @export
print.predictive_ranking <- function(x, ...) {
  cat("Predictive Treatment Ranking\n")
  cat("=============================\n\n")

  cat("Heterogeneity:\n")
  cat(sprintf("  Observed tau: %.3f\n", x$tau))
  cat(sprintf("  Predictive tau: %.3f\n", x$tau_predictive))
  cat("\n")

  # Create comparison table
  comparison <- data.frame(
    Treatment = x$treatments,
    Standard_PScore = x$standard_pscore,
    Predictive_PScore = x$predictive_pscore,
    Standard_SUCRA = x$standard_sucra,
    Predictive_SUCRA = x$predictive_sucra,
    ProbBest_Standard = x$prob_best_standard,
    ProbBest_Predictive = x$prob_best_predictive
  )

  # Sort by predictive P-score
  comparison <- comparison[order(-comparison$Predictive_PScore), ]

  cat("Treatment Rankings:\n")
  print(comparison, digits = 3, row.names = FALSE)

  cat("\nInterpretation:\n")
  cat("- Standard metrics: Performance in observed studies\n")
  cat("- Predictive metrics: Expected performance in future studies (accounts for heterogeneity)\n")

  invisible(x)
}


#' Plot method for predictive ranking
#' @export
plot.predictive_ranking <- function(x, type = c("comparison", "rankogram"), ...) {
  type <- match.arg(type)

  if (type == "comparison") {
    plot_predictive_comparison(x, ...)
  } else if (type == "rankogram") {
    plot_rankograms(x, ...)
  }
}


#' Plot standard vs predictive rankings
#' @keywords internal
plot_predictive_comparison <- function(pred_rank, ...) {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

  # Plot 1: P-scores
  treatments <- pred_rank$treatments
  n_trt <- length(treatments)

  y_pos <- 1:n_trt

  # Sort by predictive P-score
  ord <- order(pred_rank$predictive_pscore)

  plot(pred_rank$standard_pscore[ord], y_pos,
    xlim = c(0, 1),
    ylim = c(0.5, n_trt + 0.5),
    xlab = "P-score",
    ylab = "",
    yaxt = "n",
    pch = 19,
    col = "blue",
    cex = 1.5,
    main = "Standard vs Predictive P-scores"
  )

  points(pred_rank$predictive_pscore[ord], y_pos, pch = 17, col = "red", cex = 1.5)

  axis(2, at = y_pos, labels = treatments[ord], las = 1)

  legend("bottomright",
    legend = c("Standard", "Predictive"),
    pch = c(19, 17),
    col = c("blue", "red"),
    bty = "n"
  )

  # Plot 2: Probability of being best
  plot(pred_rank$prob_best_standard[ord], y_pos,
    xlim = c(0, max(pred_rank$prob_best_standard, pred_rank$prob_best_predictive)),
    ylim = c(0.5, n_trt + 0.5),
    xlab = "P(Best)",
    ylab = "",
    yaxt = "n",
    pch = 19,
    col = "blue",
    cex = 1.5,
    main = "Probability of Being Best"
  )

  points(pred_rank$prob_best_predictive[ord], y_pos, pch = 17, col = "red", cex = 1.5)

  axis(2, at = y_pos, labels = treatments[ord], las = 1)

  legend("topright",
    legend = c("Standard", "Predictive"),
    pch = c(19, 17),
    col = c("blue", "red"),
    bty = "n"
  )

  par(mfrow = c(1, 1))
}


#' Plot rankograms
#' @keywords internal
plot_rankograms <- function(pred_rank, predictive = TRUE, ...) {
  if (predictive) {
    rankogram_data <- pred_rank$rankograms_predictive
    main_title <- "Predictive Rankograms"
  } else {
    rankogram_data <- pred_rank$rankograms_standard
    main_title <- "Standard Rankograms"
  }

  treatments <- unique(rankogram_data$treatment)
  n_trt <- length(treatments)

  # Set up colors
  colors <- rainbow(n_trt)

  par(mar = c(5, 5, 4, 2))

  plot(1, 1,
    type = "n",
    xlim = c(1, n_trt),
    ylim = c(0, 1),
    xlab = "Rank",
    ylab = "Probability",
    main = main_title,
    xaxt = "n"
  )

  axis(1, at = 1:n_trt)

  for (i in 1:n_trt) {
    trt <- treatments[i]
    trt_data <- rankogram_data[rankogram_data$treatment == trt, ]

    lines(trt_data$rank, trt_data$probability,
      col = colors[i],
      lwd = 2,
      type = "b",
      pch = 19
    )
  }

  legend("topright",
    legend = treatments,
    col = colors,
    lwd = 2,
    pch = 19,
    bty = "n",
    cex = 0.8
  )
}


#' Print method for treatment hierarchy
#' @export
print.treatment_hierarchy <- function(x, ...) {
  cat("Multi-Outcome Treatment Hierarchy\n")
  cat("==================================\n\n")

  cat("Outcome weights:\n")
  for (i in seq_along(x$outcome_names)) {
    cat(sprintf("  %s: %.2f\n", x$outcome_names[i], x$outcome_weights[i]))
  }
  cat("\n")

  hierarchy <- data.frame(
    Treatment = x$treatments,
    Composite_Score = x$composite_score,
    Rank = x$composite_rank
  )

  hierarchy <- hierarchy[order(hierarchy$Rank), ]

  cat("Treatment Hierarchy (best to worst):\n")
  print(hierarchy, digits = 3, row.names = FALSE)

  invisible(x)
}
