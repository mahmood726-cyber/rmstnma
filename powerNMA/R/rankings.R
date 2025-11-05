#' Treatment Rankings and SUCRA Scores for Network Meta-Analysis
#'
#' Comprehensive ranking tools including SUCRA, probability matrices, and
#' rank-o-grams. Based on Salanti et al. (2011) and extensions.
#'
#' @name rankings
NULL

#' Calculate SUCRA Scores
#'
#' Surface Under the Cumulative Ranking curve (SUCRA) provides a single number
#' associated with each treatment, representing the probability that the
#' treatment is among the best options.
#'
#' SUCRA = 1 indicates certainty that treatment is the best
#' SUCRA = 0 indicates certainty that treatment is the worst
#'
#' @param nma_result A netmeta object from network meta-analysis
#' @param small_values How to interpret effect sizes:
#'   "good" (smaller is better, e.g., adverse events),
#'   "bad" (larger is better, e.g., response rates)
#' @param nsim Number of simulations for Bayesian approach (default: 10000)
#' @return sucra object with scores, probabilities, and rankings
#' @export
#' @references
#' Salanti G, et al. (2011). Graphical methods and numerical summaries for
#' presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. Journal of Clinical Epidemiology, 64(2):163-171.
#'
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab,
#'                         data = data, reference.group = "A", sm = "OR")
#' sucra_results <- calculate_sucra(nma, small_values = "good")
#' print(sucra_results)
#' plot(sucra_results)
#' }
calculate_sucra <- function(nma_result, small_values = c("good", "bad"), nsim = 10000) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  small_values <- match.arg(small_values)

  # Get treatment effects and variance-covariance matrix
  treatments <- rownames(nma_result$TE.random)
  n_treatments <- length(treatments)

  # Use random effects estimates
  te_estimates <- nma_result$TE.random[, nma_result$reference.group]
  vcov_matrix <- nma_result$Cov.random

  # Simulate from multivariate normal distribution
  set.seed(42)  # For reproducibility
  simulated_effects <- MASS::mvrnorm(n = nsim, mu = te_estimates, Sigma = vcov_matrix)

  # For each simulation, rank the treatments
  if (small_values == "good") {
    # Smaller effects are better
    ranks_matrix <- t(apply(simulated_effects, 1, rank))
  } else {
    # Larger effects are better
    ranks_matrix <- t(apply(simulated_effects, 1, function(x) rank(-x)))
  }

  # Calculate probability of each treatment being at each rank
  prob_matrix <- matrix(0, n_treatments, n_treatments)
  rownames(prob_matrix) <- treatments
  colnames(prob_matrix) <- paste0("Rank", 1:n_treatments)

  for (i in 1:n_treatments) {
    for (j in 1:n_treatments) {
      prob_matrix[i, j] <- sum(ranks_matrix[, i] == j) / nsim
    }
  }

  # Calculate SUCRA scores
  # SUCRA = (1/(n-1)) * sum of probabilities of being better than average rank
  sucra_scores <- rep(0, n_treatments)
  names(sucra_scores) <- treatments

  for (i in 1:n_treatments) {
    cumsum_probs <- cumsum(prob_matrix[i, ])
    sucra_scores[i] <- sum(cumsum_probs[-n_treatments]) / (n_treatments - 1)
  }

  # Probability of being best
  prob_best <- prob_matrix[, 1]

  # Mean rank
  mean_ranks <- apply(ranks_matrix, 2, mean)
  names(mean_ranks) <- treatments

  # Median rank
  median_ranks <- apply(ranks_matrix, 2, median)
  names(median_ranks) <- treatments

  # Create summary table
  summary_table <- data.frame(
    Treatment = treatments,
    SUCRA = round(sucra_scores * 100, 1),  # Convert to percentage
    ProbBest = round(prob_best * 100, 1),
    MeanRank = round(mean_ranks, 2),
    MedianRank = median_ranks,
    stringsAsFactors = FALSE
  )

  # Sort by SUCRA descending
  summary_table <- summary_table[order(-summary_table$SUCRA), ]

  output <- list(
    sucra_scores = sucra_scores,
    prob_best = prob_best,
    mean_ranks = mean_ranks,
    median_ranks = median_ranks,
    prob_matrix = prob_matrix,
    summary = summary_table,
    small_values = small_values,
    n_treatments = n_treatments,
    nsim = nsim
  )

  class(output) <- c("sucra", "list")
  output
}

#' Calculate Treatment Rankings from Bayesian NMA
#'
#' Extract treatment rankings from a Bayesian network meta-analysis object
#' (e.g., from gemtc or BUGSnet packages)
#'
#' @param bayesian_result Bayesian NMA result object
#' @param small_values Whether smaller is better ("good") or worse ("bad")
#' @return sucra object
#' @export
#' @examples
#' \dontrun{
#' # From gemtc
#' model <- gemtc::mtc.model(network)
#' results <- gemtc::mtc.run(model)
#' rankings <- rank_bayesian_nma(results)
#' }
rank_bayesian_nma <- function(bayesian_result, small_values = c("good", "bad")) {

  small_values <- match.arg(small_values)

  # Try to extract from gemtc object
  if (inherits(bayesian_result, "mtc.result")) {
    return(rank_gemtc_result(bayesian_result, small_values))
  }

  # Add support for other Bayesian packages here
  stop("Unsupported Bayesian NMA object type. Currently supports: gemtc (mtc.result)")
}

#' Rank GEMTC Result
#'
#' @param gemtc_result mtc.result object from gemtc package
#' @param small_values Whether smaller is better
#' @return sucra object
#' @keywords internal
rank_gemtc_result <- function(gemtc_result, small_values) {

  if (!requireNamespace("gemtc", quietly = TRUE)) {
    stop("Package 'gemtc' required for ranking Bayesian NMA results")
  }

  # Extract relative effects
  rel_effects <- gemtc::relative.effect(gemtc_result, preserve.extra = FALSE)

  # Get samples
  samples <- as.matrix(rel_effects$samples)

  treatments <- colnames(samples)
  n_treatments <- length(treatments)
  nsim <- nrow(samples)

  # Rank treatments for each MCMC sample
  if (small_values == "good") {
    ranks_matrix <- t(apply(samples, 1, rank))
  } else {
    ranks_matrix <- t(apply(samples, 1, function(x) rank(-x)))
  }

  # Calculate probability matrix
  prob_matrix <- matrix(0, n_treatments, n_treatments)
  rownames(prob_matrix) <- treatments
  colnames(prob_matrix) <- paste0("Rank", 1:n_treatments)

  for (i in 1:n_treatments) {
    for (j in 1:n_treatments) {
      prob_matrix[i, j] <- sum(ranks_matrix[, i] == j) / nsim
    }
  }

  # Calculate SUCRA
  sucra_scores <- rep(0, n_treatments)
  names(sucra_scores) <- treatments

  for (i in 1:n_treatments) {
    cumsum_probs <- cumsum(prob_matrix[i, ])
    sucra_scores[i] <- sum(cumsum_probs[-n_treatments]) / (n_treatments - 1)
  }

  prob_best <- prob_matrix[, 1]
  mean_ranks <- apply(ranks_matrix, 2, mean)
  median_ranks <- apply(ranks_matrix, 2, median)

  summary_table <- data.frame(
    Treatment = treatments,
    SUCRA = round(sucra_scores * 100, 1),
    ProbBest = round(prob_best * 100, 1),
    MeanRank = round(mean_ranks, 2),
    MedianRank = median_ranks,
    stringsAsFactors = FALSE
  )

  summary_table <- summary_table[order(-summary_table$SUCRA), ]

  output <- list(
    sucra_scores = sucra_scores,
    prob_best = prob_best,
    mean_ranks = mean_ranks,
    median_ranks = median_ranks,
    prob_matrix = prob_matrix,
    summary = summary_table,
    small_values = small_values,
    n_treatments = n_treatments,
    nsim = nsim
  )

  class(output) <- c("sucra", "list")
  output
}

#' Create Rank-o-gram Plot
#'
#' Visualize the probability of each treatment being at each rank position.
#' Also known as rankogram or ranking probability plot.
#'
#' @param sucra_result A sucra object from calculate_sucra()
#' @param style Plot style: "heatmap", "bar", or "line"
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' sucra_results <- calculate_sucra(nma)
#' plot_rankogram(sucra_results, style = "heatmap")
#' }
plot_rankogram <- function(sucra_result, style = c("heatmap", "bar", "line")) {

  if (!inherits(sucra_result, "sucra")) {
    stop("sucra_result must be a sucra object from calculate_sucra()")
  }

  style <- match.arg(style)

  prob_matrix <- sucra_result$prob_matrix
  treatments <- rownames(prob_matrix)

  # Convert to long format
  plot_data <- expand.grid(
    Treatment = treatments,
    Rank = 1:ncol(prob_matrix),
    stringsAsFactors = FALSE
  )

  plot_data$Probability <- as.vector(t(prob_matrix))

  # Order treatments by SUCRA
  treatment_order <- sucra_result$summary$Treatment
  plot_data$Treatment <- factor(plot_data$Treatment, levels = treatment_order)

  if (style == "heatmap") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Rank, y = Treatment, fill = Probability)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = "#f7f7f7", mid = "#67a9cf", high = "#016c59",
        midpoint = 0.5, limits = c(0, 1),
        labels = scales::percent
      ) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Probability * 100)),
                        size = 3, color = "black") +
      ggplot2::labs(
        title = "Rank-o-gram: Probability of Each Treatment at Each Rank",
        x = "Rank Position",
        y = "Treatment",
        fill = "Probability"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )

  } else if (style == "bar") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = as.factor(Rank), y = Probability, fill = Treatment)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(
        title = "Ranking Probabilities by Position",
        x = "Rank Position",
        y = "Probability",
        fill = "Treatment"
      ) +
      ggplot2::theme_minimal()

  } else {  # line
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Rank, y = Probability, color = Treatment, group = Treatment)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(
        title = "Cumulative Ranking Curves",
        x = "Rank Position",
        y = "Probability",
        color = "Treatment"
      ) +
      ggplot2::theme_minimal()
  }

  p
}

#' Plot SUCRA Scores
#'
#' Create a bar plot of SUCRA scores for all treatments
#'
#' @param sucra_result A sucra object
#' @param top_n Show only top N treatments (NULL for all)
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' sucra_results <- calculate_sucra(nma)
#' plot_sucra_scores(sucra_results)
#' }
plot_sucra_scores <- function(sucra_result, top_n = NULL) {

  if (!inherits(sucra_result, "sucra")) {
    stop("sucra_result must be a sucra object")
  }

  plot_data <- sucra_result$summary

  if (!is.null(top_n)) {
    plot_data <- plot_data[1:min(top_n, nrow(plot_data)), ]
  }

  # Ensure correct order
  plot_data$Treatment <- factor(plot_data$Treatment, levels = rev(plot_data$Treatment))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = SUCRA)) +
    ggplot2::geom_col(fill = "#2c7fb8", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", SUCRA)),
                      hjust = -0.2, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    ggplot2::labs(
      title = "SUCRA Scores (Surface Under Cumulative Ranking Curve)",
      subtitle = "Higher scores indicate higher probability of being among the best treatments",
      x = NULL,
      y = "SUCRA (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
}

#' Compare Rankings Across Multiple Outcomes
#'
#' When multiple outcomes are analyzed, compare how treatment rankings differ
#'
#' @param sucra_list Named list of sucra objects (one per outcome)
#' @return Data frame with rankings across outcomes
#' @export
#' @examples
#' \dontrun{
#' sucra_efficacy <- calculate_sucra(nma_efficacy)
#' sucra_safety <- calculate_sucra(nma_safety)
#' comparison <- compare_rankings(list(
#'   Efficacy = sucra_efficacy,
#'   Safety = sucra_safety
#' ))
#' print(comparison)
#' }
compare_rankings <- function(sucra_list) {

  if (!is.list(sucra_list) || length(sucra_list) == 0) {
    stop("sucra_list must be a named list of sucra objects")
  }

  # Check all are sucra objects
  if (!all(sapply(sucra_list, inherits, "sucra"))) {
    stop("All elements must be sucra objects")
  }

  outcome_names <- names(sucra_list)
  if (is.null(outcome_names)) {
    outcome_names <- paste0("Outcome", seq_along(sucra_list))
  }

  # Get all treatments
  all_treatments <- unique(unlist(lapply(sucra_list, function(x) x$summary$Treatment)))

  # Build comparison table
  comparison <- data.frame(Treatment = all_treatments, stringsAsFactors = FALSE)

  for (i in seq_along(sucra_list)) {
    outcome <- outcome_names[i]
    sucra_obj <- sucra_list[[i]]

    # Match treatments
    idx <- match(all_treatments, sucra_obj$summary$Treatment)

    comparison[[paste0(outcome, "_SUCRA")]] <- sucra_obj$summary$SUCRA[idx]
    comparison[[paste0(outcome, "_Rank")]] <- sucra_obj$summary$MedianRank[idx]
  }

  # Calculate mean SUCRA across outcomes
  sucra_cols <- grep("_SUCRA$", names(comparison), value = TRUE)
  comparison$Mean_SUCRA <- rowMeans(comparison[, sucra_cols], na.rm = TRUE)

  # Sort by mean SUCRA
  comparison <- comparison[order(-comparison$Mean_SUCRA), ]

  class(comparison) <- c("ranking_comparison", "data.frame")
  comparison
}

#' Print SUCRA Results
#'
#' @param x sucra object
#' @param ... Additional arguments
#' @export
print.sucra <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Treatment Rankings and SUCRA Scores\n")
  cat("══════════════════════════════════════════════════\n\n")

  cat(sprintf("Number of treatments: %d\n", x$n_treatments))
  cat(sprintf("Simulations: %d\n", x$nsim))
  cat(sprintf("Interpretation: %s values are %s\n\n",
             ifelse(x$small_values == "good", "Smaller", "Larger"),
             ifelse(x$small_values == "good", "better", "better")))

  cat("Rankings:\n")
  print(x$summary, row.names = FALSE)

  cat("\n")
  cat("SUCRA interpretation:\n")
  cat("  100%: Certainty treatment is the best\n")
  cat("    0%: Certainty treatment is the worst\n")
  cat("   50%: No preference (average)\n")

  cat("\nTop 3 treatments:\n")
  for (i in 1:min(3, nrow(x$summary))) {
    cat(sprintf("  %d. %s (SUCRA: %.1f%%, P(best): %.1f%%)\n",
               i,
               x$summary$Treatment[i],
               x$summary$SUCRA[i],
               x$summary$ProbBest[i]))
  }

  invisible(x)
}

#' Plot SUCRA Object
#'
#' @param x sucra object
#' @param type Plot type: "sucra", "rankogram", or "both"
#' @param ... Additional arguments
#' @export
plot.sucra <- function(x, type = c("sucra", "rankogram", "both"), ...) {
  type <- match.arg(type)

  if (type == "sucra") {
    plot_sucra_scores(x, ...)
  } else if (type == "rankogram") {
    plot_rankogram(x, ...)
  } else {
    # Both plots
    p1 <- plot_sucra_scores(x)
    p2 <- plot_rankogram(x, style = "heatmap")

    if (requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::wrap_plots(p1, p2, ncol = 1)
    } else {
      message("Install 'patchwork' package to view both plots together")
      print(p1)
      print(p2)
    }
  }
}

#' Print Ranking Comparison
#'
#' @param x ranking_comparison object
#' @param ... Additional arguments
#' @export
print.ranking_comparison <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Treatment Rankings Across Multiple Outcomes\n")
  cat("══════════════════════════════════════════════════\n\n")

  print(as.data.frame(x), row.names = FALSE)

  cat("\nNote: Higher SUCRA = better ranking\n")
  cat("      Lower Rank = better position (Rank 1 is best)\n")

  invisible(x)
}
