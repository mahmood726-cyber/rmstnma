#' Comprehensive Visualization for Network Meta-Analysis
#'
#' Publication-ready plots for NMA including network graphs, forest plots,
#' funnel plots, and comparison-adjusted plots.
#'
#' @name visualization
NULL

#' Plot Network Diagram
#'
#' Create a network plot showing treatments as nodes and comparisons as edges.
#' Edge thickness represents number of studies; node size represents sample size
#' or number of studies.
#'
#' @param data Pairwise NMA data with treat1, treat2, studlab
#' @param nma_result Optional netmeta object for additional information
#' @param layout Graph layout algorithm: "circle", "spring", "tree", "star"
#' @param node_size_var Variable for node sizing: "studies", "patients", "constant"
#' @param edge_width_var Variable for edge width: "studies", "precision", "constant"
#' @param color_by Color nodes by: "none", "sucra", "effect_size"
#' @param sucra_result Optional sucra object for coloring by SUCRA
#' @param label_size Size of treatment labels
#' @return ggplot2 or igraph plot object
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' plot_network(data, layout = "spring")
#'
#' # With NMA results
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' plot_network(data, nma, layout = "circle")
#' }
plot_network <- function(data,
                        nma_result = NULL,
                        layout = c("circle", "spring", "tree", "star"),
                        node_size_var = c("studies", "patients", "constant"),
                        edge_width_var = c("studies", "precision", "constant"),
                        color_by = c("none", "sucra", "effect_size"),
                        sucra_result = NULL,
                        label_size = 4) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required for network plots")
  }

  layout <- match.arg(layout)
  node_size_var <- match.arg(node_size_var)
  edge_width_var <- match.arg(edge_width_var)
  color_by <- match.arg(color_by)

  # Build network graph
  treatments <- sort(unique(c(as.character(data$treat1), as.character(data$treat2))))
  n_treatments <- length(treatments)

  # Create edge list with counts
  edge_data <- data %>%
    dplyr::group_by(treat1, treat2) %>%
    dplyr::summarise(
      n_studies = dplyr::n(),
      mean_se = mean(seTE, na.rm = TRUE),
      .groups = "drop"
    )

  # Create igraph object
  g <- igraph::graph_from_data_frame(
    d = edge_data[, c("treat1", "treat2")],
    vertices = treatments,
    directed = FALSE
  )

  # Edge attributes
  if (edge_width_var == "studies") {
    igraph::E(g)$width <- edge_data$n_studies
  } else if (edge_width_var == "precision") {
    igraph::E(g)$width <- 1 / edge_data$mean_se
  } else {
    igraph::E(g)$width <- 2
  }

  # Normalize edge widths for visualization
  widths <- igraph::E(g)$width
  igraph::E(g)$width <- 1 + 4 * (widths - min(widths)) / (max(widths) - min(widths) + 1e-6)

  # Node sizes
  if (node_size_var == "studies") {
    node_counts <- table(c(as.character(data$treat1), as.character(data$treat2)))
    igraph::V(g)$size <- as.numeric(node_counts[igraph::V(g)$name])
  } else if (node_size_var == "patients") {
    if ("n1" %in% names(data) && "n2" %in% names(data)) {
      patient_counts <- data %>%
        dplyr::group_by(treat1) %>%
        dplyr::summarise(n = sum(n1, na.rm = TRUE)) %>%
        rbind(data %>%
                dplyr::group_by(treat2 = treat2) %>%
                dplyr::summarise(n = sum(n2, na.rm = TRUE))) %>%
        dplyr::group_by(treat = c(treat1, treat2)) %>%
        dplyr::summarise(total = sum(n), .groups = "drop")

      igraph::V(g)$size <- patient_counts$total[match(igraph::V(g)$name, patient_counts$treat)]
    } else {
      warning("Patient counts not available; using constant node size")
      igraph::V(g)$size <- 15
    }
  } else {
    igraph::V(g)$size <- 15
  }

  # Normalize node sizes
  sizes <- igraph::V(g)$size
  igraph::V(g)$size <- 10 + 20 * (sizes - min(sizes)) / (max(sizes) - min(sizes) + 1e-6)

  # Node colors
  if (color_by == "sucra" && !is.null(sucra_result)) {
    sucra_values <- sucra_result$sucra_scores[igraph::V(g)$name]
    igraph::V(g)$color <- grDevices::colorRampPalette(c("#d73027", "#fee08b", "#1a9850"))(100)[
      as.integer(sucra_values * 99) + 1
    ]
  } else if (color_by == "effect_size" && !is.null(nma_result)) {
    effects <- nma_result$TE.random[igraph::V(g)$name, nma_result$reference.group]
    effect_colors <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
    effect_normalized <- (effects - min(effects, na.rm = TRUE)) /
      (max(effects, na.rm = TRUE) - min(effects, na.rm = TRUE) + 1e-6)
    igraph::V(g)$color <- effect_colors[as.integer(effect_normalized * 99) + 1]
  } else {
    igraph::V(g)$color <- "#2c7fb8"
  }

  # Layout
  if (layout == "circle") {
    layout_coords <- igraph::layout_in_circle(g)
  } else if (layout == "spring") {
    layout_coords <- igraph::layout_with_fr(g)
  } else if (layout == "tree") {
    layout_coords <- igraph::layout_as_tree(g)
  } else {
    layout_coords <- igraph::layout_as_star(g)
  }

  # Plot using base graphics (works reliably)
  par(mar = c(1, 1, 3, 1))
  plot(g,
       layout = layout_coords,
       vertex.label = igraph::V(g)$name,
       vertex.label.cex = label_size / 4,
       vertex.label.color = "black",
       vertex.label.font = 2,
       edge.color = "gray50",
       main = "Network Diagram",
       sub = sprintf("%d treatments, %d comparisons",
                    n_treatments, igraph::ecount(g)))

  # Add legend
  if (color_by == "sucra") {
    legend("bottomright",
           legend = c("High SUCRA", "Medium", "Low SUCRA"),
           fill = c("#1a9850", "#fee08b", "#d73027"),
           title = "Treatment Ranking",
           bty = "n")
  }

  invisible(g)
}

#' Forest Plot for Network Meta-Analysis
#'
#' Create a forest plot showing treatment effects vs. reference with confidence
#' intervals. Supports both frequentist and Bayesian results.
#'
#' @param nma_result A netmeta object or Bayesian NMA result
#' @param reference Reference treatment (if NULL, uses default from nma_result)
#' @param sort_by Sort treatments by: "alphabet", "effect", "sucra"
#' @param sucra_result Optional sucra object for sorting by SUCRA
#' @param show_heterogeneity Show I² and τ² statistics
#' @param xlim X-axis limits (auto-detected if NULL)
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' forest_plot(nma, sort_by = "effect")
#' }
forest_plot <- function(nma_result,
                       reference = NULL,
                       sort_by = c("alphabet", "effect", "sucra"),
                       sucra_result = NULL,
                       show_heterogeneity = TRUE,
                       xlim = NULL) {

  sort_by <- match.arg(sort_by)

  if (!inherits(nma_result, "netmeta")) {
    stop("Currently supports netmeta objects only")
  }

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  # Extract effects vs reference
  treatments <- rownames(nma_result$TE.random)
  treatments <- setdiff(treatments, reference)

  effects <- nma_result$TE.random[treatments, reference]
  se_effects <- nma_result$seTE.random[treatments, reference]
  lower_ci <- nma_result$lower.random[treatments, reference]
  upper_ci <- nma_result$upper.random[treatments, reference]

  plot_data <- data.frame(
    Treatment = treatments,
    Effect = effects,
    SE = se_effects,
    Lower = lower_ci,
    Upper = upper_ci,
    stringsAsFactors = FALSE
  )

  # Sort
  if (sort_by == "effect") {
    plot_data <- plot_data[order(plot_data$Effect), ]
  } else if (sort_by == "sucra" && !is.null(sucra_result)) {
    sucra_order <- sucra_result$summary$Treatment
    plot_data <- plot_data[match(sucra_order, plot_data$Treatment), ]
    plot_data <- plot_data[!is.na(plot_data$Treatment), ]
  } else {
    plot_data <- plot_data[order(plot_data$Treatment), ]
  }

  plot_data$Treatment <- factor(plot_data$Treatment, levels = plot_data$Treatment)

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Effect, y = Treatment)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper),
                           height = 0.2, size = 0.8, color = "#2c7fb8") +
    ggplot2::geom_point(size = 3, color = "#2c7fb8") +
    ggplot2::labs(
      title = sprintf("Forest Plot: Treatment Effects vs. %s", reference),
      subtitle = if (show_heterogeneity) {
        sprintf("Random effects model | I² = %.1f%% | τ² = %.3f",
               nma_result$I2 * 100, nma_result$tau^2)
      } else NULL,
      x = sprintf("Effect Size (%s)", nma_result$sm),
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(xlim)) {
    p <- p + ggplot2::xlim(xlim)
  }

  p
}

#' Comparison-Adjusted Funnel Plot
#'
#' Create a funnel plot adjusted for different comparisons to assess
#' publication bias in network meta-analysis.
#'
#' @param nma_result A netmeta object
#' @param method Adjustment method: "common" or "random"
#' @return ggplot2 object
#' @export
#' @references
#' Chaimani A, Salanti G. Using network meta-analysis to evaluate the existence
#' of small-study effects in a network of interventions.
#' Research Synthesis Methods. 2012;3(2):161-176.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' comparison_adjusted_funnel(nma)
#' }
comparison_adjusted_funnel <- function(nma_result, method = c("common", "random")) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  method <- match.arg(method)

  # Extract study-level data
  studlab <- nma_result$studlab
  treat1 <- nma_result$treat1
  treat2 <- nma_result$treat2
  TE <- nma_result$TE
  seTE <- nma_result$seTE

  # Get comparison-specific effects from NMA
  if (method == "random") {
    comparison_effects <- nma_result$TE.random[cbind(treat2, treat1)]
  } else {
    comparison_effects <- nma_result$TE.fixed[cbind(treat2, treat1)]
  }

  # Calculate comparison-adjusted effects
  adjusted_effects <- TE - comparison_effects

  # Precision (1/SE)
  precision <- 1 / seTE

  plot_data <- data.frame(
    Study = studlab,
    AdjustedEffect = adjusted_effects,
    Precision = precision,
    SE = seTE,
    stringsAsFactors = FALSE
  )

  # Calculate funnel plot limits
  max_precision <- max(precision, na.rm = TRUE)
  funnel_se <- seq(0, max(seTE, na.rm = TRUE), length.out = 100)
  funnel_lower <- -1.96 * funnel_se
  funnel_upper <- 1.96 * funnel_se

  funnel_lines <- data.frame(
    SE = funnel_se,
    Lower = funnel_lower,
    Upper = funnel_upper
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = AdjustedEffect, y = Precision)) +
    ggplot2::geom_line(data = funnel_lines,
                      ggplot2::aes(x = Lower, y = 1/SE),
                      linetype = "dashed", color = "gray50") +
    ggplot2::geom_line(data = funnel_lines,
                      ggplot2::aes(x = Upper, y = 1/SE),
                      linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    ggplot2::geom_point(alpha = 0.6, size = 2, color = "#2c7fb8") +
    ggplot2::labs(
      title = "Comparison-Adjusted Funnel Plot",
      subtitle = "For assessing small-study effects and publication bias",
      x = "Comparison-Adjusted Effect",
      y = "Precision (1/SE)"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Net Heat Plot
#'
#' Visualize inconsistency in the network using net heat plot.
#' Shows contribution of each design to the inconsistency in each comparison.
#'
#' @param nma_result A netmeta object with inconsistency calculations
#' @return ggplot2 heatmap
#' @export
#' @references
#' Krahn U, et al. (2013). A graphical tool for locating inconsistency in
#' network meta-analyses. BMC Medical Research Methodology, 13:35.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' net_heat_plot(nma)
#' }
net_heat_plot <- function(nma_result) {

  if (!requireNamespace("netmeta", quietly = TRUE)) {
    stop("Package 'netmeta' required")
  }

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Try to use netmeta's netheat function if available
  tryCatch({
    netmeta::netheat(nma_result, random = TRUE)
  }, error = function(e) {
    message("Net heat plot generation failed: ", e$message)
    message("This may indicate insufficient inconsistency data in the network")
  })

  invisible(NULL)
}

#' Contribution Plot
#'
#' Show the contribution of each direct comparison to each network estimate.
#' Based on the percentage contribution matrix.
#'
#' @param nma_result A netmeta object
#' @param treatment Treatment to show contributions for (NULL for all)
#' @return ggplot2 or matrix plot
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' contribution_plot(nma, treatment = "DrugA")
#' }
contribution_plot <- function(nma_result, treatment = NULL) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Get contribution matrix if available
  if (is.null(nma_result$A.matrix)) {
    stop("Contribution matrix not available in netmeta object")
  }

  # Use netmeta's netcontrib function
  contrib <- netmeta::netcontrib(nma_result)

  # Plot
  if (!is.null(treatment)) {
    # Plot for specific treatment
    netmeta::plot(contrib, comparison = treatment)
  } else {
    # Show matrix
    netmeta::plot(contrib)
  }

  invisible(contrib)
}

#' Interval Plot
#'
#' Create an interval plot showing treatment effects and prediction intervals.
#' Prediction intervals account for between-study heterogeneity.
#'
#' @param nma_result A netmeta object
#' @param reference Reference treatment
#' @param show_prediction Show prediction intervals (wider than CIs)
#' @return ggplot2 object
#' @export
#' @references
#' Riley RD, et al. (2011). Interpretation of random effects meta-analyses.
#' BMJ, 342:d549.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' interval_plot(nma, show_prediction = TRUE)
#' }
interval_plot <- function(nma_result, reference = NULL, show_prediction = TRUE) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  treatments <- rownames(nma_result$TE.random)
  treatments <- setdiff(treatments, reference)

  effects <- nma_result$TE.random[treatments, reference]
  lower_ci <- nma_result$lower.random[treatments, reference]
  upper_ci <- nma_result$upper.random[treatments, reference]

  # Calculate prediction intervals if requested
  if (show_prediction && !is.null(nma_result$tau)) {
    tau <- nma_result$tau
    se_effects <- nma_result$seTE.random[treatments, reference]

    # Prediction interval: estimate ± t * sqrt(SE² + τ²)
    # Using large-sample z-value (1.96) as approximation
    pred_se <- sqrt(se_effects^2 + tau^2)
    lower_pi <- effects - 1.96 * pred_se
    upper_pi <- effects + 1.96 * pred_se
  } else {
    lower_pi <- lower_ci
    upper_pi <- upper_ci
  }

  plot_data <- data.frame(
    Treatment = treatments,
    Effect = effects,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    Lower_PI = lower_pi,
    Upper_PI = upper_pi,
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[order(plot_data$Effect), ]
  plot_data$Treatment <- factor(plot_data$Treatment, levels = plot_data$Treatment)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Effect, y = Treatment))

  if (show_prediction) {
    # Prediction intervals (wider)
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = Lower_PI, xmax = Upper_PI),
      height = 0, size = 1.5, color = "#bdbdbd", alpha = 0.6
    )
  }

  # Confidence intervals
  p <- p +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = Lower_CI, xmax = Upper_CI),
      height = 0.2, size = 1, color = "#2c7fb8"
    ) +
    ggplot2::geom_point(size = 3, color = "#2c7fb8") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = sprintf("Treatment Effects vs. %s", reference),
      subtitle = if (show_prediction) {
        sprintf("95%% CI (colored) and 95%% Prediction Interval (gray) | τ = %.3f", nma_result$tau)
      } else {
        "95% Confidence Intervals"
      },
      x = sprintf("Effect Size (%s)", nma_result$sm),
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}

#' Generate All Standard NMA Plots
#'
#' Create a comprehensive set of publication-ready plots for NMA results
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @param sucra_result Optional sucra object
#' @param output_dir Directory to save plots (NULL for display only)
#' @return List of plot objects
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' sucra <- calculate_sucra(nma)
#' plots <- generate_nma_plots(nma, data, sucra, output_dir = "plots")
#' }
generate_nma_plots <- function(nma_result, data, sucra_result = NULL, output_dir = NULL) {

  plots <- list()

  message("Generating network diagram...")
  plots$network <- plot_network(data, nma_result)

  message("Generating forest plot...")
  plots$forest <- forest_plot(nma_result, sort_by = "effect")

  message("Generating funnel plot...")
  plots$funnel <- comparison_adjusted_funnel(nma_result)

  message("Generating interval plot...")
  plots$intervals <- interval_plot(nma_result, show_prediction = TRUE)

  if (!is.null(sucra_result)) {
    message("Generating SUCRA plots...")
    plots$sucra_scores <- plot_sucra_scores(sucra_result)
    plots$rankogram <- plot_rankogram(sucra_result, style = "heatmap")
  }

  # Save if directory provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    for (plot_name in names(plots)) {
      filename <- file.path(output_dir, paste0(plot_name, ".png"))
      message(sprintf("Saving %s...", filename))

      if (plot_name == "network") {
        grDevices::png(filename, width = 800, height = 800)
        print(plots[[plot_name]])
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename, plots[[plot_name]], width = 10, height = 8, dpi = 300)
      }
    }

    message(sprintf("\nAll plots saved to: %s", output_dir))
  }

  invisible(plots)
}
