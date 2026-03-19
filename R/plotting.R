#' Plot RMST NMA Results
#'
#' Creates various plots for RMST NMA results
#'
#' @param x An rmst_nma_fit object
#' @param type Type of plot ("forest", "rank", "rmst", "trace", "network")
#' @param tau Specific tau value for plotting (if multiple fitted)
#' @param ... Additional plotting arguments
#'
#' @return A ggplot object
#' @export
#'
#' @importFrom stats median quantile sd
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbarh geom_point geom_col
#' @importFrom ggplot2 labs theme_minimal theme element_blank scale_x_continuous annotate theme_void
#' @importFrom ggplot2 geom_line facet_wrap scale_color_manual coord_fixed geom_ribbon geom_segment geom_text
#' @importFrom dplyr bind_rows
plot.rmst_nma_fit <- function(x,
                              type = c("forest", "rank", "rmst", "trace", "network"),
                              tau = NULL,
                              ...) {

  type <- match.arg(type)

  if (is.null(tau) && length(x$tau) > 1) {
    tau <- x$tau[1]
    message(sprintf("Multiple tau values fitted. Plotting for tau = %.1f", tau))
  } else if (is.null(tau)) {
    tau <- x$tau
  }

  switch(type,
         forest = plot_forest(x, tau, ...),
         rank = plot_rankogram(x, tau, ...),
         rmst = plot_rmst_curves(x, ...),
         trace = plot_trace(x, ...),
         network = plot_network(x, ...)
  )
}

#' Forest plot of RMST contrasts
#' @keywords internal
plot_forest <- function(fit, tau, reference = NULL, ...) {

  # Check if posterior is available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for forest plots. Install it with: install.packages('posterior')")
  }

  # Get treatment names
  treatments <- fit$network$treatments
  n_trt <- length(treatments)

  if (is.null(reference)) {
    reference <- treatments[1]
  }

  ref_idx <- which(treatments == reference)
  tau_idx <- which(fit$tau == tau)

  if (length(tau_idx) == 0) {
    stop(sprintf("tau = %.1f not found in fitted model", tau))
  }

  # Extract RMST predictions from Stan generated quantities
  draws <- posterior::as_draws_df(fit$fit$draws("rmst_pred"))

  # Calculate contrasts
  contrasts <- list()

  for (i in 1:n_trt) {
    if (i != ref_idx) {
      # Get RMST draws for both treatments at specified tau
      rmst_i <- draws[[sprintf("rmst_pred[%d,%d]", i, tau_idx)]]
      rmst_ref <- draws[[sprintf("rmst_pred[%d,%d]", ref_idx, tau_idx)]]

      if (!is.null(rmst_i) && !is.null(rmst_ref)) {
        contrast_draws <- rmst_i - rmst_ref

        contrasts[[treatments[i]]] <- data.frame(
          treatment = treatments[i],
          mean = mean(contrast_draws),
          lower = quantile(contrast_draws, 0.025),
          median = median(contrast_draws),
          upper = quantile(contrast_draws, 0.975),
          prob_positive = mean(contrast_draws > 0),
          row.names = NULL
        )
      }
    }
  }

  if (length(contrasts) == 0) {
    stop("No contrasts could be calculated from the fitted model")
  }

  plot_data <- dplyr::bind_rows(contrasts)

  # Order treatments by mean difference
  plot_data$treatment <- factor(plot_data$treatment,
                                levels = plot_data$treatment[order(plot_data$mean)])

  # Create forest plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = mean, y = treatment)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = sprintf("RMST Differences vs %s at tau=%.1f", reference, tau),
      subtitle = "95% Credible Intervals",
      x = "RMST Difference",
      y = "Treatment"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

#' Rankogram plot using Stan's generated quantities
#' @keywords internal
plot_rankogram <- function(fit, tau, ...) {

  # Check if posterior is available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for rankogram plots. Install it with: install.packages('posterior')")
  }

  treatments <- fit$network$treatments
  n_trt <- length(treatments)
  tau_idx <- which(fit$tau == tau)

  if (length(tau_idx) == 0) {
    stop(sprintf("tau = %.1f not found in fitted model", tau))
  }

  # Extract rank draws from Stan generated quantities
  draws <- posterior::as_draws_df(fit$fit$draws("rank"))

  # Calculate rank probabilities
  rank_probs <- matrix(0, n_trt, n_trt,
                       dimnames = list(treatments, paste0("Rank ", 1:n_trt)))

  for (trt_idx in 1:n_trt) {
    param_name <- sprintf("rank[%d,%d]", trt_idx, tau_idx)
    if (param_name %in% names(draws)) {
      ranks <- draws[[param_name]]
      for (r in 1:n_trt) {
        rank_probs[trt_idx, r] <- mean(ranks == r)
      }
    }
  }

  # Calculate SUCRA scores
  sucra <- numeric(n_trt)
  for (i in 1:n_trt) {
    sucra[i] <- sum((n_trt - 1:n_trt) * rank_probs[i, ]) / (n_trt - 1)
  }

  # Convert to long format for plotting
  plot_data <- expand.grid(
    treatment = treatments,
    rank = 1:n_trt
  )
  plot_data$probability <- as.vector(rank_probs)

  # Order treatments by SUCRA
  treatment_order <- treatments[order(sucra, decreasing = TRUE)]
  plot_data$treatment <- factor(plot_data$treatment, levels = treatment_order)

  # Create rankogram
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = rank, y = probability,
                                    fill = treatment)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_x_continuous(breaks = 1:n_trt) +
    ggplot2::labs(
      title = sprintf("Treatment Rankings at tau=%.1f", tau),
      subtitle = sprintf("SUCRA scores: %s",
                         paste(sprintf("%s=%.2f",
                                       treatment_order,
                                       sort(sucra, decreasing = TRUE)),
                               collapse = ", ")),
      x = "Rank",
      y = "Probability",
      fill = "Treatment"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

#' Plot RMST curves across tau values
#' @keywords internal
plot_rmst_curves <- function(fit, ...) {

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required. Install it with: install.packages('posterior')")
  }

  # Extract RMST predictions
  draws <- posterior::as_draws_df(fit$fit$draws("rmst_pred"))

  treatments <- fit$network$treatments
  n_trt <- length(treatments)
  n_tau <- length(fit$tau)

  # Prepare data for plotting
  plot_data <- list()

  for (trt_idx in 1:n_trt) {
    for (tau_idx in 1:n_tau) {
      param_name <- sprintf("rmst_pred[%d,%d]", trt_idx, tau_idx)

      if (param_name %in% names(draws)) {
        rmst_draws <- draws[[param_name]]

        plot_data[[length(plot_data) + 1]] <- data.frame(
          treatment = treatments[trt_idx],
          tau = fit$tau[tau_idx],
          mean = mean(rmst_draws),
          lower = quantile(rmst_draws, 0.025),
          upper = quantile(rmst_draws, 0.975),
          lower_50 = quantile(rmst_draws, 0.25),
          upper_50 = quantile(rmst_draws, 0.75)
        )
      }
    }
  }

  plot_df <- dplyr::bind_rows(plot_data)

  # Create RMST curves plot
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = tau, y = mean,
                                    color = treatment, fill = treatment)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_50, ymax = upper_50),
                         alpha = 0.3) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::labs(
      title = "RMST Estimates Across Restriction Times",
      subtitle = "Lines show means with 50% and 95% credible intervals",
      x = "Restriction Time (tau)",
      y = "RMST",
      color = "Treatment",
      fill = "Treatment"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

#' Trace plots for MCMC diagnostics
#' @keywords internal
plot_trace <- function(fit, pars = NULL, ...) {

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required. Install it with: install.packages('posterior')")
  }

  # Default parameters to plot
  if (is.null(pars)) {
    pars <- c("d[2]", "sigma[1]", "mu[1]")
    # Filter to existing parameters
    all_pars <- fit$fit$metadata()$variables
    pars <- pars[pars %in% all_pars]

    if (length(pars) == 0) {
      # Try alternative parameter names
      pars <- grep("^d\\[", all_pars, value = TRUE)[1:min(3, length(grep("^d\\[", all_pars, value = TRUE)))]
    }
  }

  # Extract draws
  draws <- posterior::as_draws_df(fit$fit$draws(variables = pars))

  # Convert to long format for plotting
  chains <- draws$.chain
  iterations <- draws$.iteration
  draws_long <- list()

  for (par in pars) {
    if (par %in% names(draws)) {
      draws_long[[length(draws_long) + 1]] <- data.frame(
        parameter = par,
        chain = chains,
        iteration = iterations,
        value = draws[[par]]
      )
    }
  }

  plot_data <- dplyr::bind_rows(draws_long)
  plot_data$chain <- factor(plot_data$chain)

  # Create trace plots
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = iteration, y = value, color = chain)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "MCMC Trace Plots",
      x = "Iteration",
      y = "Parameter Value",
      color = "Chain"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))

  return(p)
}

#' Network plot
#' @keywords internal
plot_network <- function(fit, ...) {

  adj <- fit$network$adjacency
  treatments <- fit$network$treatments
  n_trt <- length(treatments)

  # Convert adjacency to edge list
  edges <- which(adj > 0, arr.ind = TRUE)
  edges <- edges[edges[,1] < edges[,2], , drop = FALSE]

  if (nrow(edges) > 0) {
    edge_data <- data.frame(
      from = treatments[edges[,1]],
      to = treatments[edges[,2]],
      weight = adj[edges]
    )

    # Simple circular layout for nodes
    angles <- seq(0, 2*pi, length.out = n_trt + 1)[-1]
    node_data <- data.frame(
      treatment = treatments,
      x = cos(angles),
      y = sin(angles)
    )

    # Merge coordinates for edges
    edge_coords <- merge(edge_data, node_data,
                         by.x = "from", by.y = "treatment")
    names(edge_coords)[4:5] <- c("x1", "y1")
    edge_coords <- merge(edge_coords, node_data,
                         by.x = "to", by.y = "treatment")
    names(edge_coords)[6:7] <- c("x2", "y2")

    # Create network plot
    p <- ggplot2::ggplot() +
      # Draw edges
      ggplot2::geom_segment(data = edge_coords,
                            ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2,
                                         size = weight),
                            alpha = 0.5, color = "gray50") +
      # Draw nodes
      ggplot2::geom_point(data = node_data,
                          ggplot2::aes(x = x, y = y),
                          size = 10, color = "steelblue") +
      # Add labels
      ggplot2::geom_text(data = node_data,
                         ggplot2::aes(x = x, y = y, label = treatment),
                         size = 4, color = "white") +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = sprintf("Network Structure: %d treatments, %d studies",
                                    fit$network$n_treatments,
                                    fit$network$n_studies))

  } else {
    # No edges case
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = sprintf("Network plot\n%d treatments, %d studies\nNo direct comparisons",
                                        fit$network$n_treatments,
                                        fit$network$n_studies),
                        size = 6) +
      ggplot2::theme_void()
  }

  return(p)
}
