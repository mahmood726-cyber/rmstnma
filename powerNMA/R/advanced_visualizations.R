#' Advanced Interactive Visualizations for Network Meta-Analysis
#'
#' State-of-the-art interactive visualizations using plotly, creating
#' publication-ready and interactive dashboard-quality graphics.
#'
#' @name advanced_visualizations
#' @references
#' Chaimani A, et al. (2023). Data visualisation approaches for component
#' network meta-analysis. BMC Medical Research Methodology, 23:186.
NULL

#' Create Interactive Network Dashboard
#'
#' Generate comprehensive interactive dashboard with all visualizations.
#'
#' @param nma_result NMA result object
#' @param data Original data
#' @param launch_browser Launch in browser?
#' @return Shiny app object
#' @export
create_interactive_dashboard <- function(nma_result, data, launch_browser = TRUE) {

  if (!requireNamespace("shiny", quietly = TRUE) ||
      !requireNamespace("plotly", quietly = TRUE)) {
    stop("Packages 'shiny' and 'plotly' required for interactive dashboard")
  }

  message("Creating interactive dashboard...")

  # UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("Network Meta-Analysis Interactive Dashboard"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("plot_type", "Select Visualization:",
                          choices = c("Network Graph" = "network",
                                    "Forest Plot" = "forest",
                                    "Funnel Plot" = "funnel",
                                    "Treatment Rankings" = "rankings",
                                    "Heat Map" = "heatmap",
                                    "3D Network" = "network3d")),

        shiny::sliderInput("node_size", "Node Size:", min = 5, max = 50, value = 20),
        shiny::sliderInput("edge_width", "Edge Width:", min = 1, max = 10, value = 3),

        shiny::checkboxInput("show_labels", "Show Labels", TRUE),
        shiny::checkboxInput("show_ci", "Show Confidence Intervals", TRUE),

        shiny::hr(),
        shiny::downloadButton("download_plot", "Download Plot")
      ),

      shiny::mainPanel(
        plotly::plotlyOutput("main_plot", height = "600px"),

        shiny::hr(),

        shiny::tabsetPanel(
          shiny::tabPanel("Summary",
                         shiny::verbatimTextOutput("summary_text")),
          shiny::tabPanel("Treatment Effects",
                         DT::DTOutput("effects_table")),
          shiny::tabPanel("Heterogeneity",
                         shiny::verbatimTextOutput("heterogeneity_text")),
          shiny::tabPanel("Network Metrics",
                         shiny::verbatimTextOutput("network_metrics"))
        )
      )
    )
  )

  # Server
  server <- function(input, output, session) {

    # Main plot
    output$main_plot <- plotly::renderPlotly({
      switch(input$plot_type,
             network = create_interactive_network(nma_result, data, input),
             forest = create_interactive_forest(nma_result, input),
             funnel = create_interactive_funnel(nma_result, data, input),
             rankings = create_interactive_rankings(nma_result, input),
             heatmap = create_interactive_heatmap(nma_result, input),
             network3d = create_3d_network(nma_result, data, input))
    })

    # Summary
    output$summary_text <- shiny::renderPrint({
      cat(sprintf("Treatments: %d\n", length(nma_result$trts)))
      cat(sprintf("Studies: %d\n", length(unique(nma_result$studlab))))
      cat(sprintf("Comparisons: %d\n", length(nma_result$studlab)))
      cat(sprintf("\nHeterogeneity (I²): %.1f%%\n", nma_result$I2))
      cat(sprintf("Tau: %.3f\n", nma_result$tau))
    })

    # Effects table
    output$effects_table <- DT::renderDT({
      effects_df <- data.frame(
        Treatment = rownames(nma_result$TE.random),
        Effect = nma_result$TE.random[, 1],
        SE = nma_result$seTE.random[, 1],
        Lower_95 = nma_result$lower.random[, 1],
        Upper_95 = nma_result$upper.random[, 1]
      )

      DT::datatable(effects_df,
                   options = list(pageLength = 10),
                   rownames = FALSE) %>%
        DT::formatRound(columns = 2:5, digits = 3)
    })

    # Heterogeneity
    output$heterogeneity_text <- shiny::renderPrint({
      cat("Heterogeneity Assessment\n")
      cat("========================\n\n")
      cat(sprintf("I² = %.1f%%\n", nma_result$I2))
      cat(sprintf("τ² = %.3f\n", nma_result$tau^2))
      cat(sprintf("τ = %.3f\n", nma_result$tau))

      interpretation <- if (nma_result$I2 < 25) "Low"
                       else if (nma_result$I2 < 50) "Moderate"
                       else if (nma_result$I2 < 75) "Substantial"
                       else "Considerable"

      cat(sprintf("\nInterpretation: %s heterogeneity\n", interpretation))
    })

    # Network metrics
    output$network_metrics <- shiny::renderPrint({
      geometry <- analyze_network_geometry(data)
      cat("Network Metrics\n")
      cat("===============\n\n")
      cat(sprintf("Network Type: %s\n", geometry$characteristics$network_type))
      cat(sprintf("Density: %.1f%%\n", geometry$connectivity$density * 100))
      cat(sprintf("Mean Degree: %.2f\n", geometry$characteristics$mean_degree))
      cat(sprintf("Robustness Score: %.1f%%\n", geometry$robustness$robustness_score * 100))
    })

    # Download
    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0("nma_", input$plot_type, "_", Sys.Date(), ".html")
      },
      content = function(file) {
        p <- switch(input$plot_type,
                   network = create_interactive_network(nma_result, data, input),
                   forest = create_interactive_forest(nma_result, input),
                   funnel = create_interactive_funnel(nma_result, data, input),
                   rankings = create_interactive_rankings(nma_result, input),
                   heatmap = create_interactive_heatmap(nma_result, input),
                   network3d = create_3d_network(nma_result, data, input))

        htmlwidgets::saveWidget(p, file)
      }
    )
  }

  app <- shiny::shinyApp(ui, server)

  if (launch_browser) {
    shiny::runApp(app)
  }

  invisible(app)
}

#' Create Interactive Network Plot
#' @keywords internal
create_interactive_network <- function(nma_result, data, input) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required")
  }

  # Build network
  treatments <- nma_result$trts
  adj_matrix <- matrix(0, length(treatments), length(treatments),
                      dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])
    adj_matrix[t1, t2] <- adj_matrix[t1, t2] + 1
    adj_matrix[t2, t1] <- adj_matrix[t2, t1] + 1
  }

  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

  # Layout
  layout <- igraph::layout_with_fr(g)

  # Create edge list
  edge_list <- igraph::as_edgelist(g)
  edge_trace <- list()

  for (i in seq_len(nrow(edge_list))) {
    v0 <- which(treatments == edge_list[i, 1])
    v1 <- which(treatments == edge_list[i, 2])

    edge_trace[[i]] <- list(
      x = c(layout[v0, 1], layout[v1, 1]),
      y = c(layout[v0, 2], layout[v1, 2]),
      mode = "lines",
      line = list(width = input$edge_width, color = "rgba(125,125,125,0.5)"),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }

  # Node sizes by degree
  degrees <- igraph::degree(g)

  # Create plotly
  p <- plotly::plot_ly()

  # Add edges
  for (trace in edge_trace) {
    p <- plotly::add_trace(p, x = trace$x, y = trace$y,
                          mode = trace$mode, line = trace$line,
                          hoverinfo = trace$hoverinfo,
                          showlegend = trace$showlegend,
                          type = "scatter")
  }

  # Add nodes
  p <- plotly::add_trace(p,
    x = layout[, 1],
    y = layout[, 2],
    mode = "markers+text",
    text = if (input$show_labels) treatments else NULL,
    textposition = "top center",
    marker = list(
      size = input$node_size + degrees * 3,
      color = degrees,
      colorscale = "Viridis",
      showscale = TRUE,
      colorbar = list(title = "Degree")
    ),
    hovertext = paste0("Treatment: ", treatments, "<br>Degree: ", degrees),
    hoverinfo = "text",
    type = "scatter"
  )

  p <- plotly::layout(p,
    title = "Interactive Network Graph",
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    hovermode = "closest"
  )

  p
}

#' Create Interactive Forest Plot
#' @keywords internal
create_interactive_forest <- function(nma_result, input) {

  treatments <- rownames(nma_result$TE.random)
  effects <- nma_result$TE.random[, 1]
  lower <- nma_result$lower.random[, 1]
  upper <- nma_result$upper.random[, 1]

  # Sort by effect
  order_idx <- order(effects, decreasing = TRUE)
  treatments <- treatments[order_idx]
  effects <- effects[order_idx]
  lower <- lower[order_idx]
  upper <- upper[order_idx]

  p <- plotly::plot_ly()

  # Add confidence intervals
  if (input$show_ci) {
    for (i in seq_along(treatments)) {
      p <- plotly::add_trace(p,
        x = c(lower[i], upper[i]),
        y = c(i, i),
        mode = "lines",
        line = list(color = "gray", width = 2),
        showlegend = FALSE,
        hoverinfo = "skip",
        type = "scatter"
      )
    }
  }

  # Add points
  p <- plotly::add_trace(p,
    x = effects,
    y = seq_along(treatments),
    mode = "markers",
    marker = list(size = 12, color = "#2c7fb8"),
    hovertext = paste0(
      "Treatment: ", treatments, "<br>",
      "Effect: ", round(effects, 3), "<br>",
      "95% CI: [", round(lower, 3), ", ", round(upper, 3), "]"
    ),
    hoverinfo = "text",
    showlegend = FALSE,
    type = "scatter"
  )

  # Add reference line
  p <- plotly::add_trace(p,
    x = c(0, 0),
    y = c(0.5, length(treatments) + 0.5),
    mode = "lines",
    line = list(color = "red", width = 1, dash = "dash"),
    showlegend = FALSE,
    hoverinfo = "skip",
    type = "scatter"
  )

  p <- plotly::layout(p,
    title = "Interactive Forest Plot",
    xaxis = list(title = "Treatment Effect"),
    yaxis = list(
      title = "",
      ticktext = treatments,
      tickvals = seq_along(treatments),
      tickmode = "array"
    ),
    hovermode = "closest"
  )

  p
}

#' Create Interactive Funnel Plot
#' @keywords internal
create_interactive_funnel <- function(nma_result, data, input) {

  # Extract effects and standard errors
  effects <- data$TE
  se <- data$seTE

  p <- plotly::plot_ly()

  # Add funnel lines (pseudo 95% CI)
  max_se <- max(se) * 1.1
  se_seq <- seq(0, max_se, length.out = 100)

  p <- plotly::add_trace(p,
    x = 1.96 * se_seq,
    y = se_seq,
    mode = "lines",
    line = list(color = "gray", dash = "dash"),
    showlegend = FALSE,
    hoverinfo = "skip",
    type = "scatter"
  )

  p <- plotly::add_trace(p,
    x = -1.96 * se_seq,
    y = se_seq,
    mode = "lines",
    line = list(color = "gray", dash = "dash"),
    showlegend = FALSE,
    hoverinfo = "skip",
    type = "scatter"
  )

  # Add studies
  p <- plotly::add_trace(p,
    x = effects,
    y = se,
    mode = "markers",
    marker = list(size = 8, color = "#2c7fb8", opacity = 0.6),
    hovertext = paste0(
      "Study: ", data$studlab, "<br>",
      "Effect: ", round(effects, 3), "<br>",
      "SE: ", round(se, 3)
    ),
    hoverinfo = "text",
    showlegend = FALSE,
    type = "scatter"
  )

  p <- plotly::layout(p,
    title = "Interactive Funnel Plot",
    xaxis = list(title = "Effect Size", zeroline = TRUE),
    yaxis = list(title = "Standard Error", autorange = "reversed"),
    hovermode = "closest"
  )

  p
}

#' Create Interactive Rankings Plot
#' @keywords internal
create_interactive_rankings <- function(nma_result, input) {

  sucra <- calculate_sucra(nma_result)

  treatments <- names(sucra$sucra_scores)
  scores <- sucra$sucra_scores * 100

  # Sort
  order_idx <- order(scores, decreasing = TRUE)
  treatments <- treatments[order_idx]
  scores <- scores[order_idx]

  p <- plotly::plot_ly(
    x = treatments,
    y = scores,
    type = "bar",
    marker = list(
      color = scores,
      colorscale = "RdYlGn",
      showscale = TRUE,
      colorbar = list(title = "SUCRA %")
    ),
    hovertext = paste0(
      "Treatment: ", treatments, "<br>",
      "SUCRA: ", round(scores, 1), "%<br>",
      "Rank: ", seq_along(treatments)
    ),
    hoverinfo = "text"
  )

  p <- plotly::layout(p,
    title = "Treatment Rankings (SUCRA Scores)",
    xaxis = list(title = "Treatment"),
    yaxis = list(title = "SUCRA Score (%)", range = c(0, 100)),
    hovermode = "closest"
  )

  p
}

#' Create Interactive Heatmap
#' @keywords internal
create_interactive_heatmap <- function(nma_result, input) {

  effects_matrix <- nma_result$TE.random
  treatments <- rownames(effects_matrix)

  p <- plotly::plot_ly(
    x = treatments,
    y = treatments,
    z = effects_matrix,
    type = "heatmap",
    colorscale = "RdBu",
    zmid = 0,
    hovertemplate = paste0(
      "Row: %{y}<br>",
      "Column: %{x}<br>",
      "Effect: %{z:.3f}<extra></extra>"
    )
  )

  p <- plotly::layout(p,
    title = "Treatment Effects Heatmap",
    xaxis = list(title = "Treatment"),
    yaxis = list(title = "Treatment")
  )

  p
}

#' Create 3D Network Visualization
#' @keywords internal
create_3d_network <- function(nma_result, data, input) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required")
  }

  # Build network
  treatments <- nma_result$trts
  adj_matrix <- matrix(0, length(treatments), length(treatments),
                      dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])
    adj_matrix[t1, t2] <- adj_matrix[t1, t2] + 1
    adj_matrix[t2, t1] <- adj_matrix[t2, t1] + 1
  }

  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

  # 3D layout
  layout_3d <- igraph::layout_with_fr(g, dim = 3)

  # Edges
  edge_list <- igraph::as_edgelist(g)
  edge_traces <- list()

  for (i in seq_len(nrow(edge_list))) {
    v0 <- which(treatments == edge_list[i, 1])
    v1 <- which(treatments == edge_list[i, 2])

    edge_traces[[i]] <- list(
      x = c(layout_3d[v0, 1], layout_3d[v1, 1], NA),
      y = c(layout_3d[v0, 2], layout_3d[v1, 2], NA),
      z = c(layout_3d[v0, 3], layout_3d[v1, 3], NA),
      mode = "lines",
      line = list(color = "gray", width = 2),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }

  # Create plot
  p <- plotly::plot_ly()

  # Add edges
  for (trace in edge_traces) {
    p <- plotly::add_trace(p,
      x = trace$x, y = trace$y, z = trace$z,
      mode = trace$mode,
      line = trace$line,
      hoverinfo = trace$hoverinfo,
      showlegend = trace$showlegend,
      type = "scatter3d"
    )
  }

  # Add nodes
  degrees <- igraph::degree(g)

  p <- plotly::add_trace(p,
    x = layout_3d[, 1],
    y = layout_3d[, 2],
    z = layout_3d[, 3],
    mode = "markers+text",
    text = if (input$show_labels) treatments else NULL,
    marker = list(
      size = input$node_size + degrees * 2,
      color = degrees,
      colorscale = "Viridis",
      showscale = TRUE,
      colorbar = list(title = "Degree")
    ),
    hovertext = paste0("Treatment: ", treatments, "<br>Degree: ", degrees),
    hoverinfo = "text",
    type = "scatter3d"
  )

  p <- plotly::layout(p,
    title = "3D Network Visualization",
    scene = list(
      xaxis = list(showgrid = FALSE, showticklabels = FALSE),
      yaxis = list(showgrid = FALSE, showticklabels = FALSE),
      zaxis = list(showgrid = FALSE, showticklabels = FALSE)
    )
  )

  p
}

#' Export Interactive Plot to HTML
#'
#' Save interactive plotly visualization to standalone HTML file.
#'
#' @param plot Plotly plot object
#' @param filename Output filename
#' @return File path
#' @export
export_interactive_plot <- function(plot, filename = "interactive_plot.html") {

  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("Package 'htmlwidgets' required")
  }

  htmlwidgets::saveWidget(plot, filename)

  message(sprintf("✓ Interactive plot saved: %s", filename))

  invisible(filename)
}
