#' Network Geometry and Connectivity Analysis
#'
#' Advanced tools for analyzing network meta-analysis structure, connectivity,
#' and geometry. Essential for assessing network quality and identifying
#' potential issues.
#'
#' @name network_geometry
#' @references
#' Salanti G, et al. (2008). Evaluation of networks of randomized trials.
#' Statistical Methods in Medical Research, 17(3):279-301.
#'
#' Rücker G, Schwarzer G (2014). Reduce dimension or reduce weights? Comparing
#' two approaches to multi-arm studies in network meta-analysis. Statistics in
#' Medicine, 33(25):4353-4369.
NULL

#' Comprehensive Network Geometry Analysis
#'
#' Analyze all aspects of network structure including connectivity, diversity,
#' redundancy, and potential for bias.
#'
#' @param data Pairwise comparison data
#' @param reference Reference treatment (optional)
#' @return network_geometry object
#' @export
#' @examples
#' \dontrun{
#' geometry <- analyze_network_geometry(nma_data)
#' print(geometry)
#' plot(geometry)
#' }
analyze_network_geometry <- function(data, reference = NULL) {

  message("Analyzing network geometry and connectivity...")

  # Extract treatments
  treatments <- unique(c(as.character(data$treat1), as.character(data$treat2)))
  n_treatments <- length(treatments)

  # Build adjacency matrix
  adj_matrix <- build_adjacency_matrix(data, treatments)

  # Connectivity analysis
  connectivity <- analyze_connectivity(adj_matrix, treatments)

  # Network characteristics
  characteristics <- calculate_network_characteristics(data, adj_matrix, treatments)

  # Diversity measures
  diversity <- calculate_diversity_measures(data, treatments)

  # Multi-arm studies
  multiarm <- analyze_multiarm_studies(data)

  # Graph theory metrics
  graph_metrics <- calculate_graph_metrics(adj_matrix, treatments)

  # Robustness assessment
  robustness <- assess_network_robustness(adj_matrix, data, treatments)

  # Evidence flow
  evidence_flow <- analyze_evidence_flow(adj_matrix, data, treatments, reference)

  result <- list(
    treatments = treatments,
    n_treatments = n_treatments,
    n_studies = length(unique(data$studlab)),
    n_comparisons = nrow(data),
    adjacency_matrix = adj_matrix,
    connectivity = connectivity,
    characteristics = characteristics,
    diversity = diversity,
    multiarm = multiarm,
    graph_metrics = graph_metrics,
    robustness = robustness,
    evidence_flow = evidence_flow
  )

  class(result) <- c("network_geometry", "list")
  result
}

#' Build Adjacency Matrix
#' @keywords internal
build_adjacency_matrix <- function(data, treatments) {

  n <- length(treatments)
  adj <- matrix(0, nrow = n, ncol = n, dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])

    adj[t1, t2] <- adj[t1, t2] + 1
    adj[t2, t1] <- adj[t2, t1] + 1
  }

  adj
}

#' Analyze Network Connectivity
#' @keywords internal
analyze_connectivity <- function(adj_matrix, treatments) {

  # Check if network is connected
  binary_adj <- (adj_matrix > 0) * 1

  # BFS to check connectivity
  visited <- rep(FALSE, length(treatments))
  names(visited) <- treatments

  queue <- treatments[1]
  visited[treatments[1]] <- TRUE

  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]

    neighbors <- treatments[binary_adj[current, ] > 0]
    for (neighbor in neighbors) {
      if (!visited[neighbor]) {
        visited[neighbor] <- TRUE
        queue <- c(queue, neighbor)
      }
    }
  }

  is_connected <- all(visited)

  # Identify components if disconnected
  if (!is_connected) {
    components <- identify_components(binary_adj, treatments)
  } else {
    components <- list(treatments)
  }

  # Calculate connectivity metrics
  n <- length(treatments)
  n_possible_edges <- n * (n - 1) / 2
  n_observed_edges <- sum(adj_matrix > 0) / 2

  density <- n_observed_edges / n_possible_edges

  # Minimum spanning tree
  if (is_connected) {
    mst_edges <- n - 1
    redundancy <- n_observed_edges - mst_edges
  } else {
    mst_edges <- NA
    redundancy <- NA
  }

  list(
    is_connected = is_connected,
    n_components = length(components),
    components = components,
    density = density,
    n_edges = n_observed_edges,
    n_possible_edges = n_possible_edges,
    redundancy = redundancy,
    interpretation = if (is_connected) {
      sprintf("Network is connected with %.1f%% density", density * 100)
    } else {
      sprintf("Network is DISCONNECTED with %d component(s)", length(components))
    }
  )
}

#' Identify Network Components
#' @keywords internal
identify_components <- function(adj_matrix, treatments) {

  visited <- rep(FALSE, length(treatments))
  names(visited) <- treatments
  components <- list()

  for (start_node in treatments) {
    if (!visited[start_node]) {
      # BFS from this node
      component <- character()
      queue <- start_node
      visited[start_node] <- TRUE

      while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]
        component <- c(component, current)

        neighbors <- treatments[adj_matrix[current, ] > 0]
        for (neighbor in neighbors) {
          if (!visited[neighbor]) {
            visited[neighbor] <- TRUE
            queue <- c(queue, neighbor)
          }
        }
      }

      components[[length(components) + 1]] <- component
    }
  }

  components
}

#' Calculate Network Characteristics
#' @keywords internal
calculate_network_characteristics <- function(data, adj_matrix, treatments) {

  # Degree distribution
  degrees <- rowSums(adj_matrix > 0)

  # Hub treatments (high degree)
  hub_threshold <- quantile(degrees, 0.75)
  hubs <- treatments[degrees >= hub_threshold]

  # Peripheral treatments (low degree)
  peripheral_threshold <- quantile(degrees, 0.25)
  peripheral <- treatments[degrees <= peripheral_threshold]

  # Star network detection (one central hub connected to all)
  max_degree <- max(degrees)
  is_star <- (max_degree == length(treatments) - 1) && sum(degrees == 1) >= length(treatments) - 2

  # Line network detection
  is_line <- all(degrees <= 2) && sum(degrees == 1) == 2

  # Complete network (fully connected)
  n <- length(treatments)
  n_possible <- n * (n - 1) / 2
  n_actual <- sum(adj_matrix > 0) / 2
  is_complete <- (n_actual == n_possible)

  list(
    degrees = degrees,
    mean_degree = mean(degrees),
    sd_degree = sd(degrees),
    max_degree = max_degree,
    min_degree = min(degrees),
    hubs = hubs,
    peripheral = peripheral,
    is_star = is_star,
    is_line = is_line,
    is_complete = is_complete,
    network_type = if (is_complete) "Complete" else if (is_star) "Star" else if (is_line) "Line" else "Mixed"
  )
}

#' Calculate Diversity Measures
#' @keywords internal
calculate_diversity_measures <- function(data, treatments) {

  # Treatment diversity (Shannon entropy)
  treat_counts <- table(c(as.character(data$treat1), as.character(data$treat2)))
  treat_props <- treat_counts / sum(treat_counts)
  shannon_entropy <- -sum(treat_props * log(treat_props))
  max_entropy <- log(length(treatments))
  diversity_index <- shannon_entropy / max_entropy

  # Comparison diversity
  data$comparison <- paste(
    pmin(as.character(data$treat1), as.character(data$treat2)),
    pmax(as.character(data$treat1), as.character(data$treat2)),
    sep = " vs "
  )

  comparison_counts <- table(data$comparison)

  # Gini coefficient for comparison imbalance
  sorted_counts <- sort(comparison_counts)
  n <- length(sorted_counts)
  gini <- 2 * sum((1:n) * sorted_counts) / (n * sum(sorted_counts)) - (n + 1) / n

  list(
    shannon_entropy = shannon_entropy,
    max_entropy = max_entropy,
    diversity_index = diversity_index,
    gini_coefficient = gini,
    n_unique_comparisons = length(comparison_counts),
    comparison_balance = 1 - gini,
    interpretation = sprintf(
      "Diversity: %.1f%% (%.2f) | Balance: %.1f%%",
      diversity_index * 100, shannon_entropy,
      (1 - gini) * 100
    )
  )
}

#' Analyze Multi-arm Studies
#' @keywords internal
analyze_multiarm_studies <- function(data) {

  # Identify multi-arm studies
  studies <- unique(data$studlab)

  multiarm_info <- lapply(studies, function(study) {
    study_data <- data[data$studlab == study, ]
    study_treatments <- unique(c(as.character(study_data$treat1),
                                as.character(study_data$treat2)))

    n_arms <- length(study_treatments)
    n_comparisons <- nrow(study_data)

    list(
      study = study,
      n_arms = n_arms,
      n_comparisons = n_comparisons,
      treatments = study_treatments,
      is_multiarm = n_arms > 2
    )
  })

  n_multiarm <- sum(sapply(multiarm_info, function(x) x$is_multiarm))

  # Multi-arm proportions
  all_arms <- sapply(multiarm_info, function(x) x$n_arms)

  list(
    n_studies = length(studies),
    n_multiarm = n_multiarm,
    proportion_multiarm = n_multiarm / length(studies),
    arm_distribution = table(all_arms),
    max_arms = max(all_arms),
    mean_arms = mean(all_arms),
    multiarm_studies = multiarm_info[sapply(multiarm_info, function(x) x$is_multiarm)]
  )
}

#' Calculate Graph Theory Metrics
#' @keywords internal
calculate_graph_metrics <- function(adj_matrix, treatments) {

  # Use igraph if available
  if (requireNamespace("igraph", quietly = TRUE)) {

    binary_adj <- (adj_matrix > 0) * 1
    g <- igraph::graph_from_adjacency_matrix(binary_adj, mode = "undirected")

    # Betweenness centrality
    betweenness <- igraph::betweenness(g, normalized = TRUE)
    names(betweenness) <- treatments

    # Closeness centrality
    closeness <- igraph::closeness(g, normalized = TRUE)
    names(closeness) <- treatments

    # Eigenvector centrality
    eigenvector <- igraph::eigen_centrality(g)$vector
    names(eigenvector) <- treatments

    # Clustering coefficient
    clustering <- igraph::transitivity(g, type = "local")
    names(clustering) <- treatments
    clustering[is.nan(clustering)] <- 0

    # Average path length
    if (igraph::is_connected(g)) {
      avg_path_length <- igraph::average.path.length(g)
      diameter <- igraph::diameter(g)
    } else {
      avg_path_length <- Inf
      diameter <- Inf
    }

    # Network centralization
    degree_cent <- igraph::centr_degree(g)$centralization

    list(
      betweenness = betweenness,
      closeness = closeness,
      eigenvector = eigenvector,
      clustering = clustering,
      mean_clustering = mean(clustering),
      avg_path_length = avg_path_length,
      diameter = diameter,
      centralization = degree_cent,
      most_central = names(which.max(betweenness))
    )

  } else {
    message("Install 'igraph' package for advanced graph metrics")
    list(note = "igraph not available")
  }
}

#' Assess Network Robustness
#' @keywords internal
assess_network_robustness <- function(adj_matrix, data, treatments) {

  binary_adj <- (adj_matrix > 0) * 1

  # Test robustness to node removal
  critical_nodes <- character()

  for (node in treatments) {
    # Remove node
    adj_reduced <- binary_adj[-which(treatments == node), -which(treatments == node)]
    treatments_reduced <- treatments[treatments != node]

    # Check if still connected
    if (length(treatments_reduced) > 1) {
      components <- identify_components(adj_reduced, treatments_reduced)
      if (length(components) > 1) {
        critical_nodes <- c(critical_nodes, node)
      }
    }
  }

  # Test robustness to edge removal
  edges <- which(upper.tri(adj_matrix) & adj_matrix > 0, arr.ind = TRUE)
  critical_edges <- list()

  if (nrow(edges) > 0) {
    for (i in seq_len(min(nrow(edges), 10))) {  # Test up to 10 edges
      edge <- edges[i, ]
      adj_reduced <- binary_adj
      adj_reduced[edge[1], edge[2]] <- 0
      adj_reduced[edge[2], edge[1]] <- 0

      components <- identify_components(adj_reduced, treatments)
      if (length(components) > 1) {
        critical_edges[[length(critical_edges) + 1]] <- c(
          treatments[edge[1]], treatments[edge[2]]
        )
      }
    }
  }

  list(
    n_critical_nodes = length(critical_nodes),
    critical_nodes = critical_nodes,
    n_critical_edges = length(critical_edges),
    critical_edges = critical_edges,
    robustness_score = 1 - length(critical_nodes) / length(treatments),
    interpretation = if (length(critical_nodes) == 0) {
      "Network is robust - no single node removal disconnects it"
    } else {
      sprintf("Network has %d critical node(s) - removal would disconnect network",
              length(critical_nodes))
    }
  )
}

#' Analyze Evidence Flow
#' @keywords internal
analyze_evidence_flow <- function(adj_matrix, data, treatments, reference) {

  if (is.null(reference)) {
    # Choose most connected treatment as reference
    degrees <- rowSums(adj_matrix > 0)
    reference <- treatments[which.max(degrees)]
  }

  # Calculate direct vs indirect evidence proportions
  direct_comparisons <- sum(adj_matrix > 0) / 2

  # All pairwise comparisons in connected network
  n <- length(treatments)
  all_comparisons <- n * (n - 1) / 2

  indirect_comparisons <- all_comparisons - direct_comparisons
  proportion_indirect <- indirect_comparisons / all_comparisons

  # Evidence contribution by treatment
  evidence_contribution <- rowSums(adj_matrix)
  names(evidence_contribution) <- treatments

  # Rank treatments by evidence contribution
  ranked <- sort(evidence_contribution, decreasing = TRUE)

  list(
    reference = reference,
    n_direct_comparisons = direct_comparisons,
    n_indirect_comparisons = indirect_comparisons,
    proportion_indirect = proportion_indirect,
    evidence_contribution = evidence_contribution,
    ranked_by_evidence = ranked,
    top_contributor = names(ranked)[1],
    interpretation = sprintf(
      "%.1f%% of comparisons rely on indirect evidence (reference: %s)",
      proportion_indirect * 100, reference
    )
  )
}

#' Print Network Geometry Analysis
#' @param x network_geometry object
#' @param ... Additional arguments
#' @export
print.network_geometry <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Network Geometry Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Network Size:\n"))
  cat(sprintf("  • Treatments: %d\n", x$n_treatments))
  cat(sprintf("  • Studies: %d\n", x$n_studies))
  cat(sprintf("  • Comparisons: %d\n\n", x$n_comparisons))

  cat("Connectivity:\n")
  cat(sprintf("  • Status: %s\n",
             if (x$connectivity$is_connected) "Connected ✓" else "DISCONNECTED ✗"))
  cat(sprintf("  • Density: %.1f%%\n", x$connectivity$density * 100))
  cat(sprintf("  • Edges: %d / %d possible\n",
             x$connectivity$n_edges, x$connectivity$n_possible_edges))
  if (!is.na(x$connectivity$redundancy)) {
    cat(sprintf("  • Redundant edges: %d\n", x$connectivity$redundancy))
  }
  cat("\n")

  cat("Network Structure:\n")
  cat(sprintf("  • Type: %s\n", x$characteristics$network_type))
  cat(sprintf("  • Mean degree: %.2f (SD: %.2f)\n",
             x$characteristics$mean_degree, x$characteristics$sd_degree))
  cat(sprintf("  • Hub treatments (%d): %s\n",
             length(x$characteristics$hubs),
             paste(head(x$characteristics$hubs, 3), collapse = ", ")))
  cat("\n")

  cat("Diversity & Balance:\n")
  cat(sprintf("  %s\n", x$diversity$interpretation))
  cat(sprintf("  • Unique comparisons: %d\n", x$diversity$n_unique_comparisons))
  cat("\n")

  cat("Multi-arm Studies:\n")
  cat(sprintf("  • Multi-arm studies: %d / %d (%.1f%%)\n",
             x$multiarm$n_multiarm, x$multiarm$n_studies,
             x$multiarm$proportion_multiarm * 100))
  cat(sprintf("  • Mean arms per study: %.2f\n", x$multiarm$mean_arms))
  cat(sprintf("  • Max arms: %d\n\n", x$multiarm$max_arms))

  if (!is.null(x$graph_metrics$most_central)) {
    cat("Graph Metrics:\n")
    cat(sprintf("  • Most central: %s\n", x$graph_metrics$most_central))
    cat(sprintf("  • Mean clustering: %.3f\n", x$graph_metrics$mean_clustering))
    if (is.finite(x$graph_metrics$avg_path_length)) {
      cat(sprintf("  • Avg path length: %.2f\n", x$graph_metrics$avg_path_length))
      cat(sprintf("  • Diameter: %d\n", x$graph_metrics$diameter))
    }
    cat("\n")
  }

  cat("Robustness:\n")
  cat(sprintf("  %s\n", x$robustness$interpretation))
  cat(sprintf("  • Robustness score: %.1f%%\n\n", x$robustness$robustness_score * 100))

  cat("Evidence Flow:\n")
  cat(sprintf("  %s\n", x$evidence_flow$interpretation))
  cat(sprintf("  • Top contributor: %s\n\n", x$evidence_flow$top_contributor))

  cat("Use plot() to visualize network geometry\n\n")

  invisible(x)
}

#' Plot Network Geometry
#' @param x network_geometry object
#' @param type Type of plot
#' @param ... Additional arguments
#' @export
plot.network_geometry <- function(x, type = c("network", "degree", "centrality", "robustness"), ...) {

  type <- match.arg(type)

  if (type == "network") {
    plot_network_graph(x)
  } else if (type == "degree") {
    plot_degree_distribution(x)
  } else if (type == "centrality") {
    plot_centrality_measures(x)
  } else if (type == "robustness") {
    plot_robustness_analysis(x)
  }
}

#' Plot Network Graph
#' @keywords internal
plot_network_graph <- function(geometry) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required for network plotting")
  }

  g <- igraph::graph_from_adjacency_matrix(
    geometry$adjacency_matrix,
    mode = "undirected",
    weighted = TRUE
  )

  # Node colors by degree
  degrees <- geometry$characteristics$degrees
  node_colors <- grDevices::colorRampPalette(c("#fee08b", "#d73027"))(max(degrees))
  igraph::V(g)$color <- node_colors[degrees]

  # Node sizes by evidence contribution
  igraph::V(g)$size <- sqrt(geometry$evidence_flow$evidence_contribution) * 5

  # Edge widths by number of studies
  igraph::E(g)$width <- sqrt(igraph::E(g)$weight)

  plot(g,
       layout = igraph::layout_with_fr(g),
       vertex.label.cex = 0.8,
       vertex.label.color = "black",
       main = sprintf("Network Graph (%s)", geometry$characteristics$network_type))

  legend("bottomright",
         legend = c("High degree", "Low degree"),
         fill = c("#d73027", "#fee08b"),
         cex = 0.8,
         title = "Node Degree")
}

#' Plot Degree Distribution
#' @keywords internal
plot_degree_distribution <- function(geometry) {

  degrees <- geometry$characteristics$degrees

  par(mfrow = c(1, 2))

  # Histogram
  hist(degrees,
       breaks = seq(0, max(degrees) + 1) - 0.5,
       col = "#2c7fb8",
       border = "white",
       main = "Degree Distribution",
       xlab = "Degree (Number of Connections)",
       ylab = "Frequency")

  abline(v = mean(degrees), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = sprintf("Mean: %.2f", mean(degrees)),
         col = "red", lty = 2, lwd = 2, cex = 0.8)

  # Bar plot
  sorted_degrees <- sort(degrees, decreasing = TRUE)
  barplot(sorted_degrees,
          names.arg = names(sorted_degrees),
          col = "#2c7fb8",
          border = NA,
          main = "Treatment Connectivity",
          ylab = "Degree",
          las = 2,
          cex.names = 0.7)

  par(mfrow = c(1, 1))
}

#' Plot Centrality Measures
#' @keywords internal
plot_centrality_measures <- function(geometry) {

  if (is.null(geometry$graph_metrics$betweenness)) {
    message("Graph metrics not available. Install 'igraph' package.")
    return(invisible(NULL))
  }

  # Create data frame
  cent_data <- data.frame(
    Treatment = names(geometry$graph_metrics$betweenness),
    Betweenness = geometry$graph_metrics$betweenness,
    Closeness = geometry$graph_metrics$closeness,
    Eigenvector = geometry$graph_metrics$eigenvector,
    stringsAsFactors = FALSE
  )

  # Sort by betweenness
  cent_data <- cent_data[order(-cent_data$Betweenness), ]

  par(mfrow = c(1, 3))

  # Betweenness
  barplot(cent_data$Betweenness,
          names.arg = cent_data$Treatment,
          col = "#2c7fb8",
          main = "Betweenness Centrality",
          ylab = "Normalized Betweenness",
          las = 2,
          cex.names = 0.6)

  # Closeness
  barplot(cent_data$Closeness,
          names.arg = cent_data$Treatment,
          col = "#d7191c",
          main = "Closeness Centrality",
          ylab = "Normalized Closeness",
          las = 2,
          cex.names = 0.6)

  # Eigenvector
  barplot(cent_data$Eigenvector,
          names.arg = cent_data$Treatment,
          col = "#1a9850",
          main = "Eigenvector Centrality",
          ylab = "Normalized Eigenvector",
          las = 2,
          cex.names = 0.6)

  par(mfrow = c(1, 1))
}

#' Plot Robustness Analysis
#' @keywords internal
plot_robustness_analysis <- function(geometry) {

  # Create robustness visualization
  treatments <- geometry$treatments
  n <- length(treatments)

  # Robustness scores
  robust_score <- geometry$robustness$robustness_score
  critical_nodes <- geometry$robustness$critical_nodes

  # Create color coding
  node_colors <- rep("#1a9850", n)  # Green for robust
  names(node_colors) <- treatments

  if (length(critical_nodes) > 0) {
    node_colors[critical_nodes] <- "#d73027"  # Red for critical
  }

  par(mfrow = c(1, 2))

  # Bar plot of node importance
  barplot(geometry$evidence_flow$evidence_contribution,
          col = node_colors,
          main = "Node Importance & Criticality",
          ylab = "Evidence Contribution",
          las = 2,
          cex.names = 0.7)

  legend("topright",
         legend = c("Robust", "Critical"),
         fill = c("#1a9850", "#d73027"),
         cex = 0.8)

  # Robustness score gauge
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "",
       main = "Network Robustness Score")

  rect(0.2, 0.3, 0.8, 0.5, col = "gray90", border = NA)
  rect(0.2, 0.3, 0.2 + robust_score * 0.6, 0.5,
       col = ifelse(robust_score > 0.7, "#1a9850", "#d73027"),
       border = NA)

  text(0.5, 0.6, sprintf("%.1f%%", robust_score * 100),
       cex = 2, font = 2)

  text(0.5, 0.15, geometry$robustness$interpretation,
       cex = 0.8)

  par(mfrow = c(1, 1))
}
