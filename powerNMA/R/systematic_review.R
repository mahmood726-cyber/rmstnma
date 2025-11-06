#' Generate PRISMA-NMA Checklist
#'
#' Creates a PRISMA-NMA reporting checklist for systematic reviews
#'
#' @param results powerNMA results object
#' @param protocol_info List with protocol registration, dates, etc.
#' @return Data frame with PRISMA-NMA checklist
#' @export
#' @examples
#' \dontrun{
#' protocol <- list(
#'   registration = "PROSPERO CRD42025123456",
#'   date = "2025-01-15"
#' )
#' checklist <- generate_prisma_nma_checklist(results, protocol)
#' }
generate_prisma_nma_checklist <- function(results, protocol_info = NULL) {
  checklist <- tibble::tribble(
    ~Section, ~Item, ~Description, ~Completed, ~Location,

    # Title and Abstract
    "Title", "1", "Identify as NMA", TRUE, "Title",
    "Abstract", "2", "Structured summary", NA, "Abstract",

    # Introduction
    "Rationale", "3", "Describe rationale", NA, "Introduction",
    "Objectives", "4", "Provide PICO", NA, "Introduction",

    # Methods
    "Protocol", "5", "Registration information",
    !is.null(protocol_info$registration), protocol_info$registration %||% "",

    "Eligibility", "6", "Eligibility criteria", NA, "Methods",
    "Information sources", "7", "Databases searched", NA, "Methods",
    "Search", "8", "Search strategy", NA, "Methods",

    # Selection and data
    "Study selection", "9", "Selection process", NA, "Methods",
    "Data collection", "10", "Data extraction", NA, "Methods",

    # Geometry
    "Geometry", "S1", "Specify network structure",
    !is.null(results$results$geometry), "Results",

    # Statistical methods
    "Summary measures", "11", "Effect measure specified",
    !is.null(results$config$sm), paste("SM:", results$config$sm),

    "Synthesis methods", "12", "Methods for NMA",
    !is.null(results$results$main_nma), "Statistical Analysis",

    "Inconsistency", "S2", "Methods for assessing inconsistency",
    !is.null(results$results$global_inconsistency), "Methods/Results",

    "Heterogeneity", "13", "Heterogeneity assessment",
    !is.null(results$results$main_nma$tau), "Results",

    # Risk of bias
    "Risk of bias", "14", "Risk of bias assessment", NA, "",

    # Results
    "Study selection", "15", "Study selection results", NA, "Results",
    "Study characteristics", "16", "Study characteristics", NA, "Results",

    "Network geometry", "S3", "Network diagram",
    !is.null(results$results$geometry), "Results/Figure",

    "Risk of bias summary", "17", "RoB results", NA, "",

    "Individual study results", "18", "Individual effects", NA, "Results",

    "Network estimates", "S4", "Network estimates with CI",
    !is.null(results$results$main_nma), "Results/Table",

    "Heterogeneity results", "S5", "Heterogeneity statistics",
    !is.null(results$results$main_nma$I2.random), "Results",

    "Inconsistency results", "S6", "Inconsistency assessment results",
    !is.null(results$results$global_inconsistency), "Results",

    "Treatment ranking", "S7", "Ranking with uncertainties",
    !is.null(results$results$ranking), "Results"
  )

  checklist
}

#' Create Network Diagram for Publication
#'
#' Generate publication-quality network diagram
#'
#' @param data Pairwise data
#' @param output_file Output file path (optional)
#' @return ggplot object
#' @export
create_network_diagram <- function(data, output_file = NULL) {
  if (!has_pkg("igraph")) {
    stop("Package 'igraph' required for network diagrams")
  }

  # Create edge list
  edges <- data %>%
    dplyr::mutate(
      t1 = pmin(treat1, treat2),
      t2 = pmax(treat1, treat2)
    ) %>%
    dplyr::group_by(t1, t2) %>%
    dplyr::summarise(n_studies = dplyr::n(), .groups = "drop")

  # Create igraph object
  g <- igraph::graph_from_data_frame(
    edges[, c("t1", "t2")],
    directed = FALSE
  )

  # Add edge weights
  igraph::E(g)$weight <- edges$n_studies

  # Calculate node sizes by degree
  node_sizes <- igraph::degree(g)

  # Create plot data
  layout <- igraph::layout_nicely(g)

  plot_data <- tibble::tibble(
    treatment = igraph::V(g)$name,
    x = layout[, 1],
    y = layout[, 2],
    size = node_sizes
  )

  edge_data <- as.data.frame(igraph::as_edgelist(g))
  names(edge_data) <- c("from", "to")
  edge_data$n_studies <- edges$n_studies

  edge_data <- edge_data %>%
    dplyr::left_join(plot_data, by = c("from" = "treatment")) %>%
    dplyr::rename(x_from = x, y_from = y) %>%
    dplyr::left_join(plot_data, by = c("to" = "treatment")) %>%
    dplyr::rename(x_to = x, y_to = y)

  # Create ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                  size = n_studies),
      alpha = 0.6
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = x, y = y, size = size),
      shape = 21, fill = "lightblue", color = "black"
    ) +
    ggplot2::geom_text(
      data = plot_data,
      ggplot2::aes(x = x, y = y, label = treatment),
      vjust = -1.5
    ) +
    ggplot2::scale_size_continuous(name = "Number of studies",
                                  range = c(1, 5)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "Network Diagram")

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 10, height = 8)
  }

  p
}

#' Assess Risk of Bias Integration
#'
#' Structure for integrating RoB assessments
#'
#' @param data Study data
#' @param rob_data Risk of bias assessments
#' @return Data with RoB integrated
#' @export
integrate_rob_assessment <- function(data, rob_data) {
  if (!all(c("studlab", "rob_overall") %in% names(rob_data))) {
    stop("rob_data must have 'studlab' and 'rob_overall' columns")
  }

  # Merge RoB data
  data_rob <- data %>%
    dplyr::left_join(rob_data, by = "studlab")

  # Check coverage
  missing_rob <- data_rob %>%
    dplyr::filter(is.na(rob_overall)) %>%
    dplyr::distinct(studlab) %>%
    dplyr::pull(studlab)

  if (length(missing_rob) > 0) {
    warning(sprintf(
      "Missing RoB assessment for %d studies: %s",
      length(missing_rob),
      paste(head(missing_rob, 3), collapse = ", ")
    ))
  }

  data_rob
}

#' RoB-Restricted Sensitivity Analysis
#'
#' Perform NMA restricted to low risk of bias studies
#'
#' @param data Data with RoB assessments
#' @param ref_treatment Reference treatment
#' @param sm Summary measure
#' @param rob_threshold RoB level to include ("Low" or c("Low", "Some concerns"))
#' @return NMA results
#' @export
rob_sensitivity_analysis <- function(data, ref_treatment, sm = "HR",
                                    rob_threshold = "Low") {
  if (!"rob_overall" %in% names(data)) {
    stop("Data must include 'rob_overall' column from integrate_rob_assessment()")
  }

  data_subset <- data %>%
    dplyr::filter(rob_overall %in% rob_threshold)

  msg("RoB-restricted analysis: %d/%d studies (threshold: %s)",
      length(unique(data_subset$studlab)),
      length(unique(data$studlab)),
      paste(rob_threshold, collapse = "/"))

  if (length(unique(data_subset$studlab)) < 3) {
    warning("Too few studies for RoB-restricted analysis")
    return(NULL)
  }

  robust_netmeta(data_subset, ref_treatment, sm = sm, random = TRUE)
}

#' Generate Methods Text for Manuscript
#'
#' Creates formatted methods text for systematic review manuscript
#'
#' @param results powerNMA results
#' @param format Output format ("cochrane" or "prisma")
#' @return Character string with methods text
#' @export
generate_methods_text <- function(results, format = c("cochrane", "prisma")) {
  format <- match.arg(format)

  sm_full <- switch(results$config$sm,
    "HR" = "hazard ratio",
    "OR" = "odds ratio",
    "RR" = "risk ratio",
    "MD" = "mean difference",
    "SMD" = "standardized mean difference",
    results$config$sm
  )

  text <- sprintf(
    "We performed network meta-analysis using the powerNMA package (version 1.0.0) in R (version %s). We used a random-effects model with the %s as the summary measure. Network geometry was assessed by calculating network density and connectivity. Statistical inconsistency was evaluated using the design-by-treatment interaction test. Heterogeneity was quantified using tau-squared (τ²) and the I² statistic. Treatment ranking was performed using P-scores (frequentist analogue of SUCRA).",
    paste(R.version$major, R.version$minor, sep = "."),
    sm_full
  )

  if (!is.null(results$results$global_inconsistency)) {
    text <- paste(text, "Global inconsistency was assessed using decomposition of Cochran's Q by design.")
  }

  if (!is.null(results$results$local_inconsistency)) {
    text <- paste(text, "Local inconsistency was evaluated using node-splitting.")
  }

  if (!is.null(results$results$loo)) {
    text <- paste(text, "Sensitivity to individual studies was assessed using leave-one-out analysis.")
  }

  text
}

#' Generate Evidence Profile Table
#'
#' Create GRADE-style Summary of Findings table structure
#'
#' @param results powerNMA results
#' @param comparisons Vector of treatment comparisons to include
#' @return Data frame with evidence profile structure
#' @export
generate_evidence_profile <- function(results, comparisons = NULL) {
  if (is.null(results$results$main_nma)) {
    stop("No NMA results available")
  }

  nma <- results$results$main_nma

  # Extract effect estimates
  effects <- nma$TE.random
  lower <- nma$lower.random
  upper <- nma$upper.random

  # Create table structure
  profile <- tibble::tibble(
    Comparison = rownames(effects),
    Effect = sprintf("%.2f", effects[, 1]),
    CI_95 = sprintf("[%.2f, %.2f]", lower[, 1], upper[, 1]),
    Studies = NA_integer_,  # To be filled
    Participants = NA_integer_,  # To be filled
    Certainty = "○○○○"  # Placeholder - requires GRADE assessment
  )

  profile
}

#' Protocol Specification
#'
#' Document pre-specified analysis plan
#'
#' @param research_question Research question
#' @param treatments Vector of treatment names
#' @param outcome Primary outcome
#' @param effect_measure Effect measure (HR, OR, etc.)
#' @param reference Reference treatment
#' @param subgroups Planned subgroup variables
#' @param sensitivity Planned sensitivity analyses
#' @return Protocol object
#' @export
specify_protocol <- function(research_question,
                            treatments,
                            outcome,
                            effect_measure = "HR",
                            reference = NULL,
                            subgroups = NULL,
                            sensitivity = NULL) {
  protocol <- list(
    research_question = research_question,
    treatments = treatments,
    outcome = outcome,
    effect_measure = effect_measure,
    reference = reference %||% treatments[1],
    subgroups = subgroups,
    sensitivity = sensitivity,
    date_specified = Sys.Date(),
    version = "1.0"
  )

  class(protocol) <- c("powernma_protocol", "list")
  protocol
}

#' Check Protocol Adherence
#'
#' Verify that analysis matches pre-specified protocol
#'
#' @param results Analysis results
#' @param protocol Protocol object
#' @return List of deviations
#' @export
check_protocol_adherence <- function(results, protocol) {
  deviations <- list()

  # Check effect measure
  if (results$config$sm != protocol$effect_measure) {
    deviations$effect_measure <- sprintf(
      "Protocol specified '%s', analysis used '%s'",
      protocol$effect_measure,
      results$config$sm
    )
  }

  # Check reference treatment
  if (results$ref_treatment != protocol$reference) {
    deviations$reference <- sprintf(
      "Protocol specified '%s', analysis used '%s'",
      protocol$reference,
      results$ref_treatment
    )
  }

  if (length(deviations) > 0) {
    warning("Protocol deviations detected. See check_protocol_adherence() output.")
  } else {
    msg("Analysis adheres to protocol")
  }

  deviations
}

#' Export Analysis Archive
#'
#' Create complete reproducibility archive
#'
#' @param results Analysis results
#' @param data Original data
#' @param output_dir Output directory
#' @export
archive_analysis <- function(results, data, output_dir = "analysis_archive") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save results
  saveRDS(results, file.path(output_dir, "results.rds"))

  # Save data
  readr::write_csv(data, file.path(output_dir, "data.csv"))

  # Save configuration
  config_json <- jsonlite::toJSON(results$config, pretty = TRUE)
  writeLines(config_json, file.path(output_dir, "config.json"))

  # Save session info
  session <- utils::sessionInfo()
  writeLines(utils::capture.output(session),
            file.path(output_dir, "session_info.txt"))

  # Create README
  readme <- sprintf(
    "# Analysis Archive

Generated: %s
R Version: %s
Package: powerNMA version %s

## Contents
- results.rds: Complete analysis results
- data.csv: Input data
- config.json: Analysis configuration
- session_info.txt: R session information

## Reproducibility
To reproduce this analysis:
1. Load data: data <- readr::read_csv('data.csv')
2. Load config: config <- jsonlite::fromJSON('config.json')
3. Run: results <- run_powernma(data, config = config)
",
    Sys.time(),
    paste(R.version$major, R.version$minor, sep = "."),
    "1.0.0"
  )

  writeLines(readme, file.path(output_dir, "README.md"))

  msg("Analysis archived to: %s", output_dir)
  invisible(output_dir)
}
