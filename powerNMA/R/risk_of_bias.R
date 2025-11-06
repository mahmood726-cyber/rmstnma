#' Risk of Bias Assessment Integration for NMA
#'
#' Tools for incorporating risk of bias assessments into network meta-analysis
#' following Cochrane Risk of Bias 2 (RoB 2) and ROBINS-I frameworks.
#'
#' @name risk_of_bias
#' @references
#' Sterne JAC, et al. (2019). RoB 2: A revised tool for assessing risk of bias
#' in randomised trials. BMJ, 366:l4898.
NULL

#' Create Risk of Bias Data Template
#'
#' Generate a template data frame for entering risk of bias assessments.
#'
#' @param studies Vector of study IDs
#' @param tool Assessment tool: "RoB2" (for RCTs) or "ROBINS-I" (for observational)
#' @return Data frame template for ROB assessments
#' @export
#' @examples
#' studies <- c("Study1", "Study2", "Study3")
#' rob_template <- create_rob_template(studies, tool = "RoB2")
#' # Fill in the ratings and save
#' write.csv(rob_template, "rob_assessments.csv", row.names = FALSE)
create_rob_template <- function(studies, tool = c("RoB2", "ROBINS-I")) {

  tool <- match.arg(tool)

  if (tool == "RoB2") {
    # Cochrane RoB 2 domains
    template <- data.frame(
      studlab = studies,
      domain1_randomization = NA_character_,  # Low/Some concerns/High
      domain2_deviations = NA_character_,
      domain3_missing_data = NA_character_,
      domain4_outcome_measurement = NA_character_,
      domain5_selective_reporting = NA_character_,
      overall_rob = NA_character_,
      notes = NA_character_,
      stringsAsFactors = FALSE
    )

  } else {
    # ROBINS-I domains
    template <- data.frame(
      studlab = studies,
      domain1_confounding = NA_character_,  # Low/Moderate/Serious/Critical
      domain2_selection = NA_character_,
      domain3_classification = NA_character_,
      domain4_deviations = NA_character_,
      domain5_missing_data = NA_character_,
      domain6_outcome_measurement = NA_character_,
      domain7_selective_reporting = NA_character_,
      overall_rob = NA_character_,
      notes = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  template
}

#' Summarize Risk of Bias Assessments
#'
#' Create summary statistics and visualization of ROB assessments.
#'
#' @param rob_data Data frame with risk of bias assessments
#' @param tool Assessment tool used
#' @return rob_summary object
#' @export
#' @examples
#' \dontrun{
#' rob_data <- read.csv("rob_assessments.csv")
#' rob_summary <- summarize_rob(rob_data, tool = "RoB2")
#' print(rob_summary)
#' plot(rob_summary)
#' }
summarize_rob <- function(rob_data, tool = c("RoB2", "ROBINS-I")) {

  tool <- match.arg(tool)

  # Identify domain columns
  domain_cols <- grep("^domain", names(rob_data), value = TRUE)

  # Count ratings for each domain
  domain_summaries <- lapply(domain_cols, function(col) {
    counts <- table(rob_data[[col]], useNA = "ifany")
    data.frame(
      Domain = col,
      Low = counts["Low"],
      Some_Concerns = counts["Some concerns"],
      High = counts["High"],
      Total = nrow(rob_data),
      stringsAsFactors = FALSE
    )
  })

  domain_summary_table <- do.call(rbind, domain_summaries)

  # Overall ROB summary
  overall_counts <- table(rob_data$overall_rob, useNA = "ifany")

  # Calculate percentages
  pct_high <- (overall_counts["High"] / sum(overall_counts)) * 100
  pct_some <- (overall_counts["Some concerns"] / sum(overall_counts)) * 100
  pct_low <- (overall_counts["Low"] / sum(overall_counts)) * 100

  result <- list(
    tool = tool,
    n_studies = nrow(rob_data),
    domain_summary = domain_summary_table,
    overall_counts = overall_counts,
    pct_high_rob = pct_high,
    pct_some_concerns = pct_some,
    pct_low_rob = pct_low,
    rob_data = rob_data
  )

  class(result) <- c("rob_summary", "list")
  result
}

#' Create Risk of Bias Traffic Light Plot
#'
#' Visual summary of risk of bias assessments across studies and domains.
#'
#' @param rob_data Data frame with ROB assessments
#' @param tool Assessment tool
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' rob_data <- read.csv("rob_assessments.csv")
#' plot_rob_traffic_light(rob_data, tool = "RoB2")
#' }
plot_rob_traffic_light <- function(rob_data, tool = c("RoB2", "ROBINS-I")) {

  tool <- match.arg(tool)

  # Reshape data to long format
  domain_cols <- grep("^domain", names(rob_data), value = TRUE)

  plot_data <- tidyr::pivot_longer(
    rob_data,
    cols = all_of(c(domain_cols, "overall_rob")),
    names_to = "Domain",
    values_to = "Rating"
  )

  # Clean domain names
  plot_data$Domain <- gsub("domain[0-9]+_", "", plot_data$Domain)
  plot_data$Domain <- gsub("_", " ", plot_data$Domain)
  plot_data$Domain <- tools::toTitleCase(plot_data$Domain)

  # Order studies by overall ROB
  study_order <- rob_data$studlab[order(rob_data$overall_rob)]
  plot_data$studlab <- factor(plot_data$studlab, levels = study_order)

  # Color mapping
  color_map <- c(
    "Low" = "#1a9850",
    "Some concerns" = "#fee08b",
    "High" = "#d73027",
    "Moderate" = "#fdae61",
    "Serious" = "#f46d43",
    "Critical" = "#a50026"
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Domain, y = studlab, fill = Rating)) +
    ggplot2::geom_tile(color = "white", size = 0.5) +
    ggplot2::scale_fill_manual(values = color_map, na.value = "gray90") +
    ggplot2::labs(
      title = sprintf("Risk of Bias Assessment (%s)", tool),
      x = "Domain",
      y = "Study",
      fill = "Risk of Bias"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Sensitivity Analysis by Risk of Bias
#'
#' Assess impact of excluding studies with high risk of bias.
#'
#' @param data Pairwise NMA data
#' @param rob_data Risk of bias assessments
#' @param threshold ROB threshold: "exclude_high" or "exclude_high_and_some"
#' @param reference Reference treatment
#' @param sm Summary measure
#' @return rob_sensitivity object
#' @export
#' @examples
#' \dontrun{
#' rob_sens <- rob_sensitivity_analysis(
#'   data = nma_data,
#'   rob_data = rob_assessments,
#'   threshold = "exclude_high"
#' )
#' print(rob_sens)
#' }
rob_sensitivity_analysis <- function(data,
                                     rob_data,
                                     threshold = c("exclude_high", "exclude_high_and_some"),
                                     reference = NULL,
                                     sm = "MD") {

  threshold <- match.arg(threshold)

  # Merge ROB data with NMA data
  if (!"studlab" %in% names(rob_data)) {
    stop("rob_data must have 'studlab' column")
  }

  data_with_rob <- merge(data, rob_data[, c("studlab", "overall_rob")],
                        by = "studlab", all.x = TRUE)

  # Base NMA with all studies
  base_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  # Exclude studies based on threshold
  if (threshold == "exclude_high") {
    data_filtered <- data_with_rob[data_with_rob$overall_rob != "High" |
                                   is.na(data_with_rob$overall_rob), ]
    excluded_studies <- data_with_rob$studlab[data_with_rob$overall_rob == "High"]
  } else {
    data_filtered <- data_with_rob[data_with_rob$overall_rob == "Low" |
                                   is.na(data_with_rob$overall_rob), ]
    excluded_studies <- data_with_rob$studlab[data_with_rob$overall_rob %in%
                                              c("High", "Some concerns")]
  }

  excluded_studies <- unique(excluded_studies)

  if (nrow(data_filtered) < 3) {
    stop("Too few studies remaining after ROB filtering")
  }

  # Filtered NMA
  filtered_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data_filtered,
    sm = sm,
    reference.group = reference
  )

  # Compare estimates
  treatments <- setdiff(rownames(base_nma$TE.random), reference)

  comparison_table <- data.frame(
    Treatment = treatments,
    All_Studies_Effect = base_nma$TE.random[treatments, reference],
    All_Studies_SE = base_nma$seTE.random[treatments, reference],
    Low_ROB_Effect = filtered_nma$TE.random[treatments, reference],
    Low_ROB_SE = filtered_nma$seTE.random[treatments, reference],
    Difference = abs(base_nma$TE.random[treatments, reference] -
                    filtered_nma$TE.random[treatments, reference]),
    stringsAsFactors = FALSE
  )

  result <- list(
    base_nma = base_nma,
    filtered_nma = filtered_nma,
    comparison = comparison_table,
    excluded_studies = excluded_studies,
    n_excluded = length(excluded_studies),
    threshold = threshold
  )

  class(result) <- c("rob_sensitivity", "list")
  result
}

#' Weight Studies by Risk of Bias
#'
#' Apply down-weighting to studies with higher risk of bias.
#'
#' @param data Pairwise NMA data
#' @param rob_data Risk of bias assessments
#' @param weights Named vector of weights (e.g., Low=1.0, Some=0.5, High=0.25)
#' @return Data frame with adjusted weights
#' @export
#' @examples
#' \dontrun{
#' weighted_data <- weight_by_rob(
#'   data = nma_data,
#'   rob_data = rob_assessments,
#'   weights = c("Low" = 1.0, "Some concerns" = 0.5, "High" = 0.25)
#' )
#' }
weight_by_rob <- function(data,
                          rob_data,
                          weights = c("Low" = 1.0, "Some concerns" = 0.5, "High" = 0.25)) {

  # Merge ROB
  data_with_rob <- merge(data, rob_data[, c("studlab", "overall_rob")],
                        by = "studlab", all.x = TRUE)

  # Apply weights
  data_with_rob$rob_weight <- weights[data_with_rob$overall_rob]

  # If no ROB assessment, use weight = 1.0
  data_with_rob$rob_weight[is.na(data_with_rob$rob_weight)] <- 1.0

  # Adjust standard errors (down-weight by increasing SE)
  data_with_rob$seTE_adjusted <- data_with_rob$seTE / sqrt(data_with_rob$rob_weight)

  data_with_rob
}

#' Print ROB Summary
#'
#' @param x rob_summary object
#' @param ... Additional arguments
#' @export
print.rob_summary <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  Risk of Bias Summary (%s)\n", x$tool))
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Studies assessed: %d\n\n", x$n_studies))

  cat("Overall Risk of Bias:\n")
  cat(sprintf("  • Low: %d (%.1f%%)\n",
             x$overall_counts["Low"], x$pct_low_rob))
  cat(sprintf("  • Some concerns: %d (%.1f%%)\n",
             x$overall_counts["Some concerns"], x$pct_some_concerns))
  cat(sprintf("  • High: %d (%.1f%%)\n\n",
             x$overall_counts["High"], x$pct_high_rob))

  cat("Domain Summary:\n\n")
  print(x$domain_summary, row.names = FALSE)

  cat("\n")
  cat("Use plot() to visualize traffic light plot\n\n")

  invisible(x)
}

#' Plot ROB Summary
#'
#' @param x rob_summary object
#' @param ... Additional arguments
#' @export
plot.rob_summary <- function(x, ...) {
  plot_rob_traffic_light(x$rob_data, tool = x$tool)
}

#' Print ROB Sensitivity Results
#'
#' @param x rob_sensitivity object
#' @param ... Additional arguments
#' @export
print.rob_sensitivity <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Risk of Bias Sensitivity Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Threshold: %s\n", x$threshold))
  cat(sprintf("Studies excluded: %d\n", x$n_excluded))
  if (x$n_excluded > 0 && length(x$excluded_studies) <= 10) {
    cat("  ", paste(x$excluded_studies, collapse = ", "), "\n")
  }
  cat("\n")

  cat("Comparison of Treatment Effects:\n\n")
  print(x$comparison, row.names = FALSE, digits = 3)

  cat("\n")
  max_diff <- max(x$comparison$Difference, na.rm = TRUE)
  if (max_diff > 0.2) {
    cat("⚠ WARNING: Large differences detected. Results sensitive to study quality.\n")
  } else {
    cat("✓ Results appear robust to study quality.\n")
  }

  cat("\n")

  invisible(x)
}
