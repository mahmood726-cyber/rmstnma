#' League Tables for Network Meta-Analysis
#'
#' Generate publication-ready league tables showing all pairwise treatment
#' comparisons with effect estimates and confidence intervals.
#'
#' @name league_tables
NULL

#' Generate League Table
#'
#' Create a league table (also called matrix plot) showing all pairwise
#' comparisons between treatments. The diagonal shows treatment names and
#' SUCRA scores. Lower triangle shows effect estimates, upper triangle shows
#' confidence intervals.
#'
#' @param nma_result A netmeta object
#' @param sucra_result Optional sucra object for including SUCRA scores
#' @param digits Number of decimal places (default: 2)
#' @param format Format type: "effect_ci", "effect_only", "full"
#' @param sort_by Sort treatments by: "alphabet", "sucra", "reference_first"
#' @param include_pval Include p-values for significance
#' @return league_table object (matrix with class)
#' @export
#' @references
#' Salanti G, et al. (2011). Graphical methods and numerical summaries for
#' presenting results from multiple-treatment meta-analysis.
#' Journal of Clinical Epidemiology, 64(2):163-171.
#'
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' sucra <- calculate_sucra(nma)
#' league <- create_league_table(nma, sucra, format = "effect_ci")
#' print(league)
#' write_league_table(league, "league_table.csv")
#' }
create_league_table <- function(nma_result,
                               sucra_result = NULL,
                               digits = 2,
                               format = c("effect_ci", "effect_only", "full"),
                               sort_by = c("alphabet", "sucra", "reference_first"),
                               include_pval = FALSE) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  format <- match.arg(format)
  sort_by <- match.arg(sort_by)

  # Get treatments
  treatments <- rownames(nma_result$TE.random)
  n_treatments <- length(treatments)

  # Sort treatments
  if (sort_by == "sucra" && !is.null(sucra_result)) {
    treatments <- sucra_result$summary$Treatment
  } else if (sort_by == "reference_first") {
    ref <- nma_result$reference.group
    treatments <- c(ref, setdiff(treatments, ref))
  } else {
    treatments <- sort(treatments)
  }

  # Initialize matrix
  league <- matrix("", nrow = n_treatments, ncol = n_treatments)
  rownames(league) <- treatments
  colnames(league) <- treatments

  # Fill diagonal with treatment names and SUCRA
  for (i in 1:n_treatments) {
    trt <- treatments[i]
    if (!is.null(sucra_result)) {
      sucra_val <- sucra_result$sucra_scores[trt]
      if (!is.na(sucra_val)) {
        league[i, i] <- sprintf("%s\nSUCRA: %.1f%%", trt, sucra_val * 100)
      } else {
        league[i, i] <- trt
      }
    } else {
      league[i, i] <- trt
    }
  }

  # Fill upper and lower triangles
  for (i in 1:n_treatments) {
    for (j in 1:n_treatments) {
      if (i == j) next

      row_trt <- treatments[i]
      col_trt <- treatments[j]

      # Effect of row vs col (row - col)
      effect <- nma_result$TE.random[row_trt, col_trt]
      se <- nma_result$seTE.random[row_trt, col_trt]
      lower <- nma_result$lower.random[row_trt, col_trt]
      upper <- nma_result$upper.random[row_trt, col_trt]

      # P-value
      z_score <- effect / se
      p_val <- 2 * (1 - pnorm(abs(z_score)))

      # Format effect size
      effect_str <- format_effect(effect, lower, upper, p_val,
                                  nma_result$sm, digits, format, include_pval)

      if (i < j) {
        # Upper triangle: row treatment vs column treatment
        league[i, j] <- effect_str
      } else {
        # Lower triangle: column treatment vs row treatment (flip sign)
        effect_flipped <- -effect
        lower_flipped <- -upper
        upper_flipped <- -lower

        effect_str_flipped <- format_effect(
          effect_flipped, lower_flipped, upper_flipped, p_val,
          nma_result$sm, digits, format, include_pval
        )

        league[i, j] <- effect_str_flipped
      }
    }
  }

  attr(league, "format") <- format
  attr(league, "sm") <- nma_result$sm
  attr(league, "digits") <- digits
  attr(league, "include_pval") <- include_pval

  class(league) <- c("league_table", "matrix")
  league
}

#' Format Effect Size for League Table
#'
#' @param effect Effect estimate
#' @param lower Lower CI
#' @param upper Upper CI
#' @param p_val P-value
#' @param sm Summary measure
#' @param digits Decimal places
#' @param format Format type
#' @param include_pval Include p-value
#' @return Formatted string
#' @keywords internal
format_effect <- function(effect, lower, upper, p_val, sm, digits, format, include_pval) {

  # Check if significant
  is_sig <- p_val < 0.05

  if (format == "effect_only") {
    result <- sprintf(paste0("%.", digits, "f"), effect)
  } else if (format == "effect_ci") {
    result <- sprintf(
      paste0("%.", digits, "f (", "%.", digits, "f, ", "%.", digits, "f)"),
      effect, lower, upper
    )
  } else {  # full
    result <- sprintf(
      paste0("%.", digits, "f\n[", "%.", digits, "f, ", "%.", digits, "f]"),
      effect, lower, upper
    )
  }

  # Add significance marker
  if (is_sig && include_pval) {
    result <- paste0(result, "*")
  }

  # Add p-value if requested
  if (include_pval && format == "full") {
    result <- paste0(result, sprintf("\np=%.3f", p_val))
  }

  result
}

#' Create Simple Comparison Matrix
#'
#' Create a simplified matrix showing only point estimates, useful for
#' quick reference.
#'
#' @param nma_result A netmeta object
#' @param digits Number of decimal places
#' @return Matrix of effect estimates
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' matrix <- create_comparison_matrix(nma)
#' print(matrix)
#' }
create_comparison_matrix <- function(nma_result, digits = 2) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  treatments <- rownames(nma_result$TE.random)
  n_treatments <- length(treatments)

  # Get effect estimates
  effects <- round(nma_result$TE.random, digits)

  # For symmetric display, show row vs column
  comparison_matrix <- effects

  class(comparison_matrix) <- c("comparison_matrix", "matrix")
  attr(comparison_matrix, "sm") <- nma_result$sm

  comparison_matrix
}

#' Generate League Table with Probabilities
#'
#' Create an enhanced league table that includes probability that one
#' treatment is better than another for each comparison.
#'
#' @param nma_result A netmeta object
#' @param sucra_result A sucra object (required for probabilities)
#' @param digits Number of decimal places
#' @return league_table object with probabilities
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' sucra <- calculate_sucra(nma)
#' league_prob <- create_league_table_with_probabilities(nma, sucra)
#' print(league_prob)
#' }
create_league_table_with_probabilities <- function(nma_result, sucra_result, digits = 2) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  if (is.null(sucra_result)) {
    stop("sucra_result required for probability calculations")
  }

  # Base league table
  league <- create_league_table(nma_result, sucra_result, digits = digits, format = "effect_ci")

  # Calculate probabilities
  treatments <- rownames(league)
  n_treatments <- length(treatments)

  prob_matrix <- matrix("", nrow = n_treatments, ncol = n_treatments)
  rownames(prob_matrix) <- colnames(prob_matrix) <- treatments

  for (i in 1:n_treatments) {
    for (j in 1:n_treatments) {
      if (i == j) {
        prob_matrix[i, j] <- treatments[i]
        next
      }

      # Calculate probability that row treatment is better than column treatment
      row_trt <- treatments[i]
      col_trt <- treatments[j]

      effect <- nma_result$TE.random[row_trt, col_trt]
      se <- nma_result$seTE.random[row_trt, col_trt]

      # Probability that effect > 0 (row better than column)
      prob <- pnorm(0, mean = effect, sd = se, lower.tail = FALSE)

      prob_matrix[i, j] <- sprintf("%.1f%%", prob * 100)
    }
  }

  attr(league, "probabilities") <- prob_matrix

  league
}

#' Export League Table to File
#'
#' Write league table to CSV, Excel, or formatted text file
#'
#' @param league_table A league_table object
#' @param filename Output filename (extension determines format)
#' @param include_notes Include interpretation notes
#' @return Invisible path to file
#' @export
#' @examples
#' \dontrun{
#' league <- create_league_table(nma, sucra)
#' write_league_table(league, "results/league_table.csv")
#' write_league_table(league, "results/league_table.xlsx")
#' }
write_league_table <- function(league_table, filename, include_notes = TRUE) {

  if (!inherits(league_table, "league_table")) {
    stop("Input must be a league_table object")
  }

  ext <- tools::file_ext(filename)

  if (ext == "csv") {
    # CSV export
    utils::write.csv(league_table, filename, row.names = TRUE)

    if (include_notes) {
      notes_file <- sub("\\.csv$", "_notes.txt", filename)
      write_league_notes(league_table, notes_file)
    }

  } else if (ext == "xlsx") {
    # Excel export
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' required for Excel export. Use CSV instead or install writexl.")
    }

    df <- as.data.frame(league_table)
    writexl::write_xlsx(df, filename)

  } else if (ext == "txt") {
    # Text export with formatting
    sink(filename)
    print(league_table)
    if (include_notes) {
      cat("\n\n")
      cat("═══════════════════════════════════════════════\n")
      cat("  Notes\n")
      cat("═══════════════════════════════════════════════\n\n")
      cat(generate_league_notes(league_table))
    }
    sink()

  } else {
    stop("Unsupported file format. Use .csv, .xlsx, or .txt")
  }

  message(sprintf("League table exported to: %s", filename))
  invisible(filename)
}

#' Generate Interpretation Notes for League Table
#'
#' @param league_table A league_table object
#' @return Character string with notes
#' @keywords internal
generate_league_notes <- function(league_table) {

  sm <- attr(league_table, "sm")
  format <- attr(league_table, "format")
  include_pval <- attr(league_table, "include_pval")

  notes <- sprintf("Effect measure: %s\n", sm)
  notes <- paste0(notes, sprintf("Format: %s\n\n", format))

  notes <- paste0(notes, "How to read this league table:\n")
  notes <- paste0(notes, "• Diagonal: Treatment names (and SUCRA scores if available)\n")
  notes <- paste0(notes, "• Upper triangle: Row treatment compared to column treatment\n")
  notes <- paste0(notes, "• Lower triangle: Column treatment compared to row treatment\n")
  notes <- paste0(notes, "• Values show effect estimate with 95% confidence interval\n\n")

  if (sm %in% c("OR", "RR", "HR")) {
    notes <- paste0(notes, sprintf("Interpretation of %s:\n", sm))
    notes <- paste0(notes, "• Value > 1: Row treatment has higher risk/odds than column\n")
    notes <- paste0(notes, "• Value < 1: Row treatment has lower risk/odds than column\n")
    notes <- paste0(notes, "• Value = 1: No difference between treatments\n\n")
  } else if (sm %in% c("MD", "SMD")) {
    notes <- paste0(notes, sprintf("Interpretation of %s:\n", sm))
    notes <- paste0(notes, "• Value > 0: Row treatment has higher mean than column\n")
    notes <- paste0(notes, "• Value < 0: Row treatment has lower mean than column\n")
    notes <- paste0(notes, "• Value = 0: No difference between treatments\n\n")
  }

  if (include_pval) {
    notes <- paste0(notes, "* indicates p < 0.05 (statistically significant)\n")
  }

  notes
}

#' Write League Table Notes to File
#'
#' @param league_table A league_table object
#' @param filename Output filename
#' @keywords internal
write_league_notes <- function(league_table, filename) {
  notes <- generate_league_notes(league_table)
  writeLines(notes, filename)
  message(sprintf("Notes written to: %s", filename))
}

#' Plot League Table as Heatmap
#'
#' Create a visual heatmap representation of the league table showing
#' effect sizes by color intensity.
#'
#' @param league_table A league_table object or netmeta object
#' @param nma_result If league_table is NULL, use this netmeta object
#' @param color_scale Color palette: "diverging", "sequential", "viridis"
#' @param show_values Show numerical values on heatmap
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' plot_league_heatmap(nma_result = nma)
#' }
plot_league_heatmap <- function(league_table = NULL,
                               nma_result = NULL,
                               color_scale = c("diverging", "sequential", "viridis"),
                               show_values = TRUE) {

  color_scale <- match.arg(color_scale)

  if (is.null(league_table) && is.null(nma_result)) {
    stop("Provide either league_table or nma_result")
  }

  if (!is.null(nma_result) && inherits(nma_result, "netmeta")) {
    # Use netmeta object directly
    treatments <- rownames(nma_result$TE.random)
    effects_matrix <- nma_result$TE.random
  } else if (!is.null(league_table)) {
    stop("Heatmap from league_table object not yet implemented. Use nma_result instead.")
  }

  # Convert to long format for ggplot
  plot_data <- expand.grid(
    Treatment1 = treatments,
    Treatment2 = treatments,
    stringsAsFactors = FALSE
  )

  plot_data$Effect <- as.vector(effects_matrix)

  # Remove diagonal
  plot_data <- plot_data[plot_data$Treatment1 != plot_data$Treatment2, ]

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment2, y = Treatment1, fill = Effect)) +
    ggplot2::geom_tile(color = "white", size = 0.5)

  # Color scale
  if (color_scale == "diverging") {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, name = "Effect"
    )
  } else if (color_scale == "sequential") {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#f7f7f7", high = "#2c7fb8", name = "Effect"
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(name = "Effect")
  }

  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Effect)),
      size = 3, color = "black"
    )
  }

  p <- p +
    ggplot2::labs(
      title = "Treatment Comparison Heatmap",
      subtitle = "Effect of row treatment vs. column treatment",
      x = "Treatment (Column)",
      y = "Treatment (Row)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  p
}

#' Print League Table
#'
#' @param x league_table object
#' @param ... Additional arguments
#' @export
print.league_table <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════\n")
  cat("  League Table: All Pairwise Comparisons\n")
  cat("═══════════════════════════════════════════════\n\n")

  cat(sprintf("Effect measure: %s\n", attr(x, "sm")))
  cat(sprintf("Format: %s\n", attr(x, "format")))
  cat(sprintf("Decimal places: %d\n\n", attr(x, "digits")))

  # Print matrix
  print(noquote(x))

  cat("\n")
  cat("How to read: Upper triangle shows row vs column.\n")
  cat("             Lower triangle shows column vs row.\n")
  cat("             Diagonal shows treatment names.\n")

  if (attr(x, "include_pval")) {
    cat("\n* indicates p < 0.05\n")
  }

  # Show probabilities if available
  if (!is.null(attr(x, "probabilities"))) {
    cat("\n")
    cat("═══════════════════════════════════════════════\n")
    cat("  Probability Matrix\n")
    cat("═══════════════════════════════════════════════\n")
    cat("Probability that row treatment is better than column treatment:\n\n")
    print(noquote(attr(x, "probabilities")))
  }

  cat("\n")

  invisible(x)
}

#' Print Comparison Matrix
#'
#' @param x comparison_matrix object
#' @param ... Additional arguments
#' @export
print.comparison_matrix <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════\n")
  cat("  Comparison Matrix (Point Estimates)\n")
  cat("═══════════════════════════════════════════════\n\n")

  cat(sprintf("Effect measure: %s\n\n", attr(x, "sm")))

  print(x)

  cat("\nNote: Values show row treatment vs. column treatment\n")

  invisible(x)
}
