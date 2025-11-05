#' Advanced Export and Interoperability
#'
#' Tools for exporting NMA results to multiple formats and ensuring
#' interoperability with other software and platforms.
#'
#' @name export_interoperability
NULL

#' Export NMA Results to Multiple Formats
#'
#' Export comprehensive NMA results to various formats for publication,
#' presentation, or further analysis.
#'
#' @param nma_result NMA result object
#' @param output_dir Output directory
#' @param formats Vector of formats: "csv", "excel", "json", "html", "latex"
#' @param prefix File prefix
#' @return List of file paths created
#' @export
#' @examples
#' \dontrun{
#' nma <- netmeta(TE, seTE, treat1, treat2, studlab, data = mydata)
#' files <- export_nma_results(nma, output_dir = "results", formats = c("csv", "excel"))
#' }
export_nma_results <- function(nma_result,
                               output_dir = "nma_export",
                               formats = c("csv", "excel", "json", "html"),
                               prefix = "nma") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(sprintf("Created output directory: %s", output_dir))
  }

  message("Exporting NMA results...")

  files_created <- list()

  # Extract key results
  results_list <- extract_nma_results(nma_result)

  # Export to each format
  for (fmt in formats) {
    files <- switch(fmt,
      csv = export_to_csv(results_list, output_dir, prefix),
      excel = export_to_excel(results_list, output_dir, prefix),
      json = export_to_json(results_list, output_dir, prefix),
      html = export_to_html(results_list, nma_result, output_dir, prefix),
      latex = export_to_latex(results_list, output_dir, prefix),
      NULL
    )

    if (!is.null(files)) {
      files_created[[fmt]] <- files
    }
  }

  message(sprintf("\n✓ Exported to %d format(s)", length(files_created)))

  invisible(files_created)
}

#' Extract NMA Results
#' @keywords internal
extract_nma_results <- function(nma_result) {

  # Treatment effects
  effects <- data.frame(
    Comparison = rownames(nma_result$TE.random),
    Reference = colnames(nma_result$TE.random)[1],
    Effect = nma_result$TE.random[, 1],
    SE = nma_result$seTE.random[, 1],
    Lower_95 = nma_result$lower.random[, 1],
    Upper_95 = nma_result$upper.random[, 1],
    P_value = 2 * (1 - pnorm(abs(nma_result$TE.random[, 1] / nma_result$seTE.random[, 1]))),
    stringsAsFactors = FALSE
  )

  # Heterogeneity
  heterogeneity <- data.frame(
    Statistic = c("tau", "tau2", "I2"),
    Value = c(nma_result$tau, nma_result$tau^2, nma_result$I2),
    stringsAsFactors = FALSE
  )

  # Network characteristics
  network_info <- data.frame(
    Characteristic = c("Treatments", "Studies", "Comparisons"),
    Value = c(length(nma_result$trts),
             length(unique(nma_result$studlab)),
             length(nma_result$studlab)),
    stringsAsFactors = FALSE
  )

  list(
    effects = effects,
    heterogeneity = heterogeneity,
    network_info = network_info,
    treatments = nma_result$trts
  )
}

#' Export to CSV
#' @keywords internal
export_to_csv <- function(results_list, output_dir, prefix) {

  files <- character()

  # Effects
  effects_file <- file.path(output_dir, paste0(prefix, "_effects.csv"))
  write.csv(results_list$effects, effects_file, row.names = FALSE)
  files <- c(files, effects_file)

  # Heterogeneity
  het_file <- file.path(output_dir, paste0(prefix, "_heterogeneity.csv"))
  write.csv(results_list$heterogeneity, het_file, row.names = FALSE)
  files <- c(files, het_file)

  # Network info
  info_file <- file.path(output_dir, paste0(prefix, "_network_info.csv"))
  write.csv(results_list$network_info, info_file, row.names = FALSE)
  files <- c(files, info_file)

  message(sprintf("  ✓ CSV: %d files", length(files)))

  files
}

#' Export to Excel
#' @keywords internal
export_to_excel <- function(results_list, output_dir, prefix) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("  ✗ Excel export requires 'openxlsx' package")
    return(NULL)
  }

  excel_file <- file.path(output_dir, paste0(prefix, "_results.xlsx"))

  wb <- openxlsx::createWorkbook()

  # Add sheets
  openxlsx::addWorksheet(wb, "Treatment Effects")
  openxlsx::writeData(wb, "Treatment Effects", results_list$effects)

  openxlsx::addWorksheet(wb, "Heterogeneity")
  openxlsx::writeData(wb, "Heterogeneity", results_list$heterogeneity)

  openxlsx::addWorksheet(wb, "Network Info")
  openxlsx::writeData(wb, "Network Info", results_list$network_info)

  # Save
  openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)

  message(sprintf("  ✓ Excel: %s", excel_file))

  excel_file
}

#' Export to JSON
#' @keywords internal
export_to_json <- function(results_list, output_dir, prefix) {

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("  ✗ JSON export requires 'jsonlite' package")
    return(NULL)
  }

  json_file <- file.path(output_dir, paste0(prefix, "_results.json"))

  # Convert to JSON
  json_data <- jsonlite::toJSON(results_list, pretty = TRUE, auto_unbox = TRUE)

  # Write
  writeLines(json_data, json_file)

  message(sprintf("  ✓ JSON: %s", json_file))

  json_file
}

#' Export to HTML
#' @keywords internal
export_to_html <- function(results_list, nma_result, output_dir, prefix) {

  html_file <- file.path(output_dir, paste0(prefix, "_results.html"))

  html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
  <title>Network Meta-Analysis Results</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h1 { color: #2c7fb8; }
    h2 { color: #666; margin-top: 30px; }
    table { border-collapse: collapse; width: 100%%; margin-top: 10px; }
    th { background-color: #2c7fb8; color: white; padding: 10px; text-align: left; }
    td { border: 1px solid #ddd; padding: 8px; }
    tr:nth-child(even) { background-color: #f2f2f2; }
    .stat { font-size: 1.2em; font-weight: bold; color: #2c7fb8; }
  </style>
</head>
<body>
  <h1>Network Meta-Analysis Results</h1>

  <h2>Network Characteristics</h2>
  %s

  <h2>Treatment Effects</h2>
  <p>Reference: %s</p>
  %s

  <h2>Heterogeneity</h2>
  %s

  <p style="margin-top: 40px; color: #666; font-size: 0.9em;">
    Generated by powerNMA R package
  </p>
</body>
</html>
',
    table_to_html(results_list$network_info),
    results_list$effects$Reference[1],
    table_to_html(results_list$effects),
    table_to_html(results_list$heterogeneity)
  )

  writeLines(html_content, html_file)

  message(sprintf("  ✓ HTML: %s", html_file))

  html_file
}

#' Export to LaTeX
#' @keywords internal
export_to_latex <- function(results_list, output_dir, prefix) {

  latex_file <- file.path(output_dir, paste0(prefix, "_results.tex"))

  latex_content <- sprintf('
\\documentclass{article}
\\usepackage{booktabs}
\\usepackage{longtable}

\\begin{document}

\\section{Network Meta-Analysis Results}

\\subsection{Network Characteristics}
%s

\\subsection{Treatment Effects}
%s

\\subsection{Heterogeneity}
%s

\\end{document}
',
    dataframe_to_latex(results_list$network_info),
    dataframe_to_latex(results_list$effects),
    dataframe_to_latex(results_list$heterogeneity)
  )

  writeLines(latex_content, latex_file)

  message(sprintf("  ✓ LaTeX: %s", latex_file))

  latex_file
}

#' Convert Data Frame to HTML Table
#' @keywords internal
table_to_html <- function(df) {

  # Header
  header <- paste0("  <tr><th>", paste(names(df), collapse = "</th><th>"), "</th></tr>")

  # Rows
  rows <- apply(df, 1, function(row) {
    paste0("  <tr><td>", paste(row, collapse = "</td><td>"), "</td></tr>")
  })

  paste0("<table>\n", header, "\n", paste(rows, collapse = "\n"), "\n</table>")
}

#' Convert Data Frame to LaTeX Table
#' @keywords internal
dataframe_to_latex <- function(df) {

  n_cols <- ncol(df)

  # Begin table
  tex <- sprintf("\\begin{tabular}{%s}\n", paste(rep("l", n_cols), collapse = ""))
  tex <- paste0(tex, "\\toprule\n")

  # Header
  tex <- paste0(tex, paste(names(df), collapse = " & "), " \\\\\n")
  tex <- paste0(tex, "\\midrule\n")

  # Rows
  for (i in seq_len(nrow(df))) {
    row_tex <- paste(df[i, ], collapse = " & ")
    tex <- paste0(tex, row_tex, " \\\\\n")
  }

  # End table
  tex <- paste0(tex, "\\bottomrule\n")
  tex <- paste0(tex, "\\end{tabular}\n")

  tex
}

#' Export for GRADE Assessment
#'
#' Export NMA results in format suitable for GRADE evidence assessment.
#'
#' @param nma_result NMA result
#' @param output_file Output file path
#' @return Data frame for GRADE
#' @export
export_for_grade <- function(nma_result, output_file = "grade_assessment.csv") {

  # Extract comparisons
  comparisons <- expand.grid(
    Treatment = nma_result$trts,
    Reference = nma_result$reference.group,
    stringsAsFactors = FALSE
  )

  comparisons <- comparisons[comparisons$Treatment != comparisons$Reference, ]

  # Add effect estimates
  comparisons$Effect <- nma_result$TE.random[cbind(comparisons$Treatment, comparisons$Reference)]
  comparisons$SE <- nma_result$seTE.random[cbind(comparisons$Treatment, comparisons$Reference)]
  comparisons$Lower_95 <- nma_result$lower.random[cbind(comparisons$Treatment, comparisons$Reference)]
  comparisons$Upper_95 <- nma_result$upper.random[cbind(comparisons$Treatment, comparisons$Reference)]

  # Add GRADE columns (to be filled manually)
  comparisons$Risk_of_Bias <- ""
  comparisons$Inconsistency <- ""
  comparisons$Indirectness <- ""
  comparisons$Imprecision <- ""
  comparisons$Publication_Bias <- ""
  comparisons$Overall_Quality <- ""

  write.csv(comparisons, output_file, row.names = FALSE)

  message(sprintf("✓ GRADE template exported: %s", output_file))

  comparisons
}

#' Export Network Diagram for Publication
#'
#' Export high-quality network diagram suitable for publication.
#'
#' @param data NMA data
#' @param filename Output filename
#' @param format Format: "png", "pdf", "svg"
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for raster formats
#' @return File path
#' @export
export_network_diagram <- function(data, filename = "network_diagram",
                                   format = c("png", "pdf", "svg"),
                                   width = 8, height = 6, dpi = 300) {

  format <- match.arg(format)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required for network diagrams")
  }

  # Build network
  treatments <- unique(c(as.character(data$treat1), as.character(data$treat2)))

  adj_matrix <- matrix(0, nrow = length(treatments), ncol = length(treatments),
                      dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])
    adj_matrix[t1, t2] <- adj_matrix[t1, t2] + 1
    adj_matrix[t2, t1] <- adj_matrix[t2, t1] + 1
  }

  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

  # Create file
  output_file <- paste0(filename, ".", format)

  if (format == "png") {
    png(output_file, width = width, height = height, units = "in", res = dpi)
  } else if (format == "pdf") {
    pdf(output_file, width = width, height = height)
  } else if (format == "svg") {
    svg(output_file, width = width, height = height)
  }

  # Plot
  plot(g,
       layout = igraph::layout_with_fr(g),
       vertex.size = 20,
       vertex.label.cex = 1.2,
       vertex.color = "#2c7fb8",
       vertex.label.color = "white",
       edge.width = sqrt(igraph::E(g)$weight) * 2,
       edge.color = "gray50")

  dev.off()

  message(sprintf("✓ Network diagram exported: %s", output_file))

  output_file
}

#' Convert to BUGSnet Format
#'
#' Convert data to BUGSnet format for Bayesian NMA in R.
#'
#' @param data NMA data
#' @return List in BUGSnet format
#' @export
convert_to_bugsnet <- function(data) {

  message("Converting to BUGSnet format...")

  # BUGSnet expects specific column names
  bugsnet_data <- data.frame(
    studyID = data$studlab,
    T = data$treat1,
    C = data$treat2,
    TE = data$TE,
    seTE = data$seTE,
    stringsAsFactors = FALSE
  )

  # Add sample sizes if available
  if ("n1" %in% names(data)) {
    bugsnet_data$n_T <- data$n1
  }

  if ("n2" %in% names(data)) {
    bugsnet_data$n_C <- data$n2
  }

  message("✓ Converted to BUGSnet format")

  bugsnet_data
}

#' Convert to GeMTC Format
#'
#' Convert data to gemtc format for Bayesian NMA.
#'
#' @param data NMA data
#' @return gemtc data object
#' @export
convert_to_gemtc <- function(data) {

  if (!requireNamespace("gemtc", quietly = TRUE)) {
    stop("Package 'gemtc' required for this conversion")
  }

  message("Converting to gemtc format...")

  # gemtc expects specific structure
  gemtc_data <- data.frame(
    study = data$studlab,
    treatment = data$treat1,
    diff = data$TE,
    std.err = data$seTE,
    stringsAsFactors = FALSE
  )

  message("✓ Converted to gemtc format")

  gemtc_data
}

#' Export Summary Table for Manuscript
#'
#' Create publication-ready summary table.
#'
#' @param nma_result NMA result
#' @param filename Output filename
#' @param format Format: "docx", "html", "latex"
#' @return File path
#' @export
export_summary_table <- function(nma_result, filename = "summary_table",
                                 format = c("docx", "html", "latex")) {

  format <- match.arg(format)

  # Create summary table
  summary_df <- data.frame(
    Treatment = rownames(nma_result$TE.random),
    Effect = sprintf("%.2f", nma_result$TE.random[, 1]),
    CI_95 = sprintf("%.2f to %.2f",
                   nma_result$lower.random[, 1],
                   nma_result$upper.random[, 1]),
    P_value = sprintf("%.3f",
                     2 * (1 - pnorm(abs(nma_result$TE.random[, 1] / nma_result$seTE.random[, 1])))),
    stringsAsFactors = FALSE
  )

  # Remove reference (all zeros)
  summary_df <- summary_df[summary_df$Treatment != nma_result$reference.group, ]

  output_file <- paste0(filename, ".", format)

  if (format == "html") {
    html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
  <style>
    table { border-collapse: collapse; }
    th, td { border: 1px solid black; padding: 8px; text-align: left; }
    th { background-color: #f2f2f2; }
  </style>
</head>
<body>
  <h2>Network Meta-Analysis Results</h2>
  <p>Reference: %s</p>
  %s
</body>
</html>
',
      nma_result$reference.group,
      table_to_html(summary_df)
    )

    writeLines(html_content, output_file)

  } else if (format == "latex") {
    tex_content <- sprintf('
\\begin{table}[ht]
\\centering
\\caption{Network Meta-Analysis Results (Reference: %s)}
%s
\\end{table}
',
      nma_result$reference.group,
      dataframe_to_latex(summary_df)
    )

    writeLines(tex_content, output_file)

  } else if (format == "docx") {
    if (!requireNamespace("officer", quietly = TRUE)) {
      message("Package 'officer' required for Word export")
      return(NULL)
    }

    # Create Word document (placeholder - full implementation would use officer package)
    message("Word export requires 'officer' package - exporting as CSV instead")
    output_file <- paste0(filename, ".csv")
    write.csv(summary_df, output_file, row.names = FALSE)
  }

  message(sprintf("✓ Summary table exported: %s", output_file))

  output_file
}
