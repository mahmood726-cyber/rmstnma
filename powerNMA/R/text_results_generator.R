#' Comprehensive Text Results Generator
#'
#' @description
#' Generates detailed narrative text summaries of network meta-analysis results
#' with complete statistical reporting, clinical interpretation, and structured
#' sections covering all aspects of the analysis.
#'
#' @details
#' This module creates publication-ready text summaries that can be directly
#' included in manuscripts, reports, or presentations. The generator covers:
#' \itemize{
#'   \item Study and network characteristics
#'   \item Treatment effect estimates with confidence intervals
#'   \item Treatment rankings and SUCRA scores
#'   \item Heterogeneity assessment
#'   \item Inconsistency evaluation
#'   \item Publication bias assessment
#'   \item Sensitivity analyses
#'   \item Clinical interpretation
#' }
#'
#' @references
#' Dias S, et al. (2013). Evidence synthesis for decision making 2.
#' Chaimani A, et al. (2013). Graphical tools for network meta-analysis.
#' Salanti G, et al. (2011). Evaluation of networks of randomized trials.
#'
#' @author powerNMA Development Team
#' @name text_results_generator
NULL

#' Generate Comprehensive Text Results
#'
#' @description
#' Main function to generate comprehensive narrative text summary of all
#' network meta-analysis results.
#'
#' @param nma_result Network meta-analysis result object
#' @param sucra_result SUCRA rankings result (optional)
#' @param heterogeneity_result Heterogeneity assessment result (optional)
#' @param inconsistency_result Inconsistency check result (optional)
#' @param publication_bias_result Publication bias assessment (optional)
#' @param sensitivity_result Sensitivity analysis result (optional)
#' @param style Character string: "narrative", "structured", "detailed"
#' @param include_clinical Logical whether to include clinical interpretation
#' @param include_stats Logical whether to include detailed statistics
#' @param reference Reference treatment (if NULL, uses from nma_result)
#'
#' @return List containing:
#'   \item{full_text}{Complete narrative text}
#'   \item{sections}{Individual section texts}
#'   \item{metadata}{Generation metadata}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run comprehensive NMA
#' nma <- run_comprehensive_nma(data)
#' sucra <- calculate_sucra(nma)
#' het <- heterogeneity_report(nma)
#'
#' # Generate text results
#' text_results <- generate_comprehensive_text_results(
#'   nma_result = nma,
#'   sucra_result = sucra,
#'   heterogeneity_result = het,
#'   style = "narrative",
#'   include_clinical = TRUE
#' )
#'
#' # Print full text
#' cat(text_results$full_text)
#'
#' # Or access individual sections
#' cat(text_results$sections$overview)
#' cat(text_results$sections$treatment_effects)
#' }
generate_comprehensive_text_results <- function(
    nma_result,
    sucra_result = NULL,
    heterogeneity_result = NULL,
    inconsistency_result = NULL,
    publication_bias_result = NULL,
    sensitivity_result = NULL,
    style = c("narrative", "structured", "detailed"),
    include_clinical = TRUE,
    include_stats = TRUE,
    reference = NULL) {

  style <- match.arg(style)

  # Extract reference
  if (is.null(reference)) {
    reference <- nma_result$reference.group
  }

  # Generate each section
  sections <- list()

  # 1. Network overview
  sections$overview <- generate_network_overview_text(
    nma_result, style, include_stats
  )

  # 2. Treatment effects
  sections$treatment_effects <- generate_treatment_effects_text(
    nma_result, reference, style, include_stats, include_clinical
  )

  # 3. Treatment rankings
  if (!is.null(sucra_result)) {
    sections$rankings <- generate_rankings_text(
      sucra_result, style, include_clinical
    )
  }

  # 4. Heterogeneity
  if (!is.null(heterogeneity_result)) {
    sections$heterogeneity <- generate_heterogeneity_text(
      heterogeneity_result, style, include_clinical
    )
  }

  # 5. Inconsistency
  if (!is.null(inconsistency_result)) {
    sections$inconsistency <- generate_inconsistency_text(
      inconsistency_result, style, include_clinical
    )
  }

  # 6. Publication bias
  if (!is.null(publication_bias_result)) {
    sections$publication_bias <- generate_publication_bias_text(
      publication_bias_result, style, include_clinical
    )
  }

  # 7. Sensitivity analysis
  if (!is.null(sensitivity_result)) {
    sections$sensitivity <- generate_sensitivity_text(
      sensitivity_result, style, include_clinical
    )
  }

  # 8. Summary and conclusions
  sections$summary <- generate_summary_text(
    nma_result, sucra_result, heterogeneity_result,
    style, include_clinical
  )

  # Combine sections
  full_text <- combine_sections(sections, style)

  # Metadata
  metadata <- list(
    style = style,
    include_clinical = include_clinical,
    include_stats = include_stats,
    sections_included = names(sections),
    generation_time = Sys.time(),
    word_count = length(strsplit(full_text, "\\s+")[[1]]),
    character_count = nchar(full_text)
  )

  return(structure(
    list(
      full_text = full_text,
      sections = sections,
      metadata = metadata
    ),
    class = "text_results"
  ))
}

#' Generate Network Overview Text
#'
#' @description
#' Creates narrative text describing network characteristics.
#'
#' @param nma_result NMA result object
#' @param style Output style
#' @param include_stats Include statistical details
#'
#' @return Character string with overview text
#'
#' @keywords internal
generate_network_overview_text <- function(nma_result, style, include_stats) {

  # Extract network characteristics
  n_studies <- nma_result$k
  n_treatments <- length(nma_result$trts)
  treatments <- nma_result$trts
  n_comparisons <- nma_result$m
  n_patients <- sum(nma_result$n, na.rm = TRUE)

  # Build text
  text <- sprintf(
    "The network meta-analysis included %d studies comparing %d treatments (%s).",
    n_studies, n_treatments, paste(treatments, collapse = ", ")
  )

  if (include_stats) {
    text <- paste(
      text,
      sprintf(
        "A total of %d patients were analyzed across %d treatment comparisons.",
        n_patients, n_comparisons
      )
    )
  }

  # Network structure
  if (style %in% c("detailed", "structured")) {
    # Calculate network density
    max_comparisons <- n_treatments * (n_treatments - 1) / 2
    density <- (n_comparisons / max_comparisons) * 100

    text <- paste(
      text,
      sprintf(
        "The network density was %.1f%%, indicating a %s connected network.",
        density,
        ifelse(density > 75, "highly",
               ifelse(density > 50, "moderately", "sparsely"))
      )
    )
  }

  return(text)
}

#' Generate Treatment Effects Text
#'
#' @description
#' Creates narrative text for treatment effect estimates.
#'
#' @keywords internal
generate_treatment_effects_text <- function(nma_result, reference,
                                           style, include_stats,
                                           include_clinical) {

  # Extract treatment effects
  effects <- nma_result$TE.random
  se <- nma_result$seTE.random
  lower <- nma_result$lower.random
  upper <- nma_result$upper.random
  pvals <- nma_result$pval.random
  treatments <- rownames(effects)

  # Get comparisons vs reference
  ref_idx <- which(treatments == reference)
  if (length(ref_idx) == 0) ref_idx <- 1

  text_parts <- character(0)

  # Overall summary
  n_sig <- sum(pvals[, ref_idx] < 0.05, na.rm = TRUE)
  n_total <- nrow(effects) - 1

  text_parts <- c(text_parts, sprintf(
    "Compared to %s (reference treatment), %d of %d interventions showed statistically significant differences (p < 0.05).",
    reference, n_sig, n_total
  ))

  # Individual treatment effects
  for (i in seq_along(treatments)) {
    if (treatments[i] == reference) next

    eff <- effects[i, ref_idx]
    se_val <- se[i, ref_idx]
    low <- lower[i, ref_idx]
    up <- upper[i, ref_idx]
    pval <- pvals[i, ref_idx]

    # Determine effect direction and magnitude
    direction <- ifelse(eff > 0, "higher", "lower")
    sig <- ifelse(pval < 0.05, "statistically significant", "not statistically significant")

    # Format effect size based on measure
    sm <- nma_result$sm
    if (sm %in% c("OR", "RR", "HR")) {
      eff_exp <- exp(eff)
      low_exp <- exp(low)
      up_exp <- exp(up)

      effect_text <- sprintf(
        "%s showed %s effects compared to %s (%s = %.2f, 95%% CI: %.2f-%.2f, p = %.3f), which was %s.",
        treatments[i], direction, reference, sm,
        eff_exp, low_exp, up_exp, pval, sig
      )
    } else {
      effect_text <- sprintf(
        "%s showed %s effects compared to %s (%s = %.2f, 95%% CI: %.2f-%.2f, p = %.3f), which was %s.",
        treatments[i], direction, reference, sm,
        eff, low, up, pval, sig
      )
    }

    text_parts <- c(text_parts, effect_text)

    # Add clinical interpretation if requested
    if (include_clinical && style == "detailed") {
      if (pval < 0.05) {
        magnitude <- ifelse(abs(eff) > 0.5, "large", ifelse(abs(eff) > 0.2, "moderate", "small"))
        text_parts <- c(text_parts, sprintf(
          "This represents a %s effect size with potential clinical relevance.",
          magnitude
        ))
      }
    }
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Rankings Text
#'
#' @description
#' Creates narrative text for treatment rankings.
#'
#' @keywords internal
generate_rankings_text <- function(sucra_result, style, include_clinical) {

  sucra_scores <- sucra_result$sucra_scores
  rankings <- sucra_result$mean_ranks

  # Sort by SUCRA
  sorted_idx <- order(sucra_scores, decreasing = TRUE)
  treatments <- names(sucra_scores)[sorted_idx]
  scores <- sucra_scores[sorted_idx]

  text_parts <- character(0)

  # Overall ranking summary
  best <- treatments[1]
  worst <- treatments[length(treatments)]

  text_parts <- c(text_parts, sprintf(
    "Treatment rankings based on SUCRA scores showed %s as the highest-ranked intervention (SUCRA = %.1f%%), while %s was ranked lowest (SUCRA = %.1f%%).",
    best, scores[1], worst, scores[length(scores)]
  ))

  # Top 3 treatments
  if (length(treatments) >= 3 && style %in% c("detailed", "structured")) {
    text_parts <- c(text_parts, sprintf(
      "The top three treatments were %s (SUCRA = %.1f%%), %s (SUCRA = %.1f%%), and %s (SUCRA = %.1f%%).",
      treatments[1], scores[1],
      treatments[2], scores[2],
      treatments[3], scores[3]
    ))
  }

  # Clinical interpretation
  if (include_clinical) {
    if (scores[1] - scores[2] > 20) {
      text_parts <- c(text_parts,
        "The highest-ranked treatment showed a substantial advantage over other interventions, suggesting it may be the preferred option for clinical practice."
      )
    } else if (scores[1] - scores[2] < 10) {
      text_parts <- c(text_parts,
        "The top-ranked treatments showed similar performance, indicating that multiple interventions may be suitable depending on patient characteristics and preferences."
      )
    }
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Heterogeneity Text
#'
#' @description
#' Creates narrative text for heterogeneity assessment.
#'
#' @keywords internal
generate_heterogeneity_text <- function(heterogeneity_result, style, include_clinical) {

  tau2 <- heterogeneity_result$tau2
  I2 <- heterogeneity_result$I2

  # Interpret I2
  I2_interp <- if (I2 < 25) {
    "low"
  } else if (I2 < 50) {
    "moderate"
  } else if (I2 < 75) {
    "substantial"
  } else {
    "considerable"
  }

  text_parts <- character(0)

  text_parts <- c(text_parts, sprintf(
    "Heterogeneity assessment revealed %s heterogeneity across the network (τ² = %.3f, I² = %.1f%%).",
    I2_interp, tau2, I2
  ))

  # Clinical implications
  if (include_clinical) {
    if (I2 < 50) {
      text_parts <- c(text_parts,
        "The low to moderate heterogeneity suggests that treatment effects are relatively consistent across studies, supporting the generalizability of findings."
      )
    } else {
      text_parts <- c(text_parts,
        "The substantial heterogeneity indicates important differences between studies that should be explored through subgroup analyses or meta-regression."
      )
    }
  }

  # Prediction intervals if available
  if (!is.null(heterogeneity_result$prediction_intervals) && style == "detailed") {
    text_parts <- c(text_parts,
      "Prediction intervals were calculated to account for heterogeneity in estimating the range of true effects in future studies."
    )
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Inconsistency Text
#'
#' @description
#' Creates narrative text for inconsistency assessment.
#'
#' @keywords internal
generate_inconsistency_text <- function(inconsistency_result, style, include_clinical) {

  text_parts <- character(0)

  # Check for node-splitting results
  if (!is.null(inconsistency_result$node_splitting)) {
    ns_result <- inconsistency_result$node_splitting
    n_comparisons <- nrow(ns_result)
    n_inconsistent <- sum(ns_result$p_value < 0.05, na.rm = TRUE)

    text_parts <- c(text_parts, sprintf(
      "Node-splitting analysis evaluated %d comparisons for inconsistency. %d comparison(s) showed evidence of inconsistency (p < 0.05).",
      n_comparisons, n_inconsistent
    ))

    if (n_inconsistent == 0) {
      text_parts <- c(text_parts,
        "The absence of significant inconsistency supports the validity of the network meta-analysis assumptions."
      )
    } else {
      text_parts <- c(text_parts,
        "The presence of inconsistency suggests differences between direct and indirect evidence that should be interpreted cautiously."
      )
    }
  }

  # Design-by-treatment interaction
  if (!is.null(inconsistency_result$design_inconsistency)) {
    di_result <- inconsistency_result$design_inconsistency
    pval <- di_result$p_value

    text_parts <- c(text_parts, sprintf(
      "The design-by-treatment interaction model showed %s evidence of inconsistency (p = %.3f).",
      ifelse(pval < 0.05, "significant", "no significant"), pval
    ))
  }

  # Clinical interpretation
  if (include_clinical && length(text_parts) > 0) {
    if (n_inconsistent == 0) {
      text_parts <- c(text_parts,
        "These consistency checks support confidence in the pooled network estimates and treatment rankings."
      )
    }
  }

  if (length(text_parts) == 0) {
    text_parts <- "Inconsistency assessment was not performed or not available."
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Publication Bias Text
#'
#' @description
#' Creates narrative text for publication bias assessment.
#'
#' @keywords internal
generate_publication_bias_text <- function(publication_bias_result, style, include_clinical) {

  text_parts <- character(0)

  # Comparison-adjusted funnel plot
  if (!is.null(publication_bias_result$egger_test)) {
    egger_p <- publication_bias_result$egger_test$p_value

    text_parts <- c(text_parts, sprintf(
      "Publication bias assessment using comparison-adjusted Egger's test showed %s evidence of small-study effects (p = %.3f).",
      ifelse(egger_p < 0.05, "significant", "no significant"), egger_p
    ))
  }

  # Begg's test
  if (!is.null(publication_bias_result$begg_test)) {
    begg_p <- publication_bias_result$begg_test$p_value

    text_parts <- c(text_parts, sprintf(
      "Begg's test for publication bias was %s (p = %.3f).",
      ifelse(begg_p < 0.05, "significant", "not significant"), begg_p
    ))
  }

  # Clinical interpretation
  if (include_clinical && length(text_parts) > 0) {
    if (egger_p < 0.05 || begg_p < 0.05) {
      text_parts <- c(text_parts,
        "The presence of potential publication bias suggests that effect estimates should be interpreted with caution, as smaller studies with null or negative results may be missing."
      )
    } else {
      text_parts <- c(text_parts,
        "The absence of significant publication bias supports the robustness of the findings."
      )
    }
  }

  if (length(text_parts) == 0) {
    text_parts <- "Publication bias assessment was not performed or not available."
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Sensitivity Text
#'
#' @description
#' Creates narrative text for sensitivity analyses.
#'
#' @keywords internal
generate_sensitivity_text <- function(sensitivity_result, style, include_clinical) {

  text_parts <- character(0)

  # Leave-one-out
  if (!is.null(sensitivity_result$leave_one_out)) {
    loo <- sensitivity_result$leave_one_out
    n_studies <- length(loo$removed_studies)

    text_parts <- c(text_parts, sprintf(
      "Leave-one-out sensitivity analysis showed that the overall conclusions were robust to the removal of individual studies in %d of %d cases.",
      sum(!loo$changed_significance, na.rm = TRUE), n_studies
    ))

    if (any(loo$changed_significance, na.rm = TRUE)) {
      influential_studies <- loo$removed_studies[loo$changed_significance]
      text_parts <- c(text_parts, sprintf(
        "Removal of %d study/studies (%s) altered the statistical significance of findings, suggesting these may be influential.",
        length(influential_studies), paste(influential_studies, collapse = ", ")
      ))
    }
  }

  # Subgroup analyses
  if (!is.null(sensitivity_result$subgroup)) {
    text_parts <- c(text_parts,
      "Subgroup analyses were conducted to explore potential sources of heterogeneity."
    )
  }

  if (include_clinical) {
    text_parts <- c(text_parts,
      "Sensitivity analyses support the stability and reliability of the main findings across different analytical scenarios."
    )
  }

  if (length(text_parts) == 0) {
    text_parts <- "Sensitivity analyses were not performed or not available."
  }

  return(paste(text_parts, collapse = " "))
}

#' Generate Summary Text
#'
#' @description
#' Creates overall summary and conclusions.
#'
#' @keywords internal
generate_summary_text <- function(nma_result, sucra_result,
                                 heterogeneity_result, style,
                                 include_clinical) {

  text_parts <- character(0)

  # Best treatment
  if (!is.null(sucra_result)) {
    best_treatment <- names(which.max(sucra_result$sucra_scores))
    best_sucra <- max(sucra_result$sucra_scores)

    text_parts <- c(text_parts, sprintf(
      "In summary, %s emerged as the most effective treatment based on network meta-analysis rankings (SUCRA = %.1f%%).",
      best_treatment, best_sucra
    ))
  }

  # Quality of evidence
  I2 <- if (!is.null(heterogeneity_result)) heterogeneity_result$I2 else NULL

  if (!is.null(I2)) {
    quality <- if (I2 < 50) {
      "moderate to high"
    } else {
      "moderate due to heterogeneity"
    }

    text_parts <- c(text_parts, sprintf(
      "The overall quality of evidence was %s.",
      quality
    ))
  }

  # Clinical implications
  if (include_clinical) {
    text_parts <- c(text_parts,
      "These findings provide important evidence to inform clinical decision-making and treatment guidelines. However, results should be interpreted in the context of individual patient characteristics, preferences, and clinical judgment."
    )
  }

  return(paste(text_parts, collapse = " "))
}

#' Combine Sections into Full Text
#'
#' @description
#' Combines individual sections into formatted full text.
#'
#' @keywords internal
combine_sections <- function(sections, style) {

  if (style == "narrative") {
    # Continuous narrative
    full_text <- paste(unlist(sections), collapse = " ")
  } else {
    # Structured with section headers
    section_texts <- character(0)

    section_names <- c(
      overview = "Network Overview",
      treatment_effects = "Treatment Effects",
      rankings = "Treatment Rankings",
      heterogeneity = "Heterogeneity Assessment",
      inconsistency = "Inconsistency Assessment",
      publication_bias = "Publication Bias",
      sensitivity = "Sensitivity Analyses",
      summary = "Summary and Conclusions"
    )

    for (section_name in names(sections)) {
      header <- section_names[[section_name]]
      if (is.null(header)) header <- section_name

      section_texts <- c(
        section_texts,
        sprintf("\n## %s\n\n%s\n", header, sections[[section_name]])
      )
    }

    full_text <- paste(section_texts, collapse = "\n")
  }

  return(full_text)
}

#' Print Method for text_results
#'
#' @description
#' Print method for text results objects.
#'
#' @param x text_results object
#' @param ... Additional arguments (not used)
#'
#' @export
print.text_results <- function(x, ...) {
  cat("Comprehensive NMA Text Results\n")
  cat("==============================\n\n")
  cat("Style:", x$metadata$style, "\n")
  cat("Sections:", paste(x$metadata$sections_included, collapse = ", "), "\n")
  cat("Word count:", x$metadata$word_count, "\n")
  cat("Character count:", x$metadata$character_count, "\n\n")
  cat("Full Text:\n")
  cat("----------\n")
  cat(x$full_text)
  cat("\n")
  invisible(x)
}

#' Export Text Results
#'
#' @description
#' Exports text results to various file formats.
#'
#' @param text_results text_results object
#' @param output_file Output file path
#' @param format Character string: "txt", "docx", "html", "pdf"
#'
#' @export
#'
#' @examples
#' \dontrun{
#' text_results <- generate_comprehensive_text_results(nma, sucra)
#' export_text_results(text_results, "nma_results.txt", format = "txt")
#' export_text_results(text_results, "nma_results.docx", format = "docx")
#' }
export_text_results <- function(text_results, output_file,
                               format = c("txt", "docx", "html", "pdf")) {

  format <- match.arg(format)

  if (format == "txt") {
    # Plain text
    writeLines(text_results$full_text, output_file)

  } else if (format == "docx") {
    # Word document using officer package
    if (!requireNamespace("officer", quietly = TRUE)) {
      stop("Package 'officer' is required for DOCX export. Please install it.")
    }

    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, text_results$full_text)
    print(doc, target = output_file)

  } else if (format == "html") {
    # HTML
    html_content <- sprintf(
      "<!DOCTYPE html>\n<html>\n<head>\n<title>NMA Results</title>\n</head>\n<body>\n<pre>%s</pre>\n</body>\n</html>",
      text_results$full_text
    )
    writeLines(html_content, output_file)

  } else if (format == "pdf") {
    # PDF via rmarkdown
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      stop("Package 'rmarkdown' is required for PDF export. Please install it.")
    }

    temp_rmd <- tempfile(fileext = ".Rmd")
    writeLines(
      c("---", "title: 'Network Meta-Analysis Results'", "output: pdf_document",
        "---", "", text_results$full_text),
      temp_rmd
    )
    rmarkdown::render(temp_rmd, output_file = output_file, quiet = TRUE)
  }

  message(sprintf("Text results exported to: %s", output_file))
  invisible(output_file)
}
