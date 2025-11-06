#' Comprehensive Sensitivity Analysis for Network Meta-Analysis
#'
#' Tools for assessing robustness of NMA results through various sensitivity analyses
#' including leave-one-out, subgroup exclusions, and assumption testing.
#'
#' @name sensitivity_analysis
NULL

#' Leave-One-Out Sensitivity Analysis
#'
#' Assess impact of each study by iteratively excluding one study at a time.
#'
#' @param data Pairwise NMA data
#' @param reference Reference treatment
#' @param sm Summary measure
#' @param comparison Specific comparison to assess (NULL for all)
#' @return loo_sensitivity object
#' @export
#' @examples
#' \dontrun{
#' loo_result <- leave_one_out_nma(data, reference = "Placebo", sm = "OR")
#' print(loo_result)
#' plot(loo_result)
#' }
leave_one_out_nma <- function(data, reference = NULL, sm = "MD", comparison = NULL) {

  studies <- unique(data$studlab)
  n_studies <- length(studies)

  if (n_studies < 3) {
    stop("Need at least 3 studies for leave-one-out analysis")
  }

  # Base NMA with all studies
  base_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  if (is.null(reference)) {
    reference <- base_nma$reference.group
  }

  # Get treatments
  treatments <- rownames(base_nma$TE.random)
  if (!is.null(comparison)) {
    treatments <- comparison
  } else {
    treatments <- setdiff(treatments, reference)
  }

  message(sprintf("Running leave-one-out for %d studies...", n_studies))

  # Run NMA excluding each study
  loo_results <- lapply(studies, function(study) {
    data_loo <- data[data$studlab != study, ]

    nma_loo <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = data_loo,
        sm = sm,
        reference.group = reference
      )
    }, error = function(e) {
      warning(sprintf("LOO NMA failed for study '%s': %s", study, e$message))
      NULL
    })

    if (is.null(nma_loo)) return(NULL)

    # Extract effects for each treatment
    effects <- sapply(treatments, function(trt) {
      if (trt %in% rownames(nma_loo$TE.random)) {
        nma_loo$TE.random[trt, reference]
      } else {
        NA
      }
    })

    list(
      excluded_study = study,
      nma = nma_loo,
      effects = effects
    )
  })

  # Remove failed analyses
  loo_results <- loo_results[!sapply(loo_results, is.null)]

  # Calculate influence statistics
  influence_stats <- calculate_loo_influence(base_nma, loo_results, treatments, reference)

  result <- list(
    base_nma = base_nma,
    loo_results = loo_results,
    influence_stats = influence_stats,
    studies = studies,
    treatments = treatments,
    reference = reference
  )

  class(result) <- c("loo_sensitivity", "list")
  result
}

#' Calculate Leave-One-Out Influence Statistics
#'
#' @keywords internal
calculate_loo_influence <- function(base_nma, loo_results, treatments, reference) {

  base_effects <- sapply(treatments, function(trt) {
    base_nma$TE.random[trt, reference]
  })

  influence <- lapply(treatments, function(trt) {
    # Get LOO effects for this treatment
    loo_effects <- sapply(loo_results, function(loo) loo$effects[trt])

    # Calculate influence metrics
    differences <- loo_effects - base_effects[trt]
    max_diff <- max(abs(differences), na.rm = TRUE)
    influential_study <- loo_results[[which.max(abs(differences))]]$excluded_study

    data.frame(
      Treatment = trt,
      Base_Effect = base_effects[trt],
      Mean_LOO_Effect = mean(loo_effects, na.rm = TRUE),
      Max_Difference = max_diff,
      Most_Influential_Study = influential_study,
      Range_LOO = diff(range(loo_effects, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, influence)
}

#' Cumulative Meta-Analysis
#'
#' Assess how treatment effects evolve as studies are added chronologically.
#'
#' @param data Pairwise NMA data with publication year
#' @param year_var Name of publication year variable
#' @param reference Reference treatment
#' @param sm Summary measure
#' @return cumulative_nma object
#' @export
#' @examples
#' \dontrun{
#' cumulative <- cumulative_meta_analysis(
#'   data = data,
#'   year_var = "pub_year",
#'   reference = "Placebo"
#' )
#' plot(cumulative)
#' }
cumulative_meta_analysis <- function(data, year_var, reference = NULL, sm = "MD") {

  if (!(year_var %in% names(data))) {
    stop(sprintf("Year variable '%s' not found in data", year_var))
  }

  # Sort by year
  data <- data[order(data[[year_var]]), ]

  # Get unique years
  years <- sort(unique(data[[year_var]]))

  message(sprintf("Running cumulative NMA for %d time points...", length(years)))

  cumulative_results <- lapply(seq_along(years), function(i) {
    year <- years[i]
    data_cumulative <- data[data[[year_var]] <= year, ]

    nma <- tryCatch({
      netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = studlab,
        data = data_cumulative,
        sm = sm,
        reference.group = reference
      )
    }, error = function(e) {
      NULL
    })

    list(
      year = year,
      n_studies = nrow(data_cumulative),
      nma = nma
    )
  })

  result <- list(
    cumulative_results = cumulative_results,
    years = years,
    reference = reference
  )

  class(result) <- c("cumulative_nma", "list")
  result
}

#' Fixed vs Random Effects Sensitivity
#'
#' Compare fixed and random effects models.
#'
#' @param data Pairwise NMA data
#' @param reference Reference treatment
#' @param sm Summary measure
#' @return List with both models
#' @export
#' @examples
#' \dontrun{
#' fe_re_comparison <- fixed_vs_random_nma(data, reference = "Placebo")
#' }
fixed_vs_random_nma <- function(data, reference = NULL, sm = "MD") {

  # Fixed effects
  fe_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference,
    comb.fixed = TRUE,
    comb.random = FALSE
  )

  # Random effects
  re_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference,
    comb.fixed = FALSE,
    comb.random = TRUE
  )

  # Compare estimates
  treatments <- setdiff(rownames(fe_nma$TE.fixed), reference)

  comparison_table <- data.frame(
    Treatment = treatments,
    FE_Effect = fe_nma$TE.fixed[treatments, reference],
    FE_SE = fe_nma$seTE.fixed[treatments, reference],
    RE_Effect = re_nma$TE.random[treatments, reference],
    RE_SE = re_nma$seTE.random[treatments, reference],
    Difference = abs(fe_nma$TE.fixed[treatments, reference] -
                    re_nma$TE.random[treatments, reference]),
    stringsAsFactors = FALSE
  )

  list(
    fixed_effects = fe_nma,
    random_effects = re_nma,
    comparison = comparison_table,
    tau = re_nma$tau,
    I2 = re_nma$I2
  )
}

#' Print LOO Sensitivity Results
#'
#' @param x loo_sensitivity object
#' @param ... Additional arguments
#' @export
print.loo_sensitivity <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Leave-One-Out Sensitivity Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Studies analyzed: %d\n", length(x$studies)))
  cat(sprintf("Treatments: %d\n", length(x$treatments)))
  cat(sprintf("Reference: %s\n\n", x$reference))

  cat("Influence Statistics:\n\n")
  print(x$influence_stats, row.names = FALSE, digits = 3)

  cat("\n")
  cat("Most influential study overall: ")
  max_influence_idx <- which.max(x$influence_stats$Max_Difference)
  cat(sprintf("%s (for %s, difference = %.3f)\n",
             x$influence_stats$Most_Influential_Study[max_influence_idx],
             x$influence_stats$Treatment[max_influence_idx],
             x$influence_stats$Max_Difference[max_influence_idx]))

  cat("\n")

  invisible(x)
}

#' Plot LOO Sensitivity Results
#'
#' @param x loo_sensitivity object
#' @param treatment Which treatment to plot
#' @param ... Additional arguments
#' @export
plot.loo_sensitivity <- function(x, treatment = NULL, ...) {

  if (is.null(treatment)) {
    treatment <- x$treatments[1]
  }

  # Extract LOO effects for this treatment
  loo_effects <- sapply(x$loo_results, function(loo) loo$effects[treatment])
  excluded_studies <- sapply(x$loo_results, function(loo) loo$excluded_study)
  base_effect <- x$base_nma$TE.random[treatment, x$reference]

  plot_data <- data.frame(
    Study = excluded_studies,
    Effect = loo_effects,
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[order(plot_data$Effect), ]
  plot_data$Study <- factor(plot_data$Study, levels = plot_data$Study)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Study, y = Effect)) +
    ggplot2::geom_hline(yintercept = base_effect, linetype = "dashed",
                       color = "red", size = 1) +
    ggplot2::geom_point(size = 3, color = "#2c7fb8") +
    ggplot2::labs(
      title = "Leave-One-Out Sensitivity Analysis",
      subtitle = sprintf("%s vs %s (red line = full model)", treatment, x$reference),
      x = "Excluded Study",
      y = sprintf("Treatment Effect (%s)", x$base_nma$sm)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
