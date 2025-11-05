#' Living Network Meta-Analysis Framework
#'
#' Framework for continuously updating network meta-analyses as new evidence
#' emerges. Implements sequential updating protocols, trigger rules, and
#' version control following Elliott et al. (2024).
#'
#' @name living_nma
#' @references
#' Elliott JH, et al. (2024). Living systematic reviews and network meta-analysis:
#' 3. Statistical methods for updating. Mayo Clinic Proceedings.
#'
#' Thomas J, et al. (2017). Living systematic reviews: 2. Combining human and
#' machine effort. Journal of Clinical Epidemiology, 91:31-37.
NULL

#' Initialize Living NMA Project
#'
#' Set up a living NMA project with version control and update protocols.
#'
#' @param baseline_data Initial pairwise NMA data
#' @param project_name Name of the living NMA project
#' @param update_frequency Expected update frequency ("monthly", "quarterly", "as_needed")
#' @param trigger_rules List of trigger rules for when to update
#' @param output_dir Directory to store versions and logs
#' @return living_nma_project object
#' @export
#' @examples
#' \dontrun{
#' living_project <- initialize_living_nma(
#'   baseline_data = nma_data,
#'   project_name = "Antidepressants_NMA",
#'   update_frequency = "quarterly",
#'   trigger_rules = list(
#'     new_studies = 3,         # Update when 3+ new studies
#'     relative_change = 0.15    # Update if effect changes >15%
#'   ),
#'   output_dir = "living_nma_versions"
#' )
#' }
initialize_living_nma <- function(baseline_data,
                                  project_name,
                                  update_frequency = c("monthly", "quarterly", "as_needed"),
                                  trigger_rules = NULL,
                                  output_dir = NULL) {

  update_frequency <- match.arg(update_frequency)

  # Create output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), paste0(project_name, "_living_nma"))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(sprintf("Created living NMA directory: %s", output_dir))
  }

  # Default trigger rules
  if (is.null(trigger_rules)) {
    trigger_rules <- list(
      new_studies_threshold = 5,
      relative_effect_change = 0.20,
      new_treatment = TRUE,
      time_threshold_days = 90
    )
  }

  # Run baseline NMA
  message("Running baseline network meta-analysis...")

  baseline_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = baseline_data
  )

  # Initialize version 1
  version_info <- list(
    version = 1,
    date = Sys.Date(),
    n_studies = length(unique(baseline_data$studlab)),
    n_treatments = length(baseline_nma$trts),
    study_list = unique(baseline_data$studlab),
    change_log = "Initial baseline version"
  )

  # Save baseline
  baseline_file <- file.path(output_dir, "version_001_baseline.rds")
  saveRDS(list(
    data = baseline_data,
    nma = baseline_nma,
    version_info = version_info
  ), file = baseline_file)

  message(sprintf("Baseline saved: %s", baseline_file))

  project <- list(
    project_name = project_name,
    current_version = 1,
    baseline_data = baseline_data,
    baseline_nma = baseline_nma,
    current_data = baseline_data,
    current_nma = baseline_nma,
    update_frequency = update_frequency,
    trigger_rules = trigger_rules,
    output_dir = output_dir,
    version_history = list(version_info),
    last_update = Sys.Date()
  )

  class(project) <- c("living_nma_project", "list")

  # Save project metadata
  metadata_file <- file.path(output_dir, "project_metadata.rds")
  saveRDS(project, file = metadata_file)

  message("Living NMA project initialized successfully")

  project
}

#' Update Living NMA with New Evidence
#'
#' Incorporate new studies into the living NMA and check trigger rules.
#'
#' @param project living_nma_project object
#' @param new_data New pairwise comparison data to add
#' @param force_update Force update even if triggers not met
#' @param update_note Description of what's being added
#' @return Updated living_nma_project
#' @export
#' @examples
#' \dontrun{
#' # New studies identified
#' new_studies <- read.csv("new_studies_q2.csv")
#'
#' updated_project <- update_living_nma(
#'   project = living_project,
#'   new_data = new_studies,
#'   update_note = "Q2 2025 update: 4 new RCTs"
#' )
#' }
update_living_nma <- function(project, new_data, force_update = FALSE,
                             update_note = "Regular update") {

  if (!inherits(project, "living_nma_project")) {
    stop("Input must be living_nma_project object")
  }

  message("\n")
  message("════════════════════════════════════════════════════════")
  message(sprintf("  Living NMA Update: %s", project$project_name))
  message("════════════════════════════════════════════════════════")
  message("\n")

  # Check for new studies
  existing_studies <- unique(project$current_data$studlab)
  new_studies_list <- setdiff(unique(new_data$studlab), existing_studies)
  n_new_studies <- length(new_studies_list)

  message(sprintf("Current version: %d", project$current_version))
  message(sprintf("Studies in current version: %d", length(existing_studies)))
  message(sprintf("New studies to add: %d", n_new_studies))

  if (n_new_studies == 0 && !force_update) {
    message("\nNo new studies to add. Update skipped.")
    return(project)
  }

  # Combine data
  combined_data <- rbind(project$current_data, new_data)

  # Check for new treatments
  old_treatments <- unique(c(as.character(project$current_data$treat1),
                            as.character(project$current_data$treat2)))
  new_treatments <- setdiff(
    unique(c(as.character(new_data$treat1), as.character(new_data$treat2))),
    old_treatments
  )

  if (length(new_treatments) > 0) {
    message(sprintf("\n✓ New treatments added: %s", paste(new_treatments, collapse = ", ")))
  }

  # Run updated NMA
  message("\nRunning updated network meta-analysis...")

  updated_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = combined_data
  )

  # Compare with previous version
  comparison <- compare_nma_versions(project$current_nma, updated_nma)

  message("\nEffect Size Changes:")
  print(comparison$effect_changes, row.names = FALSE, digits = 3)

  # Check trigger rules
  triggers_met <- check_trigger_rules(
    n_new_studies = n_new_studies,
    effect_changes = comparison$effect_changes,
    new_treatments = new_treatments,
    days_since_update = as.numeric(difftime(Sys.Date(), project$last_update, units = "days")),
    trigger_rules = project$trigger_rules
  )

  message("\nTrigger Rules:")
  for (rule_name in names(triggers_met)) {
    status <- if (triggers_met[[rule_name]]) "✓ MET" else "○ Not met"
    message(sprintf("  %s: %s", rule_name, status))
  }

  should_update <- any(unlist(triggers_met)) || force_update

  if (!should_update) {
    message("\n⚠ Trigger rules not met. Update not recommended.")
    message("  Use force_update = TRUE to override.")
    return(project)
  }

  # Proceed with update
  message("\n✓ Proceeding with update...")

  # Create new version
  new_version <- project$current_version + 1

  version_info <- list(
    version = new_version,
    date = Sys.Date(),
    n_studies = length(unique(combined_data$studlab)),
    n_treatments = length(updated_nma$trts),
    n_new_studies = n_new_studies,
    new_studies = new_studies_list,
    new_treatments = if (length(new_treatments) > 0) new_treatments else NULL,
    change_log = update_note,
    comparison = comparison,
    triggers_met = triggers_met
  )

  # Save new version
  version_file <- file.path(project$output_dir,
                           sprintf("version_%03d.rds", new_version))

  saveRDS(list(
    data = combined_data,
    nma = updated_nma,
    version_info = version_info
  ), file = version_file)

  message(sprintf("\nNew version saved: %s", version_file))

  # Update project
  project$current_version <- new_version
  project$current_data <- combined_data
  project$current_nma <- updated_nma
  project$last_update <- Sys.Date()
  project$version_history[[new_version]] <- version_info

  # Save updated metadata
  metadata_file <- file.path(project$output_dir, "project_metadata.rds")
  saveRDS(project, file = metadata_file)

  message("\n✓ Living NMA updated successfully")

  # Generate update report
  generate_update_report(project, version_info, project$output_dir)

  project
}

#' Compare Two NMA Versions
#'
#' Compare treatment effects between two versions of the NMA.
#'
#' @param old_nma Previous NMA result
#' @param new_nma New NMA result
#' @return Comparison results
#' @keywords internal
compare_nma_versions <- function(old_nma, new_nma) {

  # Get common treatments
  common_treatments <- intersect(old_nma$trts, new_nma$trts)

  if (length(common_treatments) < 2) {
    return(list(
      effect_changes = data.frame(),
      max_change = NA
    ))
  }

  reference <- old_nma$reference.group

  # Extract effects
  changes <- lapply(common_treatments, function(trt) {
    if (trt == reference) return(NULL)

    old_effect <- old_nma$TE.random[trt, reference]
    new_effect <- new_nma$TE.random[trt, reference]

    old_se <- old_nma$seTE.random[trt, reference]
    new_se <- new_nma$seTE.random[trt, reference]

    abs_change <- new_effect - old_effect
    rel_change <- abs_change / abs(old_effect)

    # Test for significant change
    se_diff <- sqrt(old_se^2 + new_se^2)
    z_change <- abs_change / se_diff
    p_change <- 2 * (1 - pnorm(abs(z_change)))

    data.frame(
      Treatment = trt,
      Old_Effect = old_effect,
      New_Effect = new_effect,
      Absolute_Change = abs_change,
      Relative_Change = rel_change,
      Change_Significant = p_change < 0.05,
      stringsAsFactors = FALSE
    )
  })

  changes <- do.call(rbind, changes[!sapply(changes, is.null)])

  list(
    effect_changes = changes,
    max_absolute_change = max(abs(changes$Absolute_Change), na.rm = TRUE),
    max_relative_change = max(abs(changes$Relative_Change), na.rm = TRUE)
  )
}

#' Check Trigger Rules for Update
#'
#' Determine if any trigger rules are met for updating the living NMA.
#'
#' @keywords internal
check_trigger_rules <- function(n_new_studies, effect_changes, new_treatments,
                               days_since_update, trigger_rules) {

  triggers <- list()

  # Rule 1: Number of new studies
  if (!is.null(trigger_rules$new_studies_threshold)) {
    triggers$new_studies <- n_new_studies >= trigger_rules$new_studies_threshold
  }

  # Rule 2: Relative effect change
  if (!is.null(trigger_rules$relative_effect_change) && nrow(effect_changes) > 0) {
    max_change <- max(abs(effect_changes$Relative_Change), na.rm = TRUE)
    triggers$effect_change <- max_change >= trigger_rules$relative_effect_change
  }

  # Rule 3: New treatment
  if (!is.null(trigger_rules$new_treatment)) {
    triggers$new_treatment <- length(new_treatments) > 0 && trigger_rules$new_treatment
  }

  # Rule 4: Time threshold
  if (!is.null(trigger_rules$time_threshold_days)) {
    triggers$time_threshold <- days_since_update >= trigger_rules$time_threshold_days
  }

  triggers
}

#' Generate Update Report
#'
#' Create a summary report of what changed in the update.
#'
#' @keywords internal
generate_update_report <- function(project, version_info, output_dir) {

  report_file <- file.path(output_dir,
                          sprintf("update_report_v%03d.txt", version_info$version))

  sink(report_file)

  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  Living NMA Update Report\n"))
  cat(sprintf("  Project: %s\n", project$project_name))
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Version: %d\n", version_info$version))
  cat(sprintf("Date: %s\n\n", version_info$date))

  cat("Changes:\n")
  cat(sprintf("  • New studies added: %d\n", version_info$n_new_studies))

  if (!is.null(version_info$new_treatments)) {
    cat(sprintf("  • New treatments: %s\n", paste(version_info$new_treatments, collapse = ", ")))
  }

  cat(sprintf("\nTotal studies: %d\n", version_info$n_studies))
  cat(sprintf("Total treatments: %d\n\n", version_info$n_treatments))

  if (!is.null(version_info$comparison$effect_changes)) {
    cat("Treatment Effect Changes:\n\n")
    print(version_info$comparison$effect_changes, row.names = FALSE, digits = 3)
  }

  cat("\n\nTrigger Rules Met:\n")
  for (rule in names(version_info$triggers_met)) {
    status <- if (version_info$triggers_met[[rule]]) "YES" else "NO"
    cat(sprintf("  • %s: %s\n", rule, status))
  }

  cat(sprintf("\nNotes: %s\n", version_info$change_log))

  sink()

  message(sprintf("Update report saved: %s", report_file))
}

#' Load Living NMA Project
#'
#' Load a previously saved living NMA project.
#'
#' @param project_dir Directory containing the project
#' @return living_nma_project object
#' @export
#' @examples
#' \dontrun{
#' project <- load_living_nma("Antidepressants_NMA_living_nma")
#' }
load_living_nma <- function(project_dir) {

  metadata_file <- file.path(project_dir, "project_metadata.rds")

  if (!file.exists(metadata_file)) {
    stop(sprintf("Project metadata not found in: %s", project_dir))
  }

  project <- readRDS(metadata_file)

  message(sprintf("Loaded living NMA project: %s", project$project_name))
  message(sprintf("Current version: %d", project$current_version))
  message(sprintf("Last update: %s", project$last_update))

  project
}

#' Get Version History
#'
#' Retrieve the version history of a living NMA project.
#'
#' @param project living_nma_project object
#' @return Data frame with version history
#' @export
get_version_history <- function(project) {

  if (!inherits(project, "living_nma_project")) {
    stop("Input must be living_nma_project object")
  }

  history <- lapply(project$version_history, function(v) {
    data.frame(
      Version = v$version,
      Date = as.character(v$date),
      N_Studies = v$n_studies,
      N_Treatments = v$n_treatments,
      N_New_Studies = if (is.null(v$n_new_studies)) 0 else v$n_new_studies,
      Change_Log = v$change_log,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, history)
}

#' Plot Version Timeline
#'
#' Visualize how the network has evolved over versions.
#'
#' @param project living_nma_project object
#' @return ggplot2 object
#' @export
plot_version_timeline <- function(project) {

  if (!inherits(project, "living_nma_project")) {
    stop("Input must be living_nma_project object")
  }

  history <- get_version_history(project)
  history$Date <- as.Date(history$Date)

  # Reshape for plotting
  plot_data <- data.frame(
    Version = rep(history$Version, 2),
    Date = rep(history$Date, 2),
    Count = c(history$N_Studies, history$N_Treatments),
    Type = rep(c("Studies", "Treatments"), each = nrow(history)),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Date, y = Count, color = Type, group = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(data = subset(plot_data, Type == "Studies"),
                      ggplot2::aes(label = paste0("v", Version)),
                      vjust = -1, size = 3) +
    ggplot2::scale_color_manual(values = c("Studies" = "#2c7fb8", "Treatments" = "#d7191c")) +
    ggplot2::labs(
      title = sprintf("Living NMA Timeline: %s", project$project_name),
      subtitle = sprintf("%d versions • Last update: %s", project$current_version, project$last_update),
      x = "Date",
      y = "Count",
      color = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  print(p)
  invisible(p)
}

#' Print Living NMA Project
#'
#' @param x living_nma_project object
#' @param ... Additional arguments
#' @export
print.living_nma_project <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  Living Network Meta-Analysis Project\n"))
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Project: %s\n", x$project_name))
  cat(sprintf("Current version: %d\n", x$current_version))
  cat(sprintf("Last update: %s\n", x$last_update))
  cat(sprintf("Update frequency: %s\n\n", x$update_frequency))

  cat("Current Network:\n")
  cat(sprintf("  • Studies: %d\n", length(unique(x$current_data$studlab))))
  cat(sprintf("  • Treatments: %d\n", length(x$current_nma$trts)))
  cat(sprintf("  • Comparisons: %d\n\n", nrow(x$current_data)))

  cat("Trigger Rules:\n")
  for (rule in names(x$trigger_rules)) {
    cat(sprintf("  • %s: %s\n", rule, x$trigger_rules[[rule]]))
  }

  cat(sprintf("\nOutput directory: %s\n", x$output_dir))

  cat("\n")
  cat("Use get_version_history() to view all versions\n")
  cat("Use update_living_nma() to add new evidence\n")
  cat("Use plot_version_timeline() to visualize evolution\n\n")

  invisible(x)
}

#' Compare Treatment Effect Across Versions
#'
#' Track how a specific treatment effect has changed across versions.
#'
#' @param project living_nma_project object
#' @param treatment Treatment name
#' @param reference Reference treatment (if NULL, uses project reference)
#' @return Data frame with effect over time
#' @export
track_treatment_effect <- function(project, treatment, reference = NULL) {

  if (!inherits(project, "living_nma_project")) {
    stop("Input must be living_nma_project object")
  }

  if (is.null(reference)) {
    reference <- project$current_nma$reference.group
  }

  # Load all versions
  versions <- list()

  for (v in seq_len(project$current_version)) {
    version_file <- file.path(project$output_dir, sprintf("version_%03d*.rds", v))
    files <- Sys.glob(version_file)

    if (length(files) > 0) {
      version_data <- readRDS(files[1])

      if (treatment %in% version_data$nma$trts) {
        effect <- version_data$nma$TE.random[treatment, reference]
        se <- version_data$nma$seTE.random[treatment, reference]

        versions[[v]] <- data.frame(
          Version = v,
          Date = version_data$version_info$date,
          Effect = effect,
          SE = se,
          Lower_95 = effect - 1.96 * se,
          Upper_95 = effect + 1.96 * se,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  result <- do.call(rbind, versions)
  result$Date <- as.Date(result$Date)

  # Plot
  p <- ggplot2::ggplot(result, ggplot2::aes(x = Date, y = Effect)) +
    ggplot2::geom_line(size = 1, color = "#2c7fb8") +
    ggplot2::geom_point(size = 3, color = "#2c7fb8") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
                        alpha = 0.2, fill = "#2c7fb8") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = sprintf("Treatment Effect Over Time: %s vs %s", treatment, reference),
      subtitle = sprintf("Project: %s", project$project_name),
      x = "Version Date",
      y = sprintf("Effect (95%% CI)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )

  print(p)

  invisible(result)
}
