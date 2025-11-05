#' Batch Processing and Workflow Automation
#'
#' @description
#' Advanced batch processing capabilities for running multiple NMA analyses,
#' automated workflows, parallel processing, and comprehensive result aggregation.
#'
#' @details
#' This module provides enterprise-grade batch processing features:
#' \itemize{
#'   \item Multiple dataset batch analysis
#'   \item Parallel processing with progress tracking
#'   \item Automated workflow pipelines
#'   \item Result aggregation and comparison
#'   \item Error handling and logging
#'   \item Resumable batch jobs
#' }
#'
#' @author powerNMA Development Team
#' @name batch_processing
NULL

#' Run Batch NMA Analysis
#'
#' @description
#' Processes multiple datasets in batch with comprehensive analysis pipeline.
#'
#' @param datasets List of datasets or paths to data files
#' @param analysis_config Analysis configuration (from setup_powernma())
#' @param parallel Logical whether to use parallel processing
#' @param n_cores Number of cores for parallel processing (default: auto-detect)
#' @param save_individual Logical whether to save individual results
#' @param output_dir Directory for saving results
#' @param log_file Path to log file (if NULL, uses temp file)
#' @param resume_from Resume from specific dataset index (for interrupted jobs)
#'
#' @return List with batch results and summary
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare multiple datasets
#' datasets <- list(
#'   depression = depression_data,
#'   anxiety = anxiety_data,
#'   schizophrenia = schizo_data
#' )
#'
#' # Configure analysis
#' config <- setup_powernma(
#'   sm = "SMD",
#'   use_bayesian = FALSE,
#'   run_sensitivity = TRUE
#' )
#'
#' # Run batch analysis
#' batch_results <- run_batch_nma(
#'   datasets = datasets,
#'   analysis_config = config,
#'   parallel = TRUE,
#'   n_cores = 4,
#'   output_dir = "batch_results"
#' )
#'
#' # View summary
#' print(batch_results$summary)
#'
#' # Access individual results
#' depression_results <- batch_results$results$depression
#' }
run_batch_nma <- function(datasets,
                         analysis_config = setup_powernma(),
                         parallel = FALSE,
                         n_cores = NULL,
                         save_individual = TRUE,
                         output_dir = "batch_results",
                         log_file = NULL,
                         resume_from = 1) {

  # Setup
  start_time <- Sys.time()
  n_datasets <- length(datasets)

  # Initialize logging
  if (is.null(log_file)) {
    log_file <- tempfile(pattern = "batch_nma_", fileext = ".log")
  }
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  write_log <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", timestamp, msg), file = log_con)
    message(msg)
  }

  write_log(sprintf("Starting batch NMA analysis: %d datasets", n_datasets))

  # Create output directory
  if (save_individual && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Initialize results storage
  results <- list()
  errors <- list()
  timings <- numeric(n_datasets)

  # Determine parallel settings
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    write_log(sprintf("Using parallel processing with %d cores", n_cores))
  }

  # Process function for each dataset
  process_dataset <- function(i, data_name, data) {
    write_log(sprintf("Processing dataset %d/%d: %s", i, n_datasets, data_name))

    dataset_start <- Sys.time()

    tryCatch({
      # Load data if path provided
      if (is.character(data)) {
        data <- load_data_file(data)
      }

      # Validate data
      validate_nma_input(data)

      # Run comprehensive analysis
      result <- run_ultimate_nma(
        data = data,
        sm = analysis_config$sm,
        reference = analysis_config$reference,
        assess_geometry = TRUE,
        assess_inconsistency = analysis_config$assess_inconsistency,
        calculate_rankings = TRUE,
        assess_publication_bias = analysis_config$assess_publication_bias,
        export_results = save_individual,
        export_dir = file.path(output_dir, data_name)
      )

      # Generate manuscripts if requested
      if (analysis_config$generate_manuscripts) {
        result$manuscripts <- generate_manuscripts_pipeline(result)
      }

      dataset_time <- as.numeric(difftime(Sys.time(), dataset_start, units = "secs"))

      write_log(sprintf("Completed %s in %.2f seconds", data_name, dataset_time))

      return(list(
        success = TRUE,
        result = result,
        timing = dataset_time,
        error = NULL
      ))

    }, error = function(e) {
      write_log(sprintf("ERROR in %s: %s", data_name, e$message))

      return(list(
        success = FALSE,
        result = NULL,
        timing = as.numeric(difftime(Sys.time(), dataset_start, units = "secs")),
        error = e$message
      ))
    })
  }

  # Run batch processing
  if (parallel && n_cores > 1) {
    # Parallel processing
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary functions and packages
    parallel::clusterExport(cl, c("run_ultimate_nma", "validate_nma_input"),
                           envir = environment())

    batch_output <- parallel::parLapply(
      cl,
      seq(resume_from, n_datasets),
      function(i) {
        data_name <- names(datasets)[i]
        data <- datasets[[i]]
        process_dataset(i, data_name, data)
      }
    )

  } else {
    # Sequential processing with progress bar
    pb <- txtProgressBar(min = resume_from, max = n_datasets, style = 3)

    batch_output <- lapply(seq(resume_from, n_datasets), function(i) {
      setTxtProgressBar(pb, i)
      data_name <- names(datasets)[i]
      data <- datasets[[i]]
      process_dataset(i, data_name, data)
    })

    close(pb)
  }

  # Aggregate results
  names(batch_output) <- names(datasets)[resume_from:n_datasets]

  for (name in names(batch_output)) {
    if (batch_output[[name]]$success) {
      results[[name]] <- batch_output[[name]]$result
      timings[[name]] <- batch_output[[name]]$timing
    } else {
      errors[[name]] <- batch_output[[name]]$error
    }
  }

  # Create summary
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  summary_stats <- data.frame(
    Dataset = names(datasets)[resume_from:n_datasets],
    Success = sapply(batch_output, function(x) x$success),
    Time_Seconds = sapply(batch_output, function(x) x$timing),
    Error = sapply(batch_output, function(x) x$error %||% "None"),
    stringsAsFactors = FALSE
  )

  write_log(sprintf("Batch processing complete: %d/%d successful",
                   sum(summary_stats$Success), nrow(summary_stats)))
  write_log(sprintf("Total time: %.2f seconds", total_time))

  # Aggregate statistics across successful analyses
  aggregate_stats <- if (length(results) > 0) {
    aggregate_batch_statistics(results)
  } else {
    NULL
  }

  # Create comparison table
  comparison_table <- if (length(results) > 1) {
    create_batch_comparison_table(results)
  } else {
    NULL
  }

  return(structure(
    list(
      results = results,
      summary = summary_stats,
      aggregate_statistics = aggregate_stats,
      comparison_table = comparison_table,
      errors = errors,
      total_time = total_time,
      log_file = log_file,
      output_dir = output_dir
    ),
    class = "batch_nma_results"
  ))
}

#' Aggregate Batch Statistics
#'
#' @description
#' Aggregates statistics across multiple NMA results.
#'
#' @param results List of NMA results
#'
#' @return Data frame with aggregate statistics
#'
#' @keywords internal
aggregate_batch_statistics <- function(results) {

  n_results <- length(results)

  stats <- data.frame(
    Dataset = names(results),
    N_Studies = sapply(results, function(r) r$nma_result$k),
    N_Treatments = sapply(results, function(r) length(r$nma_result$trts)),
    N_Comparisons = sapply(results, function(r) r$nma_result$m),
    Tau2 = sapply(results, function(r) r$nma_result$tau^2),
    I2 = sapply(results, function(r) {
      if (!is.null(r$heterogeneity)) r$heterogeneity$I2 else NA
    }),
    Best_Treatment = sapply(results, function(r) {
      if (!is.null(r$sucra)) {
        names(which.max(r$sucra$sucra_scores))
      } else {
        NA
      }
    }),
    Best_SUCRA = sapply(results, function(r) {
      if (!is.null(r$sucra)) {
        max(r$sucra$sucra_scores)
      } else {
        NA
      }
    }),
    Inconsistency_Detected = sapply(results, function(r) {
      if (!is.null(r$inconsistency)) {
        any(r$inconsistency$node_splitting$p_value < 0.05, na.rm = TRUE)
      } else {
        NA
      }
    }),
    stringsAsFactors = FALSE
  )

  return(stats)
}

#' Create Batch Comparison Table
#'
#' @description
#' Creates comparison table across batch results.
#'
#' @param results List of NMA results
#'
#' @return Data frame with comparisons
#'
#' @keywords internal
create_batch_comparison_table <- function(results) {

  # Get all unique treatments across datasets
  all_treatments <- unique(unlist(lapply(results, function(r) r$nma_result$trts)))

  # Create ranking matrix
  ranking_matrix <- matrix(NA, nrow = length(all_treatments), ncol = length(results))
  rownames(ranking_matrix) <- all_treatments
  colnames(ranking_matrix) <- names(results)

  for (i in seq_along(results)) {
    if (!is.null(results[[i]]$sucra)) {
      sucra_scores <- results[[i]]$sucra$sucra_scores
      for (trt in names(sucra_scores)) {
        if (trt %in% all_treatments) {
          ranking_matrix[trt, i] <- sucra_scores[[trt]]
        }
      }
    }
  }

  # Convert to data frame
  comparison_df <- as.data.frame(ranking_matrix)
  comparison_df$Treatment <- rownames(ranking_matrix)
  comparison_df$Mean_SUCRA <- rowMeans(ranking_matrix, na.rm = TRUE)
  comparison_df$SD_SUCRA <- apply(ranking_matrix, 1, sd, na.rm = TRUE)

  # Reorder columns
  comparison_df <- comparison_df[, c("Treatment", "Mean_SUCRA", "SD_SUCRA", names(results))]

  # Sort by mean SUCRA
  comparison_df <- comparison_df[order(comparison_df$Mean_SUCRA, decreasing = TRUE), ]

  return(comparison_df)
}

#' Automated Workflow Pipeline
#'
#' @description
#' Runs complete automated workflow from data to publication-ready outputs.
#'
#' @param data Dataset or path to data file
#' @param workflow_config Workflow configuration list
#' @param output_dir Output directory for all results
#'
#' @return List with all workflow outputs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define workflow
#' workflow <- list(
#'   analysis = list(sm = "OR", assess_all = TRUE),
#'   visualizations = c("network", "forest", "funnel", "rankings"),
#'   manuscripts = list(generate = TRUE, use_ai = TRUE, style = "detailed"),
#'   reports = list(format = c("docx", "pdf"), include_appendix = TRUE)
#' )
#'
#' # Run automated workflow
#' outputs <- automated_nma_workflow(
#'   data = my_data,
#'   workflow_config = workflow,
#'   output_dir = "complete_analysis"
#' )
#' }
automated_nma_workflow <- function(data,
                                  workflow_config = default_workflow_config(),
                                  output_dir = "automated_workflow") {

  message("Starting automated NMA workflow...")
  workflow_start <- Sys.time()

  # Create output directory structure
  dirs <- create_workflow_directories(output_dir)

  # Initialize results list
  outputs <- list()

  # Step 1: Data validation and preparation
  message("Step 1/7: Data validation...")
  outputs$validation <- validate_nma_input(data)

  # Step 2: Run comprehensive analysis
  message("Step 2/7: Running comprehensive analysis...")
  outputs$analysis <- run_ultimate_nma(
    data = data,
    sm = workflow_config$analysis$sm,
    assess_geometry = TRUE,
    assess_inconsistency = TRUE,
    calculate_rankings = TRUE,
    assess_publication_bias = TRUE,
    export_results = FALSE
  )

  # Step 3: Generate visualizations
  if (!is.null(workflow_config$visualizations)) {
    message("Step 3/7: Generating visualizations...")
    outputs$visualizations <- generate_all_visualizations(
      nma_result = outputs$analysis$nma_result,
      data = data,
      sucra_result = outputs$analysis$sucra,
      viz_types = workflow_config$visualizations,
      output_dir = dirs$plots
    )
  }

  # Step 4: Generate manuscripts
  if (workflow_config$manuscripts$generate) {
    message("Step 4/7: Generating manuscript sections...")
    outputs$manuscripts <- generate_manuscripts_pipeline(
      nma_result = outputs$analysis,
      use_ai = workflow_config$manuscripts$use_ai,
      style = workflow_config$manuscripts$style,
      output_dir = dirs$manuscripts
    )
  }

  # Step 5: Generate comprehensive text results
  message("Step 5/7: Generating text results...")
  outputs$text_results <- generate_comprehensive_text_results(
    nma_result = outputs$analysis$nma_result,
    sucra_result = outputs$analysis$sucra,
    heterogeneity_result = outputs$analysis$heterogeneity,
    inconsistency_result = outputs$analysis$inconsistency,
    style = "detailed",
    include_clinical = TRUE
  )

  # Step 6: Generate reports
  if (!is.null(workflow_config$reports)) {
    message("Step 6/7: Generating reports...")
    outputs$reports <- generate_automated_reports(
      analysis_results = outputs$analysis,
      formats = workflow_config$reports$format,
      include_appendix = workflow_config$reports$include_appendix,
      output_dir = dirs$reports
    )
  }

  # Step 7: Create summary dashboard
  message("Step 7/7: Creating summary dashboard...")
  outputs$dashboard <- create_workflow_summary_dashboard(
    outputs = outputs,
    output_file = file.path(output_dir, "workflow_summary.html")
  )

  workflow_time <- as.numeric(difftime(Sys.time(), workflow_start, units = "secs"))
  message(sprintf("Workflow complete! Total time: %.2f seconds", workflow_time))

  outputs$workflow_info <- list(
    start_time = workflow_start,
    end_time = Sys.time(),
    total_time = workflow_time,
    output_dir = output_dir,
    config = workflow_config
  )

  return(structure(outputs, class = "automated_workflow_results"))
}

#' Default Workflow Configuration
#'
#' @return List with default workflow settings
#' @keywords internal
default_workflow_config <- function() {
  list(
    analysis = list(
      sm = "OR",
      assess_all = TRUE
    ),
    visualizations = c("network", "forest", "funnel", "rankings", "heatmap"),
    manuscripts = list(
      generate = TRUE,
      use_ai = FALSE,
      style = "detailed"
    ),
    reports = list(
      format = c("docx", "html"),
      include_appendix = TRUE
    )
  )
}

#' Create Workflow Directory Structure
#'
#' @keywords internal
create_workflow_directories <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    plots = file.path(base_dir, "plots"),
    manuscripts = file.path(base_dir, "manuscripts"),
    reports = file.path(base_dir, "reports"),
    data = file.path(base_dir, "data"),
    tables = file.path(base_dir, "tables")
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  return(dirs)
}

#' Generate All Visualizations
#'
#' @keywords internal
generate_all_visualizations <- function(nma_result, data, sucra_result,
                                       viz_types, output_dir) {

  plots <- list()

  for (viz_type in viz_types) {
    filename <- file.path(output_dir, sprintf("%s_plot.html", viz_type))

    plot_obj <- switch(
      viz_type,
      network = create_interactive_network(data, nma_result),
      forest = create_interactive_forest(nma_result),
      funnel = create_interactive_funnel(nma_result),
      rankings = create_interactive_rankings(sucra_result),
      heatmap = create_interactive_heatmap(nma_result),
      network_3d = create_3d_network(data, nma_result)
    )

    if (!is.null(plot_obj)) {
      htmlwidgets::saveWidget(plot_obj, filename)
      plots[[viz_type]] <- filename
    }
  }

  return(plots)
}

#' Generate Manuscripts Pipeline
#'
#' @keywords internal
generate_manuscripts_pipeline <- function(nma_result, use_ai = FALSE,
                                         style = "detailed", output_dir = NULL) {

  # Generate Methods
  methods <- generate_methods_section(
    nma_result = nma_result$nma_result,
    style = style,
    use_ai = FALSE
  )

  # Generate Results
  results_section <- generate_results_section(
    nma_result = nma_result$nma_result,
    sucra_result = nma_result$sucra,
    heterogeneity_result = nma_result$heterogeneity,
    inconsistency_result = nma_result$inconsistency,
    style = style,
    use_ai = FALSE
  )

  # AI Enhancement if requested
  if (use_ai && check_ollama_available()) {
    methods <- enhance_methods_with_ollama(methods$text, enhancement_level = "moderate")
    results_section <- enhance_results_with_ollama(results_section$text,
                                                   add_clinical_context = TRUE,
                                                   enhancement_level = "moderate")
  }

  # Save to files if output_dir provided
  if (!is.null(output_dir)) {
    writeLines(methods$text %||% methods$enhanced_text,
              file.path(output_dir, "methods_section.txt"))
    writeLines(results_section$text %||% results_section$enhanced_text,
              file.path(output_dir, "results_section.txt"))
  }

  return(list(
    methods = methods,
    results = results_section
  ))
}

#' Generate Automated Reports
#'
#' @keywords internal
generate_automated_reports <- function(analysis_results, formats = c("docx"),
                                      include_appendix = TRUE, output_dir = ".") {

  report_files <- list()

  for (fmt in formats) {
    filename <- file.path(output_dir, sprintf("nma_report.%s", fmt))

    tryCatch({
      generate_nma_report(
        nma_results = analysis_results,
        output_file = filename,
        format = fmt,
        include_appendix = include_appendix
      )

      report_files[[fmt]] <- filename
    }, error = function(e) {
      warning(sprintf("Failed to generate %s report: %s", fmt, e$message))
    })
  }

  return(report_files)
}

#' Create Workflow Summary Dashboard
#'
#' @keywords internal
create_workflow_summary_dashboard <- function(outputs, output_file) {

  html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
  <title>NMA Workflow Summary</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h1 { color: #2c3e50; }
    h2 { color: #3498db; border-bottom: 2px solid #3498db; padding-bottom: 5px; }
    .section { margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 5px; }
    .success { color: #27ae60; }
    .info { color: #3498db; }
    table { border-collapse: collapse; width: 100%%; margin: 10px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #3498db; color: white; }
  </style>
</head>
<body>
  <h1>Network Meta-Analysis Workflow Summary</h1>
  <p class="info">Generated: %s</p>

  <div class="section">
    <h2>Analysis Overview</h2>
    <p><strong>Number of Studies:</strong> %d</p>
    <p><strong>Number of Treatments:</strong> %d</p>
    <p><strong>Best Treatment (SUCRA):</strong> %s (%.2f%%)</p>
    <p><strong>Heterogeneity (I²):</strong> %.1f%%</p>
  </div>

  <div class="section">
    <h2>Workflow Performance</h2>
    <p><strong>Total Time:</strong> %.2f seconds</p>
    <p><strong>Outputs Generated:</strong> %d</p>
  </div>

  <div class="section">
    <h2 class="success">✓ Workflow Complete</h2>
    <p>All requested outputs have been generated successfully.</p>
  </div>
</body>
</html>
  ',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    outputs$analysis$nma_result$k,
    length(outputs$analysis$nma_result$trts),
    names(which.max(outputs$analysis$sucra$sucra_scores)),
    max(outputs$analysis$sucra$sucra_scores),
    outputs$analysis$heterogeneity$I2 %||% 0,
    outputs$workflow_info$total_time,
    sum(!sapply(outputs, is.null))
  )

  writeLines(html_content, output_file)

  return(output_file)
}

#' Print Method for Batch Results
#'
#' @export
print.batch_nma_results <- function(x, ...) {
  cat("Batch NMA Results\n")
  cat("=================\n\n")
  cat(sprintf("Total Datasets: %d\n", nrow(x$summary)))
  cat(sprintf("Successful: %d\n", sum(x$summary$Success)))
  cat(sprintf("Failed: %d\n", sum(!x$summary$Success)))
  cat(sprintf("Total Time: %.2f seconds\n\n", x$total_time))

  cat("Summary:\n")
  print(x$summary, row.names = FALSE)

  if (length(x$errors) > 0) {
    cat("\nErrors:\n")
    for (name in names(x$errors)) {
      cat(sprintf("  %s: %s\n", name, x$errors[[name]]))
    }
  }

  invisible(x)
}

#' Print Method for Automated Workflow Results
#'
#' @export
print.automated_workflow_results <- function(x, ...) {
  cat("Automated NMA Workflow Results\n")
  cat("==============================\n\n")
  cat(sprintf("Workflow completed in %.2f seconds\n", x$workflow_info$total_time))
  cat(sprintf("Output directory: %s\n\n", x$workflow_info$output_dir))

  cat("Generated Outputs:\n")
  cat(sprintf("  ✓ Analysis results\n"))
  if (!is.null(x$visualizations)) {
    cat(sprintf("  ✓ %d visualizations\n", length(x$visualizations)))
  }
  if (!is.null(x$manuscripts)) {
    cat(sprintf("  ✓ Manuscript sections (Methods + Results)\n"))
  }
  if (!is.null(x$reports)) {
    cat(sprintf("  ✓ %d reports\n", length(x$reports)))
  }

  cat(sprintf("\n✓ Summary dashboard: %s\n", x$dashboard))

  invisible(x)
}
