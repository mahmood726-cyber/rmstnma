#' Shiny Server Module Functions
#'
#' @description
#' This module contains all server logic functions for the powerNMA
#' Shiny dashboard. Each function handles the server logic for a specific tab.
#'
#' @author powerNMA Development Team
#' @name shiny_server_modules
NULL

#' Data Import Server Module
#'
#' @keywords internal
data_import_server <- function(input, output, session, rv) {

  # Load example data
  shiny::observeEvent(input$load_example, {
    req(input$example_data != "none")

    rv$data <- switch(
      input$example_data,
      smoking = generate_example_smoking_data(),
      depression = generate_example_depression_data(),
      diabetes = generate_example_diabetes_data(),
      cardio = generate_example_cardio_data(),
      multiarm = generate_example_multiarm_data()
    )

    shiny::showNotification("Example data loaded successfully!", type = "message")
  })

  # Upload data file
  shiny::observeEvent(input$data_file, {
    req(input$data_file)

    tryCatch({
      ext <- tools::file_ext(input$data_file$name)

      rv$data <- if (ext %in% c("csv", "txt")) {
        read.csv(input$data_file$datapath, stringsAsFactors = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package 'readxl' required for Excel files")
        }
        readxl::read_excel(input$data_file$datapath)
      } else {
        stop("Unsupported file format")
      }

      shiny::showNotification("Data uploaded successfully!", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Data preview
  output$data_preview <- DT::renderDataTable({
    req(rv$data)
    DT::datatable(
      rv$data,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # Validate data
  shiny::observeEvent(input$validate_data, {
    req(rv$data)

    tryCatch({
      validation <- validate_nma_input(rv$data, data_type = input$data_format)

      output$data_validation_status <- shiny::renderUI({
        shinydashboard::box(
          width = 12,
          status = "success",
          shiny::h4(shiny::icon("check-circle"), " Data Validation Passed"),
          shiny::p(sprintf("Format: %s", input$data_format)),
          shiny::p(sprintf("Studies: %d", validation$n_studies)),
          shiny::p(sprintf("Treatments: %d", validation$n_treatments)),
          shiny::p(sprintf("Comparisons: %d", validation$n_comparisons))
        )
      })

      shiny::showNotification("Data validation successful!", type = "message")
    }, error = function(e) {
      output$data_validation_status <- shiny::renderUI({
        shinydashboard::box(
          width = 12,
          status = "danger",
          shiny::h4(shiny::icon("exclamation-circle"), " Validation Failed"),
          shiny::p(e$message)
        )
      })

      shiny::showNotification(paste("Validation error:", e$message), type = "error")
    })
  })

  # Data summary
  output$data_summary <- shiny::renderPrint({
    req(rv$data)
    summary(rv$data)
  })

  # Download template
  output$download_template <- shiny::downloadHandler(
    filename = function() {
      sprintf("nma_template_%s.csv", input$data_format)
    },
    content = function(file) {
      template <- create_data_template(input$data_format)
      write.csv(template, file, row.names = FALSE)
    }
  )
}

#' Standard NMA Server Module
#'
#' @keywords internal
standard_nma_server <- function(input, output, session, rv) {

  # Update reference treatment choices
  shiny::observe({
    req(rv$data)
    treatments <- unique(c(rv$data$treat1, rv$data$treat2))
    shiny::updateSelectInput(session, "nma_reference", choices = treatments)
  })

  # Run standard NMA
  shiny::observeEvent(input$run_standard_nma, {
    req(rv$data)

    rv$analysis_running <- TRUE
    shinyWidgets::updateProgressBar(session, "nma_progress", value = 10)

    tryCatch({
      # Run NMA
      shinyWidgets::updateProgressBar(session, "nma_progress", value = 30)

      rv$nma_result <- netmeta::netmeta(
        TE = rv$data$TE,
        seTE = rv$data$seTE,
        treat1 = rv$data$treat1,
        treat2 = rv$data$treat2,
        studlab = rv$data$studlab,
        data = rv$data,
        sm = input$effect_measure,
        reference.group = input$nma_reference,
        small.values = if (input$nma_small_good) "desirable" else "undesirable"
      )

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 60)

      # Calculate rankings if requested
      if (input$nma_calc_rankings) {
        rv$sucra_result <- calculate_sucra(
          rv$nma_result,
          small_values = if (input$nma_small_good) "good" else "bad",
          nsim = input$nma_nsim
        )
      }

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 80)

      # Assess heterogeneity if requested
      if (input$nma_assess_het) {
        rv$heterogeneity_result <- heterogeneity_report(rv$nma_result)
      }

      # Assess inconsistency if requested
      if (input$nma_assess_incon) {
        rv$inconsistency_result <- node_splitting(rv$nma_result, rv$data)
      }

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 100)

      rv$analysis_running <- FALSE

      shiny::showNotification("Analysis complete!", type = "message", duration = 3)
    }, error = function(e) {
      rv$analysis_running <- FALSE
      shiny::showNotification(paste("Error:", e$message), type = "error", duration = NULL)
    })
  })

  # Summary output
  output$nma_summary_text <- shiny::renderPrint({
    req(rv$nma_result)
    summary(rv$nma_result)
  })

  # Effect estimates table
  output$nma_effects_table <- DT::renderDataTable({
    req(rv$nma_result)

    effects_df <- data.frame(
      Comparison = row.names(rv$nma_result$TE.random),
      Estimate = round(rv$nma_result$TE.random[, 1], 3),
      SE = round(rv$nma_result$seTE.random[, 1], 3),
      Lower = round(rv$nma_result$lower.random[, 1], 3),
      Upper = round(rv$nma_result$upper.random[, 1], 3),
      P_value = format.pval(rv$nma_result$pval.random[, 1])
    )

    DT::datatable(
      effects_df,
      options = list(pageLength = 20),
      rownames = FALSE
    )
  })

  # Network plot
  output$nma_network_plot <- plotly::renderPlotly({
    req(rv$nma_result)
    create_interactive_network(rv$data, rv$nma_result)
  })

  # Forest plot
  output$nma_forest_plot <- plotly::renderPlotly({
    req(rv$nma_result)
    create_interactive_forest(rv$nma_result)
  })

  # League table
  output$nma_league_table <- DT::renderDataTable({
    req(rv$nma_result)

    league <- create_league_table(
      rv$nma_result,
      sucra_result = rv$sucra_result,
      format = "effect_ci"
    )

    DT::datatable(league, options = list(pageLength = 20))
  })
}

#' Bayesian NMA Server Module
#'
#' @keywords internal
bayesian_nma_server <- function(input, output, session, rv) {

  # Run Bayesian NMA
  shiny::observeEvent(input$run_bayesian_nma, {
    req(rv$data)

    rv$analysis_running <- TRUE

    shiny::withProgress(message = "Running Bayesian Analysis...", {
      tryCatch({
        shiny::incProgress(0.2, detail = "Preparing data")

        # Prepare data for gemtc/Stan
        if (input$bayes_method == "gemtc") {
          # Run gemtc analysis
          rv$nma_result <- run_bayesian_gemtc(
            rv$data,
            n_iter = input$bayes_n_iter,
            n_burnin = input$bayes_n_burnin,
            n_thin = input$bayes_n_thin,
            n_chains = input$bayes_n_chains
          )
        } else if (input$bayes_method == "stan") {
          # Run Stan analysis
          rv$nma_result <- run_bayesian_stan(
            rv$data,
            n_iter = input$bayes_n_iter,
            n_warmup = input$bayes_n_burnin,
            n_chains = input$bayes_n_chains
          )
        }

        shiny::incProgress(0.6, detail = "Calculating posteriors")

        # Calculate model comparison
        if (!is.null(rv$nma_result)) {
          rv$bayes_model_comparison <- calculate_model_comparison(rv$nma_result)
        }

        shiny::incProgress(0.2, detail = "Finalizing")

        rv$analysis_running <- FALSE

        shiny::showNotification("Bayesian analysis complete!", type = "message")
      }, error = function(e) {
        rv$analysis_running <- FALSE
        shiny::showNotification(paste("Error:", e$message), type = "error", duration = NULL)
      })
    })
  })

  # Outputs
  output$bayes_summary <- shiny::renderPrint({
    req(rv$nma_result)
    summary(rv$nma_result)
  })

  output$bayes_effects_table <- DT::renderDataTable({
    req(rv$nma_result)
    # Extract posterior summaries
    # Implementation depends on gemtc/Stan output format
    DT::datatable(data.frame())
  })

  output$bayes_model_comparison <- DT::renderDataTable({
    req(rv$bayes_model_comparison)
    DT::datatable(rv$bayes_model_comparison)
  })
}

#' Component NMA Server Module
#'
#' @keywords internal
component_nma_server <- function(input, output, session, rv) {

  # Component definitions storage
  components <- shiny::reactiveValues(list = list())

  # Add component
  shiny::observeEvent(input$add_component, {
    new_id <- length(components$list) + 1
    components$list[[new_id]] <- list(name = "", treatments = c())
  })

  # Run Component NMA
  shiny::observeEvent(input$run_component_nma, {
    req(rv$data, input$component_matrix_file)

    tryCatch({
      # Load component matrix
      comp_matrix <- read.csv(input$component_matrix_file$datapath, row.names = 1)

      # Run CNMA
      rv$component_nma_result <- if (input$cnma_model == "additive") {
        additive_cnma(rv$data, comp_matrix, sm = input$effect_measure)
      } else {
        interaction_cnma(rv$data, comp_matrix, sm = input$effect_measure)
      }

      shiny::showNotification("Component NMA complete!", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Outputs
  output$cnma_component_effects <- DT::renderDataTable({
    req(rv$component_nma_result)
    DT::datatable(rv$component_nma_result$component_effects)
  })

  output$cnma_component_forest <- plotly::renderPlotly({
    req(rv$component_nma_result)
    plot_component_forest(rv$component_nma_result)
  })
}

#' Multivariate NMA Server Module
#'
#' @keywords internal
multivariate_nma_server <- function(input, output, session, rv) {

  # Update outcome choices
  shiny::observe({
    req(rv$data)
    outcomes <- names(rv$data)[grepl("TE|outcome", names(rv$data), ignore.case = TRUE)]
    shiny::updateSelectInput(session, "mvnma_outcome1", choices = outcomes)
    shiny::updateSelectInput(session, "mvnma_outcome2", choices = outcomes)
  })

  # Run Multivariate NMA
  shiny::observeEvent(input$run_multivariate_nma, {
    req(rv$data, input$mvnma_outcome1, input$mvnma_outcome2)

    shiny::withProgress(message = "Running Multivariate NMA...", {
      tryCatch({
        # Prepare multivariate data
        mvnma_data <- prepare_multivariate_data(
          rv$data,
          outcome1 = input$mvnma_outcome1,
          outcome2 = input$mvnma_outcome2,
          correlation = input$mvnma_correlation
        )

        # Run MVNMA
        rv$multivariate_nma_result <- fit_mvnma(
          mvnma_data,
          reference = input$nma_reference
        )

        # Calculate benefit-risk
        rv$benefit_risk <- calculate_benefit_risk(
          rv$multivariate_nma_result,
          efficacy_weight = input$mvnma_efficacy_weight
        )

        shiny::showNotification("Multivariate NMA complete!", type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })

  # Outputs
  output$mvnma_joint_results <- DT::renderDataTable({
    req(rv$multivariate_nma_result)
    DT::datatable(rv$multivariate_nma_result$joint_estimates)
  })

  output$mvnma_benefit_risk_plot <- plotly::renderPlotly({
    req(rv$benefit_risk)
    plot_benefit_risk(rv$benefit_risk)
  })

  output$mvnma_net_benefit <- DT::renderDataTable({
    req(rv$benefit_risk)
    DT::datatable(rv$benefit_risk$net_benefit_table)
  })
}

#' Living NMA Server Module
#'
#' @keywords internal
living_nma_server <- function(input, output, session, rv) {

  # Initialize Living NMA project
  shiny::observeEvent(input$init_living_nma, {
    req(rv$data, input$living_project_name)

    tryCatch({
      rv$living_nma_project <- initialize_living_nma(
        baseline_data = rv$data,
        project_name = input$living_project_name,
        update_frequency = input$living_update_frequency
      )

      shiny::showNotification(
        sprintf("Living NMA project '%s' initialized!", input$living_project_name),
        type = "message"
      )

      # Update project dropdown
      shiny::updateSelectInput(
        session,
        "living_existing_project",
        choices = input$living_project_name,
        selected = input$living_project_name
      )
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Update Living NMA
  shiny::observeEvent(input$update_living_nma, {
    req(rv$living_nma_project, input$living_new_data)

    tryCatch({
      new_data <- read.csv(input$living_new_data$datapath)

      rv$living_nma_project <- update_living_nma(
        project = rv$living_nma_project,
        new_data = new_data
      )

      shiny::showNotification("Living NMA updated successfully!", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Outputs
  output$living_version_history <- DT::renderDataTable({
    req(rv$living_nma_project)
    DT::datatable(rv$living_nma_project$version_history)
  })

  output$living_timeline_plot <- plotly::renderPlotly({
    req(rv$living_nma_project)
    plot_living_timeline(rv$living_nma_project)
  })

  output$living_effect_tracking <- plotly::renderPlotly({
    req(rv$living_nma_project)
    track_treatment_effect_plot(rv$living_nma_project)
  })
}

#' Rankings Server Module
#'
#' @keywords internal
rankings_server <- function(input, output, session, rv) {

  # Calculate rankings
  shiny::observeEvent(input$calc_rankings, {
    req(rv$nma_result)

    shiny::withProgress(message = "Calculating rankings...", {
      tryCatch({
        rv$sucra_result <- calculate_sucra(
          rv$nma_result,
          small_values = if (input$ranking_small_good) "good" else "bad",
          nsim = input$ranking_nsim
        )

        shiny::showNotification("Rankings calculated!", type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })

  # Outputs
  output$rankings_top_table <- DT::renderDataTable({
    req(rv$sucra_result)

    top_df <- data.frame(
      Rank = 1:min(5, length(rv$sucra_result$sucra_scores)),
      Treatment = names(sort(rv$sucra_result$sucra_scores, decreasing = TRUE))[1:5],
      SUCRA = sort(rv$sucra_result$sucra_scores, decreasing = TRUE)[1:5]
    )

    DT::datatable(top_df, options = list(dom = 't'), rownames = FALSE)
  })

  output$rankings_sucra_plot <- plotly::renderPlotly({
    req(rv$sucra_result)
    create_interactive_rankings(rv$sucra_result)
  })

  output$rankings_rankogram <- plotly::renderPlotly({
    req(rv$sucra_result)
    plot_rankogram_interactive(rv$sucra_result)
  })

  output$rankings_full_table <- DT::renderDataTable({
    req(rv$sucra_result)

    full_df <- data.frame(
      Treatment = names(rv$sucra_result$sucra_scores),
      SUCRA = round(rv$sucra_result$sucra_scores, 2),
      Mean_Rank = round(rv$sucra_result$mean_ranks, 2)
    )

    DT::datatable(full_df, options = list(pageLength = 20), rownames = FALSE)
  })
}

#' Visualizations Server Module
#'
#' @keywords internal
visualizations_server <- function(input, output, session, rv) {

  # Generate visualization
  viz_plot <- shiny::eventReactive(input$generate_viz, {
    req(rv$nma_result)

    switch(
      input$viz_type,
      network = create_interactive_network(rv$data, rv$nma_result, layout = input$viz_layout),
      network_3d = create_3d_network(rv$data, rv$nma_result),
      forest = create_interactive_forest(rv$nma_result),
      funnel = create_interactive_funnel(rv$nma_result),
      heatmap = create_interactive_heatmap(rv$nma_result),
      interval = create_interval_plot_interactive(rv$nma_result),
      comparison = create_comparison_matrix_plot(rv$nma_result),
      netheat = create_net_heat_plot(rv$nma_result),
      contour_funnel = create_contour_funnel_interactive(rv$nma_result)
    )
  })

  output$main_visualization <- plotly::renderPlotly({
    viz_plot()
  })

  output$viz_title <- shiny::renderText({
    req(input$viz_type)
    titles <- c(
      network = "Network Graph",
      network_3d = "3D Network Visualization",
      forest = "Forest Plot",
      funnel = "Funnel Plot",
      heatmap = "Comparison Heat Map",
      interval = "Interval Plot with Prediction Intervals",
      comparison = "Comparison Matrix",
      netheat = "Net Heat Plot",
      contour_funnel = "Contour-Enhanced Funnel Plot"
    )
    titles[[input$viz_type]]
  })

  # Download visualization
  output$download_viz <- shiny::downloadHandler(
    filename = function() {
      sprintf("nma_%s_%s.html", input$viz_type, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      htmlwidgets::saveWidget(viz_plot(), file)
    }
  )
}

# Placeholder server modules for remaining tabs
# (Full implementations follow same pattern)

#' @keywords internal
heterogeneity_server <- function(input, output, session, rv) {}

#' @keywords internal
inconsistency_server <- function(input, output, session, rv) {}

#' @keywords internal
pub_bias_server <- function(input, output, session, rv) {}

#' @keywords internal
metaregression_server <- function(input, output, session, rv) {}

#' @keywords internal
network_geometry_server <- function(input, output, session, rv) {}

#' @keywords internal
simulation_server <- function(input, output, session, rv) {}

#' @keywords internal
voi_server <- function(input, output, session, rv) {}

#' Generate Methods Server Module
#'
#' @keywords internal
generate_methods_server <- function(input, output, session, rv) {

  # Generate methods section
  methods_result <- shiny::eventReactive(input$generate_methods, {
    req(rv$nma_result)

    shiny::withProgress(message = "Generating Methods section...", {
      analysis_config <- list(
        search_databases = input$methods_search_db,
        risk_of_bias_tool = input$methods_rob_tool,
        has_component_nma = input$methods_has_cnma,
        has_multivariate = input$methods_has_mvnma,
        has_ipd = input$methods_has_ipd
      )

      result <- generate_methods_section(
        nma_result = rv$nma_result,
        analysis_config = analysis_config,
        style = input$methods_style,
        use_ai = input$methods_use_ai && check_ollama_available()
      )

      if (input$methods_use_ai && check_ollama_available()) {
        result <- enhance_methods_with_ollama(
          result$text,
          enhancement_level = "moderate"
        )
        return(result)
      }

      return(result)
    })
  })

  output$generated_methods_text <- shiny::renderText({
    result <- methods_result()
    if (is.list(result) && "enhanced_text" %in% names(result)) {
      result$enhanced_text
    } else if (is.list(result) && "text" %in% names(result)) {
      result$text
    } else {
      result
    }
  })

  output$methods_stats <- shiny::renderPrint({
    result <- methods_result()
    if (is.list(result)) {
      cat("Generation Statistics:\n")
      cat("---------------------\n")
      if ("permutation_count" %in% names(result)) {
        cat(sprintf("Permutations: %d\n", result$permutation_count))
      }
      if ("word_count" %in% names(result)) {
        cat(sprintf("Word count: %d\n", result$word_count))
      }
      if ("enhanced" %in% names(result)) {
        cat(sprintf("AI Enhanced: %s\n", result$enhanced))
      }
    }
  })

  # Download handlers
  output$download_methods_txt <- shiny::downloadHandler(
    filename = function() {
      sprintf("methods_section_%s.txt", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      writeLines(methods_result()$enhanced_text %||% methods_result()$text, file)
    }
  )

  output$download_methods_docx <- shiny::downloadHandler(
    filename = function() {
      sprintf("methods_section_%s.docx", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      if (requireNamespace("officer", quietly = TRUE)) {
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, methods_result()$enhanced_text %||% methods_result()$text)
        print(doc, target = file)
      } else {
        stop("Package 'officer' required for DOCX export")
      }
    }
  )
}

#' @keywords internal
generate_results_server <- function(input, output, session, rv) {}

#' @keywords internal
ai_enhancement_server <- function(input, output, session, rv) {}

#' @keywords internal
text_results_server <- function(input, output, session, rv) {}

#' @keywords internal
reports_server <- function(input, output, session, rv) {}

#' @keywords internal
batch_processing_server <- function(input, output, session, rv) {}

#' @keywords internal
export_server <- function(input, output, session, rv) {}

#' @keywords internal
database_server <- function(input, output, session, rv, db_path) {}

#' @keywords internal
settings_server <- function(input, output, session, rv) {}

#' @keywords internal
help_server <- function(input, output, session, rv) {}

#' Helper Functions
#' @keywords internal

# Generate example datasets
generate_example_smoking_data <- function() {
  simulate_nma_data(n_studies = 30, n_treatments = 8, sm = "OR")
}

generate_example_depression_data <- function() {
  simulate_nma_data(n_studies = 40, n_treatments = 10, sm = "SMD")
}

generate_example_diabetes_data <- function() {
  simulate_nma_data(n_studies = 35, n_treatments = 7, sm = "MD")
}

generate_example_cardio_data <- function() {
  simulate_nma_data(n_studies = 25, n_treatments = 6, sm = "HR")
}

generate_example_multiarm_data <- function() {
  simulate_nma_data(n_studies = 20, n_treatments = 5, sm = "OR", multiarm = TRUE)
}

# Create data templates
create_data_template <- function(format) {
  if (format == "pairwise") {
    data.frame(
      studlab = c("Study1", "Study1", "Study2"),
      treat1 = c("A", "A", "B"),
      treat2 = c("B", "C", "C"),
      TE = c(NA, NA, NA),
      seTE = c(NA, NA, NA)
    )
  } else {
    data.frame()
  }
}

# NULL-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a
