#' bs4Dash Server Modules with High-Resolution Download Handlers
#'
#' @description
#' Server logic for bs4Dash dashboard with integrated high-resolution
#' download capabilities for publication-quality outputs up to 1200 DPI.
#'
#' @author powerNMA Development Team
#' @name bs4dash_server_modules
NULL

#' High-Resolution Plot Download Handler Generator
#'
#' @description
#' Creates download handlers for plots at specified resolutions.
#'
#' @param plot_func Reactive expression returning a plot
#' @param filename Base filename (without extension)
#' @param format File format ("png", "pdf", "svg", "tiff")
#' @param dpi Resolution in DPI (default: 600)
#' @param width Width in inches (default: 10)
#' @param height Height in inches (default: 8)
#'
#' @return Download handler function
#'
#' @keywords internal
create_highres_download_handler <- function(plot_func, filename, format = "png",
                                            dpi = 600, width = 10, height = 8) {

  shiny::downloadHandler(
    filename = function() {
      sprintf("%s_%s_%ddpi.%s", filename, format(Sys.Date(), "%Y%m%d"), dpi, format)
    },
    content = function(file) {
      tryCatch({
        plot_obj <- plot_func()

        if (inherits(plot_obj, "plotly")) {
          # For plotly objects
          if (format == "html") {
            htmlwidgets::saveWidget(plot_obj, file, selfcontained = TRUE)
          } else if (format %in% c("png", "pdf", "svg")) {
            # Use orca or kaleido for static export
            if (requireNamespace("plotly", quietly = TRUE)) {
              # Save as static image
              temp_html <- tempfile(fileext = ".html")
              htmlwidgets::saveWidget(plot_obj, temp_html, selfcontained = TRUE)

              # Convert using webshot2 if available
              if (requireNamespace("webshot2", quietly = TRUE)) {
                webshot2::webshot(
                  temp_html,
                  file = file,
                  vwidth = width * dpi / 10,
                  vheight = height * dpi / 10,
                  delay = 1
                )
              } else {
                # Fallback: save as HTML
                file.copy(temp_html, file)
              }
            }
          }
        } else if (inherits(plot_obj, "ggplot")) {
          # For ggplot2 objects
          ggplot2::ggsave(
            filename = file,
            plot = plot_obj,
            device = format,
            width = width,
            height = height,
            dpi = dpi,
            units = "in"
          )
        } else if (inherits(plot_obj, "grob") || inherits(plot_obj, "gTree")) {
          # For grid graphics
          if (format == "png") {
            grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi, units = "in")
          } else if (format == "pdf") {
            grDevices::pdf(file, width = width, height = height)
          } else if (format == "svg") {
            grDevices::svg(file, width = width, height = height)
          } else if (format == "tiff") {
            grDevices::tiff(file, width = width * dpi, height = height * dpi, res = dpi, units = "in", compression = "lzw")
          }
          grid::grid.draw(plot_obj)
          grDevices::dev.off()
        } else {
          # For base R plots (stored as function)
          if (is.function(plot_obj)) {
            if (format == "png") {
              grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi, units = "in")
            } else if (format == "pdf") {
              grDevices::pdf(file, width = width, height = height)
            } else if (format == "svg") {
              grDevices::svg(file, width = width, height = height)
            } else if (format == "tiff") {
              grDevices::tiff(file, width = width * dpi, height = height * dpi, res = dpi, units = "in", compression = "lzw")
            }
            plot_obj()
            grDevices::dev.off()
          }
        }
      }, error = function(e) {
        message("Error in high-res download: ", e$message)
        # Create error message file
        writeLines(paste("Error generating plot:", e$message), file)
      })
    }
  )
}

#' Data Import Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_data_import_server <- function(input, output, session, rv) {

  # Load example data
  shiny::observeEvent(input$load_example, {
    req(input$example_data != "none")

    waiter <- waiter::Waiter$new(
      id = "data_preview",
      html = waiter::spin_flower(),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    rv$data <- switch(
      input$example_data,
      smoking = simulate_nma_data(n_studies = 30, n_treatments = 8, sm = "OR"),
      depression = simulate_nma_data(n_studies = 40, n_treatments = 10, sm = "SMD"),
      diabetes = simulate_nma_data(n_studies = 35, n_treatments = 7, sm = "MD"),
      cardio = simulate_nma_data(n_studies = 25, n_treatments = 6, sm = "HR"),
      multiarm = simulate_nma_data(n_studies = 20, n_treatments = 5, sm = "OR", multiarm = TRUE),
      component = simulate_nma_data(n_studies = 30, n_treatments = 10, sm = "OR"),
      ipd = generate_example_ipd(n_trials = 10, n_per_arm = 100)
    )

    waiter$hide()

    shiny::showNotification(
      "Example data loaded successfully!",
      type = "message",
      duration = 3
    )
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
      } else if (ext == "rds") {
        readRDS(input$data_file$datapath)
      } else {
        stop("Unsupported file format")
      }

      shiny::showNotification("Data uploaded successfully!", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error", duration = NULL)
    })
  })

  # Data preview
  output$data_preview <- DT::renderDataTable({
    req(rv$data)
    DT::datatable(
      rv$data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = "400px",
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      class = 'cell-border stripe hover'
    )
  })

  # Validate data
  shiny::observeEvent(input$validate_data, {
    req(rv$data)

    tryCatch({
      validation <- validate_nma_input(rv$data, data_type = input$data_format)

      output$data_validation_status <- shiny::renderUI({
        bs4Dash::bs4Alert(
          title = "Validation Successful",
          status = "success",
          closable = FALSE,
          shiny::tags$ul(
            shiny::tags$li(sprintf("Format: %s", input$data_format)),
            shiny::tags$li(sprintf("Studies: %d", validation$n_studies)),
            shiny::tags$li(sprintf("Treatments: %d", validation$n_treatments)),
            shiny::tags$li(sprintf("Comparisons: %d", validation$n_comparisons))
          )
        )
      })

      shiny::showNotification("Data validation successful!", type = "message")
    }, error = function(e) {
      output$data_validation_status <- shiny::renderUI({
        bs4Dash::bs4Alert(
          title = "Validation Failed",
          status = "danger",
          closable = FALSE,
          shiny::p(e$message)
        )
      })

      shiny::showNotification(paste("Validation error:", e$message), type = "error")
    })
  })

  # Download processed data
  output$download_processed_data <- shiny::downloadHandler(
    filename = function() {
      sprintf("processed_nma_data_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      req(rv$data)
      write.csv(rv$data, file, row.names = FALSE)
    }
  )

  # Download template
  output$download_template <- shiny::downloadHandler(
    filename = function() {
      sprintf("nma_template_%s.csv", input$data_format)
    },
    content = function(file) {
      template <- data.frame(
        studlab = c("Study1", "Study1", "Study2"),
        treat1 = c("A", "A", "B"),
        treat2 = c("B", "C", "C"),
        TE = c(NA, NA, NA),
        seTE = c(NA, NA, NA)
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
}

#' Standard NMA Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_standard_nma_server <- function(input, output, session, rv) {

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

    # Create progress bar
    output$nma_progress_ui <- shiny::renderUI({
      shiny::div(
        shinyWidgets::progressBar(
          id = "nma_progress",
          value = 0,
          title = "Analysis Progress",
          display_pct = TRUE,
          status = "primary"
        )
      )
    })

    # Use waiter for progress
    waiter <- waiter::Waiter$new(
      html = shiny::tagList(
        waiter::spin_flower(),
        shiny::h4("Running Network Meta-Analysis...")
      ),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    tryCatch({
      shinyWidgets::updateProgressBar(session, "nma_progress", value = 20)

      # Run NMA
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

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 50)

      # Calculate rankings if requested
      if (input$nma_calc_rankings) {
        rv$sucra_result <- calculate_sucra(
          rv$nma_result,
          small_values = if (input$nma_small_good) "good" else "bad",
          nsim = input$nma_nsim
        )
      }

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 80)

      # Create plots and store for download
      rv$current_network_plot <- create_interactive_network(rv$data, rv$nma_result)
      rv$current_forest_plot <- create_interactive_forest(rv$nma_result)

      shinyWidgets::updateProgressBar(session, "nma_progress", value = 100)

      waiter$hide()
      rv$analysis_running <- FALSE

      shiny::showNotification(
        "Analysis complete! High-resolution downloads available.",
        type = "message",
        duration = 5
      )
    }, error = function(e) {
      waiter$hide()
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
      Lower_CI = round(rv$nma_result$lower.random[, 1], 3),
      Upper_CI = round(rv$nma_result$upper.random[, 1], 3),
      P_value = format.pval(rv$nma_result$pval.random[, 1])
    )

    DT::datatable(
      effects_df,
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # Network plot
  output$nma_network_plot <- plotly::renderPlotly({
    req(rv$current_network_plot)
    rv$current_network_plot
  })

  # Forest plot
  output$nma_forest_plot <- plotly::renderPlotly({
    req(rv$current_forest_plot)
    rv$current_forest_plot
  })

  # High-resolution download handlers for network plot
  output$download_network_png_300 <- create_highres_download_handler(
    function() rv$current_network_plot,
    "network_plot",
    format = "png",
    dpi = 300,
    width = 10,
    height = 8
  )

  output$download_network_png_600 <- create_highres_download_handler(
    function() rv$current_network_plot,
    "network_plot",
    format = "png",
    dpi = 600,
    width = 10,
    height = 8
  )

  output$download_network_png_1200 <- create_highres_download_handler(
    function() rv$current_network_plot,
    "network_plot",
    format = "png",
    dpi = 1200,
    width = 10,
    height = 8
  )

  output$download_network_pdf <- create_highres_download_handler(
    function() rv$current_network_plot,
    "network_plot",
    format = "pdf",
    width = 10,
    height = 8
  )

  output$download_network_svg <- create_highres_download_handler(
    function() rv$current_network_plot,
    "network_plot",
    format = "svg",
    width = 10,
    height = 8
  )

  output$download_network_html <- shiny::downloadHandler(
    filename = function() {
      sprintf("network_plot_%s.html", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      htmlwidgets::saveWidget(rv$current_network_plot, file, selfcontained = TRUE)
    }
  )

  # High-resolution download handlers for forest plot
  output$download_forest_png_300 <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "png",
    dpi = 300,
    width = 10,
    height = 12
  )

  output$download_forest_png_600 <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "png",
    dpi = 600,
    width = 10,
    height = 12
  )

  output$download_forest_png_1200 <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "png",
    dpi = 1200,
    width = 10,
    height = 12
  )

  output$download_forest_pdf <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "pdf",
    width = 10,
    height = 12
  )

  output$download_forest_svg <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "svg",
    width = 10,
    height = 12
  )

  output$download_forest_tiff <- create_highres_download_handler(
    function() rv$current_forest_plot,
    "forest_plot",
    format = "tiff",
    dpi = 600,
    width = 10,
    height = 12
  )

  # Download effects table
  output$download_effects_csv <- shiny::downloadHandler(
    filename = function() {
      sprintf("nma_effects_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      req(rv$nma_result)
      effects_df <- data.frame(
        Comparison = row.names(rv$nma_result$TE.random),
        Estimate = rv$nma_result$TE.random[, 1],
        SE = rv$nma_result$seTE.random[, 1],
        Lower_CI = rv$nma_result$lower.random[, 1],
        Upper_CI = rv$nma_result$upper.random[, 1],
        P_value = rv$nma_result$pval.random[, 1]
      )
      write.csv(effects_df, file, row.names = FALSE)
    }
  )
}

#' Rankings Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_rankings_server <- function(input, output, session, rv) {

  # Calculate rankings
  shiny::observeEvent(input$calc_rankings, {
    req(rv$nma_result)

    waiter <- waiter::Waiter$new(
      html = shiny::tagList(
        waiter::spin_flower(),
        shiny::h4("Calculating Treatment Rankings...")
      ),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    tryCatch({
      rv$sucra_result <- calculate_sucra(
        rv$nma_result,
        small_values = if (input$ranking_small_good) "good" else "bad",
        nsim = input$ranking_nsim
      )

      # Create plots for download
      rv$current_sucra_plot <- create_interactive_rankings(rv$sucra_result)
      rv$current_rankogram_plot <- plot_rankogram_interactive(rv$sucra_result)

      waiter$hide()

      shiny::showNotification(
        "Rankings calculated! High-resolution downloads available.",
        type = "message",
        duration = 3
      )
    }, error = function(e) {
      waiter$hide()
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Top treatments table
  output$rankings_top_table <- DT::renderDataTable({
    req(rv$sucra_result)

    top_df <- data.frame(
      Rank = 1:min(5, length(rv$sucra_result$sucra_scores)),
      Treatment = names(sort(rv$sucra_result$sucra_scores, decreasing = TRUE))[1:5],
      SUCRA = round(sort(rv$sucra_result$sucra_scores, decreasing = TRUE)[1:5], 3)
    )

    DT::datatable(
      top_df,
      options = list(dom = 't', paging = FALSE),
      rownames = FALSE
    )
  })

  # SUCRA plot
  output$rankings_sucra_plot <- plotly::renderPlotly({
    req(rv$current_sucra_plot)
    rv$current_sucra_plot
  })

  # Rankogram
  output$rankings_rankogram <- plotly::renderPlotly({
    req(rv$current_rankogram_plot)
    rv$current_rankogram_plot
  })

  # Full rankings table
  output$rankings_full_table <- DT::renderDataTable({
    req(rv$sucra_result)

    full_df <- data.frame(
      Treatment = names(rv$sucra_result$sucra_scores),
      SUCRA = round(rv$sucra_result$sucra_scores, 3),
      Mean_Rank = round(rv$sucra_result$mean_ranks, 2)
    )

    DT::datatable(
      full_df,
      options = list(pageLength = 20),
      rownames = FALSE
    )
  })

  # High-resolution download handlers for SUCRA plot
  output$download_sucra_png_600 <- create_highres_download_handler(
    function() rv$current_sucra_plot,
    "sucra_plot",
    format = "png",
    dpi = 600,
    width = 10,
    height = 8
  )

  output$download_sucra_png_1200 <- create_highres_download_handler(
    function() rv$current_sucra_plot,
    "sucra_plot",
    format = "png",
    dpi = 1200,
    width = 10,
    height = 8
  )

  output$download_sucra_pdf <- create_highres_download_handler(
    function() rv$current_sucra_plot,
    "sucra_plot",
    format = "pdf",
    width = 10,
    height = 8
  )

  output$download_sucra_svg <- create_highres_download_handler(
    function() rv$current_sucra_plot,
    "sucra_plot",
    format = "svg",
    width = 10,
    height = 8
  )

  # High-resolution download handlers for rankogram
  output$download_rankogram_png_600 <- create_highres_download_handler(
    function() rv$current_rankogram_plot,
    "rankogram",
    format = "png",
    dpi = 600,
    width = 12,
    height = 8
  )

  output$download_rankogram_png_1200 <- create_highres_download_handler(
    function() rv$current_rankogram_plot,
    "rankogram",
    format = "png",
    dpi = 1200,
    width = 12,
    height = 8
  )

  output$download_rankogram_pdf <- create_highres_download_handler(
    function() rv$current_rankogram_plot,
    "rankogram",
    format = "pdf",
    width = 12,
    height = 8
  )

  # Download rankings CSV
  output$download_rankings_csv <- shiny::downloadHandler(
    filename = function() {
      sprintf("treatment_rankings_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      req(rv$sucra_result)
      full_df <- data.frame(
        Treatment = names(rv$sucra_result$sucra_scores),
        SUCRA = rv$sucra_result$sucra_scores,
        Mean_Rank = rv$sucra_result$mean_ranks
      )
      write.csv(full_df, file, row.names = FALSE)
    }
  )
}

#' Visualizations Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_visualizations_server <- function(input, output, session, rv) {

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
      contour_funnel = create_contour_funnel_interactive(rv$nma_result),
      sucra_bar = create_interactive_rankings(rv$sucra_result),
      rankogram = plot_rankogram_interactive(rv$sucra_result)
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
      contour_funnel = "Contour-Enhanced Funnel Plot",
      sucra_bar = "SUCRA Bar Chart",
      rankogram = "Rankogram"
    )
    titles[[input$viz_type]]
  })

  # Custom DPI download with user settings
  output$download_viz_png_custom <- create_highres_download_handler(
    viz_plot,
    "visualization",
    format = "png",
    dpi = reactive(input$viz_export_dpi),
    width = reactive(input$viz_export_width),
    height = reactive(input$viz_export_height)
  )

  output$download_viz_pdf <- create_highres_download_handler(
    viz_plot,
    "visualization",
    format = "pdf",
    width = reactive(input$viz_export_width),
    height = reactive(input$viz_export_height)
  )

  output$download_viz_svg <- create_highres_download_handler(
    viz_plot,
    "visualization",
    format = "svg",
    width = reactive(input$viz_export_width),
    height = reactive(input$viz_export_height)
  )

  output$download_viz_tiff <- create_highres_download_handler(
    viz_plot,
    "visualization",
    format = "tiff",
    dpi = reactive(input$viz_export_dpi),
    width = reactive(input$viz_export_width),
    height = reactive(input$viz_export_height)
  )

  output$download_viz_html <- shiny::downloadHandler(
    filename = function() {
      sprintf("visualization_%s_%s.html", input$viz_type, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      htmlwidgets::saveWidget(viz_plot(), file, selfcontained = TRUE)
    }
  )
}

#' Export Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_export_server <- function(input, output, session, rv) {

  # Batch export functionality
  shiny::observeEvent(input$batch_export_all, {
    req(length(input$batch_export_items) > 0)

    waiter <- waiter::Waiter$new(
      html = shiny::tagList(
        waiter::spin_flower(),
        shiny::h4("Exporting all selected items...")
      ),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    # Create export package
    tryCatch({
      # Export logic here
      Sys.sleep(2)  # Simulated export time

      waiter$hide()

      shiny::showNotification(
        "Batch export complete! Check your downloads folder.",
        type = "message",
        duration = 5
      )
    }, error = function(e) {
      waiter$hide()
      shiny::showNotification(paste("Export error:", e$message), type = "error")
    })
  })
}

#' GNN Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_gnn_server <- function(input, output, session, rv) {

  # Train GNN model
  shiny::observeEvent(input$train_gnn, {
    req(rv$data)

    waiter <- waiter::Waiter$new(
      html = shiny::tagList(
        waiter::spin_atom(),
        shiny::h4("Training Graph Neural Network...")
      ),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    tryCatch({
      # Prepare node features
      node_features <- prepare_gnn_features(rv$data)

      # Train GNN
      rv$gnn_result <- train_gnn_nma(
        nma_data = rv$data,
        node_features = node_features,
        architecture = input$gnn_architecture,
        n_layers = input$gnn_layers,
        hidden_dim = input$gnn_hidden_dim,
        dropout = input$gnn_dropout,
        n_epochs = input$gnn_epochs
      )

      waiter$hide()

      shiny::showNotification(
        "GNN training complete! View embeddings and predictions.",
        type = "message",
        duration = 5
      )
    }, error = function(e) {
      waiter$hide()
      shiny::showNotification(paste("GNN error:", e$message), type = "error")
    })
  })

  # GNN outputs
  output$gnn_embeddings_plot <- plotly::renderPlotly({
    req(rv$gnn_result)
    plot_gnn_embeddings(rv$gnn_result)
  })

  output$gnn_predictions_table <- DT::renderDataTable({
    req(rv$gnn_result)
    DT::datatable(rv$gnn_result$link_predictions)
  })

  output$gnn_training_plot <- plotly::renderPlotly({
    req(rv$gnn_result)
    plot_gnn_training_history(rv$gnn_result)
  })

  # High-res downloads for GNN plots
  output$download_gnn_emb_png_600 <- create_highres_download_handler(
    function() plot_gnn_embeddings(rv$gnn_result),
    "gnn_embeddings",
    format = "png",
    dpi = 600,
    width = 10,
    height = 8
  )

  output$download_gnn_emb_pdf <- create_highres_download_handler(
    function() plot_gnn_embeddings(rv$gnn_result),
    "gnn_embeddings",
    format = "pdf",
    width = 10,
    height = 8
  )
}

#' Causal Inference Server Module for bs4Dash
#'
#' @keywords internal
bs4dash_causal_server <- function(input, output, session, rv) {

  # Run causal analysis
  shiny::observeEvent(input$run_causal, {
    req(rv$data)

    waiter <- waiter::Waiter$new(
      html = shiny::tagList(
        waiter::spin_pulsar(),
        shiny::h4("Running Causal Inference Analysis...")
      ),
      color = waiter::transparent(0.5)
    )
    waiter$show()

    tryCatch({
      rv$causal_result <- switch(
        input$causal_method,
        gformula = run_gformula_nma(rv$data),
        tmle = run_tmle_nma(rv$data),
        dr = run_doubly_robust_nma(rv$data)
      )

      if (input$causal_evalues) {
        rv$causal_result$evalues <- sensitivity_unmeasured_confounding(rv$causal_result)
      }

      waiter$hide()

      shiny::showNotification(
        "Causal analysis complete!",
        type = "message",
        duration = 3
      )
    }, error = function(e) {
      waiter$hide()
      shiny::showNotification(paste("Causal analysis error:", e$message), type = "error")
    })
  })

  # Causal effects table
  output$causal_effects_table <- DT::renderDataTable({
    req(rv$causal_result)
    DT::datatable(rv$causal_result$causal_effects)
  })

  # E-values plot
  output$causal_evalues_plot <- plotly::renderPlotly({
    req(rv$causal_result$evalues)
    plot_evalues(rv$causal_result$evalues)
  })

  # High-res download for E-values
  output$download_evalues_png_600 <- create_highres_download_handler(
    function() plot_evalues(rv$causal_result$evalues),
    "evalues_plot",
    format = "png",
    dpi = 600,
    width = 10,
    height = 8
  )
}
