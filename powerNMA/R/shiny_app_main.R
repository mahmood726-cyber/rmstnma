#' Production-Ready Shiny Dashboard for Network Meta-Analysis
#'
#' @description
#' Comprehensive, production-ready Shiny application for conducting network
#' meta-analysis with advanced features including:
#' \itemize{
#'   \item Interactive data import wizard
#'   \item Real-time analysis with progress tracking
#'   \item Advanced visualization gallery
#'   \item Batch processing capabilities
#'   \item Automated report generation
#'   \item Session management and result caching
#'   \item Export capabilities in multiple formats
#' }
#'
#' @details
#' This module provides a complete Shiny dashboard built with shinydashboard
#' and modern UI/UX patterns. It integrates ALL powerNMA functionality
#' (Phases 5-9) into an intuitive interface.
#'
#' @references
#' RStudio Shiny: https://shiny.rstudio.com/
#' shinydashboard: https://rstudio.github.io/shinydashboard/
#'
#' @author powerNMA Development Team
#' @name shiny_app_main
NULL

#' Launch powerNMA Shiny Dashboard
#'
#' @description
#' Main function to launch the comprehensive powerNMA Shiny dashboard.
#'
#' @param launch_browser Logical whether to launch in browser (default: TRUE)
#' @param port Port number for the app (default: auto)
#' @param host Host address (default: "127.0.0.1")
#' @param theme Dashboard theme: "default", "dark", "light" (default: "default")
#' @param max_upload_mb Maximum upload file size in MB (default: 50)
#' @param enable_db Logical whether to enable database features (default: FALSE)
#' @param db_path Path to SQLite database (if enable_db = TRUE)
#'
#' @return Runs Shiny app (does not return value)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch with default settings
#' launch_powernma_app()
#'
#' # Launch with dark theme
#' launch_powernma_app(theme = "dark")
#'
#' # Launch with database enabled
#' launch_powernma_app(enable_db = TRUE, db_path = "nma_results.db")
#' }
launch_powernma_app <- function(launch_browser = TRUE,
                                port = NULL,
                                host = "127.0.0.1",
                                theme = c("default", "dark", "light"),
                                max_upload_mb = 50,
                                enable_db = FALSE,
                                db_path = NULL) {

  theme <- match.arg(theme)

  # Check required packages
  required_pkgs <- c("shiny", "shinydashboard", "shinyWidgets",
                     "DT", "plotly", "shinyjs", "shinyBS")

  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Required packages missing: %s\nInstall with: install.packages(c(%s))",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
    ))
  }

  # Set max upload size
  options(shiny.maxRequestSize = max_upload_mb * 1024^2)

  # Load UI and Server
  ui <- create_powernma_ui(theme, enable_db)
  server <- create_powernma_server(enable_db, db_path)

  # Launch app
  shiny::shinyApp(
    ui = ui,
    server = server,
    options = list(
      port = port,
      host = host,
      launch.browser = launch_browser
    )
  )
}

#' Create powerNMA UI
#'
#' @description
#' Creates the user interface for the powerNMA dashboard.
#'
#' @param theme Dashboard theme
#' @param enable_db Enable database features
#'
#' @return Shiny UI object
#'
#' @keywords internal
create_powernma_ui <- function(theme = "default", enable_db = FALSE) {

  # Define header
  header <- shinydashboard::dashboardHeader(
    title = "powerNMA: Ultimate Network Meta-Analysis",
    titleWidth = 400,
    shinydashboard::dropdownMenu(
      type = "notifications",
      icon = shiny::icon("info-circle"),
      badgeStatus = "info",
      shinydashboard::notificationItem(
        text = "Phase 5-9 Complete",
        icon = shiny::icon("check"),
        status = "success"
      ),
      shinydashboard::notificationItem(
        text = "500+ Methods Rules",
        icon = shiny::icon("file-alt"),
        status = "info"
      ),
      shinydashboard::notificationItem(
        text = "AI Enhancement Available",
        icon = shiny::icon("robot"),
        status = "warning"
      )
    )
  )

  # Define sidebar
  sidebar <- shinydashboard::dashboardSidebar(
    width = 250,
    shinydashboard::sidebarMenu(
      id = "sidebar_menu",
      shinydashboard::menuItem("Home", tabName = "home", icon = shiny::icon("home")),
      shinydashboard::menuItem("Data Import", tabName = "data_import", icon = shiny::icon("upload")),
      shinydashboard::menuItem("Analysis", icon = shiny::icon("calculator"),
        shinydashboard::menuSubItem("Standard NMA", tabName = "standard_nma"),
        shinydashboard::menuSubItem("Bayesian NMA", tabName = "bayesian_nma"),
        shinydashboard::menuSubItem("Component NMA", tabName = "component_nma"),
        shinydashboard::menuSubItem("Multivariate NMA", tabName = "multivariate_nma"),
        shinydashboard::menuSubItem("Living NMA", tabName = "living_nma")
      ),
      shinydashboard::menuItem("Rankings", tabName = "rankings", icon = shiny::icon("trophy")),
      shinydashboard::menuItem("Visualizations", tabName = "visualizations", icon = shiny::icon("chart-bar")),
      shinydashboard::menuItem("Diagnostics", icon = shiny::icon("stethoscope"),
        shinydashboard::menuSubItem("Heterogeneity", tabName = "heterogeneity"),
        shinydashboard::menuSubItem("Inconsistency", tabName = "inconsistency"),
        shinydashboard::menuSubItem("Publication Bias", tabName = "pub_bias")
      ),
      shinydashboard::menuItem("Advanced", icon = shiny::icon("cogs"),
        shinydashboard::menuSubItem("Meta-Regression", tabName = "metareg"),
        shinydashboard::menuSubItem("Network Geometry", tabName = "network_geom"),
        shinydashboard::menuSubItem("Simulation", tabName = "simulation"),
        shinydashboard::menuSubItem("Value of Info", tabName = "voi")
      ),
      shinydashboard::menuItem("Manuscripts", icon = shiny::icon("file-alt"),
        shinydashboard::menuSubItem("Generate Methods", tabName = "gen_methods"),
        shinydashboard::menuSubItem("Generate Results", tabName = "gen_results"),
        shinydashboard::menuSubItem("AI Enhancement", tabName = "ai_enhance"),
        shinydashboard::menuSubItem("Text Results", tabName = "text_results")
      ),
      shinydashboard::menuItem("Reports", tabName = "reports", icon = shiny::icon("file-pdf")),
      shinydashboard::menuItem("Batch Processing", tabName = "batch", icon = shiny::icon("tasks")),
      shinydashboard::menuItem("Export", tabName = "export", icon = shiny::icon("download")),
      if (enable_db) {
        shinydashboard::menuItem("Database", tabName = "database", icon = shiny::icon("database"))
      },
      shinydashboard::menuItem("Settings", tabName = "settings", icon = shiny::icon("cog")),
      shinydashboard::menuItem("Help", tabName = "help", icon = shiny::icon("question-circle"))
    )
  )

  # Define body
  body <- shinydashboard::dashboardBody(
    shinyjs::useShinyjs(),

    # Custom CSS
    shiny::tags$head(
      shiny::tags$style(shiny::HTML(get_custom_css(theme)))
    ),

    shinydashboard::tabItems(
      # Home tab
      shinydashboard::tabItem(
        tabName = "home",
        create_home_tab()
      ),

      # Data Import tab
      shinydashboard::tabItem(
        tabName = "data_import",
        create_data_import_tab()
      ),

      # Standard NMA tab
      shinydashboard::tabItem(
        tabName = "standard_nma",
        create_standard_nma_tab()
      ),

      # Bayesian NMA tab
      shinydashboard::tabItem(
        tabName = "bayesian_nma",
        create_bayesian_nma_tab()
      ),

      # Component NMA tab
      shinydashboard::tabItem(
        tabName = "component_nma",
        create_component_nma_tab()
      ),

      # Multivariate NMA tab
      shinydashboard::tabItem(
        tabName = "multivariate_nma",
        create_multivariate_nma_tab()
      ),

      # Living NMA tab
      shinydashboard::tabItem(
        tabName = "living_nma",
        create_living_nma_tab()
      ),

      # Rankings tab
      shinydashboard::tabItem(
        tabName = "rankings",
        create_rankings_tab()
      ),

      # Visualizations tab
      shinydashboard::tabItem(
        tabName = "visualizations",
        create_visualizations_tab()
      ),

      # Heterogeneity tab
      shinydashboard::tabItem(
        tabName = "heterogeneity",
        create_heterogeneity_tab()
      ),

      # Inconsistency tab
      shinydashboard::tabItem(
        tabName = "inconsistency",
        create_inconsistency_tab()
      ),

      # Publication Bias tab
      shinydashboard::tabItem(
        tabName = "pub_bias",
        create_pub_bias_tab()
      ),

      # Meta-Regression tab
      shinydashboard::tabItem(
        tabName = "metareg",
        create_metaregression_tab()
      ),

      # Network Geometry tab
      shinydashboard::tabItem(
        tabName = "network_geom",
        create_network_geometry_tab()
      ),

      # Simulation tab
      shinydashboard::tabItem(
        tabName = "simulation",
        create_simulation_tab()
      ),

      # Value of Information tab
      shinydashboard::tabItem(
        tabName = "voi",
        create_voi_tab()
      ),

      # Generate Methods tab
      shinydashboard::tabItem(
        tabName = "gen_methods",
        create_generate_methods_tab()
      ),

      # Generate Results tab
      shinydashboard::tabItem(
        tabName = "gen_results",
        create_generate_results_tab()
      ),

      # AI Enhancement tab
      shinydashboard::tabItem(
        tabName = "ai_enhance",
        create_ai_enhancement_tab()
      ),

      # Text Results tab
      shinydashboard::tabItem(
        tabName = "text_results",
        create_text_results_tab()
      ),

      # Reports tab
      shinydashboard::tabItem(
        tabName = "reports",
        create_reports_tab()
      ),

      # Batch Processing tab
      shinydashboard::tabItem(
        tabName = "batch",
        create_batch_processing_tab()
      ),

      # Export tab
      shinydashboard::tabItem(
        tabName = "export",
        create_export_tab()
      ),

      # Database tab (if enabled)
      if (enable_db) {
        shinydashboard::tabItem(
          tabName = "database",
          create_database_tab()
        )
      },

      # Settings tab
      shinydashboard::tabItem(
        tabName = "settings",
        create_settings_tab()
      ),

      # Help tab
      shinydashboard::tabItem(
        tabName = "help",
        create_help_tab()
      )
    )
  )

  # Combine into dashboard
  shinydashboard::dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body,
    skin = ifelse(theme == "dark", "black", "blue")
  )
}

#' Get Custom CSS
#'
#' @description
#' Returns custom CSS for dashboard styling.
#'
#' @param theme Theme name
#'
#' @return Character string with CSS
#'
#' @keywords internal
get_custom_css <- function(theme) {

  base_css <- "
    .main-header .logo {
      font-weight: bold;
      font-size: 18px;
    }

    .box {
      border-top: 3px solid #3c8dbc;
    }

    .btn-primary {
      background-color: #3c8dbc;
      border-color: #3c8dbc;
    }

    .btn-success {
      background-color: #00a65a;
      border-color: #00a65a;
    }

    .progress-bar {
      background-color: #3c8dbc;
    }

    .info-box-icon {
      font-size: 50px;
    }

    .shiny-notification {
      position: fixed;
      top: 60px;
      right: 10px;
    }

    .dataTables_wrapper {
      margin-top: 20px;
    }
  "

  if (theme == "dark") {
    base_css <- paste0(base_css, "
      .content-wrapper {
        background-color: #222d32;
        color: #b8c7ce;
      }

      .box {
        background-color: #2c3b41;
        border: 1px solid #1a2226;
        color: #b8c7ce;
      }
    ")
  }

  return(base_css)
}

#' Create Home Tab
#'
#' @description
#' Creates the home/welcome tab.
#'
#' @return Shiny UI elements
#'
#' @keywords internal
create_home_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 12,
      shinydashboard::box(
        title = "Welcome to powerNMA",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        shiny::h3("The Ultimate Network Meta-Analysis Platform"),
        shiny::p("powerNMA provides cutting-edge tools for conducting comprehensive network meta-analyses with:"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Phase 5-9 Features:"), " Advanced rankings, inconsistency assessment, comprehensive visualizations, and more"),
          shiny::tags$li(shiny::strong("Breakthrough 2024-2025 Methods:"), " Component NMA, Multivariate NMA, Living NMA, Personalized Predictions, Value of Information"),
          shiny::tags$li(shiny::strong("Network Analysis:"), " Graph theory metrics, simulation, power analysis"),
          shiny::tags$li(shiny::strong("Interactive Visualizations:"), " Real-time plotly/shiny dashboards with 3D networks"),
          shiny::tags$li(shiny::strong("AI-Powered Manuscripts:"), " 500+ rules per section, 10,000+ permutations, local Ollama integration")
        )
      )
    ),

    # Info boxes
    shiny::column(
      width = 3,
      shinydashboard::infoBox(
        "Features",
        "300+",
        "Functions Available",
        icon = shiny::icon("cogs"),
        color = "blue",
        fill = TRUE
      )
    ),
    shiny::column(
      width = 3,
      shinydashboard::infoBox(
        "Methods",
        "20+",
        "NMA Approaches",
        icon = shiny::icon("calculator"),
        color = "green",
        fill = TRUE
      )
    ),
    shiny::column(
      width = 3,
      shinydashboard::infoBox(
        "Visualizations",
        "30+",
        "Plot Types",
        icon = shiny::icon("chart-line"),
        color = "yellow",
        fill = TRUE
      )
    ),
    shiny::column(
      width = 3,
      shinydashboard::infoBox(
        "Exports",
        "10+",
        "File Formats",
        icon = shiny::icon("download"),
        color = "red",
        fill = TRUE
      )
    ),

    # Quick Start
    shiny::column(
      width = 12,
      shinydashboard::box(
        title = "Quick Start Guide",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        shiny::h4("Getting Started:"),
        shiny::tags$ol(
          shiny::tags$li("Go to ", shiny::strong("Data Import"), " to upload your data or use example datasets"),
          shiny::tags$li("Choose an ", shiny::strong("Analysis"), " type (Standard, Bayesian, Component, Multivariate, or Living NMA)"),
          shiny::tags$li("Explore ", shiny::strong("Rankings"), " and ", shiny::strong("Visualizations")),
          shiny::tags$li("Check ", shiny::strong("Diagnostics"), " (heterogeneity, inconsistency, publication bias)"),
          shiny::tags$li("Generate ", shiny::strong("Manuscripts"), " with AI-powered text"),
          shiny::tags$li("Create ", shiny::strong("Reports"), " and ", shiny::strong("Export"), " your results")
        )
      )
    ),

    # Recent updates
    shiny::column(
      width = 12,
      shinydashboard::box(
        title = "Latest Updates",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Phase 9 Complete:"), " Interactive visualizations with plotly/shiny + AI-powered manuscript generation"),
          shiny::tags$li(shiny::strong("500+ Rules:"), " Methods and Results sections with 10,000+ permutations"),
          shiny::tags$li(shiny::strong("Ollama Integration:"), " Local AI enhancement with 6 modes (readability, clinical, polish, expand, concise, academic)"),
          shiny::tags$li(shiny::strong("3D Visualizations:"), " Network graphs in 3D space"),
          shiny::tags$li(shiny::strong("Comprehensive Text Results:"), " Publication-ready summaries in multiple formats")
        )
      )
    )
  )
}

#' Create powerNMA Server
#'
#' @description
#' Creates the server logic for the powerNMA dashboard.
#'
#' @param enable_db Enable database features
#' @param db_path Path to database
#'
#' @return Shiny server function
#'
#' @keywords internal
create_powernma_server <- function(enable_db = FALSE, db_path = NULL) {

  function(input, output, session) {

    # Reactive values for storing results
    rv <- shiny::reactiveValues(
      data = NULL,
      nma_result = NULL,
      sucra_result = NULL,
      heterogeneity_result = NULL,
      inconsistency_result = NULL,
      publication_bias_result = NULL,
      metareg_result = NULL,
      network_geom_result = NULL,
      component_nma_result = NULL,
      multivariate_nma_result = NULL,
      living_nma_project = NULL,
      voi_result = NULL,
      simulation_result = NULL,
      methods_text = NULL,
      results_text = NULL,
      text_results = NULL,
      analysis_running = FALSE,
      analysis_progress = 0
    )

    # Initialize session
    shiny::observe({
      shiny::showNotification(
        "Welcome to powerNMA! Ready to start analysis.",
        type = "message",
        duration = 3
      )
    })

    # Load tab-specific server logic
    data_import_server(input, output, session, rv)
    standard_nma_server(input, output, session, rv)
    bayesian_nma_server(input, output, session, rv)
    component_nma_server(input, output, session, rv)
    multivariate_nma_server(input, output, session, rv)
    living_nma_server(input, output, session, rv)
    rankings_server(input, output, session, rv)
    visualizations_server(input, output, session, rv)
    heterogeneity_server(input, output, session, rv)
    inconsistency_server(input, output, session, rv)
    pub_bias_server(input, output, session, rv)
    metaregression_server(input, output, session, rv)
    network_geometry_server(input, output, session, rv)
    simulation_server(input, output, session, rv)
    voi_server(input, output, session, rv)
    generate_methods_server(input, output, session, rv)
    generate_results_server(input, output, session, rv)
    ai_enhancement_server(input, output, session, rv)
    text_results_server(input, output, session, rv)
    reports_server(input, output, session, rv)
    batch_processing_server(input, output, session, rv)
    export_server(input, output, session, rv)
    if (enable_db) {
      database_server(input, output, session, rv, db_path)
    }
    settings_server(input, output, session, rv)
    help_server(input, output, session, rv)

    # Session cleanup
    session$onSessionEnded(function() {
      message("powerNMA session ended")
    })
  }
}
