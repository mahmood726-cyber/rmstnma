#' bs4Dash UI Tab Components with High-Resolution Export
#'
#' @description
#' Modern Bootstrap 4 UI components for powerNMA dashboard tabs with
#' integrated high-resolution download capabilities for publication-quality outputs.
#'
#' @author powerNMA Development Team
#' @name bs4dash_ui_tabs
NULL

#' Create bs4Dash Home Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_home_tab <- function() {
  shiny::fluidRow(
    # Welcome banner
    shiny::column(
      width = 12,
      bs4Dash::bs4Jumbotron(
        title = "Welcome to powerNMA",
        lead = "The Ultimate Network Meta-Analysis Platform - Phase 14 Complete",
        status = "primary",
        btnName = "Get Started",
        href = "#shiny-tab-data_import"
      )
    ),

    # Value boxes for key metrics
    shiny::column(
      width = 3,
      bs4Dash::bs4InfoBox(
        title = "NMA Methods",
        value = "25+",
        icon = shiny::icon("calculator"),
        status = "primary",
        elevation = 2,
        width = 12
      )
    ),
    shiny::column(
      width = 3,
      bs4Dash::bs4InfoBox(
        title = "Visualizations",
        value = "40+",
        icon = shiny::icon("chart-line"),
        status = "success",
        elevation = 2,
        width = 12
      )
    ),
    shiny::column(
      width = 3,
      bs4Dash::bs4InfoBox(
        title = "Export Formats",
        value = "15+",
        icon = shiny::icon("download"),
        status = "warning",
        elevation = 2,
        width = 12
      )
    ),
    shiny::column(
      width = 3,
      bs4Dash::bs4InfoBox(
        title = "High-Res DPI",
        value = "Up to 1200",
        icon = shiny::icon("image"),
        status = "danger",
        elevation = 2,
        width = 12
      )
    ),

    # Features overview
    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Core Features",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        shiny::tags$ul(
          class = "list-unstyled",
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Standard & Bayesian NMA"),
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Component & Multivariate Analysis"),
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Living NMA with Version Control"),
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Publication Bias Detection"),
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Meta-Regression Analysis"),
          shiny::tags$li(shiny::icon("check-circle", class = "text-success"), " Interactive Visualizations")
        )
      )
    ),

    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Advanced Methods (Phase 12-14)",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        shiny::tags$ul(
          class = "list-unstyled",
          shiny::tags$li(shiny::icon("brain", class = "text-primary"), " Graph Neural Networks for NMA"),
          shiny::tags$li(shiny::icon("project-diagram", class = "text-primary"), " Causal Inference Framework"),
          shiny::tags$li(shiny::icon("chart-area", class = "text-primary"), " Distributional NMA"),
          shiny::tags$li(shiny::icon("infinity", class = "text-primary"), " Stan-based Bayesian Analysis"),
          shiny::tags$li(shiny::icon("user-md", class = "text-primary"), " IPD Network Meta-Analysis"),
          shiny::tags$li(shiny::icon("random", class = "text-primary"), " Competing Risks Models")
        )
      )
    ),

    # Quick start guide
    shiny::column(
      width = 12,
      bs4Dash::box(
        title = "Quick Start Guide",
        status = "success",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        width = 12,
        elevation = 2,
        bs4Dash::bs4Timeline(
          width = 12,
          bs4Dash::bs4TimelineEnd(color = "primary"),
          bs4Dash::bs4TimelineLabel("Getting Started", color = "teal"),
          bs4Dash::bs4TimelineItem(
            title = "1. Import Data",
            icon = "upload",
            color = "primary",
            time = "Step 1",
            footer = "Upload your data or use example datasets",
            "Navigate to Data Import tab and load your NMA data in pairwise or IPD format."
          ),
          bs4Dash::bs4TimelineItem(
            title = "2. Run Analysis",
            icon = "calculator",
            color = "success",
            time = "Step 2",
            footer = "Choose from 25+ NMA methods",
            "Select your preferred analysis method (Standard, Bayesian, Component, etc.)."
          ),
          bs4Dash::bs4TimelineItem(
            title = "3. Visualize Results",
            icon = "chart-bar",
            color = "warning",
            time = "Step 3",
            footer = "40+ interactive visualizations",
            "Explore network plots, forest plots, rankings, and more."
          ),
          bs4Dash::bs4TimelineItem(
            title = "4. Export High-Resolution",
            icon = "download",
            color = "danger",
            time = "Step 4",
            footer = "Up to 1200 DPI publication-quality",
            "Download plots and results in multiple formats with highest resolution."
          ),
          bs4Dash::bs4TimelineEnd(color = "success")
        )
      )
    ),

    # Recent updates
    shiny::column(
      width = 12,
      bs4Dash::box(
        title = "Latest Updates - Phase 14",
        status = "danger",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        width = 12,
        elevation = 2,
        bs4Dash::bs4Alert(
          title = "Deep Learning Integration",
          status = "success",
          closable = FALSE,
          "Graph Neural Networks (GNN/GCN/GAT) now available for network structure learning!"
        ),
        bs4Dash::bs4Alert(
          title = "Causal Inference Framework",
          status = "info",
          closable = FALSE,
          "G-formula, TMLE, and E-values for rigorous causal effect estimation."
        ),
        bs4Dash::bs4Alert(
          title = "High-Resolution Exports",
          status = "warning",
          closable = FALSE,
          "Publication-quality downloads up to 1200 DPI in PNG, PDF, SVG, and TIFF formats."
        )
      )
    )
  )
}

#' Create bs4Dash Data Import Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_data_import_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Upload Your Data",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        shiny::fileInput(
          "data_file",
          "Choose CSV/Excel File",
          accept = c(".csv", ".xlsx", ".xls", ".txt", ".rds"),
          width = "100%",
          buttonLabel = "Browse...",
          placeholder = "No file selected"
        ),
        shiny::radioButtons(
          "data_format",
          "Data Format:",
          choices = c(
            "Pairwise (long format)" = "pairwise",
            "Contrast-based" = "contrast",
            "Arm-based" = "arm",
            "Individual Patient Data (IPD)" = "ipd"
          ),
          selected = "pairwise",
          inline = FALSE
        ),
        shiny::selectInput(
          "effect_measure",
          "Effect Measure:",
          choices = c("Odds Ratio (OR)" = "OR",
                     "Risk Ratio (RR)" = "RR",
                     "Hazard Ratio (HR)" = "HR",
                     "Mean Difference (MD)" = "MD",
                     "Standardized MD (SMD)" = "SMD",
                     "Risk Difference (RD)" = "RD"),
          selected = "OR"
        ),
        shiny::actionButton(
          "process_data",
          "Process Data",
          icon = shiny::icon("cogs"),
          class = "btn btn-primary btn-block"
        )
      )
    ),

    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Example Datasets",
        status = "success",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        shiny::selectInput(
          "example_data",
          "Select Example Dataset:",
          choices = c(
            "None" = "none",
            "Smoking cessation (OR)" = "smoking",
            "Depression trials (SMD)" = "depression",
            "Diabetes drugs (MD)" = "diabetes",
            "Cardiovascular (HR)" = "cardio",
            "Multi-arm trial example" = "multiarm",
            "Component NMA example" = "component",
            "IPD example" = "ipd"
          ),
          width = "100%"
        ),
        shiny::actionButton(
          "load_example",
          "Load Example Data",
          icon = shiny::icon("upload"),
          class = "btn btn-success btn-block"
        ),
        shiny::hr(),
        bs4Dash::descriptionBlock(
          number = "0",
          numberColor = "primary",
          numberIcon = shiny::icon("table"),
          header = "Studies",
          text = "Loaded",
          rightBorder = TRUE,
          marginBottom = FALSE
        ),
        bs4Dash::descriptionBlock(
          number = "0",
          numberColor = "success",
          numberIcon = shiny::icon("pills"),
          header = "Treatments",
          text = "In Network",
          rightBorder = FALSE,
          marginBottom = FALSE
        )
      )
    ),

    # Data preview with validation
    shiny::column(
      width = 12,
      bs4Dash::box(
        title = "Data Preview & Validation",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        maximizable = TRUE,
        shiny::uiOutput("data_validation_status"),
        DT::dataTableOutput("data_preview"),
        shiny::br(),
        shiny::div(
          class = "btn-group-highres",
          shiny::actionButton(
            "validate_data",
            "Validate Data",
            icon = shiny::icon("check"),
            class = "btn btn-success"
          ),
          shiny::downloadButton(
            "download_template",
            "Download Template",
            class = "btn btn-secondary"
          ),
          shiny::downloadButton(
            "download_processed_data",
            "Download Processed Data",
            class = "btn btn-info"
          )
        )
      )
    )
  )
}

#' Create bs4Dash Standard NMA Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_standard_nma_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 4,
      bs4Dash::box(
        title = "Analysis Configuration",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        elevation = 2,
        shiny::selectInput(
          "nma_method",
          "NMA Method:",
          choices = c("Frequentist (netmeta)" = "frequentist",
                     "Bayesian (gemtc)" = "bayesian",
                     "Stan (multinma)" = "stan")
        ),
        shiny::selectInput(
          "nma_model",
          "Model:",
          choices = c("Random Effects" = "random",
                     "Fixed Effects" = "fixed",
                     "Common Effects" = "common")
        ),
        shiny::selectInput(
          "nma_reference",
          "Reference Treatment:",
          choices = NULL
        ),
        shiny::checkboxInput("nma_small_good", "Lower is better", FALSE),
        shiny::hr(),
        shiny::actionButton(
          "run_standard_nma",
          "Run Analysis",
          icon = shiny::icon("play"),
          class = "btn btn-success btn-lg btn-block"
        ),
        shiny::br(),
        shiny::uiOutput("nma_progress_ui")
      ),

      bs4Dash::box(
        title = "Advanced Options",
        status = "warning",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        width = 12,
        elevation = 2,
        shiny::checkboxInput("nma_assess_incon", "Assess Inconsistency", TRUE),
        shiny::checkboxInput("nma_assess_het", "Assess Heterogeneity", TRUE),
        shiny::checkboxInput("nma_assess_pub", "Assess Publication Bias", TRUE),
        shiny::checkboxInput("nma_calc_rankings", "Calculate Rankings", TRUE),
        shiny::numericInput("nma_nsim", "N Simulations:", value = 10000, min = 1000, max = 100000)
      )
    ),

    shiny::column(
      width = 8,
      bs4Dash::tabBox(
        title = "Analysis Results",
        id = "nma_results_tabs",
        width = 12,
        elevation = 2,
        side = "right",
        type = "tabs",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        maximizable = TRUE,

        shiny::tabPanel(
          "Summary",
          icon = shiny::icon("info-circle"),
          shiny::verbatimTextOutput("nma_summary_text"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_summary_pdf", "PDF", class = "btn btn-sm btn-primary"),
            shiny::downloadButton("download_summary_docx", "DOCX", class = "btn btn-sm btn-info"),
            shiny::downloadButton("download_summary_txt", "TXT", class = "btn btn-sm btn-secondary")
          )
        ),

        shiny::tabPanel(
          "Effect Estimates",
          icon = shiny::icon("table"),
          DT::dataTableOutput("nma_effects_table"),
          shiny::hr(),
          shiny::downloadButton("download_effects_csv", "Download CSV", class = "btn btn-sm btn-success")
        ),

        shiny::tabPanel(
          "Network Plot",
          icon = shiny::icon("project-diagram"),
          plotly::plotlyOutput("nma_network_plot", height = "600px"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_network_png_300", "PNG 300 DPI", class = "btn btn-sm btn-primary"),
            shiny::downloadButton("download_network_png_600", "PNG 600 DPI", class = "btn btn-sm btn-success"),
            shiny::downloadButton("download_network_png_1200", "PNG 1200 DPI", class = "btn btn-sm btn-danger"),
            shiny::downloadButton("download_network_pdf", "PDF Vector", class = "btn btn-sm btn-info"),
            shiny::downloadButton("download_network_svg", "SVG", class = "btn btn-sm btn-warning"),
            shiny::downloadButton("download_network_html", "Interactive HTML", class = "btn btn-sm btn-dark")
          )
        ),

        shiny::tabPanel(
          "Forest Plot",
          icon = shiny::icon("chart-bar"),
          plotly::plotlyOutput("nma_forest_plot", height = "700px"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_forest_png_300", "PNG 300 DPI"),
            shiny::downloadButton("download_forest_png_600", "PNG 600 DPI"),
            shiny::downloadButton("download_forest_png_1200", "PNG 1200 DPI"),
            shiny::downloadButton("download_forest_pdf", "PDF Vector"),
            shiny::downloadButton("download_forest_svg", "SVG"),
            shiny::downloadButton("download_forest_tiff", "TIFF")
          )
        ),

        shiny::tabPanel(
          "League Table",
          icon = shiny::icon("th"),
          DT::dataTableOutput("nma_league_table"),
          shiny::hr(),
          shiny::downloadButton("download_league_csv", "Download CSV")
        )
      )
    )
  )
}

#' Create bs4Dash Rankings Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_rankings_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 4,
      bs4Dash::box(
        title = "Ranking Configuration",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        shiny::radioButtons(
          "ranking_method",
          "Ranking Method:",
          choices = c("SUCRA" = "sucra",
                     "P-score" = "pscore",
                     "Mean Rank" = "meanrank"),
          selected = "sucra"
        ),
        shiny::checkboxInput("ranking_small_good", "Lower is better", FALSE),
        shiny::numericInput("ranking_nsim", "N Simulations:", value = 10000, min = 1000),
        shiny::actionButton(
          "calc_rankings",
          "Calculate Rankings",
          icon = shiny::icon("trophy"),
          class = "btn btn-success btn-block"
        )
      ),

      bs4Dash::bs4Card(
        title = "Top 5 Treatments",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        collapsible = FALSE,
        DT::dataTableOutput("rankings_top_table")
      )
    ),

    shiny::column(
      width = 8,
      bs4Dash::tabBox(
        title = "Treatment Rankings",
        width = 12,
        elevation = 2,
        type = "tabs",
        side = "right",

        shiny::tabPanel(
          "SUCRA Plot",
          plotly::plotlyOutput("rankings_sucra_plot", height = "600px"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_sucra_png_600", "PNG 600 DPI"),
            shiny::downloadButton("download_sucra_png_1200", "PNG 1200 DPI"),
            shiny::downloadButton("download_sucra_pdf", "PDF"),
            shiny::downloadButton("download_sucra_svg", "SVG")
          )
        ),

        shiny::tabPanel(
          "Rankogram",
          plotly::plotlyOutput("rankings_rankogram", height = "600px"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_rankogram_png_600", "PNG 600 DPI"),
            shiny::downloadButton("download_rankogram_png_1200", "PNG 1200 DPI"),
            shiny::downloadButton("download_rankogram_pdf", "PDF")
          )
        ),

        shiny::tabPanel(
          "Full Rankings Table",
          DT::dataTableOutput("rankings_full_table"),
          shiny::downloadButton("download_rankings_csv", "Download CSV")
        )
      )
    )
  )
}

#' Create bs4Dash Visualizations Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_visualizations_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 3,
      bs4Dash::box(
        title = "Visualization Gallery",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        shiny::selectInput(
          "viz_type",
          "Select Visualization:",
          choices = c(
            "Network Graph" = "network",
            "3D Network" = "network_3d",
            "Forest Plot" = "forest",
            "Funnel Plot" = "funnel",
            "Heat Map" = "heatmap",
            "Interval Plot" = "interval",
            "Comparison Matrix" = "comparison",
            "Net Heat Plot" = "netheat",
            "Contour Funnel" = "contour_funnel",
            "SUCRA Bar Chart" = "sucra_bar",
            "Rankogram" = "rankogram"
          )
        ),
        shiny::hr(),
        shiny::h5("Customization"),
        shiny::selectInput(
          "viz_layout",
          "Network Layout:",
          choices = c("Circle" = "circle",
                     "Spring" = "spring",
                     "Tree" = "tree",
                     "Star" = "star")
        ),
        shiny::sliderInput(
          "viz_size",
          "Node Size:",
          min = 10,
          max = 100,
          value = 40
        ),
        shiny::checkboxInput("viz_labels", "Show Labels", TRUE),
        shiny::checkboxInput("viz_weights", "Weight by Evidence", TRUE),
        shiny::actionButton(
          "generate_viz",
          "Generate Visualization",
          icon = shiny::icon("chart-bar"),
          class = "btn btn-primary btn-block"
        )
      ),

      bs4Dash::box(
        title = "High-Resolution Export",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        collapsible = TRUE,
        shiny::sliderInput(
          "viz_export_dpi",
          "DPI (Resolution):",
          min = 300,
          max = 1200,
          value = 600,
          step = 100
        ),
        shiny::sliderInput(
          "viz_export_width",
          "Width (inches):",
          min = 5,
          max = 20,
          value = 10,
          step = 1
        ),
        shiny::sliderInput(
          "viz_export_height",
          "Height (inches):",
          min = 5,
          max = 20,
          value = 8,
          step = 1
        ),
        shiny::div(
          class = "btn-group-highres",
          style = "display: flex; flex-direction: column; gap: 5px;",
          shiny::downloadButton("download_viz_png_custom", "PNG (Custom DPI)", class = "btn btn-sm btn-primary btn-block"),
          shiny::downloadButton("download_viz_pdf", "PDF (Vector)", class = "btn btn-sm btn-success btn-block"),
          shiny::downloadButton("download_viz_svg", "SVG (Vector)", class = "btn btn-sm btn-info btn-block"),
          shiny::downloadButton("download_viz_tiff", "TIFF (Lossless)", class = "btn btn-sm btn-warning btn-block"),
          shiny::downloadButton("download_viz_html", "HTML (Interactive)", class = "btn btn-sm btn-dark btn-block")
        )
      )
    ),

    shiny::column(
      width = 9,
      bs4Dash::box(
        title = shiny::textOutput("viz_title"),
        status = "info",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        maximizable = TRUE,
        collapsible = TRUE,
        plotly::plotlyOutput("main_visualization", height = "750px")
      )
    )
  )
}

#' Create bs4Dash GNN Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_gnn_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 12,
      bs4Dash::bs4Jumbotron(
        title = "Graph Neural Networks for NMA",
        lead = "Revolutionary deep learning approach for network structure learning",
        status = "danger",
        btnName = NULL
      )
    ),

    shiny::column(
      width = 4,
      bs4Dash::box(
        title = "GNN Configuration",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        shiny::selectInput(
          "gnn_architecture",
          "Architecture:",
          choices = c("GCN (Graph Convolutional)" = "GCN",
                     "GAT (Graph Attention)" = "GAT",
                     "GraphSAGE" = "GraphSAGE",
                     "GIN (Graph Isomorphism)" = "GIN")
        ),
        shiny::numericInput("gnn_layers", "Number of Layers:", value = 3, min = 1, max = 10),
        shiny::numericInput("gnn_hidden_dim", "Hidden Dimension:", value = 64, min = 16, max = 256, step = 16),
        shiny::sliderInput("gnn_dropout", "Dropout Rate:", min = 0, max = 0.9, value = 0.5, step = 0.1),
        shiny::numericInput("gnn_epochs", "Training Epochs:", value = 200, min = 50, max = 1000),
        shiny::actionButton(
          "train_gnn",
          "Train GNN Model",
          icon = shiny::icon("brain"),
          class = "btn btn-danger btn-block"
        )
      )
    ),

    shiny::column(
      width = 8,
      bs4Dash::tabBox(
        title = "GNN Results",
        width = 12,
        shiny::tabPanel(
          "Treatment Embeddings",
          plotly::plotlyOutput("gnn_embeddings_plot", height = "600px"),
          shiny::hr(),
          shiny::div(
            class = "btn-group-highres",
            shiny::downloadButton("download_gnn_emb_png_600", "PNG 600 DPI"),
            shiny::downloadButton("download_gnn_emb_pdf", "PDF")
          )
        ),
        shiny::tabPanel(
          "Link Predictions",
          DT::dataTableOutput("gnn_predictions_table")
        ),
        shiny::tabPanel(
          "Training History",
          plotly::plotlyOutput("gnn_training_plot")
        )
      )
    )
  )
}

#' Create bs4Dash Causal Inference Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_causal_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 12,
      bs4Dash::bs4Jumbotron(
        title = "Causal Inference Framework",
        lead = "Rigorous causal effect estimation beyond associations",
        status = "info",
        btnName = NULL
      )
    ),

    shiny::column(
      width = 4,
      bs4Dash::box(
        title = "Causal Methods",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        shiny::selectInput(
          "causal_method",
          "Method:",
          choices = c("G-formula" = "gformula",
                     "TMLE" = "tmle",
                     "Doubly Robust" = "dr")
        ),
        shiny::checkboxInput("causal_sensitivity", "Run Sensitivity Analysis", TRUE),
        shiny::checkboxInput("causal_evalues", "Calculate E-values", TRUE),
        shiny::actionButton(
          "run_causal",
          "Run Causal Analysis",
          icon = shiny::icon("project-diagram"),
          class = "btn btn-info btn-block"
        )
      )
    ),

    shiny::column(
      width = 8,
      bs4Dash::tabBox(
        title = "Causal Results",
        width = 12,
        shiny::tabPanel(
          "Causal Effects",
          DT::dataTableOutput("causal_effects_table")
        ),
        shiny::tabPanel(
          "E-values Plot",
          plotly::plotlyOutput("causal_evalues_plot", height = "600px"),
          shiny::hr(),
          shiny::downloadButton("download_evalues_png_600", "PNG 600 DPI")
        )
      )
    )
  )
}

#' Create bs4Dash Export Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_export_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 12,
      bs4Dash::box(
        title = "High-Resolution Export Center",
        status = "success",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        shiny::h4("Publication-Quality Outputs"),
        shiny::p("Export all results and visualizations at the highest resolution for publication, presentation, or archival purposes.")
      )
    ),

    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Export Configuration",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        shiny::sliderInput(
          "global_export_dpi",
          "Global DPI Setting:",
          min = 300,
          max = 1200,
          value = 600,
          step = 100
        ),
        shiny::selectInput(
          "global_export_format",
          "Preferred Format:",
          choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg",
                     "TIFF" = "tiff", "HTML" = "html")
        ),
        shiny::checkboxInput("export_include_data", "Include Raw Data", TRUE),
        shiny::checkboxInput("export_include_code", "Include R Code", FALSE)
      )
    ),

    shiny::column(
      width = 6,
      bs4Dash::box(
        title = "Batch Export",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        shiny::checkboxGroupInput(
          "batch_export_items",
          "Select Items to Export:",
          choices = c("All Plots" = "plots",
                     "All Tables" = "tables",
                     "Summary Statistics" = "summary",
                     "Raw Data" = "data",
                     "Analysis Code" = "code")
        ),
        shiny::actionButton(
          "batch_export_all",
          "Export All Selected",
          icon = shiny::icon("download"),
          class = "btn btn-warning btn-block"
        )
      )
    )
  )
}

#' Create bs4Dash Help Tab
#'
#' @return Shiny UI elements
#' @keywords internal
create_bs4dash_help_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      width = 12,
      bs4Dash::box(
        title = "powerNMA Help & Documentation",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        elevation = 2,
        bs4Dash::bs4Accordion(
          id = "help_accordion",
          bs4Dash::bs4AccordionItem(
            title = "Getting Started",
            status = "primary",
            collapsed = FALSE,
            "Follow the Quick Start Guide on the Home dashboard to begin your analysis."
          ),
          bs4Dash::bs4AccordionItem(
            title = "High-Resolution Exports",
            status = "success",
            collapsed = TRUE,
            "All plots can be exported at resolutions up to 1200 DPI. Use the download buttons below each visualization or configure global settings in the Export tab."
          ),
          bs4Dash::bs4AccordionItem(
            title = "Advanced Methods",
            status = "warning",
            collapsed = TRUE,
            "Phase 12-14 features include Graph Neural Networks, Causal Inference, Stan-based Bayesian NMA, and more. Access these from the Advanced Methods menu."
          ),
          bs4Dash::bs4AccordionItem(
            title = "Support",
            status = "danger",
            collapsed = TRUE,
            "For questions or issues, visit our GitHub repository or consult the package documentation with ?powerNMA"
          )
        )
      )
    )
  )
}
