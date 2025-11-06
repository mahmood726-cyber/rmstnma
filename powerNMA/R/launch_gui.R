# ============================================================================
# Launch powerNMA Shiny GUI
# ============================================================================

#' Launch powerNMA Interactive GUI
#'
#' Launches an interactive Shiny dashboard for powerNMA using bs4dash.
#' Provides easy access to all four analysis pathways through a modern
#' web-based interface.
#'
#' @param port The TCP port that the application should listen on (default: random)
#' @param launch.browser If TRUE, launch system's default web browser (default: TRUE)
#' @param host The IPv4 address that the application should listen on (default: "127.0.0.1")
#'
#' @details
#' The powerNMA GUI provides:
#'
#' \strong{Features:}
#' \itemize{
#'   \item Welcome page with pathway overview
#'   \item Data upload (CSV, Excel, example datasets)
#'   \item 4 analysis pathways:
#'     \itemize{
#'       \item Pathway 1: Manual Standard (full control, validated methods)
#'       \item Pathway 2: Manual Experimental (full control, cutting-edge)
#'       \item Pathway 3: Auto Standard (automatic, validated)
#'       \item Pathway 4: Auto Experimental (automatic, cutting-edge)
#'     }
#'   \item Interactive results visualization
#'   \item Export reports (HTML, PDF, CSV)
#'   \item Help and documentation
#' }
#'
#' \strong{Requirements:}
#' The following packages are required:
#' \itemize{
#'   \item shiny
#'   \item bs4Dash
#'   \item DT
#'   \item plotly
#'   \item shinyWidgets
#' }
#'
#' @return Launches the Shiny application. Returns nothing.
#'
#' @examples
#' \dontrun{
#' # Launch the GUI
#' launch_powernma_gui()
#'
#' # Launch on specific port
#' launch_powernma_gui(port = 8080)
#'
#' # Launch without opening browser
#' launch_powernma_gui(launch.browser = FALSE)
#' }
#'
#' @export
launch_powernma_gui <- function(port = NULL,
                                launch.browser = TRUE,
                                host = "127.0.0.1") {

  # Check if required packages are installed
  required_packages <- c("shiny", "bs4Dash", "DT", "plotly", "shinyWidgets")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(paste0(
      "The following packages are required to run the GUI but are not installed:\n",
      paste(missing_packages, collapse = ", "), "\n\n",
      "Install them with:\n",
      "install.packages(c('", paste(missing_packages, collapse = "', '"), "'))"
    ))
  }

  # Get path to Shiny app
  app_dir <- system.file("shiny", package = "powerNMA")

  if (app_dir == "") {
    stop("Could not find Shiny app directory. Please reinstall powerNMA.")
  }

  # Welcome message
  message("=============================================================")
  message("Launching powerNMA Interactive GUI")
  message("=============================================================")
  message("")
  message("Welcome to powerNMA v3.0!")
  message("")
  message("The GUI provides access to:")
  message("  - 4 analysis pathways (Manual + Auto, Standard + Experimental)")
  message("  - 11 network meta-analysis methods")
  message("  - Interactive data upload and visualization")
  message("  - Report generation and export")
  message("")
  message("To stop the GUI, press Ctrl+C (or Esc in RStudio)")
  message("=============================================================")
  message("")

  # Launch Shiny app
  if (is.null(port)) {
    shiny::runApp(
      appDir = app_dir,
      launch.browser = launch.browser,
      host = host
    )
  } else {
    shiny::runApp(
      appDir = app_dir,
      port = port,
      launch.browser = launch.browser,
      host = host
    )
  }
}


#' Alias for launch_powernma_gui
#'
#' @inheritParams launch_powernma_gui
#' @export
powernma_gui <- function(port = NULL,
                         launch.browser = TRUE,
                         host = "127.0.0.1") {
  launch_powernma_gui(port = port, launch.browser = launch.browser, host = host)
}


#' Check GUI Dependencies
#'
#' Checks if all required packages for the Shiny GUI are installed.
#'
#' @return Invisible TRUE if all dependencies are installed, otherwise prints
#'   missing packages and returns FALSE.
#'
#' @examples
#' \dontrun{
#' # Check if GUI can be launched
#' check_gui_dependencies()
#' }
#'
#' @export
check_gui_dependencies <- function() {

  required_packages <- c("shiny", "bs4Dash", "DT", "plotly", "shinyWidgets")

  installed <- sapply(required_packages, requireNamespace, quietly = TRUE)

  if (all(installed)) {
    message("✓ All GUI dependencies are installed!")
    message("  You can launch the GUI with: launch_powernma_gui()")
    return(invisible(TRUE))
  } else {
    missing <- required_packages[!installed]

    message("✗ Missing GUI dependencies:")
    message("  ", paste(missing, collapse = ", "))
    message("")
    message("Install them with:")
    message("  install.packages(c('", paste(missing, collapse = "', '"), "'))")

    return(invisible(FALSE))
  }
}


#' Install GUI Dependencies
#'
#' Installs all required packages for the Shiny GUI.
#'
#' @param repos CRAN mirror to use (default: "https://cloud.r-project.org")
#'
#' @return Invisible NULL
#'
#' @examples
#' \dontrun{
#' # Install all GUI dependencies
#' install_gui_dependencies()
#' }
#'
#' @export
install_gui_dependencies <- function(repos = "https://cloud.r-project.org") {

  required_packages <- c("shiny", "bs4Dash", "DT", "plotly", "shinyWidgets")

  message("Installing GUI dependencies...")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("  Installing ", pkg, "...")
      utils::install.packages(pkg, repos = repos)
    } else {
      message("  ", pkg, " already installed")
    }
  }

  message("")
  message("✓ All GUI dependencies installed!")
  message("  Launch the GUI with: launch_powernma_gui()")

  invisible(NULL)
}
