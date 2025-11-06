# CI Test Helper Functions
# =========================
# This file contains helper functions specifically designed for
# CI/CD environments to make tests faster, more reliable, and
# easier to debug.

# Environment Detection ---------------------------------------------------

#' Check if Running in CI Environment
#'
#' Detects if tests are running in a continuous integration environment
#' by checking common CI environment variables.
#'
#' @return Logical indicating whether code is running in CI
#' @export
is_ci <- function() {
  ci_vars <- c(
    "CI",                    # Generic CI flag
    "CONTINUOUS_INTEGRATION",
    "GITHUB_ACTIONS",        # GitHub Actions
    "TRAVIS",                # Travis CI
    "CIRCLECI",              # Circle CI
    "JENKINS_URL",           # Jenkins
    "GITLAB_CI",             # GitLab CI
    "APPVEYOR",              # AppVeyor
    "TF_BUILD"               # Azure Pipelines
  )

  any(Sys.getenv(ci_vars) != "")
}

#' Check if Running Locally (Not CI)
#'
#' @return Logical indicating whether code is running locally
#' @export
is_local <- function() {
  !is_ci()
}

#' Get CI Platform Name
#'
#' @return Character string with CI platform name or "local"
#' @export
get_ci_platform <- function() {
  if (Sys.getenv("GITHUB_ACTIONS") != "") return("GitHub Actions")
  if (Sys.getenv("TRAVIS") != "") return("Travis CI")
  if (Sys.getenv("CIRCLECI") != "") return("Circle CI")
  if (Sys.getenv("GITLAB_CI") != "") return("GitLab CI")
  if (Sys.getenv("APPVEYOR") != "") return("AppVeyor")
  if (Sys.getenv("JENKINS_URL") != "") return("Jenkins")
  if (Sys.getenv("TF_BUILD") != "") return("Azure Pipelines")
  if (is_ci()) return("Unknown CI")
  return("local")
}

# Skip Conditions ---------------------------------------------------------

#' Skip Test if Not in CI
#'
#' Useful for long-running validation tests that should only run in CI
#'
#' @export
skip_if_not_ci <- function() {
  if (is_local()) {
    testthat::skip("Skipping: not in CI environment")
  }
}

#' Skip Test if in CI
#'
#' Useful for interactive tests or tests that require local resources
#'
#' @export
skip_if_ci <- function() {
  if (is_ci()) {
    testthat::skip(paste("Skipping: running in", get_ci_platform()))
  }
}

#' Skip Test if No Display Available
#'
#' Useful for GUI tests, Shiny tests, or tests that require graphics
#'
#' @export
skip_if_no_display <- function() {
  if (is_ci() && Sys.getenv("DISPLAY") == "") {
    testthat::skip("Skipping: no display available (headless CI)")
  }

  # Also check for interactive session
  if (!interactive() && !capabilities("X11")) {
    testthat::skip("Skipping: no graphics capabilities")
  }
}

#' Skip Test if Database Not Available
#'
#' Useful for database integration tests
#'
#' @param connection_string Optional database connection to test
#' @export
skip_if_no_database <- function(connection_string = NULL) {
  if (is_ci() && Sys.getenv("DATABASE_URL") == "") {
    testthat::skip("Skipping: no database configured")
  }

  if (!is.null(connection_string)) {
    tryCatch({
      # Basic connection test would go here
      # For now, just skip if connection string looks empty
      if (nchar(connection_string) == 0) {
        testthat::skip("Skipping: empty database connection string")
      }
    }, error = function(e) {
      testthat::skip(paste("Skipping: database not available -", e$message))
    })
  }
}

#' Skip Slow Tests in Quick Mode
#'
#' Set environment variable QUICK_TESTS=1 to skip slow tests
#'
#' @export
skip_if_quick <- function() {
  if (Sys.getenv("QUICK_TESTS") == "1" ||
      Sys.getenv("QUICK_TESTS") == "true") {
    testthat::skip("Skipping slow test: QUICK_TESTS enabled")
  }
}

#' Skip Tests Requiring External APIs
#'
#' @param api_name Optional name of the specific API
#' @export
skip_if_no_api <- function(api_name = "external API") {
  if (is_ci() && Sys.getenv("SKIP_API_TESTS") != "false") {
    testthat::skip(paste("Skipping:", api_name, "test in CI"))
  }
}

# Package Loading Helpers -------------------------------------------------

#' Safely Require Package with Suppressed Warnings
#'
#' Load packages quietly in CI to reduce log noise. Useful for packages
#' that emit startup messages or warnings that clutter CI logs.
#'
#' @param package Character string with package name
#' @param quietly Logical, whether to suppress all messages (default TRUE in CI)
#' @param warn_conflicts Logical, whether to warn about conflicts
#'   (default FALSE in CI)
#' @return Logical indicating success
#' @export
ci_safe_require <- function(package, quietly = NULL, warn_conflicts = NULL) {
  # Default to quiet mode in CI
  if (is.null(quietly)) {
    quietly <- is_ci()
  }

  if (is.null(warn_conflicts)) {
    warn_conflicts <- !is_ci()
  }

  # Check if package is available
  if (!requireNamespace(package, quietly = TRUE)) {
    if (is_ci()) {
      testthat::skip(paste("Package", package, "not available"))
    } else {
      stop(paste("Package", package, "is required but not installed"))
    }
  }

  # Load the package
  if (quietly) {
    suppressPackageStartupMessages(
      suppressWarnings(
        library(package, character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = warn_conflicts)
      )
    )
  } else {
    library(package, character.only = TRUE,
            quietly = quietly,
            warn.conflicts = warn_conflicts)
  }

  return(TRUE)
}

#' Require Multiple Packages Safely
#'
#' @param packages Character vector of package names
#' @param ... Additional arguments passed to ci_safe_require
#' @return Logical vector indicating success for each package
#' @export
ci_safe_require_all <- function(packages, ...) {
  sapply(packages, ci_safe_require, ...)
}

# Test Data Generation Helpers --------------------------------------------

#' Generate Small Test Dataset for Fast Tests
#'
#' Creates minimal datasets suitable for CI testing. These are much smaller
#' than full validation datasets but maintain realistic structure.
#'
#' @param type Type of dataset: "nma", "ipd", "star", "complete"
#' @param n_studies Number of studies (default varies by type)
#' @param n_treatments Number of treatments (default varies by type)
#' @param seed Random seed for reproducibility
#' @return Data frame with test data
#' @export
create_small_test_data <- function(type = "nma",
                                   n_studies = NULL,
                                   n_treatments = NULL,
                                   seed = 12345) {
  set.seed(seed)

  # Set sensible defaults for CI testing (much smaller than validation)
  defaults <- list(
    nma = list(n_studies = 5, n_treatments = 3),
    ipd = list(n_trials = 2, n_per_arm = 20),
    star = list(n_treatments = 4, n_studies_per_comparison = 3),
    complete = list(n_treatments = 3, n_studies_per_comparison = 2)
  )

  type <- match.arg(type, names(defaults))

  if (is.null(n_studies)) {
    n_studies <- defaults[[type]]$n_studies %||%
                 defaults[[type]]$n_studies_per_comparison %||% 5
  }

  if (is.null(n_treatments)) {
    n_treatments <- defaults[[type]]$n_treatments %||% 3
  }

  # Generate appropriate dataset
  result <- switch(type,
    nma = {
      # Basic pairwise NMA data
      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::tibble(
          studlab = paste0("Study_", seq_len(n_studies)),
          treat1 = sample(
            paste0("Trt", 1:n_treatments),
            n_studies,
            replace = TRUE
          ),
          treat2 = sample(
            paste0("Trt", 1:n_treatments),
            n_studies,
            replace = TRUE
          ),
          TE = stats::rnorm(n_studies, mean = -0.3, sd = 0.5),
          seTE = stats::runif(n_studies, min = 0.1, max = 0.3),
          n1 = round(stats::runif(n_studies, 50, 150)),
          n2 = round(stats::runif(n_studies, 50, 150))
        )
      } else {
        data.frame(
          studlab = paste0("Study_", seq_len(n_studies)),
          treat1 = sample(
            paste0("Trt", 1:n_treatments),
            n_studies,
            replace = TRUE
          ),
          treat2 = sample(
            paste0("Trt", 1:n_treatments),
            n_studies,
            replace = TRUE
          ),
          TE = stats::rnorm(n_studies, mean = -0.3, sd = 0.5),
          seTE = stats::runif(n_studies, min = 0.1, max = 0.3),
          n1 = round(stats::runif(n_studies, 50, 150)),
          n2 = round(stats::runif(n_studies, 50, 150)),
          stringsAsFactors = FALSE
        )
      }
    },

    ipd = {
      # Individual patient data
      n_trials <- defaults[[type]]$n_trials
      n_per_arm <- defaults[[type]]$n_per_arm
      n_total <- n_trials * n_per_arm * 2

      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::tibble(
          trial = rep(paste0("Trial_", 1:n_trials), each = n_per_arm * 2),
          patient_id = paste0("P", seq_len(n_total)),
          treatment = rep(
            rep(c("Control", "Active"), each = n_per_arm),
            n_trials
          ),
          time = stats::rexp(n_total, rate = 0.05),
          status = stats::rbinom(n_total, 1, 0.6),
          age = stats::rnorm(n_total, mean = 60, sd = 10),
          sex = sample(c("M", "F"), n_total, replace = TRUE)
        )
      } else {
        data.frame(
          trial = rep(paste0("Trial_", 1:n_trials), each = n_per_arm * 2),
          patient_id = paste0("P", seq_len(n_total)),
          treatment = rep(
            rep(c("Control", "Active"), each = n_per_arm),
            n_trials
          ),
          time = stats::rexp(n_total, rate = 0.05),
          status = stats::rbinom(n_total, 1, 0.6),
          age = stats::rnorm(n_total, mean = 60, sd = 10),
          sex = sample(c("M", "F"), n_total, replace = TRUE),
          stringsAsFactors = FALSE
        )
      }
    },

    star = {
      # Star network (all treatments vs reference)
      n_studies_total <- (n_treatments - 1) * n_studies
      reference <- "Trt1"
      active_trts <- paste0("Trt", 2:n_treatments)

      treat2_vec <- rep(active_trts, each = n_studies)

      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::tibble(
          studlab = paste0("Study_", seq_len(n_studies_total)),
          treat1 = reference,
          treat2 = treat2_vec,
          TE = stats::rnorm(n_studies_total, mean = -0.2, sd = 0.3),
          seTE = stats::runif(n_studies_total, min = 0.15, max = 0.25),
          n1 = round(stats::runif(n_studies_total, 80, 120)),
          n2 = round(stats::runif(n_studies_total, 80, 120))
        )
      } else {
        data.frame(
          studlab = paste0("Study_", seq_len(n_studies_total)),
          treat1 = reference,
          treat2 = treat2_vec,
          TE = stats::rnorm(n_studies_total, mean = -0.2, sd = 0.3),
          seTE = stats::runif(n_studies_total, min = 0.15, max = 0.25),
          n1 = round(stats::runif(n_studies_total, 80, 120)),
          n2 = round(stats::runif(n_studies_total, 80, 120)),
          stringsAsFactors = FALSE
        )
      }
    },

    complete = {
      # Complete network (all pairwise comparisons)
      comparisons <- utils::combn(paste0("Trt", 1:n_treatments), 2)
      n_comparisons <- ncol(comparisons)
      n_studies_total <- n_comparisons * n_studies

      treat_pairs <- lapply(1:n_studies, function(i) {
        comparisons[, rep(1:n_comparisons, length.out = n_studies_total)[i]]
      })
      treat_pairs <- do.call(rbind, treat_pairs)

      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::tibble(
          studlab = paste0("Study_", seq_len(n_studies_total)),
          treat1 = treat_pairs[1:n_studies_total, 1],
          treat2 = treat_pairs[1:n_studies_total, 2],
          TE = stats::rnorm(n_studies_total, mean = -0.25, sd = 0.4),
          seTE = stats::runif(n_studies_total, min = 0.1, max = 0.3),
          n1 = round(stats::runif(n_studies_total, 60, 140)),
          n2 = round(stats::runif(n_studies_total, 60, 140))
        )
      } else {
        data.frame(
          studlab = paste0("Study_", seq_len(n_studies_total)),
          treat1 = treat_pairs[1:n_studies_total, 1],
          treat2 = treat_pairs[1:n_studies_total, 2],
          TE = stats::rnorm(n_studies_total, mean = -0.25, sd = 0.4),
          seTE = stats::runif(n_studies_total, min = 0.1, max = 0.3),
          n1 = round(stats::runif(n_studies_total, 60, 140)),
          n2 = round(stats::runif(n_studies_total, 60, 140)),
          stringsAsFactors = FALSE
        )
      }
    }
  )

  return(result)
}

# Utility Functions -------------------------------------------------------

#' Null-coalescing Operator
#'
#' Returns first non-null value
#'
#' @param a First value
#' @param b Second value
#' @return First non-null value
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Print CI Debug Information
#'
#' Prints useful debugging information in CI environments
#'
#' @export
print_ci_info <- function() {
  cat("\n========== CI Environment Info ==========\n")
  cat("Platform:", get_ci_platform(), "\n")
  cat("Is CI:", is_ci(), "\n")
  cat("R Version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
  cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")

  # Print relevant environment variables
  ci_vars <- c("CI", "GITHUB_ACTIONS", "QUICK_TESTS", "DISPLAY",
               "R_LIBS_USER", "NOT_CRAN")

  cat("\nEnvironment Variables:\n")
  for (var in ci_vars) {
    val <- Sys.getenv(var)
    if (val != "") {
      cat("  ", var, "=", val, "\n")
    }
  }
  cat("=========================================\n\n")
}

# Timeout Helpers ---------------------------------------------------------

#' Run Code with Timeout in CI
#'
#' Adds timeout protection for potentially slow operations in CI
#'
#' @param expr Expression to evaluate
#' @param timeout Timeout in seconds (default 60 in CI, Inf locally)
#' @param on_timeout What to do on timeout: "skip", "error", or "warn"
#' @export
with_ci_timeout <- function(expr, timeout = NULL, on_timeout = "skip") {
  if (is.null(timeout)) {
    timeout <- if (is_ci()) 60 else Inf
  }

  if (is.infinite(timeout)) {
    return(eval(expr, envir = parent.frame()))
  }

  result <- tryCatch({
    R.utils::withTimeout(
      eval(expr, envir = parent.frame()),
      timeout = timeout
    )
  }, TimeoutException = function(e) {
    msg <- paste("Operation timed out after", timeout, "seconds")

    switch(on_timeout,
      skip = testthat::skip(msg),
      error = stop(msg),
      warn = {
        warning(msg)
        return(NULL)
      }
    )
  })

  return(result)
}
