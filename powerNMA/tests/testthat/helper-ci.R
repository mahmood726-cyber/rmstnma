# Test Helper Functions for CI Environment
# Provides utilities for conditional test execution

#' Skip test if running on CI
#' @keywords internal
skip_on_ci <- function() {
  if (identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping on CI")
  }
}

#' Skip test if package is not installed
#' @param pkg Package name
#' @keywords internal
skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    testthat::skip(paste("Package", pkg, "not installed"))
  }
}

#' Skip test if running on specific OS
#' @param os Operating system name ("windows", "mac", "linux")
#' @keywords internal
skip_on_os <- function(os) {
  os_match <- switch(
    tolower(os),
    windows = .Platform$OS.type == "windows",
    mac = Sys.info()[["sysname"]] == "Darwin",
    linux = Sys.info()[["sysname"]] == "Linux",
    FALSE
  )

  if (os_match) {
    testthat::skip(paste("Skipping on", os))
  }
}

#' Skip long-running tests unless explicitly requested
#' @keywords internal
skip_if_quick <- function() {
  if (identical(Sys.getenv("RUN_SLOW_TESTS"), "false") ||
      identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping slow test")
  }
}

#' Skip tests requiring graphical display
#' @keywords internal
skip_if_no_display <- function() {
  if (identical(Sys.getenv("DISPLAY"), "") ||
      identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping test requiring display")
  }
}

#' Skip tests requiring database connection
#' @keywords internal
skip_if_no_database <- function() {
  if (identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping database test on CI")
  }
}

#' Skip tests requiring Stan
#' @keywords internal
skip_if_no_stan <- function() {
  if (!requireNamespace("rstan", quietly = TRUE) ||
      identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping Stan test")
  }
}

#' Skip tests requiring network/API access
#' @keywords internal
skip_if_offline <- function() {
  if (identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping online test on CI")
  }
}

#' Create small test dataset for CI
#' @keywords internal
create_small_test_data <- function() {
  # Smaller dataset for faster CI tests
  data.frame(
    study = rep(1:3, each = 2),
    treatment = rep(c("A", "B"), 3),
    mean = rnorm(6),
    se = runif(6, 0.1, 0.3),
    n = rep(50, 6)
  )
}

#' Suppress warnings from specific packages in CI
#' @keywords internal
ci_safe_require <- function(pkg) {
  suppressWarnings(
    suppressMessages(
      requireNamespace(pkg, quietly = TRUE)
    )
  )
}
