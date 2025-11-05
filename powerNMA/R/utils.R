#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Check if package is available
#' @param pkg Package name
#' @keywords internal
has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Message with timestamp
#' @param ... Message components
#' @keywords internal
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(...), "\n", sep = "")
}

#' Safe numeric clipping
#' @param x Numeric vector
#' @param lo Lower bound
#' @param hi Upper bound
#' @keywords internal
safe_clip <- function(x, lo, hi) {
  pmin(hi, pmax(lo, x))
}

#' Coalesce NA values
#' @param x Vector
#' @param val Replacement value
#' @keywords internal
vcoalesce <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}

#' Safe try-catch wrapper
#' @param expr Expression to evaluate
#' @param context Context description
#' @param silent Suppress error messages
#' @keywords internal
.safe_try <- function(expr, context = "", silent = TRUE) {
  out <- try(expr, silent = silent)
  if (inherits(out, "try-error")) {
    msg("ERROR in %s: %s", context,
        as.character(attr(out, "condition")$message %||% out))
  }
  out
}

#' Validate that argument is not NULL
#' @param x Value to check
#' @param arg_name Name of the argument
#' @param func_name Name of the calling function
#' @keywords internal
assert_not_null <- function(x, arg_name, func_name = "") {
  if (is.null(x)) {
    context <- if (nchar(func_name) > 0) sprintf("[%s] ", func_name) else ""
    stop(sprintf("%sArgument '%s' cannot be NULL", context, arg_name),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate that argument is a data frame
#' @param x Value to check
#' @param arg_name Name of the argument
#' @param func_name Name of the calling function
#' @param allow_empty Allow empty data frames
#' @keywords internal
assert_data_frame <- function(x, arg_name, func_name = "", allow_empty = FALSE) {
  context <- if (nchar(func_name) > 0) sprintf("[%s] ", func_name) else ""

  if (!is.data.frame(x)) {
    stop(sprintf("%sArgument '%s' must be a data.frame, got %s",
                 context, arg_name, class(x)[1]),
         call. = FALSE)
  }

  if (!allow_empty && nrow(x) == 0) {
    stop(sprintf("%sArgument '%s' cannot be an empty data frame",
                 context, arg_name),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate that argument is numeric
#' @param x Value to check
#' @param arg_name Name of the argument
#' @param func_name Name of the calling function
#' @param min_value Minimum allowed value (inclusive)
#' @param max_value Maximum allowed value (inclusive)
#' @keywords internal
assert_numeric <- function(x, arg_name, func_name = "",
                          min_value = -Inf, max_value = Inf) {
  context <- if (nchar(func_name) > 0) sprintf("[%s] ", func_name) else ""

  if (!is.numeric(x)) {
    stop(sprintf("%sArgument '%s' must be numeric, got %s",
                 context, arg_name, class(x)[1]),
         call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(sprintf("%sArgument '%s' contains non-finite values",
                 context, arg_name),
         call. = FALSE)
  }

  if (!all(x >= min_value & x <= max_value)) {
    stop(sprintf("%sArgument '%s' must be in range [%g, %g]",
                 context, arg_name, min_value, max_value),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate that argument is a positive integer
#' @param x Value to check
#' @param arg_name Name of the argument
#' @param func_name Name of the calling function
#' @keywords internal
assert_positive_integer <- function(x, arg_name, func_name = "") {
  context <- if (nchar(func_name) > 0) sprintf("[%s] ", func_name) else ""

  if (!is.numeric(x) || length(x) != 1) {
    stop(sprintf("%sArgument '%s' must be a single numeric value",
                 context, arg_name),
         call. = FALSE)
  }

  if (!is.finite(x) || x < 1 || x != floor(x)) {
    stop(sprintf("%sArgument '%s' must be a positive integer, got %g",
                 context, arg_name, x),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate required columns exist in data frame
#' @param data Data frame to check
#' @param required_cols Vector of required column names
#' @param func_name Name of the calling function
#' @keywords internal
assert_columns_exist <- function(data, required_cols, func_name = "") {
  context <- if (nchar(func_name) > 0) sprintf("[%s] ", func_name) else ""

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("%sMissing required columns: %s",
                 context, paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Simple memoization cache
#' @keywords internal
.powernma_cache <- new.env(parent = emptyenv())

#' Memoize function results
#' @param key Cache key
#' @param expr Expression to evaluate
#' @param enable_cache Enable caching
#' @keywords internal
memoize <- function(key, expr, enable_cache = FALSE) {
  if (!enable_cache) return(eval.parent(substitute(expr)))

  if (exists(key, envir = .powernma_cache, inherits = FALSE)) {
    return(get(key, envir = .powernma_cache))
  }

  val <- eval.parent(substitute(expr))
  assign(key, val, envir = .powernma_cache)
  val
}

#' Clear memoization cache
#' @export
clear_powernma_cache <- function() {
  rm(list = ls(envir = .powernma_cache), envir = .powernma_cache)
  msg("Cache cleared")
}

#' Validate IPD structure
#'
#' Checks if an IPD data frame has the required columns and data types.
#'
#' @param ipd The IPD data frame to validate
#' @return `TRUE` if validation passes, otherwise throws an error
#' @export
#' @examples
#' \dontrun{
#' ipd <- generate_example_ipd()
#' validate_ipd(ipd)
#' }
validate_ipd <- function(ipd) {
  # Enhanced validation with better error messages
  assert_not_null(ipd, "ipd", "validate_ipd")
  assert_data_frame(ipd, "ipd", "validate_ipd", allow_empty = FALSE)

  required_cols <- c("trial", "treatment", "time", "status")
  assert_columns_exist(ipd, required_cols, "validate_ipd")

  # Type validation
  if (!is.numeric(ipd$time)) {
    stop("[validate_ipd] Column 'time' must be numeric, got ",
         class(ipd$time)[1], call. = FALSE)
  }

  # Status column validation
  if (!all(ipd$status %in% c(0, 1, NA))) {
    invalid_vals <- unique(ipd$status[!ipd$status %in% c(0, 1, NA)])
    stop(sprintf("[validate_ipd] Column 'status' must contain only 0 and 1, found: %s",
                 paste(invalid_vals, collapse = ", ")),
         call. = FALSE)
  }

  # Negative time check
  if (any(ipd$time < 0, na.rm = TRUE)) {
    n_negative <- sum(ipd$time < 0, na.rm = TRUE)
    warning(sprintf("[validate_ipd] Detected %d negative time values", n_negative),
            call. = FALSE)
  }

  # NA check
  if (any(is.na(ipd$time))) {
    n_na <- sum(is.na(ipd$time))
    warning(sprintf("[validate_ipd] Detected %d NA values in time column", n_na),
            call. = FALSE)
  }

  TRUE
}

#' Validate NMA pairwise data
#'
#' Validates pairwise data for network meta-analysis
#'
#' @param data Data frame with NMA data
#' @return `TRUE` if validation passes
#' @export
validate_nma_input <- function(data) {
  # Enhanced validation with better error messages
  assert_not_null(data, "data", "validate_nma_input")
  assert_data_frame(data, "data", "validate_nma_input", allow_empty = FALSE)

  # Check required columns
  required_cols <- c("studlab", "treat1", "treat2", "TE", "seTE")
  assert_columns_exist(data, required_cols, "validate_nma_input")

  # Validate TE column
  if (any(!is.finite(data$TE))) {
    n_bad <- sum(!is.finite(data$TE))
    stop(sprintf("[validate_nma_input] Column 'TE' contains %d non-finite values (Inf, -Inf, or NA)",
                 n_bad),
         call. = FALSE)
  }

  # Validate seTE column
  if (any(!is.finite(data$seTE))) {
    n_bad <- sum(!is.finite(data$seTE))
    stop(sprintf("[validate_nma_input] Column 'seTE' contains %d non-finite values (Inf, -Inf, or NA)",
                 n_bad),
         call. = FALSE)
  }

  if (any(data$seTE <= 0, na.rm = TRUE)) {
    n_bad <- sum(data$seTE <= 0, na.rm = TRUE)
    min_val <- min(data$seTE[data$seTE <= 0], na.rm = TRUE)
    stop(sprintf("[validate_nma_input] Column 'seTE' must be strictly positive. Found %d values <= 0 (min: %g)",
                 n_bad, min_val),
         call. = FALSE)
  }

  # Check for sufficient data
  if (nrow(data) < 2) {
    stop("[validate_nma_input] Need at least 2 comparisons for network meta-analysis",
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Clean NMA data
#'
#' Cleans and standardizes NMA data
#'
#' @param data Input data frame
#' @return Cleaned data frame
#' @export
clean_nma_data <- function(data) {
  data <- data %>%
    dplyr::mutate(
      studlab = as.character(studlab),
      treat1 = as.character(treat1),
      treat2 = as.character(treat2)
    )

  data <- data %>%
    dplyr::filter(is.finite(TE), is.finite(seTE), seTE > 0)

  data
}

#' Print method for RMST NMA results
#' @param x An `rmst_nma` object
#' @param ... Additional arguments
#' @export
#' @method print rmst_nma
print.rmst_nma <- function(x, ...) {
  cat("RMST Network Meta-Analysis Results\n")
  cat("===================================\n")
  cat("Time points analyzed:",
      paste(gsub("tau_", "", names(x)), collapse = ", "),
      "days\n\n")

  for (tau_name in names(x)) {
    cat("--- Analysis at", gsub("tau_", "", tau_name), "days ---\n")
    print(x[[tau_name]], ...)
    cat("\n")
  }
}

#' Print method for Milestone NMA results
#' @param x A `milestone_nma` object
#' @param ... Additional arguments
#' @export
#' @method print milestone_nma
print.milestone_nma <- function(x, ...) {
  cat("Milestone Survival Network Meta-Analysis Results\n")
  cat("=================================================\n")
  cat("Time points analyzed:",
      paste(gsub("day_", "", names(x)), collapse = ", "),
      "days\n\n")

  for (time_name in names(x)) {
    cat("--- Analysis at", gsub("day_", "", time_name), "days ---\n")
    print(x[[time_name]], ...)
    cat("\n")
  }
}

#' Print method for powerNMA results
#' @param x A `powernma_result` object
#' @param ... Additional arguments
#' @export
#' @method print powernma_result
print.powernma_result <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  powerNMA Analysis Results\n")
  cat("═══════════════════════════════════════════════════\n\n")

  # Display mode information
  mode <- x$mode %||% "unknown"
  cat("Mode:", toupper(mode), "\n")

  if (mode == "standard") {
    cat("Status: \u2713 VALIDATED - Suitable for systematic reviews\n")
    cat("Based on:", x$based_on %||% "netmeta, gemtc", "\n")
  } else if (mode == "experimental") {
    cat("Status: \u26a0 EXPERIMENTAL - Research use only\n")
    cat("Validation: See VALIDATION_PLAN.md\n")
  }

  cat("\n")
  cat("Reference treatment:", x$ref_treatment %||% "Not specified", "\n")

  # Show method-specific information
  if (mode == "standard") {
    # Standard mode: show standard NMA info
    if (!is.null(x$network)) {
      cat(sprintf("Heterogeneity \u03c4: %.4f\n", x$network$tau))
      cat(sprintf("I\u00b2: %.1f%%\n", x$network$I2 * 100))
    }
  } else {
    # Experimental mode: show warnings if present
    if (!is.null(x$warnings) && length(x$warnings) > 0) {
      cat("\nWarnings:\n")
      for (w in x$warnings[1:min(3, length(x$warnings))]) {
        cat("  \u2022", w, "\n")
      }
      if (length(x$warnings) > 3) {
        cat("  ... and", length(x$warnings) - 3, "more\n")
      }
    }
  }

  # Show completed analyses
  cat("\nAnalyses completed:\n")

  if (mode == "standard") {
    components <- list(
      network = "Frequentist NMA (netmeta)",
      bayesian = "Bayesian NMA (gemtc)",
      sensitivity = "Leave-one-out sensitivity",
      geometry = "Network geometry",
      inconsistency = "Inconsistency assessment"
    )
  } else {
    components <- list(
      network = "Standard NMA",
      rmst = "Time-varying RMST NMA",
      milestone = "Milestone Survival NMA",
      transportability = "Transportability weighting",
      publication_bias = "PET-PEESE bias assessment",
      loto = "Leave-one-treatment-out",
      multiverse = "Multiverse robustness"
    )
  }

  for (comp_name in names(components)) {
    status <- if (!is.null(x[[comp_name]])) "\u2713" else "\u2717"
    cat(sprintf("  %s %s\n", status, components[[comp_name]]))
  }

  cat("\n")
  if (mode == "experimental") {
    cat("\u26a0 REMINDER: Experimental methods not validated for clinical use\n")
  }

  cat("\nUse summary(x) for detailed results\n")
  invisible(x)
}

#' Summary method for powerNMA results
#' @param object A `powernma_result` object
#' @param ... Additional arguments
#' @export
#' @method summary powernma_result
summary.powernma_result <- function(object, ...) {
  print(object, ...)

  if (!is.null(object$results$summary)) {
    cat("\n--- Top 5 Treatments (P-score) ---\n")
    print(object$results$summary$top5_pscore)
  }

  invisible(object)
}
