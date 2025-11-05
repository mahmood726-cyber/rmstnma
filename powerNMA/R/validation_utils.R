#' Comprehensive Input Validation Utilities
#'
#' @description
#' Production-grade validation functions to ensure data quality and
#' provide clear, actionable error messages with function context.
#'
#' These utilities implement defensive programming practices identified
#' in code review as critical improvements for robustness.
#'
#' @author powerNMA Development Team
#' @name validation_utils
NULL

#' Validate Not Null
#'
#' @description
#' Checks that an object is not NULL and provides informative error message.
#'
#' @param x Object to check
#' @param name Name of the parameter (for error message)
#' @param func Function name calling this validation (optional)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_not_null(data, "data", "run_nma")
#' }
assert_not_null <- function(x, name, func = NULL) {
  if (is.null(x)) {
    func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""
    stop(sprintf("%sParameter '%s' cannot be NULL", func_context, name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate Data Frame
#'
#' @description
#' Validates that input is a data frame with required properties.
#'
#' @param df Object to check
#' @param name Parameter name (for error message)
#' @param func Function name (optional)
#' @param min_rows Minimum number of rows required (default: 1)
#' @param required_cols Character vector of required column names (optional)
#' @param allow_empty Logical, allow empty data frame (default: FALSE)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_data_frame(data, "data", "run_nma",
#'                   required_cols = c("studlab", "treat1", "treat2"))
#' }
assert_data_frame <- function(df, name, func = NULL, min_rows = 1,
                               required_cols = NULL, allow_empty = FALSE) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check if data frame
  if (!is.data.frame(df)) {
    stop(sprintf("%sParameter '%s' must be a data.frame, got %s",
                 func_context, name, paste(class(df), collapse = ", ")),
         call. = FALSE)
  }

  # Check if empty
  if (!allow_empty && nrow(df) == 0) {
    stop(sprintf("%sParameter '%s' is an empty data frame (0 rows)",
                 func_context, name),
         call. = FALSE)
  }

  # Check minimum rows
  if (nrow(df) < min_rows) {
    stop(sprintf("%sParameter '%s' must have at least %d rows, got %d",
                 func_context, name, min_rows, nrow(df)),
         call. = FALSE)
  }

  # Check required columns
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(sprintf("%sParameter '%s' missing required columns: %s",
                   func_context, name, paste(missing_cols, collapse = ", ")),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate Numeric
#'
#' @description
#' Validates numeric inputs with range checking.
#'
#' @param x Object to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param min Minimum allowed value (optional)
#' @param max Maximum allowed value (optional)
#' @param allow_na Allow NA values (default: FALSE)
#' @param allow_inf Allow infinite values (default: FALSE)
#' @param integer_only Require integer values (default: FALSE)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_numeric(n_studies, "n_studies", "simulate_nma_data",
#'                min = 1, integer_only = TRUE)
#' }
assert_numeric <- function(x, name, func = NULL, min = NULL, max = NULL,
                           allow_na = FALSE, allow_inf = FALSE,
                           integer_only = FALSE) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check if numeric
  if (!is.numeric(x)) {
    stop(sprintf("%sParameter '%s' must be numeric, got %s",
                 func_context, name, paste(class(x), collapse = ", ")),
         call. = FALSE)
  }

  # Check NA
  if (!allow_na && any(is.na(x))) {
    stop(sprintf("%sParameter '%s' contains NA values (not allowed)",
                 func_context, name),
         call. = FALSE)
  }

  # Check infinite
  if (!allow_inf && any(is.infinite(x[!is.na(x)]))) {
    stop(sprintf("%sParameter '%s' contains infinite values (not allowed)",
                 func_context, name),
         call. = FALSE)
  }

  # Check min
  if (!is.null(min)) {
    if (any(x[!is.na(x)] < min)) {
      stop(sprintf("%sParameter '%s' must be >= %s, got min value %s",
                   func_context, name, format(min), format(min(x, na.rm = TRUE))),
           call. = FALSE)
    }
  }

  # Check max
  if (!is.null(max)) {
    if (any(x[!is.na(x)] > max)) {
      stop(sprintf("%sParameter '%s' must be <= %s, got max value %s",
                   func_context, name, format(max), format(max(x, na.rm = TRUE))),
           call. = FALSE)
    }
  }

  # Check integer
  if (integer_only) {
    if (!all(x[!is.na(x)] == as.integer(x[!is.na(x)]))) {
      stop(sprintf("%sParameter '%s' must be integer values",
                   func_context, name),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate Character
#'
#' @description
#' Validates character inputs with choices checking.
#'
#' @param x Object to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param choices Allowed values (optional)
#' @param min_length Minimum string length (default: 1)
#' @param allow_empty Allow empty strings (default: FALSE)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_character(model_type, "model_type", "run_nma",
#'                  choices = c("random", "fixed", "common"))
#' }
assert_character <- function(x, name, func = NULL, choices = NULL,
                             min_length = 1, allow_empty = FALSE) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check if character
  if (!is.character(x)) {
    stop(sprintf("%sParameter '%s' must be character, got %s",
                 func_context, name, paste(class(x), collapse = ", ")),
         call. = FALSE)
  }

  # Check empty
  if (!allow_empty && any(nchar(x) == 0)) {
    stop(sprintf("%sParameter '%s' contains empty strings (not allowed)",
                 func_context, name),
         call. = FALSE)
  }

  # Check minimum length
  if (any(nchar(x) < min_length)) {
    stop(sprintf("%sParameter '%s' must have at least %d characters",
                 func_context, name, min_length),
         call. = FALSE)
  }

  # Check choices
  if (!is.null(choices)) {
    invalid <- setdiff(x, choices)
    if (length(invalid) > 0) {
      stop(sprintf("%sParameter '%s' has invalid values: %s. Allowed: %s",
                   func_context, name,
                   paste(invalid, collapse = ", "),
                   paste(choices, collapse = ", ")),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate Logical
#'
#' @description
#' Validates logical inputs.
#'
#' @param x Object to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param length Expected length (optional)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_logical(use_bayesian, "use_bayesian", "run_nma", length = 1)
#' }
assert_logical <- function(x, name, func = NULL, length = NULL) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check if logical
  if (!is.logical(x)) {
    stop(sprintf("%sParameter '%s' must be logical (TRUE/FALSE), got %s",
                 func_context, name, paste(class(x), collapse = ", ")),
         call. = FALSE)
  }

  # Check length
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("%sParameter '%s' must have length %d, got %d",
                 func_context, name, length, length(x)),
           call. = FALSE)
  }

  # Check NA
  if (any(is.na(x))) {
    stop(sprintf("%sParameter '%s' contains NA values (not allowed for logical)",
                 func_context, name),
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate NMA Object
#'
#' @description
#' Validates that input is a valid NMA object.
#'
#' @param obj Object to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param require_class Required class(es)
#' @param require_components Required list components
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_nma_object(nma_result, "nma_result", "calculate_sucra",
#'                   require_class = c("netmeta", "gemtc"))
#' }
assert_nma_object <- function(obj, name, func = NULL,
                               require_class = c("netmeta", "gemtc", "auto_nma"),
                               require_components = NULL) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check class
  if (!inherits(obj, require_class)) {
    stop(sprintf("%sParameter '%s' must inherit from one of: %s. Got: %s",
                 func_context, name,
                 paste(require_class, collapse = ", "),
                 paste(class(obj), collapse = ", ")),
         call. = FALSE)
  }

  # Check required components
  if (!is.null(require_components)) {
    if (!is.list(obj)) {
      stop(sprintf("%sParameter '%s' must be a list with components",
                   func_context, name),
           call. = FALSE)
    }

    missing_components <- setdiff(require_components, names(obj))
    if (length(missing_components) > 0) {
      stop(sprintf("%sParameter '%s' missing required components: %s",
                   func_context, name, paste(missing_components, collapse = ", ")),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate Positive Integer
#'
#' @description
#' Convenience function for validating positive integers (common case).
#'
#' @param x Value to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param min Minimum value (default: 1)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_positive_integer(n_studies, "n_studies", "simulate_nma_data")
#' }
assert_positive_integer <- function(x, name, func = NULL, min = 1) {
  assert_not_null(x, name, func)
  assert_numeric(x, name, func, min = min, integer_only = TRUE, allow_na = FALSE)
  invisible(TRUE)
}

#' Validate Probability
#'
#' @description
#' Validates values are probabilities (between 0 and 1).
#'
#' @param x Value(s) to check
#' @param name Parameter name
#' @param func Function name (optional)
#' @param allow_zero Allow exactly 0 (default: TRUE)
#' @param allow_one Allow exactly 1 (default: TRUE)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_probability(alpha, "alpha", "hypothesis_test")
#' }
assert_probability <- function(x, name, func = NULL,
                               allow_zero = TRUE, allow_one = TRUE) {
  min_val <- if (allow_zero) 0 else .Machine$double.eps
  max_val <- if (allow_one) 1 else (1 - .Machine$double.eps)

  assert_numeric(x, name, func, min = min_val, max = max_val, allow_na = FALSE)
  invisible(TRUE)
}

#' Validate Seed
#'
#' @description
#' Validates random seed parameter.
#'
#' @param seed Seed value (numeric or NULL)
#' @param name Parameter name (default: "seed")
#' @param func Function name (optional)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_seed(seed, func = "simulate_nma_data")
#' }
assert_seed <- function(seed, name = "seed", func = NULL) {
  if (!is.null(seed)) {
    assert_numeric(seed, name, func, integer_only = TRUE, allow_na = FALSE)
  }
  invisible(TRUE)
}

#' Safe Try-Catch Wrapper
#'
#' @description
#' Executes an expression with error handling and context.
#'
#' @param expr Expression to evaluate
#' @param func Function name for error context
#' @param error_handler Function to handle errors (optional)
#' @param warning_handler Function to handle warnings (optional)
#' @param silent Suppress messages (default: FALSE)
#'
#' @return Result of expr, or error object if failed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' result <- safe_try({
#'   risky_operation()
#' }, func = "my_function")
#' }
safe_try <- function(expr, func = NULL, error_handler = NULL,
                     warning_handler = NULL, silent = FALSE) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  tryCatch(
    {
      if (!silent) {
        withCallingHandlers(
          expr,
          warning = function(w) {
            if (!is.null(warning_handler)) {
              warning_handler(w)
            } else {
              message(sprintf("%sWarning: %s", func_context, conditionMessage(w)))
            }
            invokeRestart("muffleWarning")
          }
        )
      } else {
        suppressWarnings(expr)
      }
    },
    error = function(e) {
      if (!is.null(error_handler)) {
        error_handler(e)
      } else {
        stop(sprintf("%sError: %s", func_context, conditionMessage(e)),
             call. = FALSE)
      }
    }
  )
}

#' Validate Network Structure
#'
#' @description
#' Validates basic network structure properties.
#'
#' @param data Network data (data frame)
#' @param treat1_col Name of first treatment column
#' @param treat2_col Name of second treatment column
#' @param studlab_col Name of study label column
#' @param func Function name (optional)
#' @param min_treatments Minimum number of unique treatments (default: 2)
#' @param min_studies Minimum number of unique studies (default: 2)
#'
#' @return Invisible TRUE if valid
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' assert_network_structure(data, "treat1", "treat2", "studlab")
#' }
assert_network_structure <- function(data, treat1_col = "treat1", treat2_col = "treat2",
                                     studlab_col = "studlab", func = NULL,
                                     min_treatments = 2, min_studies = 2) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""

  # Check data frame
  assert_data_frame(data, "data", func,
                    required_cols = c(treat1_col, treat2_col, studlab_col))

  # Check number of unique treatments
  treatments <- unique(c(data[[treat1_col]], data[[treat2_col]]))
  n_treatments <- length(treatments)

  if (n_treatments < min_treatments) {
    stop(sprintf("%sNetwork must have at least %d treatments, got %d",
                 func_context, min_treatments, n_treatments),
         call. = FALSE)
  }

  # Check number of unique studies
  studies <- unique(data[[studlab_col]])
  n_studies <- length(studies)

  if (n_studies < min_studies) {
    stop(sprintf("%sNetwork must have at least %d studies, got %d",
                 func_context, min_studies, n_studies),
         call. = FALSE)
  }

  # Check for self-loops (treatment compared to itself)
  self_loops <- data[[treat1_col]] == data[[treat2_col]]
  if (any(self_loops, na.rm = TRUE)) {
    n_loops <- sum(self_loops, na.rm = TRUE)
    warning(sprintf("%sFound %d self-loops (treatment compared to itself)",
                    func_context, n_loops),
            call. = FALSE)
  }

  invisible(TRUE)
}

#' Create Standardized Error
#'
#' @description
#' Creates error with standard format including function context.
#'
#' @param message Error message
#' @param func Function name
#' @param class Error class (default: "validation_error")
#'
#' @return Error condition
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' stop(create_error("Invalid input", "my_function"))
#' }
create_error <- function(message, func = NULL,
                        class = c("validation_error", "error", "condition")) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""
  full_message <- paste0(func_context, message)

  structure(
    list(message = full_message, call = sys.call(-1)),
    class = class
  )
}

#' Create Standardized Warning
#'
#' @description
#' Creates warning with standard format including function context.
#'
#' @param message Warning message
#' @param func Function name
#'
#' @return Warning condition
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' warning(create_warning("Unusual value detected", "my_function"))
#' }
create_warning <- function(message, func = NULL) {
  func_context <- if (!is.null(func)) sprintf("[%s] ", func) else ""
  full_message <- paste0(func_context, message)

  structure(
    list(message = full_message, call = sys.call(-1)),
    class = c("validation_warning", "warning", "condition")
  )
}
