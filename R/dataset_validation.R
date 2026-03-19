#' NMA Dataset Schema and Validation
#'
#' Utilities to standardize and validate survival NMA datasets so they share a
#' common structure across studies and networks.
#'
#' The standard (minimal) column set is:
#' - `study_id` (character)
#' - `treatment` (character)
#' - `time` (numeric, months)
#' - `survival` (numeric, in [0, 1])
#'
#' Optional but recommended columns:
#' - `source` (character, e.g. "km", "ipdfromkm")
#' - `arm_id` (character or integer, per-study arm label)
#' - `n_risk` (numeric, number at risk at time)
#'
#' @name nma_dataset_validation
NULL

#' Return the canonical NMA dataset schema
#'
#' @return A list with `required` and `optional` named vectors of types
#' @export
nma_schema <- function() {
  list(
    required = c(
      study_id = "character",
      treatment = "character",
      time = "numeric",
      survival = "numeric"
    ),
    optional = c(
      source = "character",
      arm_id = "character",
      n_risk = "numeric"
    )
  )
}

#' Standardize common column names to canonical names
#'
#' Attempts to map a variety of common alternative column names to the
#' standardized set expected by this package.
#'
#' @param x A data.frame or tibble
#' @param mappings Optional named list overriding default synonyms
#' @return The input with columns renamed where possible
#' @export
standardize_nma_columns <- function(x, mappings = NULL) {
  stopifnot(is.data.frame(x))

  default <- list(
    study_id = c("study_id", "study", "trial", "studyid"),
    treatment = c("treatment", "trt", "arm", "tx", "treat"),
    time = c("time", "t", "month", "months", "time_months"),
    survival = c("survival", "S", "surv", "prob_surv", "km")
  )

  if (!is.null(mappings)) {
    for (nm in names(mappings)) {
      default[[nm]] <- unique(c(mappings[[nm]], default[[nm]]))
    }
  }

  # Build rename vector: first match wins
  rename_map <- c()
  current <- names(x)
  for (target in names(default)) {
    alts <- default[[target]]
    hit <- intersect(alts, current)
    if (length(hit) > 0 && target %in% current) {
      next
    }
    if (length(hit) > 0) {
      rename_map[hit[1]] <- target
    }
  }

  dplyr::rename(x, !!!rlang::syms(names(rename_map)) := !!!rlang::syms(unname(rename_map)))
}

#' Validate a survival NMA dataset against the schema
#'
#' @param x A data.frame with (at least) required columns
#' @param strict If TRUE, apply stricter checks (monotonic survival, time 0==1)
#' @param quiet If TRUE, suppress message output
#' @return A list with `ok` (logical), `issues` (character vector), and
#'   `summary` (list of quick stats)
#' @export
validate_nma_dataset <- function(x, strict = FALSE, quiet = FALSE) {
  issues <- character()
  schema <- nma_schema()

  if (!is.data.frame(x)) {
    return(list(ok = FALSE, issues = "Input is not a data.frame", summary = NULL))
  }

  # Presence
  miss <- setdiff(names(schema$required), names(x))
  if (length(miss) > 0) issues <- c(issues, paste0("Missing columns: ", paste(miss, collapse = ", ")))

  # Types (best-effort check)
  for (nm in intersect(names(schema$required), names(x))) {
    want <- schema$required[[nm]]
    got <- class(x[[nm]])[1]
    if (want == "numeric" && !is.numeric(x[[nm]])) issues <- c(issues, paste0("Column `", nm, "` must be numeric"))
    if (want == "character" && !is.character(x[[nm]])) issues <- c(issues, paste0("Column `", nm, "` must be character"))
  }

  # Basic value constraints
  if ("time" %in% names(x)) {
    if (any(!is.finite(x$time))) issues <- c(issues, "`time` contains non-finite values")
    if (any(x$time < 0, na.rm = TRUE)) issues <- c(issues, "`time` contains negative values")
  }
  if ("survival" %in% names(x)) {
    if (any(!is.finite(x$survival))) issues <- c(issues, "`survival` contains non-finite values")
    if (any(x$survival < 0 | x$survival > 1, na.rm = TRUE)) issues <- c(issues, "`survival` not in [0,1]")
  }

  # Uniqueness of (study, treatment, time)
  if (all(c("study_id", "treatment", "time") %in% names(x))) {
    key <- paste(x$study_id, x$treatment, x$time, sep = "\r")
    dups <- any(duplicated(key))
    if (dups) issues <- c(issues, "Duplicate (study_id, treatment, time) rows present")
  }

  # Network connectivity sanity: at least two treatments per study
  if (all(c("study_id", "treatment") %in% names(x))) {
    trt_per_study <- tapply(x$treatment, x$study_id, function(v) length(unique(v)))
    if (any(trt_per_study < 2)) issues <- c(issues, "Some studies have <2 treatments")
  }

  # Strict checks: monotone non-increasing survival per study/arm, time 0==1
  if (strict && all(c("study_id", "treatment", "time", "survival") %in% names(x))) {
    .chk <- function(df) {
      o <- order(df$time, na.last = NA)
      s <- df$survival[o]
      bad <- any(diff(s) > 1e-6, na.rm = TRUE)
      bad
    }
    bad_groups <- by(x, interaction(x$study_id, x$treatment, drop = TRUE), .chk)
    if (any(unlist(bad_groups))) issues <- c(issues, "Survival increases over time in some arms (strict)")

    # time zero survival ~ 1
    .t0 <- by(x, interaction(x$study_id, x$treatment, drop = TRUE), function(df) {
      i <- which.min(abs(df$time))
      if (length(i) == 0) return(TRUE)
      abs(df$time[i]) <= 1e-8 && abs(df$survival[i] - 1) <= 1e-2
    })
    if (!all(unlist(.t0))) issues <- c(issues, "Some arms lack survival≈1 at time≈0 (strict)")
  }

  ok <- length(issues) == 0

  # Quick summary
  smry <- NULL
  if (all(c("study_id", "treatment", "time") %in% names(x))) {
    smry <- list(
      n_rows = nrow(x),
      n_studies = length(unique(x$study_id)),
      n_treatments = length(unique(x$treatment)),
      time_range = if ("time" %in% names(x)) range(x$time, na.rm = TRUE) else NULL
    )
  }

  if (!quiet) {
    if (ok) message("NMA dataset validation: OK") else message("NMA dataset validation: issues found")
  }

  list(ok = ok, issues = issues, summary = smry)
}

#' List available packaged NMA datasets
#'
#' Lists `.rda` datasets shipped in the package `data/` directory that match
#' the expected format.
#'
#' @return A character vector of dataset object names
#' @export
list_nma_datasets <- function() {
  pkg_path <- system.file(package = utils::packageName())
  # During development, fall back to working directory
  data_dir <- if (nzchar(pkg_path)) file.path(pkg_path, "data") else "data"
  if (!dir.exists(data_dir)) return(character())
  files <- list.files(data_dir, pattern = "\\.rda$", full.names = FALSE)
  sub("\\.rda$", "", files)
}

