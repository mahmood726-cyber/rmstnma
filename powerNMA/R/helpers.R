# ============================================================================
# Helper Functions for Common Data Transformations
# ============================================================================

#' Convert Wide Format to Pairwise Format
#'
#' Transforms data from wide format (one row per study) to pairwise format
#' for network meta-analysis.
#'
#' @param data Data frame in wide format
#' @param study_col Name of study identifier column
#' @param treat_cols Vector of treatment column names
#' @param outcome_cols Vector of outcome column names (corresponding to treat_cols)
#' @param se_cols Vector of standard error column names
#' @return Data frame in pairwise format
#' @export
#' @examples
#' \dontrun{
#' # Wide format data
#' wide_data <- data.frame(
#'   study = c("S1", "S2"),
#'   treat_a = c("Placebo", "Placebo"),
#'   outcome_a = c(0.5, 0.6),
#'   se_a = c(0.1, 0.12),
#'   treat_b = c("DrugA", "DrugB"),
#'   outcome_b = c(0.3, 0.4),
#'   se_b = c(0.11, 0.13)
#' )
#'
#' pairwise_data <- wide_to_pairwise(
#'   wide_data,
#'   study_col = "study",
#'   treat_cols = c("treat_a", "treat_b"),
#'   outcome_cols = c("outcome_a", "outcome_b"),
#'   se_cols = c("se_a", "se_b")
#' )
#' }
wide_to_pairwise <- function(data, study_col, treat_cols, outcome_cols, se_cols) {
  # Validate inputs
  assert_data_frame(data, "data", "wide_to_pairwise")
  assert_columns_exist(data, c(study_col, treat_cols, outcome_cols, se_cols),
                      "wide_to_pairwise")

  if (length(treat_cols) != length(outcome_cols) ||
      length(treat_cols) != length(se_cols)) {
    stop("[wide_to_pairwise] treat_cols, outcome_cols, and se_cols must have same length",
         call. = FALSE)
  }

  pairwise_list <- list()

  for (i in seq_len(nrow(data))) {
    study_id <- data[[study_col]][i]
    n_arms <- length(treat_cols)

    # Create all pairwise combinations
    for (j in 1:(n_arms - 1)) {
      for (k in (j + 1):n_arms) {
        treat1 <- data[[treat_cols[j]]][i]
        treat2 <- data[[treat_cols[k]]][i]
        outcome1 <- data[[outcome_cols[j]]][i]
        outcome2 <- data[[outcome_cols[k]]][i]
        se1 <- data[[se_cols[j]]][i]
        se2 <- data[[se_cols[k]]][i]

        # Calculate treatment effect and SE
        TE <- outcome2 - outcome1
        seTE <- sqrt(se1^2 + se2^2)

        pairwise_list[[length(pairwise_list) + 1]] <- data.frame(
          studlab = study_id,
          treat1 = treat1,
          treat2 = treat2,
          TE = TE,
          seTE = seTE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, pairwise_list)
}

#' Add Study-Level Covariates
#'
#' Merges study-level covariates with pairwise NMA data.
#'
#' @param pairwise_data Pairwise NMA data
#' @param covariate_data Data frame with study-level covariates
#' @param study_col Name of study identifier column (default: "studlab")
#' @return Merged data frame
#' @export
add_study_covariates <- function(pairwise_data, covariate_data, study_col = "studlab") {
  assert_data_frame(pairwise_data, "pairwise_data", "add_study_covariates")
  assert_data_frame(covariate_data, "covariate_data", "add_study_covariates")

  if (!study_col %in% names(pairwise_data)) {
    stop(sprintf("[add_study_covariates] Column '%s' not found in pairwise_data",
                 study_col), call. = FALSE)
  }

  if (!study_col %in% names(covariate_data)) {
    stop(sprintf("[add_study_covariates] Column '%s' not found in covariate_data",
                 study_col), call. = FALSE)
  }

  merged <- dplyr::left_join(pairwise_data, covariate_data, by = study_col)

  # Check for unmatched studies
  unmatched <- sum(is.na(merged[[names(covariate_data)[2]]]))
  if (unmatched > 0) {
    warning(sprintf("[add_study_covariates] %d comparison(s) have no matching covariates",
                   unmatched), call. = FALSE)
  }

  merged
}

#' Filter Network by Treatment
#'
#' Removes specific treatments from the network.
#'
#' @param data Pairwise NMA data
#' @param exclude_treatments Vector of treatment names to exclude
#' @return Filtered data frame
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 20)
#' filtered <- filter_treatments(data, exclude_treatments = c("DrugD"))
filter_treatments <- function(data, exclude_treatments) {
  assert_data_frame(data, "data", "filter_treatments")
  assert_not_null(exclude_treatments, "exclude_treatments", "filter_treatments")

  original_n <- nrow(data)

  filtered <- data %>%
    dplyr::filter(!(treat1 %in% exclude_treatments | treat2 %in% exclude_treatments))

  removed_n <- original_n - nrow(filtered)

  if (removed_n > 0) {
    msg("Removed %d comparison(s) involving: %s",
        removed_n, paste(exclude_treatments, collapse = ", "))
  }

  # Check if network is still viable
  remaining_treatments <- length(unique(c(filtered$treat1, filtered$treat2)))
  if (remaining_treatments < 2) {
    warning("[filter_treatments] Less than 2 treatments remaining after filtering",
            call. = FALSE)
  }

  filtered
}

#' Standardize Treatment Names
#'
#' Cleans and standardizes treatment names for consistency.
#'
#' @param data Pairwise NMA data
#' @param mapping Optional named vector for custom mappings
#' @param to_lower Convert to lowercase
#' @param trim_whitespace Remove leading/trailing whitespace
#' @return Data frame with standardized treatment names
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 10)
#' data$treat1 <- paste0(" ", data$treat1, " ")  # Add whitespace
#' clean_data <- standardize_treatment_names(data)
standardize_treatment_names <- function(data, mapping = NULL,
                                       to_lower = FALSE, trim_whitespace = TRUE) {
  assert_data_frame(data, "data", "standardize_treatment_names")

  result <- data

  if (trim_whitespace) {
    result$treat1 <- trimws(result$treat1)
    result$treat2 <- trimws(result$treat2)
  }

  if (to_lower) {
    result$treat1 <- tolower(result$treat1)
    result$treat2 <- tolower(result$treat2)
  }

  if (!is.null(mapping)) {
    # Apply custom mappings
    for (old_name in names(mapping)) {
      new_name <- mapping[old_name]
      result$treat1[result$treat1 == old_name] <- new_name
      result$treat2[result$treat2 == old_name] <- new_name
    }
    msg("Applied %d treatment name mapping(s)", length(mapping))
  }

  result
}

#' Split Data by Subgroup
#'
#' Splits NMA data into subgroups for subgroup analysis.
#'
#' @param data Pairwise NMA data with subgroup column
#' @param subgroup_col Name of subgroup column
#' @return Named list of data frames (one per subgroup)
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 30)
#' # Assuming data has a 'study_design' column
#' subgroups <- split_by_subgroup(data, "study_design")
split_by_subgroup <- function(data, subgroup_col) {
  assert_data_frame(data, "data", "split_by_subgroup")
  assert_columns_exist(data, subgroup_col, "split_by_subgroup")

  subgroups <- split(data, data[[subgroup_col]])

  msg("Split data into %d subgroup(s)", length(subgroups))
  for (sg in names(subgroups)) {
    n_studies <- length(unique(subgroups[[sg]]$studlab))
    n_comparisons <- nrow(subgroups[[sg]])
    msg("  • %s: %d studies, %d comparisons", sg, n_studies, n_comparisons)
  }

  subgroups
}

#' Create Comparison Matrix
#'
#' Creates a matrix showing which treatments are directly compared.
#'
#' @param data Pairwise NMA data
#' @return Matrix with direct comparison counts
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 20)
#' matrix <- create_comparison_matrix(data)
#' print(matrix)
create_comparison_matrix <- function(data) {
  assert_data_frame(data, "data", "create_comparison_matrix")

  treatments <- sort(unique(c(data$treat1, data$treat2)))
  n_treat <- length(treatments)

  mat <- matrix(0, nrow = n_treat, ncol = n_treat,
               dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- data$treat1[i]
    t2 <- data$treat2[i]
    mat[t1, t2] <- mat[t1, t2] + 1
    mat[t2, t1] <- mat[t2, t1] + 1  # Symmetric
  }

  mat
}

#' Get Treatment Pairs
#'
#' Extracts all unique treatment pairs from the network.
#'
#' @param data Pairwise NMA data
#' @param include_indirect Logical; if TRUE, includes indirect comparisons
#' @return Data frame with treatment pairs
#' @export
get_treatment_pairs <- function(data, include_indirect = FALSE) {
  assert_data_frame(data, "data", "get_treatment_pairs")

  # Direct comparisons
  direct <- data %>%
    dplyr::select(treat1, treat2) %>%
    dplyr::distinct() %>%
    dplyr::mutate(comparison_type = "direct")

  if (!include_indirect) {
    return(direct)
  }

  # All possible pairs (including indirect)
  treatments <- unique(c(data$treat1, data$treat2))
  all_pairs <- expand.grid(treat1 = treatments, treat2 = treatments,
                          stringsAsFactors = FALSE) %>%
    dplyr::filter(treat1 < treat2)  # Remove duplicates

  # Mark direct vs indirect
  all_pairs$comparison_type <- "indirect"
  for (i in seq_len(nrow(all_pairs))) {
    t1 <- all_pairs$treat1[i]
    t2 <- all_pairs$treat2[i]

    if (any((data$treat1 == t1 & data$treat2 == t2) |
            (data$treat1 == t2 & data$treat2 == t1))) {
      all_pairs$comparison_type[i] <- "direct"
    }
  }

  all_pairs
}
