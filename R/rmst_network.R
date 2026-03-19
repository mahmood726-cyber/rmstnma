#' Create RMST Network Object
#'
#' Builds a network meta-analysis dataset for RMST analysis from survival data
#'
#' @param data Data frame containing survival data
#' @param study Column name for study identifier
#' @param trt Column name for treatment identifier
#' @param time Column name for time points
#' @param surv Column name for survival probabilities
#' @param events Column name for number of events (optional)
#' @param at_risk Column name for number at risk (optional)
#' @param source Character vector specifying data source type per arm
#' @param ipd_data Individual patient data (optional)
#'
#' @return An object of class 'rmst_network'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename arrange
#' @importFrom rlang sym
#' @importFrom checkmate assert_data_frame assert_character
#'
#' @examples
#' data(example_network)
#' net <- rmst_network(
#'   data = example_network,
#'   study = "study_id",
#'   trt = "treatment",
#'   time = "time",
#'   surv = "survival"
#' )
rmst_network <- function(data,
                         study,
                         trt,
                         time,
                         surv,
                         events = NULL,
                         at_risk = NULL,
                         source = "km",
                         ipd_data = NULL) {

  # Input validation using checkmate
  checkmate::assert_data_frame(data)
  checkmate::assert_character(study)
  checkmate::assert_character(trt)
  checkmate::assert_character(time)
  checkmate::assert_character(surv)

  # Convert to tidyeval symbols
  study_sym <- rlang::sym(study)
  trt_sym <- rlang::sym(trt)
  time_sym <- rlang::sym(time)
  surv_sym <- rlang::sym(surv)

  # Process the network data
  network_data <- data %>%
    dplyr::rename(
      study = !!study_sym,
      treatment = !!trt_sym,
      time = !!time_sym,
      survival = !!surv_sym
    ) %>%
    dplyr::arrange(study, treatment, time)

  # Add optional columns if provided
  if (!is.null(events)) {
    network_data$events <- data[[events]]
  }

  if (!is.null(at_risk)) {
    network_data$at_risk <- data[[at_risk]]
  }

  # Tag source type
  if (length(source) == 1) {
    network_data$source <- source
  } else {
    checkmate::assert_character(source, len = nrow(network_data))
    network_data$source <- source
  }

  # Calculate network properties
  studies <- unique(network_data$study)
  treatments <- unique(network_data$treatment)
  n_studies <- length(studies)
  n_treatments <- length(treatments)

  # Create adjacency matrix
  adj_matrix <- matrix(0, n_treatments, n_treatments,
                       dimnames = list(treatments, treatments))

  for (s in studies) {
    study_trts <- unique(network_data$treatment[network_data$study == s])
    if (length(study_trts) > 1) {
      for (i in 1:(length(study_trts)-1)) {
        for (j in (i+1):length(study_trts)) {
          adj_matrix[study_trts[i], study_trts[j]] <-
            adj_matrix[study_trts[i], study_trts[j]] + 1
          adj_matrix[study_trts[j], study_trts[i]] <-
            adj_matrix[study_trts[j], study_trts[i]] + 1
        }
      }
    }
  }

  # Create network object
  network <- list(
    data = network_data,
    studies = studies,
    treatments = treatments,
    n_studies = n_studies,
    n_treatments = n_treatments,
    adjacency = adj_matrix,
    ipd_data = ipd_data,
    call = match.call()
  )

  class(network) <- c("rmst_network", "list")

  # Print summary
  message("RMST Network created:")
  message(sprintf("  - %d studies", n_studies))
  message(sprintf("  - %d treatments", n_treatments))
  message(sprintf("  - Data sources: %s",
                  paste(unique(network_data$source), collapse = ", ")))

  return(network)
}

#' Print method for rmst_network
#'
#' @param x An rmst_network object
#' @param ... Additional arguments (ignored)
#' @export
print.rmst_network <- function(x, ...) {
  cat("RMST Network Meta-Analysis Dataset\n")
  cat("-----------------------------------\n")
  cat(sprintf("Studies: %d\n", x$n_studies))
  cat(sprintf("Treatments: %d\n", x$n_treatments))
  cat(sprintf("Total observations: %d\n", nrow(x$data)))
  cat("\nTreatments:\n")
  cat(paste(" -", x$treatments), sep = "\n")
  invisible(x)
}
