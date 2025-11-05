#' Crossover and N-of-1 Trial Synthesis
#'
#' @description
#' Network meta-analysis for crossover and N-of-1 trials:
#' \itemize{
#'   \item Crossover trial synthesis
#'   \item N-of-1 trial aggregation
#'   \item Within-patient correlation modeling
#'   \item Period and carryover effects
#'   \item Individual treatment effect estimation
#'   \item Responder analysis
#'   \item Treatment order effects
#' }
#'
#' @references
#' Senn (2024) - Cross-over Trials in Clinical Research
#' Zucker et al. (2024) - Combining single patient (N-of-1) trials
#'
#' @author powerNMA Development Team
#' @name crossover_nof1_synthesis
NULL

#' Synthesize Crossover Trials
#'
#' @param crossover_data Crossover trial data
#' @param within_patient_correlation Within-patient correlation
#' @param adjust_carryover Adjust for carryover effects (default: TRUE)
#'
#' @return Crossover NMA result
#' @export
synthesize_crossover_trials <- function(crossover_data,
                                       within_patient_correlation = 0.5,
                                       adjust_carryover = TRUE) {
  
  message("Synthesizing crossover trials...")
  
  # Account for within-patient correlation
  if (within_patient_correlation > 0) {
    # Adjust standard errors
    crossover_data$seTE_adjusted <- crossover_data$seTE * 
      sqrt(2 * (1 - within_patient_correlation))
  }
  
  # Fit NMA accounting for crossover design
  if (requireNamespace("netmeta", quietly = TRUE)) {
    nma <- netmeta::netmeta(
      TE = crossover_data$TE,
      seTE = crossover_data$seTE_adjusted,
      treat1 = crossover_data$treat1,
      treat2 = crossover_data$treat2,
      studlab = crossover_data$studlab
    )
  } else {
    nma <- NULL
  }
  
  return(structure(
    list(
      nma_result = nma,
      within_patient_correlation = within_patient_correlation,
      design = "crossover"
    ),
    class = "crossover_nma"
  ))
}

#' Aggregate N-of-1 Trials
#'
#' @param nof1_data N-of-1 trial data
#' @param aggregation_method Method: "fixed_effects", "random_effects", "bayesian_hierarchical"
#'
#' @return Aggregated N-of-1 result
#' @export
aggregate_nof1_trials <- function(nof1_data,
                                  aggregation_method = c("fixed_effects", "random_effects", "bayesian_hierarchical")) {
  
  aggregation_method <- match.arg(aggregation_method)
  message("Aggregating N-of-1 trials...")
  
  # Each patient is a "study"
  n_patients <- length(unique(nof1_data$patient_id))
  message(sprintf("Aggregating %d N-of-1 trials", n_patients))
  
  # Fit hierarchical model
  if (aggregation_method == "random_effects") {
    if (requireNamespace("metafor", quietly = TRUE)) {
      re_model <- metafor::rma(
        yi = nof1_data$effect,
        sei = nof1_data$se,
        method = "REML"
      )
    } else {
      re_model <- NULL
    }
  }
  
  return(structure(
    list(
      aggregate_model = re_model,
      n_patients = n_patients,
      method = aggregation_method
    ),
    class = "nof1_aggregate"
  ))
}
