# R/data.R
# Documentation for package datasets

#' Example Network Meta-Analysis Dataset
#'
#' A synthetic dataset containing Kaplan-Meier survival data from three studies
#' comparing three treatments. This dataset demonstrates the expected format
#' for RMST network meta-analysis.
#'
#' @format A data frame with 113 rows and 5 columns:
#' \describe{
#'   \item{study_id}{Character. Study identifier (STUDY001, STUDY002, STUDY003)}
#'   \item{treatment}{Character. Treatment name (Control, Treatment A, Treatment B)}
#'   \item{time}{Numeric. Time points in months (0 to 36)}
#'   \item{survival}{Numeric. Survival probability at each time point (0 to 1)}
#'   \item{source}{Character. Data source type ("km" for Kaplan-Meier)}
#' }
#'
#' @details
#' The dataset contains:
#' \itemize{
#'   \item STUDY001: Control vs Treatment A
#'   \item STUDY002: Control vs Treatment B
#'   \item STUDY003: Treatment A vs Treatment B
#' }
#'
#' Treatment effects simulate realistic scenarios:
#' \itemize{
#'   \item Control: Baseline survival
#'   \item Treatment A: 30\% hazard reduction vs Control
#'   \item Treatment B: 50\% hazard reduction vs Control
#' }
#'
#' @source Generated synthetically for package examples
#'
#' @examples
#' data(example_network)
#'
#' # View structure
#' str(example_network)
#'
#' # Summary by study and treatment
#' table(example_network$study_id, example_network$treatment)
#'
#' # Plot survival curves (if ggplot2 available)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(example_network, aes(x = time, y = survival,
#'                               color = treatment, linetype = study_id)) +
#'     geom_line() +
#'     labs(title = "Example Survival Curves",
#'          x = "Time (months)",
#'          y = "Survival Probability") +
#'     theme_minimal()
#' }
#'
#' # Create RMST network
#' net <- rmst_network(
#'   data = example_network,
#'   study = "study_id",
#'   trt = "treatment",
#'   time = "time",
#'   surv = "survival"
#' )
#' print(net)
"example_network"

