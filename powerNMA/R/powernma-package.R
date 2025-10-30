#' powerNMA: Comprehensive Network Meta-Analysis Suite
#'
#' @description
#' A comprehensive framework for network meta-analysis combining time-varying
#' methods (RMST and milestone survival) with advanced NMA techniques including
#' transportability weighting, Bayesian inference, meta-regression with splines,
#' publication bias assessment, sensitivity analyses, and machine learning for
#' heterogeneity exploration.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{run_powernma}}: Execute comprehensive NMA pipeline
#'   \item \code{\link{rmst_nma}}: Time-varying RMST analysis
#'   \item \code{\link{milestone_nma}}: Milestone survival analysis
#'   \item \code{\link{setup_powernma}}: Configure analysis parameters
#' }
#'
#' @section Data Functions:
#' \itemize{
#'   \item \code{\link{simulate_nma_data}}: Generate example pairwise data
#'   \item \code{\link{generate_example_ipd}}: Generate example IPD
#'   \item \code{\link{reconstruct_ipd}}: Reconstruct IPD from KM curves
#'   \item \code{\link{validate_ipd}}: Validate IPD structure
#'   \item \code{\link{validate_nma_input}}: Validate pairwise data
#' }
#'
#' @section Statistical Methods:
#'
#' **Standard Network Meta-Analysis**
#'
#' Uses the \pkg{netmeta} package for frequentist random-effects and
#' fixed-effects models. Supports multiple effect measures (HR, OR, RR, MD, SMD).
#'
#' **Time-Varying Methods**
#'
#' \describe{
#'   \item{RMST}{Restricted Mean Survival Time analysis compares average
#'   event-free time up to specific time horizons. Provides clinically
#'   interpretable measures of treatment benefit.}
#'   \item{Milestone}{Analyzes odds of survival at specific time points
#'   (e.g., 90-day, 1-year survival). Includes continuity correction for
#'   zero-event handling.}
#' }
#'
#' **Bayesian Inference**
#'
#' Optional Bayesian NMA using \pkg{gemtc} with JAGS backend. Provides
#' posterior distributions, credible intervals, and probabilistic statements
#' about treatment rankings.
#'
#' **Sensitivity Analyses**
#'
#' \itemize{
#'   \item Leave-one-study-out (LOO)
#'   \item Leave-one-treatment-out (LOTO)
#'   \item Global inconsistency (design-by-treatment)
#'   \item Local inconsistency (node-splitting)
#' }
#'
#' @section Advanced Features:
#'
#' **Transportability Weighting**
#'
#' Adjusts study weights to reflect a target population using covariate
#' distances. Supports Mahalanobis and Euclidean metrics with Gaussian
#' or tricube kernels.
#'
#' **Meta-Regression**
#'
#' Network meta-regression with:
#' \itemize{
#'   \item Continuous and categorical moderators
#'   \item Spline smoothing with cross-validation
#'   \item Cluster-robust standard errors (CR2)
#' }
#'
#' **Publication Bias**
#'
#' Multiple approaches:
#' \itemize{
#'   \item PET-PEESE regression
#'   \item Selection models (weightr)
#'   \item Copas sensitivity analysis
#'   \item Trim-and-fill
#' }
#'
#' @section Workflow:
#'
#' \preformatted{
#' # 1. Prepare data
#' data <- simulate_nma_data(n_studies = 40)
#' # or
#' ipd <- generate_example_ipd(n_trials = 10)
#'
#' # 2. Configure analysis
#' config <- setup_powernma(
#'   sm = "HR",
#'   use_bayesian = TRUE,
#'   use_timevarying = TRUE,  # for IPD
#'   export_results = TRUE
#' )
#'
#' # 3. Run analysis
#' results <- run_powernma(
#'   data = data,
#'   data_type = "pairwise",  # or "ipd"
#'   config = config
#' )
#'
#' # 4. Explore results
#' print(results)
#' summary(results)
#' }
#'
#' @section References:
#'
#' **Time-Varying Methods:**
#' \itemize{
#'   \item Guyot P, Ades AE, Ouwens MJ, Welton NJ (2012). Enhanced secondary
#'   analysis of survival data: reconstructing the data from published
#'   Kaplan-Meier survival curves. \emph{BMC Med Res Methodol} 12:9.
#'
#'   \item Wei Y, Royston P (2017). Reconstructing time-to-event data from
#'   published Kaplan-Meier curves. \emph{Stata J} 17(4):786-802.
#' }
#'
#' **Network Meta-Analysis:**
#' \itemize{
#'   \item Rücker G, Schwarzer G (2015). Ranking treatments in frequentist
#'   network meta-analysis works without resampling methods.
#'   \emph{BMC Med Res Methodol} 15:58.
#'
#'   \item Dias S, Welton NJ, Caldwell DM, Ades AE (2010). Checking consistency
#'   in mixed treatment comparison meta-analysis.
#'   \emph{Stat Med} 29(7-8):932-944.
#' }
#'
#' **Transportability:**
#' \itemize{
#'   \item Dahabreh IJ, Robertson SE, Hernán MA (2020). Extending inferences
#'   from a randomized trial to a target population.
#'   \emph{Eur J Epidemiol} 35:719-722.
#' }
#'
#' @docType package
#' @name powerNMA-package
#' @aliases powerNMA
#' @keywords package
NULL

# Suppress R CMD check notes
#' @importFrom stats approx rbinom rexp rnorm runif
#' @importFrom utils combn
NULL
