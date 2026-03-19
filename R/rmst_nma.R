#'#' Fit RMST Network Meta-Analysis
#'
#' Fits a Bayesian network meta-analysis model for RMST outcomes
#'
#' @param network An rmst_network object
#' @param tau Restriction time(s) for RMST calculation
#' @param baseline Baseline hazard model ("rp_spline" or "piecewise")
#' @param random_effects Random effects distribution ("normal" or "student_t")
#' @param inconsistency Include inconsistency parameters
#' @param propagate_recon_uncertainty Propagate reconstruction uncertainty
#' @param recon Reconstruction options
#' @param chains Number of MCMC chains
#' @param iter Number of iterations per chain
#' @param warmup Number of warmup iterations
#' @param adapt_delta Target acceptance probability
#' @param ... Additional arguments to cmdstanr
#'
#' @return An object of class 'rmst_nma_fit'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise pull n_distinct n
#' @importFrom checkmate assert_class assert_numeric assert_choice assert_logical
rmst_nma <- function(network,
                     tau = 12,
                     baseline = "rp_spline",
                     random_effects = "normal",
                     inconsistency = FALSE,
                     propagate_recon_uncertainty = FALSE,
                     recon = NULL,
                     chains = 4,
                     iter = 2000,
                     warmup = NULL,
                     adapt_delta = 0.9,
                     ...) {

  # Check for CmdStan dependencies
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("This function requires the 'cmdstanr' package.\n",
         "Install it with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))")
  }

  # Check CmdStan installation
  tryCatch({
    cmdstan_ver <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
    if (is.null(cmdstan_ver) || is.na(cmdstan_ver)) {
      stop("CmdStan is not installed. Install it with: cmdstanr::install_cmdstan()")
    }
  }, error = function(e) {
    stop("CmdStan is not installed. Install it with: cmdstanr::install_cmdstan()")
  })

  # Validate inputs
  checkmate::assert_class(network, "rmst_network")
  checkmate::assert_numeric(tau, lower = 0, min.len = 1)
  checkmate::assert_choice(baseline, c("rp_spline", "piecewise"))
  checkmate::assert_choice(random_effects, c("normal", "student_t", "none"))
  checkmate::assert_logical(inconsistency, len = 1)
  checkmate::assert_logical(propagate_recon_uncertainty, len = 1)

  if (is.null(warmup)) {
    warmup <- floor(iter / 2)
  }

  # Prepare Stan data
  stan_data <- prepare_stan_data(
    network = network,
    tau = tau,
    baseline = baseline,
    random_effects = random_effects,
    inconsistency = inconsistency
  )

  # Handle reconstruction uncertainty if requested
  if (propagate_recon_uncertainty && !is.null(recon)) {
    message("Propagating reconstruction uncertainty...")
    stan_data <- add_reconstruction_uncertainty(stan_data, network, recon)
  }

  # Compile Stan model
  model_file <- system.file("stan", "rmst_arm_based.stan",
                            package = "rmstnma")

  if (!file.exists(model_file)) {
    stop("Stan model file not found. Please ensure package is properly installed.")
  }

  model <- cmdstanr::cmdstan_model(model_file, compile = TRUE)

  # Fit the model
  message("Fitting RMST NMA model...")
  message(sprintf("  - Baseline: %s", baseline))
  message(sprintf("  - Random effects: %s", random_effects))
  message(sprintf("  - Tau: %s", paste(tau, collapse = ", ")))

  fit <- model$sample(
    data = stan_data,
    chains = chains,
    iter_sampling = iter - warmup,
    iter_warmup = warmup,
    adapt_delta = adapt_delta,
    ...
  )

  # Create output object
  result <- list(
    fit = fit,
    network = network,
    tau = tau,
    baseline = baseline,
    random_effects = random_effects,
    inconsistency = inconsistency,
    stan_data = stan_data,
    call = match.call()
  )

  class(result) <- c("rmst_nma_fit", "list")

  return(result)
}

#' Prepare Stan data
#' @keywords internal
prepare_stan_data <- function(network, tau, baseline,
                              random_effects, inconsistency) {

  # Extract network data
  data <- network$data

  # Create study and treatment indices
  study_idx <- as.numeric(factor(data$study, levels = network$studies))
  trt_idx <- as.numeric(factor(data$treatment, levels = network$treatments))

  # Count arms per study
  arms_per_study <- data %>%
    dplyr::group_by(study) %>%
    dplyr::summarise(n_arms = dplyr::n_distinct(treatment), .groups = "drop") %>%
    dplyr::pull(n_arms)

  # Get max times across all arms
  max_times <- data %>%
    dplyr::group_by(study, treatment) %>%
    dplyr::summarise(n_times = dplyr::n(), .groups = "drop") %>%
    dplyr::pull(n_times) %>%
    max()

  # Create arm-level data
  arm_data <- data %>%
    dplyr::group_by(study, treatment) %>%
    dplyr::summarise(
      study_idx = dplyr::first(as.numeric(factor(study, levels = network$studies))),
      trt_idx = dplyr::first(as.numeric(factor(treatment, levels = network$treatments))),
      .groups = "drop"
    )

  n_arms <- nrow(arm_data)

  # Initialize matrices for times and survival
  times_matrix <- matrix(0, n_arms, max_times)
  surv_matrix <- matrix(0, n_arms, max_times)
  n_times_vec <- integer(n_arms)
  study_vec <- integer(n_arms)
  treatment_vec <- integer(n_arms)

  # Fill matrices
  for (i in 1:n_arms) {
    arm_study <- arm_data$study[i]
    arm_trt <- arm_data$treatment[i]

    arm_subset <- data[data$study == arm_study & data$treatment == arm_trt, ]
    n_obs <- nrow(arm_subset)

    n_times_vec[i] <- n_obs
    study_vec[i] <- arm_data$study_idx[i]
    treatment_vec[i] <- arm_data$trt_idx[i]

    if (n_obs > 0) {
      times_matrix[i, 1:n_obs] <- arm_subset$time
      surv_matrix[i, 1:n_obs] <- arm_subset$survival
    }
  }

  # Create Stan data list
  stan_data <- list(
    # Dimensions
    N_studies = network$n_studies,
    N_treatments = network$n_treatments,
    N_arms = n_arms,
    N_tau = length(tau),

    # Survival data
    max_times = max_times,
    times = times_matrix,
    surv = surv_matrix,
    n_times = n_times_vec,

    # Restriction times
    tau = tau,

    # Study/treatment mappings
    study = study_vec,
    treatment = treatment_vec,
    n_arms_by_study = arms_per_study,

    # Model settings
    use_re = as.integer(random_effects != "none"),
    re_dist = ifelse(random_effects == "student_t", 1L, 0L),
    baseline_model = ifelse(baseline == "rp_spline", 1L, 0L),

    # Priors
    prior_mean = 0.0,
    prior_sd = 2.5
  )

  return(stan_data)
}
