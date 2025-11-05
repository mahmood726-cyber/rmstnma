#' Generate Example IPD Data
#'
#' Creates simulated IPD for testing and examples
#'
#' @param n_trials Number of trials (positive integer)
#' @param n_per_arm Sample size per arm (positive integer)
#' @param seed Random seed (integer or NULL for no seed)
#' @return Data frame with simulated IPD
#' @export
#' @examples
#' ipd <- generate_example_ipd(n_trials = 3, n_per_arm = 50)
#' head(ipd)
generate_example_ipd <- function(n_trials = 5, n_per_arm = 100, seed = 42) {
  # Input validation
  assert_positive_integer(n_trials, "n_trials", "generate_example_ipd")
  assert_positive_integer(n_per_arm, "n_per_arm", "generate_example_ipd")

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("[generate_example_ipd] Argument 'seed' must be a single numeric value or NULL",
           call. = FALSE)
    }
    set.seed(seed)
  }

  treatments <- c("Control", "DrugA", "DrugB", "DrugC", "DrugD")
  event_rates <- get_event_rates()

  designs <- list(
    c("Control", "DrugA"), c("Control", "DrugB"), c("Control", "DrugC"),
    c("DrugA", "DrugB"), c("DrugA", "DrugC"), c("DrugB", "DrugD"),
    c("Control", "DrugA", "DrugB"), c("Control", "DrugC", "DrugD")
  )

  ipd <- data.frame()

  for (i in seq_len(n_trials)) {
    design <- designs[[(i - 1) %% length(designs) + 1]]

    for (trt in design) {
      rate <- event_rates[trt]
      if (is.na(rate)) rate <- .EVENT_RATE_CONTROL  # Fallback for unknown treatments

      event_times <- stats::rexp(n_per_arm, rate)
      censor_times <- stats::runif(n_per_arm, 0, .DEFAULT_CENSOR_MAX)

      obs_times <- pmin(event_times, censor_times)
      obs_status <- as.integer(event_times <= censor_times)

      trial_data <- data.frame(
        trial = sprintf("Trial_%02d", i),
        treatment = trt,
        time = obs_times,
        status = obs_status,
        stringsAsFactors = FALSE
      )

      ipd <- rbind(ipd, trial_data)
    }
  }

  ipd
}

#' Generate Simulated NMA Data
#'
#' Creates realistic network meta-analysis data with covariates
#'
#' @param n_studies Number of studies (positive integer)
#' @param seed Random seed (integer or NULL for no seed)
#' @return Data frame with pairwise comparisons and covariates
#' @export
#' @examples
#' data <- simulate_nma_data(n_studies = 40)
#' head(data)
simulate_nma_data <- function(n_studies = 40, seed = 42) {
  # Input validation
  assert_positive_integer(n_studies, "n_studies", "simulate_nma_data")

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("[simulate_nma_data] Argument 'seed' must be a single numeric value or NULL",
           call. = FALSE)
    }
    set.seed(seed)
  }

  treatments <- c("Placebo", "DrugA", "DrugB", "DrugC", "DrugD")
  true_eff <- get_true_effects()

  designs <- list(
    c("Placebo", "DrugA"), c("Placebo", "DrugB"), c("Placebo", "DrugC"),
    c("DrugA", "DrugB"), c("DrugA", "DrugC"), c("DrugB", "DrugD"),
    c("Placebo", "DrugA", "DrugB"), c("Placebo", "DrugC", "DrugD")
  )

  study_list <- sample(designs, n_studies, replace = TRUE)

  pw <- purrr::map_dfr(seq_along(study_list), function(i) {
    stud <- sprintf("Study_%02d", i)
    tr <- study_list[[i]]
    pairs <- utils::combn(tr, 2, simplify = FALSE)

    purrr::map_dfr(pairs, function(p) {
      t1 <- p[1]
      t2 <- p[2]
      true_diff <- true_eff[t1] - true_eff[t2]
      study_re <- stats::rnorm(1, 0, .DEFAULT_STUDY_RE_SD)
      se <- stats::runif(1, .DEFAULT_SE_MIN, .DEFAULT_SE_MAX)
      TE <- stats::rnorm(1, true_diff + study_re, se)

      tibble::tibble(
        studlab = stud,
        treat1 = t1,
        treat2 = t2,
        TE = TE,
        seTE = se
      )
    })
  })

  # Add study-level covariates
  info <- pw %>%
    dplyr::distinct(studlab) %>%
    dplyr::mutate(
      year = sample(.STUDY_YEAR_MIN:.STUDY_YEAR_MAX, dplyr::n(), replace = TRUE),
      is_rct = stats::rbinom(dplyr::n(), 1, .DEFAULT_RCT_PROB),
      study_design = factor(ifelse(is_rct == 1, "RCT", "Non-RCT")),
      grade = factor(
        sample(names(.GRADE_PROBS), dplyr::n(), replace = TRUE, prob = .GRADE_PROBS),
        levels = names(.GRADE_PROBS)
      ),
      age_mean = round(stats::rnorm(dplyr::n(), .DEFAULT_AGE_MEAN, .DEFAULT_AGE_SD), 1),
      female_pct = pmin(0.8, pmax(0.2, stats::rnorm(dplyr::n(), .DEFAULT_FEMALE_PCT, .DEFAULT_FEMALE_SD))),
      bmi_mean = round(stats::rnorm(dplyr::n(), .DEFAULT_BMI_MEAN, .DEFAULT_BMI_SD), 1),
      charlson = round(pmax(0, stats::rnorm(dplyr::n(), .DEFAULT_CHARLSON_MEAN, .DEFAULT_CHARLSON_SD)), 1)
    )

  dplyr::left_join(pw, info, by = "studlab")
}

#' Convert Arm-Level to Pairwise Data
#'
#' Converts arm-based data to pairwise format for NMA
#'
#' @param df Data frame with arm-level data
#' @param outcome Type: "binary" or "continuous"
#' @param sm Summary measure (OR, RR, MD, etc.)
#' @return Pairwise data frame
#' @export
make_pairwise_from_arms <- function(df,
                                    outcome = c("binary", "continuous"),
                                    sm = "OR") {
  outcome <- match.arg(outcome)

  if (!has_pkg("netmeta")) {
    stop("Package 'netmeta' required for pairwise conversion")
  }

  if (outcome == "binary") {
    if (!all(c("treat", "event", "n", "studlab") %in% names(df))) {
      stop("Binary data needs: studlab, treat, event, n")
    }
    pw <- netmeta::pairwise(
      treat = treat, event = event, n = n,
      studlab = studlab, data = df, sm = sm
    )
  } else {
    if (!all(c("treat", "mean", "sd", "n", "studlab") %in% names(df))) {
      stop("Continuous data needs: studlab, treat, mean, sd, n")
    }
    pw <- netmeta::pairwise(
      treat = treat, mean = mean, sd = sd, n = n,
      studlab = studlab, data = df, sm = sm
    )
  }

  tibble::as_tibble(pw)
}

#' Read example data from package
#'
#' Load built-in example datasets
#'
#' @param dataset Name of dataset ("ipd", "nma", or "km")
#' @return Data frame
#' @export
#' @examples
#' \dontrun{
#' data <- read_example_data("nma")
#' }
read_example_data <- function(dataset = c("ipd", "nma", "km")) {
  dataset <- match.arg(dataset)

  path <- system.file("extdata",
                     paste0("example_", dataset, ".csv"),
                     package = "powerNMA")

  if (path == "" || !file.exists(path)) {
    warning("Example data not found, generating simulated data")
    return(switch(dataset,
      "ipd" = generate_example_ipd(),
      "nma" = simulate_nma_data(),
      "km" = NULL
    ))
  }

  readr::read_csv(path, col_types = readr::cols(), show_col_types = FALSE)
}
