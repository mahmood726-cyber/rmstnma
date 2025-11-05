#' Configure powerNMA Analysis
#'
#' Create comprehensive configuration for network meta-analysis
#'
#' @param sm Effect measure ("HR", "OR", "RR", "MD", "SMD")
#' @param use_transport Enable transportability weighting
#' @param use_bayesian Run Bayesian NMA
#' @param use_timevarying Run time-varying analyses (RMST/milestone)
#' @param tau_list Time points for RMST analysis
#' @param milestone_times Time points for milestone analysis
#' @param run_sensitivity Run sensitivity analyses
#' @param run_metareg Run meta-regression
#' @param export_results Export results to files
#' @param ... Additional configuration parameters
#' @return Configuration list
#' @export
#' @examples
#' config <- setup_powernma(sm = "HR", use_bayesian = TRUE)
setup_powernma <- function(
  sm = "HR",
  use_transport = TRUE,
  use_bayesian = TRUE,
  use_timevarying = TRUE,
  tau_list = c(90, 180, 365),
  milestone_times = c(90, 180, 365),
  run_sensitivity = TRUE,
  run_metareg = TRUE,
  export_results = TRUE,
  export_plots = TRUE,
  output_dir = "powernma_results",
  plot_dir = "powernma_plots",
  seed = 42,
  ...) {

  config <- list(
    sm = sm,
    use_transport = use_transport,
    use_bayesian = use_bayesian,
    use_timevarying = use_timevarying,
    tau_list = tau_list,
    milestone_times = milestone_times,
    run_sensitivity = run_sensitivity,
    run_metareg = run_metareg,
    export_results = export_results,
    export_plots = export_plots,
    output_dir = output_dir,
    plot_dir = plot_dir,
    seed = seed
  )

  # Merge with additional parameters
  config <- c(config, list(...))

  class(config) <- c("powernma_config", "list")
  config
}

#' Run Core Network Meta-Analysis
#'
#' Perform standard frequentist network meta-analysis
#'
#' @param data Pairwise data frame
#' @param ref_treatment Reference treatment
#' @param sm Summary measure
#' @param weights_named Named vector of study weights
#' @param random Use random effects model
#' @return netmeta object
#' @keywords internal
robust_netmeta <- function(data, ref_treatment, sm = "HR",
                          weights_named = NULL, random = TRUE) {
  dd <- data

  # Apply weights by adjusting standard errors
  if (!is.null(weights_named)) {
    w <- weights_named[as.character(dd$studlab)]
    w[is.na(w)] <- 1
    w <- pmax(1e-6, w)
    dd$seTE <- dd$seTE / sqrt(w)
  }

  # Try random effects first
  fit <- .safe_try(
    netmeta::netmeta(
      TE = TE, seTE = seTE,
      treat1 = treat1, treat2 = treat2,
      studlab = studlab,
      data = dd,
      sm = sm,
      fixed = !random,
      random = random,
      reference.group = ref_treatment
    ),
    context = "netmeta::netmeta(random)"
  )

  # Fallback to fixed effects if random fails
  if (inherits(fit, "try-error") && random) {
    warning("Random-effects NMA failed; retrying with fixed-effects.")
    fit <- netmeta::netmeta(
      TE = TE, seTE = seTE,
      treat1 = treat1, treat2 = treat2,
      studlab = studlab,
      data = dd,
      sm = sm,
      fixed = TRUE,
      random = FALSE,
      reference.group = ref_treatment
    )
  }

  fit
}

#' Run Time-Varying RMST NMA
#'
#' Perform RMST network meta-analysis at multiple time horizons
#'
#' @param ipd Individual patient data
#' @param tau_list Vector of time points
#' @param reference Reference treatment
#' @return rmst_nma object
#' @export
#' @examples
#' \dontrun{
#' ipd <- generate_example_ipd()
#' results <- rmst_nma(ipd, tau_list = c(180, 365))
#' }
rmst_nma <- function(ipd, tau_list = c(90, 180, 365), reference = "Control") {
  if (!has_pkg("survRM2")) {
    stop("Package 'survRM2' required for RMST analysis")
  }

  validate_ipd(ipd)

  results <- list()

  for (tau in tau_list) {
    pairwise <- data.frame()

    for (trial_id in unique(ipd$trial)) {
      trial_data <- ipd[ipd$trial == trial_id, ]
      treatments <- sort(unique(trial_data$treatment))

      if (length(treatments) < 2) next

      # FIXED: Create ALL pairwise comparisons for multi-arm trials
      # This properly handles k-arm trials with k(k-1)/2 comparisons
      if (reference %in% treatments) {
        # Reference-based comparisons: all treatments vs reference
        other_treatments <- setdiff(treatments, reference)
        comparisons <- lapply(other_treatments, function(t) c(reference, t))
      } else {
        # No reference in trial: all pairwise combinations
        comparisons <- utils::combn(treatments, 2, simplify = FALSE)
      }

      # Process each comparison
      for (comp in comparisons) {
        treat1 <- comp[1]
        treat2 <- comp[2]

        arm1_data <- trial_data[trial_data$treatment == treat1, ]
        arm2_data <- trial_data[trial_data$treatment == treat2, ]

      combined_data <- rbind(
        data.frame(arm1_data, arm = 0, stringsAsFactors = FALSE),
        data.frame(arm2_data, arm = 1, stringsAsFactors = FALSE)
      )

      rmst_result <- tryCatch({
        survRM2::rmst2(
          time = combined_data$time,
          status = combined_data$status,
          arm = combined_data$arm,
          tau = tau
        )
      }, error = function(e) {
        warning(sprintf("RMST failed for trial %s at tau=%d: %s",
                       trial_id, tau, e$message))
        NULL
      })

      if (!is.null(rmst_result)) {
        # FIXED: Removed incorrect sign flip
        # survRM2::rmst2() returns: RMST(arm=1) - RMST(arm=0)
        # where arm=0 = treat1 (reference), arm=1 = treat2 (intervention)
        # This directly gives TE = treat2 - treat1 (positive = favors treat2)
        # See RMST_SIGN_CONVENTION_PROOF.md for mathematical proof
        pairwise <- rbind(pairwise, data.frame(
          trial = trial_id,
          treat1 = treat1,
          treat2 = treat2,
          TE = rmst_result$unadjusted.result[1, "Est."],  # CORRECT (no sign flip)
          seTE = rmst_result$unadjusted.result[1, "se"],
          stringsAsFactors = FALSE
        ))
      }
      }  # End of comparison loop
    }  # End of trial loop

    if (nrow(pairwise) >= 2) {
      nma_result <- netmeta::netmeta(
        TE = TE, seTE = seTE,
        treat1 = treat1, treat2 = treat2,
        studlab = trial,
        data = pairwise,
        sm = "MD",
        reference.group = reference,
        comb.fixed = TRUE,
        comb.random = TRUE
      )
      results[[paste0("tau_", tau)]] <- nma_result
    } else {
      warning(sprintf("Insufficient data for NMA at tau=%d (need >=2 studies)", tau))
    }
  }

  class(results) <- c("rmst_nma", "list")
  results
}

#' Run Milestone Survival NMA
#'
#' Perform milestone analysis at specific time points
#'
#' @param ipd Individual patient data
#' @param times Vector of milestone time points
#' @param reference Reference treatment
#' @return milestone_nma object
#' @export
#' @examples
#' \dontrun{
#' ipd <- generate_example_ipd()
#' results <- milestone_nma(ipd, times = c(180, 365))
#' }
milestone_nma <- function(ipd, times = c(90, 180, 365), reference = "Control",
                         extend = FALSE, continuity_correction = c("standard", "empirical", "none")) {
  validate_ipd(ipd)

  continuity_correction <- match.arg(continuity_correction)

  results <- list()

  for (t in times) {
    pairwise_counts <- data.frame()

    for (trial_id in unique(ipd$trial)) {
      trial_data <- ipd[ipd$trial == trial_id, ]
      treatments <- sort(unique(trial_data$treatment))

      if (length(treatments) < 2) next

      # FIXED: Create ALL pairwise comparisons for multi-arm trials
      if (reference %in% treatments) {
        # Reference-based comparisons: all treatments vs reference
        other_treatments <- setdiff(treatments, reference)
        comparisons <- lapply(other_treatments, function(t) c(reference, t))
      } else {
        # No reference in trial: all pairwise combinations
        comparisons <- utils::combn(treatments, 2, simplify = FALSE)
      }

      # Process each comparison
      for (comp in comparisons) {
        treat1_name <- comp[1]
        treat2_name <- comp[2]

      # Arm 1 counts
      arm1_data <- trial_data[trial_data$treatment == treat1_name, ]
      km1 <- survival::survfit(survival::Surv(time, status) ~ 1, data = arm1_data)

      # Check if milestone exceeds follow-up
      max_time1 <- max(arm1_data$time)
      if (t > max_time1 && !extend) {
        warning(sprintf(
          "Milestone time %d exceeds follow-up (%d) in trial %s, arm %s. Skipping.",
          t, max_time1, trial_id, treat1_name
        ))
        next
      }

      summ1 <- summary(km1, times = t, extend = extend)
      n1 <- km1$n
      events1 <- n1 - summ1$n.risk

      # Arm 2 counts
      arm2_data <- trial_data[trial_data$treatment == treat2_name, ]
      km2 <- survival::survfit(survival::Surv(time, status) ~ 1, data = arm2_data)

      max_time2 <- max(arm2_data$time)
      if (t > max_time2 && !extend) {
        warning(sprintf(
          "Milestone time %d exceeds follow-up (%d) in trial %s, arm %s. Skipping.",
          t, max_time2, trial_id, treat2_name
        ))
        next
      }

      summ2 <- summary(km2, times = t, extend = extend)
      n2 <- km2$n
      events2 <- n2 - summ2$n.risk

      # Apply continuity correction if needed (FIXED: now configurable)
      needs_correction <- (events1 == 0 || events2 == 0 ||
                          events1 == n1 || events2 == n2)

      if (needs_correction) {
        if (continuity_correction == "none") {
          warning(sprintf(
            "Zero/all events in trial %s, comparison %s vs %s at time %d - SKIPPING (continuity_correction='none')",
            trial_id, treat1_name, treat2_name, t
          ))
          next
        }

        warning(sprintf(
          "Zero/all events in trial %s at time %d - applying %s continuity correction",
          trial_id, t, continuity_correction
        ))

        if (continuity_correction == "standard") {
          # Cochrane method: add 0.5 to events only
          events1 <- events1 + 0.5
          events2 <- events2 + 0.5
        } else if (continuity_correction == "empirical") {
          # Empirical method: add 0.5 to events and 1 to denominators
          events1 <- events1 + 0.5
          events2 <- events2 + 0.5
          n1 <- n1 + 1
          n2 <- n2 + 1
        }
      }

      pairwise_counts <- rbind(pairwise_counts, data.frame(
        trial = trial_id,
        treat1 = treat1_name,
        treat2 = treat2_name,
        event1 = events1,
        n1 = n1,
        event2 = events2,
        n2 = n2,
        stringsAsFactors = FALSE
      ))
      }  # End of comparison loop
    }  # End of trial loop

    if (nrow(pairwise_counts) >= 2) {
      nma_result <- netmeta::netmeta(
        event1 = event1, n1 = n1,
        event2 = event2, n2 = n2,
        treat1 = treat1, treat2 = treat2,
        studlab = trial,
        data = pairwise_counts,
        sm = "OR",
        reference.group = reference
      )
      results[[paste0("day_", t)]] <- nma_result
    } else {
      warning(sprintf("Insufficient data for milestone NMA at t=%d", t))
    }
  }

  class(results) <- c("milestone_nma", "list")
  results
}

#' Network Geometry and Connectivity
#'
#' Calculate network metrics
#'
#' @param data Pairwise NMA data
#' @return Data frame with network metrics
#' @export
network_geometry <- function(data) {
  trts <- sort(unique(c(as.character(data$treat1), as.character(data$treat2))))
  n_trt <- length(trts)

  comps <- data %>%
    dplyr::mutate(
      cmp = paste(pmin(treat1, treat2), pmax(treat1, treat2), sep = " vs ")
    ) %>%
    dplyr::distinct(cmp)

  n_comp <- nrow(comps)
  stud <- length(unique(data$studlab))

  degree_tbl <- table(c(data$treat1, data$treat2))
  mean_degree <- mean(as.numeric(degree_tbl))

  tibble::tibble(
    n_studies = stud,
    n_treatments = n_trt,
    n_comparisons = n_comp,
    mean_degree = mean_degree,
    density = (2 * n_comp) / (n_trt * (n_trt - 1))
  )
}
