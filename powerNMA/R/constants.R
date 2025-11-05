# ============================================================================
# Package-Level Constants
# ============================================================================
#
# This file defines constants used throughout the powerNMA package.
# Using named constants improves code readability and maintainability.
#

#' @keywords internal

# Simulation defaults
.DEFAULT_SE_MIN <- 0.08
.DEFAULT_SE_MAX <- 0.25
.DEFAULT_STUDY_RE_SD <- 0.10
.DEFAULT_CENSOR_MAX <- 500

# Numerical constants
.MIN_WEIGHT <- 1e-6
.MIN_STUDIES_FOR_LOO <- 5
.MIN_COMPARISONS_FOR_NMA <- 2
.NUMERIC_TOLERANCE <- 1e-10

# Time horizons (days)
.DEFAULT_TAU_LIST <- c(90, 180, 365)
.DEFAULT_MILESTONE_TIMES <- c(90, 180, 365)

# Event rates by treatment (for simulation)
.EVENT_RATE_CONTROL <- 0.003
.EVENT_RATE_DRUGA <- 0.002
.EVENT_RATE_DRUGB <- 0.0015
.EVENT_RATE_DRUGC <- 0.0025
.EVENT_RATE_DRUGD <- 0.0018

# True treatment effects (for simulation)
.TRUE_HR_DRUGA <- log(0.85)
.TRUE_HR_DRUGB <- log(0.75)
.TRUE_HR_DRUGC <- log(0.90)
.TRUE_HR_DRUGD <- log(0.78)

# Continuity correction values
.CONTINUITY_CORRECTION_STANDARD <- 0.5
.CONTINUITY_CORRECTION_EMPIRICAL_EVENTS <- 0.5
.CONTINUITY_CORRECTION_EMPIRICAL_N <- 1

# MCMC defaults
.DEFAULT_N_ADAPT <- 5000
.DEFAULT_N_ITER <- 10000
.DEFAULT_N_CHAINS <- 3

# Covariate defaults (for simulation)
.DEFAULT_AGE_MEAN <- 65
.DEFAULT_AGE_SD <- 5
.DEFAULT_FEMALE_PCT <- 0.45
.DEFAULT_FEMALE_SD <- 0.10
.DEFAULT_BMI_MEAN <- 28
.DEFAULT_BMI_SD <- 2
.DEFAULT_CHARLSON_MEAN <- 1.5
.DEFAULT_CHARLSON_SD <- 0.5
.DEFAULT_RCT_PROB <- 0.8

# GRADE quality probabilities
.GRADE_PROBS <- c(
  "High" = 0.4,
  "Moderate" = 0.4,
  "Low" = 0.15,
  "Very low" = 0.05
)

# Study year range
.STUDY_YEAR_MIN <- 2010
.STUDY_YEAR_MAX <- 2025

# Export constants to package namespace
# (These are internal and not exported to users)

#' Get simulation default constants
#' @keywords internal
get_simulation_defaults <- function() {
  list(
    se_min = .DEFAULT_SE_MIN,
    se_max = .DEFAULT_SE_MAX,
    study_re_sd = .DEFAULT_STUDY_RE_SD,
    censor_max = .DEFAULT_CENSOR_MAX,
    age_mean = .DEFAULT_AGE_MEAN,
    age_sd = .DEFAULT_AGE_SD,
    female_pct = .DEFAULT_FEMALE_PCT,
    female_sd = .DEFAULT_FEMALE_SD,
    bmi_mean = .DEFAULT_BMI_MEAN,
    bmi_sd = .DEFAULT_BMI_SD,
    charlson_mean = .DEFAULT_CHARLSON_MEAN,
    charlson_sd = .DEFAULT_CHARLSON_SD,
    rct_prob = .DEFAULT_RCT_PROB
  )
}

#' Get event rate constants
#' @keywords internal
get_event_rates <- function() {
  c(
    Control = .EVENT_RATE_CONTROL,
    DrugA = .EVENT_RATE_DRUGA,
    DrugB = .EVENT_RATE_DRUGB,
    DrugC = .EVENT_RATE_DRUGC,
    DrugD = .EVENT_RATE_DRUGD
  )
}

#' Get true treatment effect constants
#' @keywords internal
get_true_effects <- function() {
  c(
    Control = 0,
    Placebo = 0,
    DrugA = .TRUE_HR_DRUGA,
    DrugB = .TRUE_HR_DRUGB,
    DrugC = .TRUE_HR_DRUGC,
    DrugD = .TRUE_HR_DRUGD
  )
}
