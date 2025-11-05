#' Effect Size Conversions for Network Meta-Analysis
#'
#' Convert between different effect size metrics (OR, RR, HR, MD, SMD, etc.)
#' and transform effect sizes for different scales.
#'
#' @name effect_size_conversions
NULL

#' Convert Odds Ratio to Risk Ratio
#'
#' Convert odds ratios to risk ratios given the baseline risk in the control group.
#'
#' @param OR Odds ratio
#' @param p0 Baseline risk in control group (probability, 0-1)
#' @return Risk ratio
#' @export
#' @references
#' Zhang J, Yu KF. What's the relative risk? A method of correcting the odds
#' ratio in cohort studies of common outcomes. JAMA. 1998;280(19):1690-1.
#'
#' @examples
#' # OR = 2.0, baseline risk = 0.20
#' or_to_rr(2.0, 0.20)  # Returns RR ≈ 1.67
or_to_rr <- function(OR, p0) {
  if (any(p0 < 0 | p0 > 1)) {
    stop("Baseline risk p0 must be between 0 and 1")
  }

  if (any(OR <= 0)) {
    stop("OR must be positive")
  }

  # Formula: RR = OR / (1 - p0 + (p0 * OR))
  RR <- OR / (1 - p0 + (p0 * OR))
  RR
}

#' Convert Risk Ratio to Odds Ratio
#'
#' Convert risk ratios to odds ratios given the baseline risk.
#'
#' @param RR Risk ratio
#' @param p0 Baseline risk in control group
#' @return Odds ratio
#' @export
#' @examples
#' rr_to_or(1.5, 0.20)
rr_to_or <- function(RR, p0) {
  if (any(p0 < 0 | p0 > 1)) {
    stop("Baseline risk p0 must be between 0 and 1")
  }

  if (any(RR <= 0)) {
    stop("RR must be positive")
  }

  # Formula: OR = RR * (1 - p0) / (1 - RR * p0)
  OR <- RR * (1 - p0) / (1 - RR * p0)
  OR
}

#' Convert Log Odds Ratio to Log Risk Ratio
#'
#' @param logOR Log odds ratio
#' @param p0 Baseline risk
#' @return Log risk ratio
#' @export
#' @examples
#' logor_to_logrr(0.693, 0.20)  # log(2.0) to log(RR)
logor_to_logrr <- function(logOR, p0) {
  OR <- exp(logOR)
  RR <- or_to_rr(OR, p0)
  log(RR)
}

#' Convert Cohen's d to Hedges' g
#'
#' Apply small-sample bias correction to Cohen's d to get Hedges' g.
#'
#' @param d Cohen's d (standardized mean difference)
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @return Hedges' g (bias-corrected SMD)
#' @export
#' @references
#' Hedges LV. Distribution theory for Glass's estimator of effect size and
#' related estimators. Journal of Educational Statistics. 1981;6(2):107-128.
#'
#' @examples
#' cohens_d_to_hedges_g(0.5, n1 = 30, n2 = 30)
cohens_d_to_hedges_g <- function(d, n1, n2) {
  df <- n1 + n2 - 2

  # Correction factor J
  J <- 1 - (3 / (4 * df - 1))

  g <- d * J
  g
}

#' Convert Mean Difference to Standardized Mean Difference
#'
#' Standardize a mean difference using the pooled standard deviation.
#'
#' @param MD Mean difference
#' @param sd1 Standard deviation in group 1
#' @param sd2 Standard deviation in group 2
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @return Standardized mean difference (Cohen's d)
#' @export
#' @examples
#' md_to_smd(10, sd1 = 15, sd2 = 18, n1 = 50, n2 = 50)
md_to_smd <- function(MD, sd1, sd2, n1, n2) {
  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

  # Cohen's d
  d <- MD / pooled_sd
  d
}

#' Convert Standardized Mean Difference to Mean Difference
#'
#' Convert SMD back to MD using a specified standard deviation.
#'
#' @param SMD Standardized mean difference
#' @param sd Standard deviation for scaling
#' @return Mean difference
#' @export
#' @examples
#' smd_to_md(0.5, sd = 15)  # Returns MD = 7.5
smd_to_md <- function(SMD, sd) {
  MD <- SMD * sd
  MD
}

#' Convert Correlation to Fisher's Z
#'
#' Fisher's z-transformation for correlation coefficients.
#'
#' @param r Correlation coefficient (-1 to 1)
#' @return Fisher's z
#' @export
#' @examples
#' r_to_fishers_z(0.5)
r_to_fishers_z <- function(r) {
  if (any(abs(r) >= 1)) {
    stop("Correlation r must be between -1 and 1 (exclusive)")
  }

  z <- 0.5 * log((1 + r) / (1 - r))
  z
}

#' Convert Fisher's Z to Correlation
#'
#' @param z Fisher's z
#' @return Correlation coefficient
#' @export
#' @examples
#' fishers_z_to_r(0.5493)  # Returns r ≈ 0.5
fishers_z_to_r <- function(z) {
  r <- (exp(2 * z) - 1) / (exp(2 * z) + 1)
  r
}

#' Transform Effect Sizes in NMA Data
#'
#' Apply transformations to effect sizes in a pairwise NMA dataset.
#' Useful for sensitivity analyses or combining different metrics.
#'
#' @param data Pairwise data with TE (treatment effect) column
#' @param from Original effect size type
#' @param to Target effect size type
#' @param baseline_risk Baseline risk for OR/RR conversions
#' @param scale_sd Standard deviation for SMD/MD conversions
#' @return Data frame with transformed effects
#' @export
#' @examples
#' \dontrun{
#' # Convert log odds ratios to log risk ratios
#' data_rr <- transform_effect_sizes(data, from = "logOR", to = "logRR",
#'                                   baseline_risk = 0.15)
#' }
transform_effect_sizes <- function(data,
                                  from = c("logOR", "logRR", "MD", "SMD"),
                                  to = c("logRR", "logOR", "SMD", "MD"),
                                  baseline_risk = NULL,
                                  scale_sd = NULL) {

  from <- match.arg(from)
  to <- match.arg(to)

  if (!("TE" %in% names(data))) {
    stop("Data must contain 'TE' column for treatment effects")
  }

  # Make a copy
  transformed_data <- data

  if (from == "logOR" && to == "logRR") {
    if (is.null(baseline_risk)) {
      stop("baseline_risk required for OR to RR conversion")
    }
    transformed_data$TE <- logor_to_logrr(data$TE, baseline_risk)
    transformed_data$sm <- "RR"

  } else if (from == "logRR" && to == "logOR") {
    if (is.null(baseline_risk)) {
      stop("baseline_risk required for RR to OR conversion")
    }
    RR <- exp(data$TE)
    OR <- rr_to_or(RR, baseline_risk)
    transformed_data$TE <- log(OR)
    transformed_data$sm <- "OR"

  } else if (from == "MD" && to == "SMD") {
    if (is.null(scale_sd)) {
      stop("scale_sd required for MD to SMD conversion")
    }
    transformed_data$TE <- data$TE / scale_sd
    # SE also needs to be scaled
    if ("seTE" %in% names(data)) {
      transformed_data$seTE <- data$seTE / scale_sd
    }
    transformed_data$sm <- "SMD"

  } else if (from == "SMD" && to == "MD") {
    if (is.null(scale_sd)) {
      stop("scale_sd required for SMD to MD conversion")
    }
    transformed_data$TE <- data$TE * scale_sd
    if ("seTE" %in% names(data)) {
      transformed_data$seTE <- data$seTE * scale_sd
    }
    transformed_data$sm <- "MD"

  } else {
    stop(sprintf("Conversion from %s to %s not supported", from, to))
  }

  message(sprintf("Transformed %d effect sizes from %s to %s", nrow(data), from, to))

  transformed_data
}

#' Calculate Number Needed to Treat (NNT)
#'
#' Calculate NNT from odds ratio, risk ratio, or absolute risk difference.
#'
#' @param effect Effect size (OR, RR, or ARD)
#' @param type Effect type: "OR", "RR", or "ARD"
#' @param baseline_risk Baseline risk in control group (for OR/RR)
#' @return Number needed to treat
#' @export
#' @references
#' Altman DG. Confidence intervals for the number needed to treat.
#' BMJ. 1998;317(7168):1309-1312.
#'
#' @examples
#' # From OR
#' calculate_nnt(effect = 0.5, type = "OR", baseline_risk = 0.20)
#'
#' # From absolute risk difference
#' calculate_nnt(effect = 0.10, type = "ARD")
calculate_nnt <- function(effect, type = c("OR", "RR", "ARD"), baseline_risk = NULL) {

  type <- match.arg(type)

  if (type == "ARD") {
    # Direct from absolute risk difference
    NNT <- 1 / abs(effect)

  } else if (type == "OR") {
    if (is.null(baseline_risk)) {
      stop("baseline_risk required for OR to NNT conversion")
    }

    # Convert OR to ARD
    # p1 (treatment) = (OR * p0) / (1 - p0 + OR * p0)
    p1 <- (effect * baseline_risk) / (1 - baseline_risk + effect * baseline_risk)
    ARD <- p1 - baseline_risk
    NNT <- 1 / abs(ARD)

  } else if (type == "RR") {
    if (is.null(baseline_risk)) {
      stop("baseline_risk required for RR to NNT conversion")
    }

    # p1 = RR * p0
    p1 <- effect * baseline_risk
    ARD <- p1 - baseline_risk
    NNT <- 1 / abs(ARD)

  } else {
    stop("Invalid type")
  }

  NNT
}

#' Calculate NNT with Confidence Interval
#'
#' Calculate NNT with confidence intervals from effect estimate and CI.
#'
#' @param effect Effect estimate
#' @param lower Lower CI
#' @param upper Upper CI
#' @param type Effect type
#' @param baseline_risk Baseline risk
#' @return List with NNT and CI
#' @export
#' @examples
#' \dontrun{
#' nnt_with_ci(effect = 0.5, lower = 0.3, upper = 0.8,
#'            type = "OR", baseline_risk = 0.20)
#' }
nnt_with_ci <- function(effect, lower, upper, type = c("OR", "RR", "ARD"),
                       baseline_risk = NULL) {

  type <- match.arg(type)

  NNT <- calculate_nnt(effect, type, baseline_risk)
  NNT_lower <- calculate_nnt(upper, type, baseline_risk)  # Reversed for CI
  NNT_upper <- calculate_nnt(lower, type, baseline_risk)

  list(
    NNT = NNT,
    CI_lower = NNT_lower,
    CI_upper = NNT_upper,
    interpretation = sprintf(
      "NNT = %.1f (95%% CI: %.1f to %.1f)",
      NNT, NNT_lower, NNT_upper
    )
  )
}

#' Batch Convert Effect Sizes
#'
#' Convert multiple effect sizes at once with proper error handling.
#'
#' @param effects Vector of effect sizes
#' @param conversion_func Conversion function to apply
#' @param ... Additional arguments passed to conversion function
#' @return Vector of converted effects
#' @export
#' @examples
#' \dontrun{
#' ORs <- c(1.5, 2.0, 2.5)
#' RRs <- batch_convert(ORs, or_to_rr, p0 = 0.20)
#' }
batch_convert <- function(effects, conversion_func, ...) {

  results <- sapply(effects, function(e) {
    tryCatch({
      conversion_func(e, ...)
    }, error = function(err) {
      warning(sprintf("Conversion failed for effect = %.3f: %s", e, err$message))
      NA
    })
  })

  results
}

#' Get Approximate Conversion
#'
#' Provide rule-of-thumb approximations for effect size conversions when
#' exact information is unavailable.
#'
#' @param from Original metric
#' @param to Target metric
#' @param value Effect size value
#' @return Approximate converted value with warning
#' @export
#' @examples
#' # Approximate OR to RR (assumes moderate baseline risk ~0.20)
#' approximate_conversion("OR", "RR", 2.0)
approximate_conversion <- function(from, to, value) {

  if (from == "OR" && to == "RR") {
    # Rough approximation assuming baseline risk ~0.15-0.20
    warning("Using approximate conversion OR to RR (assumes baseline risk ~0.20)")
    RR_approx <- or_to_rr(value, 0.20)
    return(RR_approx)
  }

  if (from == "RR" && to == "OR") {
    warning("Using approximate conversion RR to OR (assumes baseline risk ~0.20)")
    OR_approx <- rr_to_or(value, 0.20)
    return(OR_approx)
  }

  if (from == "d" && to == "g") {
    # Assumes df ~50
    warning("Using approximate correction factor J ~0.98 (assumes df ~50)")
    return(value * 0.98)
  }

  stop(sprintf("No approximation available for %s to %s", from, to))
}

#' Create Effect Size Conversion Table
#'
#' Generate a reference table showing conversions for a range of effect sizes.
#'
#' @param values Vector of effect sizes to convert
#' @param from Original metric
#' @param to Target metric
#' @param ... Additional parameters (baseline_risk, scale_sd, etc.)
#' @return Data frame with conversions
#' @export
#' @examples
#' create_conversion_table(values = seq(0.5, 3.0, 0.5),
#'                        from = "OR", to = "RR", baseline_risk = 0.20)
create_conversion_table <- function(values, from, to, ...) {

  conversion_func <- switch(
    paste(from, to, sep = "_to_"),
    "OR_to_RR" = or_to_rr,
    "RR_to_OR" = rr_to_or,
    "MD_to_SMD" = md_to_smd,
    "SMD_to_MD" = smd_to_md,
    stop(sprintf("Conversion %s to %s not available", from, to))
  )

  converted <- batch_convert(values, conversion_func, ...)

  table <- data.frame(
    Original = values,
    Converted = converted
  )

  names(table) <- c(from, to)

  table
}
