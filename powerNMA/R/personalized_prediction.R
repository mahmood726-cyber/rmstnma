#' Personalized Treatment Prediction
#'
#' Framework for providing personalized treatment recommendations using
#' individual patient data (IPD) network meta-analysis and effect modifiers.
#'
#' @name personalized_prediction
#' @references
#' Debray TPA, et al. (2015). Individual participant data meta-analysis for a
#' binary outcome: one-stage or two-stage? PLOS ONE, 10(8):e0133865.
#'
#' Riley RD, et al. (2020). Individual participant data meta-analysis to
#' examine interactions between treatment effect and participant-level
#' covariates. Statistics in Medicine, 39(11):1693-1712.
NULL

#' Fit IPD Network Meta-Regression
#'
#' Fit network meta-regression using individual patient data to identify
#' effect modifiers and treatment-covariate interactions.
#'
#' @param ipd_data Individual patient data with columns: study, patient_id,
#'   treatment, outcome, and patient characteristics
#' @param outcome_var Name of outcome variable
#' @param treatment_var Name of treatment variable
#' @param study_var Name of study variable
#' @param effect_modifiers Character vector of effect modifier variable names
#' @param reference Reference treatment
#' @param family Outcome family ("gaussian", "binomial", "poisson")
#' @return ipd_nma object
#' @export
#' @examples
#' \dontrun{
#' # IPD with patient characteristics
#' ipd <- data.frame(
#'   study = rep(1:10, each = 100),
#'   patient_id = 1:1000,
#'   treatment = sample(c("A", "B", "C"), 1000, replace = TRUE),
#'   outcome = rnorm(1000),
#'   age = rnorm(1000, 50, 10),
#'   sex = sample(c("M", "F"), 1000, replace = TRUE),
#'   baseline_severity = rnorm(1000, 0, 1)
#' )
#'
#' ipd_model <- fit_ipd_metaregression(
#'   ipd_data = ipd,
#'   outcome_var = "outcome",
#'   treatment_var = "treatment",
#'   study_var = "study",
#'   effect_modifiers = c("age", "baseline_severity"),
#'   reference = "A"
#' )
#' }
fit_ipd_metaregression <- function(ipd_data, outcome_var, treatment_var,
                                   study_var, effect_modifiers,
                                   reference = NULL, family = "gaussian") {

  # Validate data
  required_vars <- c(outcome_var, treatment_var, study_var, effect_modifiers)
  missing_vars <- setdiff(required_vars, names(ipd_data))

  if (length(missing_vars) > 0) {
    stop(sprintf("Missing variables in IPD data: %s",
                paste(missing_vars, collapse = ", ")))
  }

  message("Fitting IPD network meta-regression...")
  message(sprintf("  Outcome: %s", outcome_var))
  message(sprintf("  Treatments: %d", length(unique(ipd_data[[treatment_var]]))))
  message(sprintf("  Studies: %d", length(unique(ipd_data[[study_var]]))))
  message(sprintf("  Patients: %d", nrow(ipd_data)))
  message(sprintf("  Effect modifiers: %s", paste(effect_modifiers, collapse = ", ")))

  # Set reference
  if (is.null(reference)) {
    reference <- sort(unique(ipd_data[[treatment_var]]))[1]
    message(sprintf("\nUsing '%s' as reference treatment", reference))
  }

  # Relevel treatment factor
  ipd_data[[treatment_var]] <- relevel(factor(ipd_data[[treatment_var]]), ref = reference)

  # Build formula with interactions
  # Main effects: treatment + covariates
  # Interactions: treatment * each covariate

  formula_parts <- c(
    outcome_var,
    " ~ ",
    treatment_var,
    " + ",
    paste(effect_modifiers, collapse = " + ")
  )

  # Add treatment-covariate interactions
  interactions <- paste(paste(treatment_var, effect_modifiers, sep = ":"), collapse = " + ")
  formula_parts <- c(formula_parts, " + ", interactions)

  # Add random effect for study
  formula_parts <- c(formula_parts, " + (1|", study_var, ")")

  formula_str <- paste(formula_parts, collapse = "")
  model_formula <- as.formula(formula_str)

  message("\nModel formula:")
  message(sprintf("  %s", formula_str))

  # Fit mixed-effects model
  if (requireNamespace("lme4", quietly = TRUE)) {

    if (family == "gaussian") {
      model <- lme4::lmer(model_formula, data = ipd_data)
    } else if (family == "binomial") {
      model <- lme4::glmer(model_formula, data = ipd_data, family = binomial())
    } else if (family == "poisson") {
      model <- lme4::glmer(model_formula, data = ipd_data, family = poisson())
    } else {
      stop("Family must be 'gaussian', 'binomial', or 'poisson'")
    }

  } else {
    stop("Package 'lme4' required for IPD meta-regression. Please install it.")
  }

  # Extract treatment effects and interactions
  coef_summary <- summary(model)$coefficients

  # Identify treatment effects (main and interaction terms)
  treatment_levels <- levels(ipd_data[[treatment_var]])
  treatment_effects <- coef_summary[grepl(treatment_var, rownames(coef_summary)), , drop = FALSE]

  # Identify effect modifier interactions
  modifier_interactions <- list()

  for (modifier in effect_modifiers) {
    interaction_pattern <- paste0(treatment_var, ".*:", modifier, "|", modifier, ":", treatment_var)
    interactions <- coef_summary[grepl(interaction_pattern, rownames(coef_summary)), , drop = FALSE]

    if (nrow(interactions) > 0) {
      modifier_interactions[[modifier]] <- interactions
    }
  }

  result <- list(
    model = model,
    formula = model_formula,
    treatment_effects = treatment_effects,
    modifier_interactions = modifier_interactions,
    effect_modifiers = effect_modifiers,
    reference = reference,
    family = family,
    ipd_data = ipd_data,
    outcome_var = outcome_var,
    treatment_var = treatment_var,
    study_var = study_var
  )

  class(result) <- c("ipd_nma", "list")

  message("\n✓ IPD meta-regression completed")

  result
}

#' Predict Personalized Treatment Effects
#'
#' Generate personalized treatment effect predictions for new patients.
#'
#' @param ipd_model ipd_nma object from fit_ipd_metaregression()
#' @param newdata Data frame with new patient characteristics
#' @param treatments Treatments to compare (if NULL, all treatments)
#' @param reference Reference treatment for comparisons
#' @return Data frame with personalized predictions
#' @export
#' @examples
#' \dontrun{
#' # New patient
#' new_patient <- data.frame(
#'   age = 65,
#'   baseline_severity = 1.5
#' )
#'
#' predictions <- predict_personalized_effects(
#'   ipd_model = ipd_model,
#'   newdata = new_patient
#' )
#' print(predictions)
#' }
predict_personalized_effects <- function(ipd_model, newdata, treatments = NULL,
                                        reference = NULL) {

  if (!inherits(ipd_model, "ipd_nma")) {
    stop("Input must be ipd_nma object from fit_ipd_metaregression()")
  }

  # Get all treatments
  all_treatments <- levels(ipd_model$ipd_data[[ipd_model$treatment_var]])

  if (is.null(treatments)) {
    treatments <- all_treatments
  }

  if (is.null(reference)) {
    reference <- ipd_model$reference
  }

  # Validate newdata has required covariates
  missing_covars <- setdiff(ipd_model$effect_modifiers, names(newdata))
  if (length(missing_covars) > 0) {
    stop(sprintf("New data missing covariates: %s",
                paste(missing_covars, collapse = ", ")))
  }

  message(sprintf("Generating personalized predictions for %d patient(s)...", nrow(newdata)))

  # Generate predictions for each treatment
  predictions_list <- list()

  for (i in seq_len(nrow(newdata))) {
    patient_data <- newdata[i, , drop = FALSE]

    treatment_predictions <- lapply(treatments, function(trt) {
      # Create prediction data
      pred_data <- patient_data
      pred_data[[ipd_model$study_var]] <- NA  # Average over studies
      pred_data[[ipd_model$treatment_var]] <- factor(trt, levels = all_treatments)

      # Predict
      pred_value <- predict(ipd_model$model, newdata = pred_data, re.form = NA)

      # Get prediction for reference
      pred_data_ref <- patient_data
      pred_data_ref[[ipd_model$study_var]] <- NA
      pred_data_ref[[ipd_model$treatment_var]] <- factor(reference, levels = all_treatments)
      pred_ref <- predict(ipd_model$model, newdata = pred_data_ref, re.form = NA)

      # Treatment effect = difference from reference
      effect <- pred_value - pred_ref

      data.frame(
        Patient = i,
        Treatment = trt,
        Predicted_Outcome = pred_value,
        Effect_vs_Reference = effect,
        stringsAsFactors = FALSE
      )
    })

    predictions_list[[i]] <- do.call(rbind, treatment_predictions)
  }

  predictions <- do.call(rbind, predictions_list)

  # Add patient characteristics
  patient_chars <- newdata[rep(seq_len(nrow(newdata)), each = length(treatments)), , drop = FALSE]
  predictions <- cbind(predictions, patient_chars)

  # Rank treatments for each patient
  predictions$Rank <- ave(predictions$Effect_vs_Reference,
                         predictions$Patient,
                         FUN = function(x) rank(-x))

  predictions
}

#' Identify Best Treatment for Patient
#'
#' Recommend the optimal treatment based on patient characteristics.
#'
#' @param ipd_model ipd_nma object
#' @param newdata Patient characteristics
#' @param outcome_direction Direction for optimization ("maximize" or "minimize")
#' @return Data frame with treatment recommendations
#' @export
identify_best_treatment <- function(ipd_model, newdata,
                                   outcome_direction = c("maximize", "minimize")) {

  outcome_direction <- match.arg(outcome_direction)

  predictions <- predict_personalized_effects(ipd_model, newdata)

  # Find best treatment for each patient
  best_treatments <- lapply(unique(predictions$Patient), function(pid) {
    patient_preds <- predictions[predictions$Patient == pid, ]

    if (outcome_direction == "maximize") {
      best_idx <- which.max(patient_preds$Effect_vs_Reference)
    } else {
      best_idx <- which.min(patient_preds$Effect_vs_Reference)
    }

    best <- patient_preds[best_idx, ]

    # Get top 3
    top3_idx <- order(patient_preds$Effect_vs_Reference,
                     decreasing = (outcome_direction == "maximize"))[1:min(3, nrow(patient_preds))]
    top3 <- patient_preds$Treatment[top3_idx]

    data.frame(
      Patient = pid,
      Recommended_Treatment = best$Treatment,
      Expected_Effect = best$Effect_vs_Reference,
      Alternative_1 = if (length(top3) > 1) top3[2] else NA,
      Alternative_2 = if (length(top3) > 2) top3[3] else NA,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, best_treatments)

  # Add patient characteristics
  result <- merge(result, newdata, by.x = "Patient", by.y = "row.names", all.x = TRUE)

  result
}

#' Visualize Treatment by Subgroup
#'
#' Plot how treatment effects vary across levels of an effect modifier.
#'
#' @param ipd_model ipd_nma object
#' @param effect_modifier Name of effect modifier to visualize
#' @param treatments Treatments to include
#' @return ggplot2 object
#' @export
plot_treatment_by_subgroup <- function(ipd_model, effect_modifier, treatments = NULL) {

  if (!effect_modifier %in% ipd_model$effect_modifiers) {
    stop(sprintf("'%s' not an effect modifier in model", effect_modifier))
  }

  # Get all treatments
  all_treatments <- levels(ipd_model$ipd_data[[ipd_model$treatment_var]])

  if (is.null(treatments)) {
    treatments <- setdiff(all_treatments, ipd_model$reference)
  }

  # Determine range of effect modifier
  modifier_values <- ipd_model$ipd_data[[effect_modifier]]

  if (is.numeric(modifier_values)) {
    # Create grid of values
    modifier_range <- seq(min(modifier_values, na.rm = TRUE),
                         max(modifier_values, na.rm = TRUE),
                         length.out = 50)

    # Create newdata
    newdata <- data.frame(
      modifier = modifier_range
    )
    names(newdata) <- effect_modifier

    # Add mean values for other modifiers
    for (mod in setdiff(ipd_model$effect_modifiers, effect_modifier)) {
      if (is.numeric(ipd_model$ipd_data[[mod]])) {
        newdata[[mod]] <- mean(ipd_model$ipd_data[[mod]], na.rm = TRUE)
      } else {
        newdata[[mod]] <- names(sort(table(ipd_model$ipd_data[[mod]]), decreasing = TRUE))[1]
      }
    }

    # Predict for each treatment
    predictions <- predict_personalized_effects(ipd_model, newdata, treatments = treatments)

    # Plot
    predictions[[effect_modifier]] <- newdata[[effect_modifier]][predictions$Patient]

    p <- ggplot2::ggplot(predictions, ggplot2::aes_string(x = effect_modifier,
                                                          y = "Effect_vs_Reference",
                                                          color = "Treatment")) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = sprintf("Treatment Effects by %s", effect_modifier),
        subtitle = sprintf("Reference: %s", ipd_model$reference),
        x = effect_modifier,
        y = "Treatment Effect vs Reference",
        color = "Treatment"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "bottom"
      )

  } else {
    # Categorical modifier
    modifier_levels <- unique(modifier_values)

    predictions_list <- list()

    for (level in modifier_levels) {
      newdata <- data.frame(level = level)
      names(newdata) <- effect_modifier

      # Add defaults for other modifiers
      for (mod in setdiff(ipd_model$effect_modifiers, effect_modifier)) {
        if (is.numeric(ipd_model$ipd_data[[mod]])) {
          newdata[[mod]] <- mean(ipd_model$ipd_data[[mod]], na.rm = TRUE)
        } else {
          newdata[[mod]] <- names(sort(table(ipd_model$ipd_data[[mod]]), decreasing = TRUE))[1]
        }
      }

      preds <- predict_personalized_effects(ipd_model, newdata, treatments = treatments)
      preds[[effect_modifier]] <- level

      predictions_list[[level]] <- preds
    }

    predictions <- do.call(rbind, predictions_list)

    p <- ggplot2::ggplot(predictions, ggplot2::aes_string(x = effect_modifier,
                                                          y = "Effect_vs_Reference",
                                                          fill = "Treatment")) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = sprintf("Treatment Effects by %s", effect_modifier),
        subtitle = sprintf("Reference: %s", ipd_model$reference),
        x = effect_modifier,
        y = "Treatment Effect vs Reference",
        fill = "Treatment"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
  }

  print(p)
  invisible(p)
}

#' Test for Effect Modification
#'
#' Test whether effect modifiers significantly interact with treatment.
#'
#' @param ipd_model ipd_nma object
#' @return Data frame with interaction tests
#' @export
test_effect_modification <- function(ipd_model) {

  if (!inherits(ipd_model, "ipd_nma")) {
    stop("Input must be ipd_nma object")
  }

  message("Testing for treatment-covariate interactions...")

  # Extract interaction terms
  coef_summary <- summary(ipd_model$model)$coefficients

  interaction_tests <- list()

  for (modifier in ipd_model$effect_modifiers) {
    interactions <- ipd_model$modifier_interactions[[modifier]]

    if (!is.null(interactions) && nrow(interactions) > 0) {

      # Extract p-values
      p_values <- interactions[, "Pr(>|t|)"]  # For lmer
      if (is.null(p_values)) {
        p_values <- interactions[, "Pr(>|z|)"]  # For glmer
      }

      min_p <- min(p_values, na.rm = TRUE)
      any_significant <- any(p_values < 0.05, na.rm = TRUE)

      interaction_tests[[modifier]] <- data.frame(
        Effect_Modifier = modifier,
        N_Interactions = nrow(interactions),
        Min_P_Value = min_p,
        Any_Significant = any_significant,
        Interpretation = if (any_significant) {
          "Treatment effects vary by this modifier"
        } else {
          "No evidence of effect modification"
        },
        stringsAsFactors = FALSE
      )
    } else {
      interaction_tests[[modifier]] <- data.frame(
        Effect_Modifier = modifier,
        N_Interactions = 0,
        Min_P_Value = NA,
        Any_Significant = FALSE,
        Interpretation = "No interactions in model",
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, interaction_tests)
  rownames(result) <- NULL

  result
}

#' Print IPD NMA Results
#'
#' @param x ipd_nma object
#' @param ... Additional arguments
#' @export
print.ipd_nma <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Individual Patient Data Network Meta-Regression\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Outcome: %s (%s)\n", x$outcome_var, x$family))
  cat(sprintf("Reference treatment: %s\n", x$reference))
  cat(sprintf("Number of patients: %d\n", nrow(x$ipd_data)))
  cat(sprintf("Number of studies: %d\n", length(unique(x$ipd_data[[x$study_var]]))))
  cat(sprintf("Number of treatments: %d\n\n", length(unique(x$ipd_data[[x$treatment_var]]))))

  cat("Effect Modifiers:\n")
  cat(sprintf("  %s\n", paste(x$effect_modifiers, collapse = ", ")))
  cat("\n")

  cat("Treatment Effects (Main Effects):\n\n")
  print(x$treatment_effects, digits = 3)

  cat("\n")
  cat("Treatment-Covariate Interactions:\n\n")

  if (length(x$modifier_interactions) == 0) {
    cat("  No significant interactions detected\n")
  } else {
    for (modifier in names(x$modifier_interactions)) {
      cat(sprintf("\n  %s:\n", modifier))
      print(x$modifier_interactions[[modifier]], digits = 3)
    }
  }

  cat("\n")
  cat("Use predict_personalized_effects() for new patients\n")
  cat("Use identify_best_treatment() for treatment recommendations\n")
  cat("Use plot_treatment_by_subgroup() to visualize effect modification\n\n")

  invisible(x)
}

#' Create Treatment Selection Heatmap
#'
#' Visualize which treatment is best for different patient profiles.
#'
#' @param ipd_model ipd_nma object
#' @param modifier1 First effect modifier
#' @param modifier2 Second effect modifier
#' @param outcome_direction Direction for optimization
#' @return ggplot2 object
#' @export
plot_treatment_selection_heatmap <- function(ipd_model, modifier1, modifier2,
                                            outcome_direction = c("maximize", "minimize")) {

  outcome_direction <- match.arg(outcome_direction)

  if (!all(c(modifier1, modifier2) %in% ipd_model$effect_modifiers)) {
    stop("Both modifiers must be in the model")
  }

  # Check both are numeric
  if (!is.numeric(ipd_model$ipd_data[[modifier1]]) ||
     !is.numeric(ipd_model$ipd_data[[modifier2]])) {
    stop("Both modifiers must be numeric for heatmap")
  }

  # Create grid
  m1_seq <- seq(min(ipd_model$ipd_data[[modifier1]], na.rm = TRUE),
               max(ipd_model$ipd_data[[modifier1]], na.rm = TRUE),
               length.out = 20)

  m2_seq <- seq(min(ipd_model$ipd_data[[modifier2]], na.rm = TRUE),
               max(ipd_model$ipd_data[[modifier2]], na.rm = TRUE),
               length.out = 20)

  grid <- expand.grid(m1 = m1_seq, m2 = m2_seq)
  names(grid) <- c(modifier1, modifier2)

  # Add defaults for other modifiers
  for (mod in setdiff(ipd_model$effect_modifiers, c(modifier1, modifier2))) {
    if (is.numeric(ipd_model$ipd_data[[mod]])) {
      grid[[mod]] <- mean(ipd_model$ipd_data[[mod]], na.rm = TRUE)
    } else {
      grid[[mod]] <- names(sort(table(ipd_model$ipd_data[[mod]]), decreasing = TRUE))[1]
    }
  }

  # Get best treatment for each grid point
  best <- identify_best_treatment(ipd_model, grid, outcome_direction)

  grid$Best_Treatment <- best$Recommended_Treatment

  p <- ggplot2::ggplot(grid, ggplot2::aes_string(x = modifier1, y = modifier2,
                                                 fill = "Best_Treatment")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_brewer(palette = "Set3") +
    ggplot2::labs(
      title = "Personalized Treatment Selection",
      subtitle = sprintf("Optimal treatment by patient characteristics (to %s outcome)",
                        outcome_direction),
      x = modifier1,
      y = modifier2,
      fill = "Best Treatment"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  print(p)
  invisible(p)
}
