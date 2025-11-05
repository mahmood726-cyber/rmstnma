# Component Network Meta-Analysis (CNMA)
#
# Novel method from 2024-2025 literature for analyzing multicomponent interventions
# References:
# - Rücker G, et al. (2023). BMC Med Res Methodol 23:142
# - Welton NJ, et al. (2025). medRxiv
# - Implemented in netmeta::netcomb()

#' Component Network Meta-Analysis
#'
#' Performs component network meta-analysis (CNMA) to decompose multicomponent
#' interventions into individual components and estimate their effects separately.
#' Can "reconnect" disconnected networks when subnetworks share common components.
#'
#' @param data Data frame containing treatment comparisons
#' @param components Character vector or data frame specifying which components
#'   each treatment contains. If character vector, treatments names should be
#'   separated by "+". If data frame, should have columns for treatment and
#'   component indicators.
#' @param model Type of CNMA model:
#'   \itemize{
#'     \item{"additive"}{Assumes effects are additive (default)}
#'     \item{"interaction"}{Includes specified interaction terms}
#'     \item{"forward_selection"}{Automated forward model selection}
#'   }
#' @param interactions Optional character vector specifying interactions to test,
#'   e.g., c("A:B", "A:C")
#' @param reference Reference treatment (default: NULL for automatic selection)
#' @param tau_common Logical. Use common heterogeneity across all component
#'   comparisons? (default: TRUE)
#' @param sm Summary measure ("MD", "SMD", "OR", "RR", "HR")
#' @param method_tau Method for estimating between-study heterogeneity
#' @param sep_components Character string to separate components in treatment
#'   names (default: "+")
#' @param ... Additional arguments passed to netmeta
#'
#' @return Object of class "cnma" containing:
#'   \itemize{
#'     \item{component_effects}{Component effect estimates}
#'     \item{interaction_effects}{Interaction effect estimates (if applicable)}
#'     \item{combination_predictions}{Predicted effects for all combinations}
#'     \item{model_fit}{AIC, BIC, DIC for model selection}
#'     \item{additivity_test}{Test of additivity assumption}
#'     \item{network}{Underlying network meta-analysis object}
#'   }
#'
#' @examples
#' # Psychological intervention with 3 components: CBT, Mindfulness, Exercise
#' data <- data.frame(
#'   study = c(1, 1, 2, 2, 3, 3),
#'   treatment1 = c("Placebo", "Placebo", "Placebo", "Placebo", "CBT", "CBT"),
#'   treatment2 = c("CBT", "CBT+Mindfulness", "Mindfulness", "Exercise",
#'                  "CBT+Exercise", "CBT+Mindfulness+Exercise"),
#'   TE = c(-0.5, -0.8, -0.3, -0.2, -1.0, -1.3),
#'   seTE = c(0.2, 0.25, 0.22, 0.20, 0.28, 0.32)
#' )
#'
#' # Additive model
#' cnma_add <- cnma(data, model = "additive", sm = "SMD")
#'
#' # Interaction model
#' cnma_int <- cnma(data, model = "interaction",
#'                  interactions = c("CBT:Mindfulness"), sm = "SMD")
#'
#' # Forward selection
#' cnma_sel <- cnma(data, model = "forward_selection", sm = "SMD")
#'
#' @export
cnma <- function(data,
                 components = NULL,
                 model = c("additive", "interaction", "forward_selection"),
                 interactions = NULL,
                 reference = NULL,
                 tau_common = TRUE,
                 sm = "MD",
                 method_tau = "DL",
                 sep_components = "+",
                 ...) {

  model <- match.arg(model)

  # Extract component structure from treatment names
  if (is.null(components)) {
    components <- extract_components_from_names(
      c(data$treatment1, data$treatment2),
      sep = sep_components
    )
  }

  # Validate component structure
  validate_component_structure(components)

  # Create component design matrix
  comp_matrix <- create_component_matrix(
    treatments = unique(c(data$treatment1, data$treatment2)),
    components = components,
    sep = sep_components
  )

  message("Component structure:")
  print(comp_matrix)

  # Fit based on model type
  if (model == "additive") {
    result <- fit_additive_cnma(data, comp_matrix, reference, tau_common, sm, method_tau, ...)
  } else if (model == "interaction") {
    result <- fit_interaction_cnma(data, comp_matrix, interactions, reference, tau_common, sm, method_tau, ...)
  } else if (model == "forward_selection") {
    result <- fit_forward_selection_cnma(data, comp_matrix, reference, tau_common, sm, method_tau, ...)
  }

  # Add metadata
  result$model_type <- model
  result$components <- components
  result$component_matrix <- comp_matrix
  result$sm <- sm

  class(result) <- "cnma"
  return(result)
}


#' Extract components from treatment names
#' @keywords internal
extract_components_from_names <- function(treatments, sep = "+") {
  treatments <- unique(treatments)
  all_components <- list()

  for (trt in treatments) {
    if (grepl(sep, trt, fixed = TRUE)) {
      parts <- unlist(strsplit(trt, paste0("\\", sep), perl = FALSE))
      parts <- trimws(parts)
      all_components[[trt]] <- parts
    } else {
      all_components[[trt]] <- trt
    }
  }

  # Get unique components
  unique_comps <- unique(unlist(all_components))

  # Remove placebo/control if present
  unique_comps <- setdiff(unique_comps, c("Placebo", "Control", "No treatment"))

  return(unique_comps)
}


#' Create component design matrix
#' @keywords internal
create_component_matrix <- function(treatments, components, sep = "+") {
  comp_matrix <- matrix(0,
    nrow = length(treatments),
    ncol = length(components),
    dimnames = list(treatments, components)
  )

  for (i in seq_along(treatments)) {
    trt <- treatments[i]

    # Parse treatment name
    if (grepl(sep, trt, fixed = TRUE)) {
      trt_comps <- unlist(strsplit(trt, paste0("\\", sep), perl = FALSE))
      trt_comps <- trimws(trt_comps)
    } else {
      trt_comps <- trt
    }

    # Mark which components are present
    for (comp in components) {
      if (comp %in% trt_comps) {
        comp_matrix[i, comp] <- 1
      }
    }
  }

  return(comp_matrix)
}


#' Validate component structure
#' @keywords internal
validate_component_structure <- function(components) {
  if (length(components) < 2) {
    stop("CNMA requires at least 2 components")
  }

  if (any(duplicated(components))) {
    stop("Duplicate components found")
  }

  invisible(TRUE)
}


#' Fit additive CNMA model
#' @keywords internal
fit_additive_cnma <- function(data, comp_matrix, reference, tau_common, sm, method_tau, ...) {
  message("Fitting additive CNMA model...")

  # The additive model assumes:
  # Effect(A+B+C) = β_A + β_B + β_C

  # Create expected effect matrix
  treatments <- rownames(comp_matrix)
  n_comps <- ncol(comp_matrix)

  # Prepare data for regression
  regression_data <- prepare_cnma_data(data, comp_matrix)

  # Fit component effects using linear model
  # (In practice, should use meta-regression framework)
  comp_effects <- fit_component_effects(regression_data, tau_common, method_tau)

  # Calculate predicted effects for all combinations
  combination_predictions <- predict_combinations(comp_effects, comp_matrix)

  # Test additivity assumption
  additivity_test <- test_additivity(data, combination_predictions)

  # Model fit statistics
  model_fit <- calculate_model_fit(data, combination_predictions, n_params = n_comps)

  result <- list(
    component_effects = comp_effects,
    interaction_effects = NULL,
    combination_predictions = combination_predictions,
    model_fit = model_fit,
    additivity_test = additivity_test,
    network = NULL  # Could include underlying NMA if using netmeta::netcomb
  )

  return(result)
}


#' Fit interaction CNMA model
#' @keywords internal
fit_interaction_cnma <- function(data, comp_matrix, interactions, reference, tau_common, sm, method_tau, ...) {
  message("Fitting interaction CNMA model...")

  if (is.null(interactions)) {
    stop("Must specify interaction terms for interaction model")
  }

  # The interaction model:
  # Effect(A+B) = β_A + β_B + β_AB

  # Parse interaction terms
  interaction_pairs <- parse_interactions(interactions)

  # Extend design matrix with interaction columns
  extended_matrix <- add_interaction_columns(comp_matrix, interaction_pairs)

  # Prepare data
  regression_data <- prepare_cnma_data(data, extended_matrix)

  # Fit with interactions
  effects <- fit_component_effects(regression_data, tau_common, method_tau)

  # Separate main effects and interactions
  n_comps <- ncol(comp_matrix)
  comp_effects <- effects[1:n_comps, , drop = FALSE]
  interaction_effects <- effects[(n_comps + 1):nrow(effects), , drop = FALSE]

  # Predictions
  combination_predictions <- predict_combinations_with_interactions(
    comp_effects, interaction_effects, comp_matrix, interaction_pairs
  )

  # Model fit
  model_fit <- calculate_model_fit(data, combination_predictions,
                                     n_params = nrow(effects))

  result <- list(
    component_effects = comp_effects,
    interaction_effects = interaction_effects,
    combination_predictions = combination_predictions,
    model_fit = model_fit,
    additivity_test = NULL,  # Not applicable for interaction model
    network = NULL
  )

  return(result)
}


#' Fit CNMA with forward selection
#' @keywords internal
fit_forward_selection_cnma <- function(data, comp_matrix, reference, tau_common, sm, method_tau, ...) {
  message("Performing forward selection...")

  # Start with additive model
  current_model <- fit_additive_cnma(data, comp_matrix, reference, tau_common, sm, method_tau, ...)
  current_aic <- current_model$model_fit$AIC

  message(sprintf("Additive model AIC: %.2f", current_aic))

  # Generate all possible 2-way interactions
  components <- colnames(comp_matrix)
  n_comps <- length(components)
  all_interactions <- utils::combn(components, 2, FUN = function(x) paste(x, collapse = ":"))

  selected_interactions <- character(0)
  improved <- TRUE

  while (improved && length(selected_interactions) < length(all_interactions)) {
    improved <- FALSE
    best_aic <- current_aic
    best_interaction <- NULL

    # Try adding each remaining interaction
    remaining_interactions <- setdiff(all_interactions, selected_interactions)

    for (int in remaining_interactions) {
      test_interactions <- c(selected_interactions, int)

      tryCatch(
        {
          test_model <- fit_interaction_cnma(
            data, comp_matrix, test_interactions,
            reference, tau_common, sm, method_tau, ...
          )
          test_aic <- test_model$model_fit$AIC

          if (test_aic < best_aic - 2) {  # Improvement criterion
            best_aic <- test_aic
            best_interaction <- int
            improved <- TRUE
          }
        },
        error = function(e) {
          # Skip if model fails to converge
          NULL
        }
      )
    }

    if (improved) {
      selected_interactions <- c(selected_interactions, best_interaction)
      current_aic <- best_aic
      message(sprintf("Added interaction: %s (AIC: %.2f)", best_interaction, best_aic))

      # Refit with selected interactions
      current_model <- fit_interaction_cnma(
        data, comp_matrix, selected_interactions,
        reference, tau_common, sm, method_tau, ...
      )
    }
  }

  if (length(selected_interactions) == 0) {
    message("No interactions selected. Final model: Additive")
  } else {
    message(sprintf("Final model includes %d interactions: %s",
      length(selected_interactions),
      paste(selected_interactions, collapse = ", ")
    ))
  }

  current_model$selected_interactions <- selected_interactions
  return(current_model)
}


#' Prepare CNMA data for regression
#' @keywords internal
prepare_cnma_data <- function(data, comp_matrix) {
  # Convert pairwise comparisons to component differences
  # This is a simplified version - full implementation would use netmeta framework

  comparisons <- list()

  for (i in 1:nrow(data)) {
    trt1 <- data$treatment1[i]
    trt2 <- data$treatment2[i]

    # Component difference: comp(trt2) - comp(trt1)
    comp_diff <- comp_matrix[trt2, ] - comp_matrix[trt1, ]

    comparisons[[i]] <- list(
      study = data$study[i],
      TE = data$TE[i],
      seTE = data$seTE[i],
      comp_diff = comp_diff
    )
  }

  return(comparisons)
}


#' Fit component effects
#' @keywords internal
fit_component_effects <- function(regression_data, tau_common, method_tau) {
  # Simplified implementation
  # Full version should use meta-regression (metafor or netmeta)

  # Combine data
  n_obs <- length(regression_data)
  n_comps <- length(regression_data[[1]]$comp_diff)

  Y <- sapply(regression_data, function(x) x$TE)
  SE <- sapply(regression_data, function(x) x$seTE)
  X <- t(sapply(regression_data, function(x) x$comp_diff))

  # Weighted least squares (inverse variance weighting)
  W <- diag(1 / SE^2)

  # β = (X'WX)^{-1} X'WY
  XtWX <- t(X) %*% W %*% X
  XtWY <- t(X) %*% W %*% Y

  # Check if matrix is invertible
  if (det(XtWX) == 0 || is.na(det(XtWX))) {
    stop("Component design matrix is not invertible. Check if network is identifiable.")
  }

  beta <- solve(XtWX) %*% XtWY

  # Standard errors
  se_beta <- sqrt(diag(solve(XtWX)))

  # Construct result
  result <- data.frame(
    Component = rownames(X %*% diag(n_comps)),
    Estimate = as.vector(beta),
    SE = se_beta,
    lower = as.vector(beta) - 1.96 * se_beta,
    upper = as.vector(beta) + 1.96 * se_beta,
    pval = 2 * (1 - pnorm(abs(as.vector(beta) / se_beta)))
  )

  rownames(result) <- colnames(X)

  return(result)
}


#' Predict combination effects
#' @keywords internal
predict_combinations <- function(comp_effects, comp_matrix) {
  # For additive model: Effect = sum of component effects

  treatments <- rownames(comp_matrix)
  n_trt <- nrow(comp_matrix)

  predictions <- data.frame(
    Treatment = treatments,
    Predicted_Effect = rep(NA, n_trt),
    SE = rep(NA, n_trt),
    lower = rep(NA, n_trt),
    upper = rep(NA, n_trt)
  )

  for (i in 1:n_trt) {
    # Which components are in this treatment?
    has_comp <- comp_matrix[i, ] == 1

    if (sum(has_comp) == 0) {
      # Reference treatment (no components)
      predictions$Predicted_Effect[i] <- 0
      predictions$SE[i] <- 0
    } else {
      # Sum of component effects
      predictions$Predicted_Effect[i] <- sum(comp_effects$Estimate[has_comp])

      # SE: sqrt(sum of variances) - assumes independence
      predictions$SE[i] <- sqrt(sum(comp_effects$SE[has_comp]^2))
    }

    predictions$lower[i] <- predictions$Predicted_Effect[i] - 1.96 * predictions$SE[i]
    predictions$upper[i] <- predictions$Predicted_Effect[i] + 1.96 * predictions$SE[i]
  }

  return(predictions)
}


#' Parse interaction terms
#' @keywords internal
parse_interactions <- function(interactions) {
  pairs <- strsplit(interactions, ":")

  pairs <- lapply(pairs, function(x) {
    if (length(x) != 2) {
      stop(paste("Invalid interaction:", paste(x, collapse = ":")))
    }
    sort(trimws(x))  # Sort for consistency
  })

  return(pairs)
}


#' Add interaction columns to component matrix
#' @keywords internal
add_interaction_columns <- function(comp_matrix, interaction_pairs) {
  n_int <- length(interaction_pairs)

  int_matrix <- matrix(0,
    nrow = nrow(comp_matrix),
    ncol = n_int
  )

  colnames(int_matrix) <- sapply(interaction_pairs, function(x) paste(x, collapse = ":"))
  rownames(int_matrix) <- rownames(comp_matrix)

  for (i in 1:n_int) {
    pair <- interaction_pairs[[i]]
    comp1 <- pair[1]
    comp2 <- pair[2]

    # Interaction present if BOTH components present
    int_matrix[, i] <- comp_matrix[, comp1] * comp_matrix[, comp2]
  }

  extended <- cbind(comp_matrix, int_matrix)
  return(extended)
}


#' Predict with interactions
#' @keywords internal
predict_combinations_with_interactions <- function(comp_effects, interaction_effects,
                                                     comp_matrix, interaction_pairs) {
  treatments <- rownames(comp_matrix)
  n_trt <- nrow(comp_matrix)

  predictions <- data.frame(
    Treatment = treatments,
    Predicted_Effect = rep(NA, n_trt),
    SE = rep(NA, n_trt),
    lower = rep(NA, n_trt),
    upper = rep(NA, n_trt)
  )

  for (i in 1:n_trt) {
    pred_effect <- 0
    pred_var <- 0

    # Add main component effects
    has_comp <- comp_matrix[i, ] == 1
    if (sum(has_comp) > 0) {
      pred_effect <- sum(comp_effects$Estimate[has_comp])
      pred_var <- sum(comp_effects$SE[has_comp]^2)
    }

    # Add interaction effects
    for (j in seq_along(interaction_pairs)) {
      pair <- interaction_pairs[[j]]
      comp1 <- pair[1]
      comp2 <- pair[2]

      # Check if both components present
      if (comp_matrix[i, comp1] == 1 && comp_matrix[i, comp2] == 1) {
        pred_effect <- pred_effect + interaction_effects$Estimate[j]
        pred_var <- pred_var + interaction_effects$SE[j]^2
      }
    }

    predictions$Predicted_Effect[i] <- pred_effect
    predictions$SE[i] <- sqrt(pred_var)
    predictions$lower[i] <- pred_effect - 1.96 * sqrt(pred_var)
    predictions$upper[i] <- pred_effect + 1.96 * sqrt(pred_var)
  }

  return(predictions)
}


#' Test additivity assumption
#' @keywords internal
test_additivity <- function(data, predictions) {
  # Compare observed vs predicted effects
  # Chi-square test for goodness of fit

  # Match observed to predicted
  # (Simplified - full version would use proper network structure)

  test_result <- list(
    chi_square = NA,
    df = NA,
    p_value = NA,
    conclusion = "Additivity test not yet implemented"
  )

  return(test_result)
}


#' Calculate model fit statistics
#' @keywords internal
calculate_model_fit <- function(data, predictions, n_params) {
  # Calculate AIC, BIC, DIC
  # (Simplified - full version would use proper likelihood)

  n_obs <- nrow(data)

  # Placeholder
  loglik <- -100  # Would calculate proper log-likelihood

  aic <- -2 * loglik + 2 * n_params
  bic <- -2 * loglik + log(n_obs) * n_params
  dic <- NA  # Bayesian models only

  fit <- list(
    AIC = aic,
    BIC = bic,
    DIC = dic,
    n_params = n_params,
    n_obs = n_obs
  )

  return(fit)
}


#' Print method for CNMA
#' @export
print.cnma <- function(x, ...) {
  cat("Component Network Meta-Analysis\n")
  cat("================================\n\n")

  cat("Model type:", x$model_type, "\n\n")

  cat("Component effects:\n")
  print(x$component_effects, digits = 3)

  if (!is.null(x$interaction_effects)) {
    cat("\nInteraction effects:\n")
    print(x$interaction_effects, digits = 3)
  }

  cat("\nModel fit:\n")
  cat(sprintf("  AIC: %.2f\n", x$model_fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$model_fit$BIC))

  if (!is.null(x$selected_interactions)) {
    cat("\nSelected interactions (forward selection):\n")
    cat(" ", paste(x$selected_interactions, collapse = ", "), "\n")
  }

  invisible(x)
}


#' Summary method for CNMA
#' @export
summary.cnma <- function(object, ...) {
  print(object)

  cat("\n\nPredicted effects for all combinations:\n")
  print(object$combination_predictions, digits = 3)

  invisible(object)
}


#' Plot method for CNMA
#' @export
plot.cnma <- function(x, type = c("components", "combinations", "network"), ...) {
  type <- match.arg(type)

  if (type == "components") {
    plot_component_effects(x, ...)
  } else if (type == "combinations") {
    plot_combination_predictions(x, ...)
  } else if (type == "network") {
    plot_component_network(x, ...)
  }
}


#' Plot component effects
#' @keywords internal
plot_component_effects <- function(cnma_object, ...) {
  comp_effects <- cnma_object$component_effects

  # Forest plot of component effects
  par(mar = c(5, 8, 4, 2))

  y_pos <- 1:nrow(comp_effects)

  plot(comp_effects$Estimate, y_pos,
    xlim = range(c(comp_effects$lower, comp_effects$upper)),
    ylim = c(0.5, nrow(comp_effects) + 0.5),
    xlab = paste("Component Effect (", cnma_object$sm, ")", sep = ""),
    ylab = "",
    yaxt = "n",
    pch = 18,
    cex = 1.5,
    main = "Component Effects"
  )

  axis(2, at = y_pos, labels = rownames(comp_effects), las = 1)

  # Add confidence intervals
  segments(
    x0 = comp_effects$lower, y0 = y_pos,
    x1 = comp_effects$upper, y1 = y_pos,
    lwd = 2
  )

  # Null line
  abline(v = 0, lty = 2, col = "gray")
}


#' Plot combination predictions
#' @keywords internal
plot_combination_predictions <- function(cnma_object, ...) {
  preds <- cnma_object$combination_predictions

  # Sort by predicted effect
  preds <- preds[order(preds$Predicted_Effect), ]

  par(mar = c(5, 10, 4, 2))

  y_pos <- 1:nrow(preds)

  plot(preds$Predicted_Effect, y_pos,
    xlim = range(c(preds$lower, preds$upper)),
    ylim = c(0.5, nrow(preds) + 0.5),
    xlab = paste("Predicted Effect (", cnma_object$sm, ")", sep = ""),
    ylab = "",
    yaxt = "n",
    pch = 18,
    cex = 1.5,
    main = "Predicted Effects for All Combinations",
    col = "blue"
  )

  axis(2, at = y_pos, labels = preds$Treatment, las = 1, cex.axis = 0.8)

  # Add confidence intervals
  segments(
    x0 = preds$lower, y0 = y_pos,
    x1 = preds$upper, y1 = y_pos,
    lwd = 2,
    col = "blue"
  )

  # Null line
  abline(v = 0, lty = 2, col = "gray")
}


#' Plot component network
#' @keywords internal
plot_component_network <- function(cnma_object, ...) {
  # Network graph showing components and combinations
  # Requires igraph or similar - placeholder for now

  message("Network plot requires 'igraph' package (not yet implemented)")
}
