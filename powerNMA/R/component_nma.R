#' Component Network Meta-Analysis (CNMA)
#'
#' Analyze complex interventions by decomposing treatments into their individual
#' components and estimating component-level effects. Implements both additive
#' and interaction models following Veroniki et al. (2025).
#'
#' @name component_nma
#' @references
#' Veroniki AA, et al. (2025). Analysing complex interventions using component
#' network meta-analysis. medRxiv preprint.
#'
#' Welton NJ, et al. (2009). Models for potentially biased evidence in
#' meta-analysis using empirically based priors. Journal of the Royal
#' Statistical Society: Series A, 172(1):119-136.
NULL

#' Define Component Structure for Treatments
#'
#' Create a matrix specifying which components are present in each treatment.
#' This is the key input for component network meta-analysis.
#'
#' @param treatments Character vector of treatment names
#' @param components Character vector of component names
#' @param interactive Logical; if TRUE, prompts user to specify components
#' @return Matrix of 0s and 1s (treatments × components)
#' @export
#' @examples
#' treatments <- c("Placebo", "Drug_A", "Drug_B", "Drug_A+B")
#' components <- c("Component_A", "Component_B")
#' comp_matrix <- create_component_matrix(
#'   treatments = treatments,
#'   components = components
#' )
#' # Manually specify:
#' comp_matrix["Placebo", ] <- c(0, 0)
#' comp_matrix["Drug_A", ] <- c(1, 0)
#' comp_matrix["Drug_B", ] <- c(0, 1)
#' comp_matrix["Drug_A+B", ] <- c(1, 1)
create_component_matrix <- function(treatments, components, interactive = FALSE) {

  # Initialize matrix
  comp_matrix <- matrix(0, nrow = length(treatments), ncol = length(components),
                       dimnames = list(treatments, components))

  if (interactive) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  Component Matrix Setup\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    cat("For each treatment, specify which components are present (1) or absent (0).\n\n")

    for (i in seq_along(treatments)) {
      cat(sprintf("Treatment: %s\n", treatments[i]))
      for (j in seq_along(components)) {
        response <- readline(sprintf("  Contains '%s'? (1/0): ", components[j]))
        comp_matrix[i, j] <- as.numeric(response)
      }
      cat("\n")
    }
  }

  comp_matrix
}

#' Validate Component Matrix
#'
#' Check that component matrix is properly specified and network is connected.
#'
#' @param comp_matrix Component matrix (treatments × components)
#' @param data NMA data with treatment comparisons
#' @return List with validation results
#' @keywords internal
validate_component_matrix <- function(comp_matrix, data) {

  # Check for all-zero treatments (placebo should exist)
  all_zero <- which(rowSums(comp_matrix) == 0)
  if (length(all_zero) == 0) {
    warning("No placebo treatment (all components = 0) found. Consider adding one.")
  }

  # Check for duplicate component profiles
  duplicates <- duplicated(comp_matrix)
  if (any(duplicates)) {
    dup_treatments <- rownames(comp_matrix)[duplicates]
    stop(sprintf("Duplicate component profiles found: %s",
                paste(dup_treatments, collapse = ", ")))
  }

  # Check that all treatments in data are in component matrix
  treatments_in_data <- unique(c(as.character(data$treat1), as.character(data$treat2)))
  missing_treatments <- setdiff(treatments_in_data, rownames(comp_matrix))
  if (length(missing_treatments) > 0) {
    stop(sprintf("Treatments in data not found in component matrix: %s",
                paste(missing_treatments, collapse = ", ")))
  }

  # Check network connectivity
  # Build adjacency matrix
  treatments <- rownames(comp_matrix)
  adj_matrix <- matrix(0, nrow = length(treatments), ncol = length(treatments),
                      dimnames = list(treatments, treatments))

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])
    adj_matrix[t1, t2] <- 1
    adj_matrix[t2, t1] <- 1
  }

  # Use BFS to check connectivity
  visited <- rep(FALSE, length(treatments))
  names(visited) <- treatments
  queue <- treatments[1]
  visited[treatments[1]] <- TRUE

  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    neighbors <- treatments[adj_matrix[current, ] == 1]
    for (neighbor in neighbors) {
      if (!visited[neighbor]) {
        visited[neighbor] <- TRUE
        queue <- c(queue, neighbor)
      }
    }
  }

  if (!all(visited)) {
    disconnected <- treatments[!visited]
    warning(sprintf("Network is disconnected. Isolated treatments: %s",
                   paste(disconnected, collapse = ", ")))
  }

  list(
    valid = TRUE,
    has_placebo = length(all_zero) > 0,
    n_treatments = nrow(comp_matrix),
    n_components = ncol(comp_matrix),
    connected = all(visited)
  )
}

#' Additive Component Network Meta-Analysis
#'
#' Fit component NMA assuming additive effects: the effect of a multi-component
#' treatment equals the sum of its component effects.
#'
#' @param data Pairwise NMA data
#' @param comp_matrix Component matrix (treatments × components)
#' @param sm Summary measure
#' @param reference Reference treatment (typically placebo)
#' @param method Estimation method: "frequentist" or "bayesian"
#' @return cnma_additive object
#' @export
#' @examples
#' \dontrun{
#' # Define component structure
#' comp_matrix <- create_component_matrix(
#'   treatments = c("Placebo", "A", "B", "A+B"),
#'   components = c("Component_A", "Component_B")
#' )
#' comp_matrix["Placebo", ] <- c(0, 0)
#' comp_matrix["A", ] <- c(1, 0)
#' comp_matrix["B", ] <- c(0, 1)
#' comp_matrix["A+B", ] <- c(1, 1)
#'
#' # Fit additive CNMA
#' cnma_result <- additive_cnma(
#'   data = nma_data,
#'   comp_matrix = comp_matrix,
#'   reference = "Placebo"
#' )
#' print(cnma_result)
#' plot(cnma_result)
#' }
additive_cnma <- function(data, comp_matrix, sm = "MD", reference = NULL,
                         method = c("frequentist", "bayesian")) {

  method <- match.arg(method)

  # Validate inputs
  validation <- validate_component_matrix(comp_matrix, data)

  if (!validation$connected) {
    stop("Network is disconnected. Cannot perform CNMA.")
  }

  # Set reference to placebo if not specified
  if (is.null(reference)) {
    all_zero <- which(rowSums(comp_matrix) == 0)
    if (length(all_zero) > 0) {
      reference <- rownames(comp_matrix)[all_zero[1]]
      message(sprintf("Using '%s' as reference (placebo)", reference))
    } else {
      stop("No reference specified and no placebo found in component matrix")
    }
  }

  message("Fitting additive component network meta-analysis...")

  # Run standard NMA first
  base_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  # Extract treatment effects vs reference
  treatments <- rownames(base_nma$TE.random)
  treatment_effects <- base_nma$TE.random[, reference]
  treatment_se <- base_nma$seTE.random[, reference]

  # Build design matrix for component effects
  # Each treatment effect = sum of its component effects
  X <- comp_matrix[treatments, , drop = FALSE]

  # Weighted least squares to estimate component effects
  # y = X * beta, where y = treatment effects, beta = component effects
  y <- treatment_effects
  weights <- 1 / (treatment_se^2)

  # Weighted regression
  W <- diag(weights)
  XtWX <- t(X) %*% W %*% X
  XtWy <- t(X) %*% W %*% y

  # Check invertibility
  if (det(XtWX) < 1e-10) {
    stop("Component effects are not identifiable (singular design matrix). ",
         "This can occur if:\n",
         "  1. Network lacks sufficient diversity in component combinations\n",
         "  2. Some components always appear together\n",
         "  3. Network is disconnected")
  }

  component_effects <- solve(XtWX) %*% XtWy
  component_se <- sqrt(diag(solve(XtWX)))

  # Calculate confidence intervals
  component_lower <- component_effects - 1.96 * component_se
  component_upper <- component_effects + 1.96 * component_se

  # P-values
  component_z <- component_effects / component_se
  component_p <- 2 * (1 - pnorm(abs(component_z)))

  # Create results table
  component_results <- data.frame(
    Component = colnames(comp_matrix),
    Effect = as.vector(component_effects),
    SE = as.vector(component_se),
    Lower_95 = as.vector(component_lower),
    Upper_95 = as.vector(component_upper),
    Z_value = as.vector(component_z),
    P_value = as.vector(component_p),
    Significant = as.vector(component_p) < 0.05,
    stringsAsFactors = FALSE
  )

  # Predict treatment effects from component model
  predicted_effects <- X %*% component_effects
  residuals <- y - predicted_effects

  # Model fit statistics
  rss <- sum(weights * residuals^2)  # Weighted residual sum of squares
  total_ss <- sum(weights * (y - weighted.mean(y, weights))^2)
  r_squared <- 1 - (rss / total_ss)

  # Calculate consistency between standard NMA and component model
  consistency_check <- data.frame(
    Treatment = treatments,
    NMA_Effect = treatment_effects,
    CNMA_Predicted = as.vector(predicted_effects),
    Difference = as.vector(residuals),
    stringsAsFactors = FALSE
  )

  result <- list(
    method = "additive",
    component_effects = component_results,
    component_matrix = comp_matrix,
    base_nma = base_nma,
    consistency_check = consistency_check,
    model_fit = list(
      R_squared = r_squared,
      RSS = rss,
      n_components = ncol(comp_matrix),
      n_treatments = length(treatments)
    ),
    sm = sm,
    reference = reference
  )

  class(result) <- c("cnma_additive", "cnma", "list")
  result
}

#' Interaction Component Network Meta-Analysis
#'
#' Fit component NMA allowing for interaction effects between components.
#' This allows for synergistic or antagonistic effects.
#'
#' @param data Pairwise NMA data
#' @param comp_matrix Component matrix (treatments × components)
#' @param interactions Character vector of interaction terms (e.g., "A:B")
#' @param sm Summary measure
#' @param reference Reference treatment
#' @return cnma_interaction object
#' @export
#' @examples
#' \dontrun{
#' # Fit interaction CNMA
#' cnma_int <- interaction_cnma(
#'   data = nma_data,
#'   comp_matrix = comp_matrix,
#'   interactions = c("Component_A:Component_B"),
#'   reference = "Placebo"
#' )
#' }
interaction_cnma <- function(data, comp_matrix, interactions = NULL,
                            sm = "MD", reference = NULL) {

  # Validate inputs
  validation <- validate_component_matrix(comp_matrix, data)

  if (!validation$connected) {
    stop("Network is disconnected. Cannot perform CNMA.")
  }

  # Set reference
  if (is.null(reference)) {
    all_zero <- which(rowSums(comp_matrix) == 0)
    if (length(all_zero) > 0) {
      reference <- rownames(comp_matrix)[all_zero[1]]
    } else {
      stop("No reference specified and no placebo found")
    }
  }

  message("Fitting interaction component network meta-analysis...")

  # Run standard NMA
  base_nma <- netmeta::netmeta(
    TE = TE, seTE = seTE,
    treat1 = treat1, treat2 = treat2,
    studlab = studlab,
    data = data,
    sm = sm,
    reference.group = reference
  )

  treatments <- rownames(base_nma$TE.random)
  y <- base_nma$TE.random[, reference]
  weights <- 1 / (base_nma$seTE.random[, reference]^2)

  # Build design matrix with main effects
  X <- comp_matrix[treatments, , drop = FALSE]

  # Add interaction terms
  if (!is.null(interactions) && length(interactions) > 0) {

    interaction_cols <- list()
    interaction_names <- character()

    for (int_term in interactions) {
      # Parse interaction (e.g., "A:B")
      components_in_int <- strsplit(int_term, ":")[[1]]
      components_in_int <- trimws(components_in_int)

      if (length(components_in_int) != 2) {
        stop(sprintf("Invalid interaction term: '%s'. Use format 'ComponentA:ComponentB'",
                    int_term))
      }

      comp1 <- components_in_int[1]
      comp2 <- components_in_int[2]

      if (!(comp1 %in% colnames(comp_matrix)) || !(comp2 %in% colnames(comp_matrix))) {
        stop(sprintf("Components in interaction '%s' not found in component matrix", int_term))
      }

      # Interaction column = product of component indicators
      int_col <- comp_matrix[treatments, comp1] * comp_matrix[treatments, comp2]
      interaction_cols[[length(interaction_cols) + 1]] <- int_col
      interaction_names <- c(interaction_names, sprintf("%s:%s", comp1, comp2))
    }

    # Combine main effects and interactions
    X_interaction <- do.call(cbind, interaction_cols)
    colnames(X_interaction) <- interaction_names
    X <- cbind(X, X_interaction)
  }

  # Weighted least squares
  W <- diag(weights)
  XtWX <- t(X) %*% W %*% X
  XtWy <- t(X) %*% W %*% y

  if (det(XtWX) < 1e-10) {
    stop("Component and interaction effects are not identifiable. ",
         "Reduce model complexity or add more data.")
  }

  effects <- solve(XtWX) %*% XtWy
  effects_se <- sqrt(diag(solve(XtWX)))

  # Separate main effects and interaction effects
  n_main <- ncol(comp_matrix)
  main_effects <- effects[1:n_main, , drop = FALSE]
  main_se <- effects_se[1:n_main]

  main_results <- data.frame(
    Component = colnames(comp_matrix),
    Effect = as.vector(main_effects),
    SE = main_se,
    Lower_95 = main_effects - 1.96 * main_se,
    Upper_95 = main_effects + 1.96 * main_se,
    P_value = 2 * (1 - pnorm(abs(main_effects / main_se))),
    stringsAsFactors = FALSE
  )

  # Interaction effects
  if (length(effects) > n_main) {
    int_effects <- effects[(n_main + 1):length(effects), , drop = FALSE]
    int_se <- effects_se[(n_main + 1):length(effects_se)]

    interaction_results <- data.frame(
      Interaction = rownames(X)[(n_main + 1):ncol(X)],
      Effect = as.vector(int_effects),
      SE = int_se,
      Lower_95 = int_effects - 1.96 * int_se,
      Upper_95 = int_effects + 1.96 * int_se,
      P_value = 2 * (1 - pnorm(abs(int_effects / int_se))),
      Significant = 2 * (1 - pnorm(abs(int_effects / int_se))) < 0.05,
      stringsAsFactors = FALSE
    )
  } else {
    interaction_results <- NULL
  }

  # Predicted effects
  predicted_effects <- X %*% effects
  residuals <- y - predicted_effects

  # Model fit
  rss <- sum(weights * residuals^2)
  total_ss <- sum(weights * (y - weighted.mean(y, weights))^2)
  r_squared <- 1 - (rss / total_ss)

  # Consistency check
  consistency_check <- data.frame(
    Treatment = treatments,
    NMA_Effect = y,
    CNMA_Predicted = as.vector(predicted_effects),
    Difference = as.vector(residuals),
    stringsAsFactors = FALSE
  )

  result <- list(
    method = "interaction",
    main_effects = main_results,
    interaction_effects = interaction_results,
    component_matrix = comp_matrix,
    base_nma = base_nma,
    consistency_check = consistency_check,
    model_fit = list(
      R_squared = r_squared,
      RSS = rss,
      n_main = n_main,
      n_interactions = if (is.null(interaction_results)) 0 else nrow(interaction_results)
    ),
    sm = sm,
    reference = reference
  )

  class(result) <- c("cnma_interaction", "cnma", "list")
  result
}

#' Compare Additive and Interaction CNMA Models
#'
#' Fit both models and compare to determine if interaction effects improve fit.
#'
#' @param data Pairwise NMA data
#' @param comp_matrix Component matrix
#' @param interactions Interaction terms to test
#' @param sm Summary measure
#' @param reference Reference treatment
#' @return cnma_comparison object
#' @export
#' @examples
#' \dontrun{
#' comparison <- compare_cnma_models(
#'   data = nma_data,
#'   comp_matrix = comp_matrix,
#'   interactions = c("A:B"),
#'   reference = "Placebo"
#' )
#' print(comparison)
#' }
compare_cnma_models <- function(data, comp_matrix, interactions,
                               sm = "MD", reference = NULL) {

  message("Fitting additive model...")
  additive_model <- additive_cnma(data, comp_matrix, sm, reference)

  message("Fitting interaction model...")
  interaction_model <- interaction_cnma(data, comp_matrix, interactions, sm, reference)

  # Compare model fits
  aic_additive <- 2 * additive_model$model_fit$n_components +
                 additive_model$model_fit$RSS

  aic_interaction <- 2 * (interaction_model$model_fit$n_main +
                         interaction_model$model_fit$n_interactions) +
                    interaction_model$model_fit$RSS

  # Likelihood ratio test (approximate)
  df_diff <- interaction_model$model_fit$n_interactions
  rss_diff <- additive_model$model_fit$RSS - interaction_model$model_fit$RSS

  # Calculate F-statistic
  n <- nrow(interaction_model$consistency_check)
  f_stat <- (rss_diff / df_diff) / (interaction_model$model_fit$RSS / (n - interaction_model$model_fit$n_main - interaction_model$model_fit$n_interactions))
  p_value <- 1 - pf(f_stat, df_diff, n - interaction_model$model_fit$n_main - interaction_model$model_fit$n_interactions)

  comparison <- data.frame(
    Model = c("Additive", "Interaction"),
    R_squared = c(additive_model$model_fit$R_squared,
                 interaction_model$model_fit$R_squared),
    RSS = c(additive_model$model_fit$RSS,
           interaction_model$model_fit$RSS),
    N_parameters = c(additive_model$model_fit$n_components,
                    interaction_model$model_fit$n_main + interaction_model$model_fit$n_interactions),
    AIC = c(aic_additive, aic_interaction),
    stringsAsFactors = FALSE
  )

  result <- list(
    additive_model = additive_model,
    interaction_model = interaction_model,
    comparison_table = comparison,
    f_test = list(
      F_statistic = f_stat,
      df1 = df_diff,
      df2 = n - interaction_model$model_fit$n_main - interaction_model$model_fit$n_interactions,
      p_value = p_value
    ),
    preferred_model = if (p_value < 0.05) "interaction" else "additive"
  )

  class(result) <- c("cnma_comparison", "list")
  result
}

#' Print Additive CNMA Results
#'
#' @param x cnma_additive object
#' @param ... Additional arguments
#' @export
print.cnma_additive <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Additive Component Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Summary measure: %s\n", x$sm))
  cat(sprintf("Reference treatment: %s\n", x$reference))
  cat(sprintf("Number of components: %d\n", x$model_fit$n_components))
  cat(sprintf("Number of treatments: %d\n", x$model_fit$n_treatments))
  cat(sprintf("Model R²: %.3f\n\n", x$model_fit$R_squared))

  cat("Component Effects:\n\n")
  print(x$component_effects, row.names = FALSE, digits = 3)

  cat("\n")
  cat("Model Assumption: Treatment effect = Sum of component effects\n")
  cat("\n")
  cat("Consistency Check (NMA vs CNMA predictions):\n\n")
  print(x$consistency_check, row.names = FALSE, digits = 3)

  cat("\n")
  max_diff <- max(abs(x$consistency_check$Difference))
  if (max_diff > 0.5) {
    cat("⚠ WARNING: Large discrepancies detected. Additive model may be inadequate.\n")
    cat("  Consider fitting interaction model with interaction_cnma().\n")
  } else {
    cat("✓ Additive model fits reasonably well.\n")
  }

  cat("\n")

  invisible(x)
}

#' Print Interaction CNMA Results
#'
#' @param x cnma_interaction object
#' @param ... Additional arguments
#' @export
print.cnma_interaction <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Interaction Component Network Meta-Analysis\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Summary measure: %s\n", x$sm))
  cat(sprintf("Reference treatment: %s\n", x$reference))
  cat(sprintf("Number of main effects: %d\n", x$model_fit$n_main))
  cat(sprintf("Number of interaction terms: %d\n", x$model_fit$n_interactions))
  cat(sprintf("Model R²: %.3f\n\n", x$model_fit$R_squared))

  cat("Main Component Effects:\n\n")
  print(x$main_effects, row.names = FALSE, digits = 3)

  if (!is.null(x$interaction_effects)) {
    cat("\n")
    cat("Interaction Effects:\n\n")
    print(x$interaction_effects, row.names = FALSE, digits = 3)

    sig_interactions <- sum(x$interaction_effects$Significant)
    if (sig_interactions > 0) {
      cat(sprintf("\n✓ %d significant interaction(s) detected.\n", sig_interactions))
      cat("  Components show synergistic or antagonistic effects.\n")
    } else {
      cat("\nNo significant interactions detected.\n")
      cat("  Additive model may be sufficient.\n")
    }
  }

  cat("\n")
  cat("Consistency Check:\n\n")
  print(x$consistency_check, row.names = FALSE, digits = 3)

  cat("\n")

  invisible(x)
}

#' Plot Component Effects
#'
#' @param x cnma object
#' @param type Plot type: "forest" or "consistency"
#' @param ... Additional arguments
#' @export
plot.cnma_additive <- function(x, type = c("forest", "consistency"), ...) {

  type <- match.arg(type)

  if (type == "forest") {
    # Forest plot of component effects
    plot_data <- x$component_effects
    plot_data$Component <- factor(plot_data$Component,
                                  levels = plot_data$Component[order(plot_data$Effect)])

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Component, y = Effect)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
                            width = 0.3, size = 1) +
      ggplot2::geom_point(size = 4, color = "#2c7fb8") +
      ggplot2::labs(
        title = "Component Effects (Additive Model)",
        subtitle = sprintf("Reference: %s", x$reference),
        x = "Component",
        y = sprintf("Effect (%s)", x$sm)
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(size = 11)
      )

  } else {
    # Consistency plot: NMA vs CNMA predictions
    plot_data <- x$consistency_check

    max_val <- max(abs(c(plot_data$NMA_Effect, plot_data$CNMA_Predicted)))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = NMA_Effect, y = CNMA_Predicted)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "red", size = 1) +
      ggplot2::geom_point(size = 4, color = "#2c7fb8", alpha = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = Treatment), vjust = -1, size = 3) +
      ggplot2::labs(
        title = "Consistency: Standard NMA vs Component NMA",
        subtitle = sprintf("R² = %.3f", x$model_fit$R_squared),
        x = "Standard NMA Effect",
        y = "Component NMA Predicted Effect"
      ) +
      ggplot2::xlim(-max_val, max_val) +
      ggplot2::ylim(-max_val, max_val) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")
      )
  }

  print(p)
  invisible(p)
}

#' @export
plot.cnma_interaction <- function(x, type = c("forest", "consistency", "interactions"), ...) {

  type <- match.arg(type)

  if (type == "forest") {
    # Forest plot of main effects
    plot_data <- x$main_effects
    plot_data$Component <- factor(plot_data$Component,
                                  levels = plot_data$Component[order(plot_data$Effect)])

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Component, y = Effect)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
                            width = 0.3, size = 1) +
      ggplot2::geom_point(size = 4, color = "#2c7fb8") +
      ggplot2::labs(
        title = "Main Component Effects (Interaction Model)",
        x = "Component",
        y = sprintf("Effect (%s)", x$sm)
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal()

  } else if (type == "interactions") {
    # Plot interaction effects
    if (is.null(x$interaction_effects)) {
      stop("No interaction effects in model")
    }

    plot_data <- x$interaction_effects
    plot_data$Interaction <- factor(plot_data$Interaction,
                                    levels = plot_data$Interaction[order(plot_data$Effect)])
    plot_data$Color <- ifelse(plot_data$Significant, "Significant", "Not significant")

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Interaction, y = Effect, color = Color)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
                            width = 0.3, size = 1) +
      ggplot2::geom_point(size = 4) +
      ggplot2::scale_color_manual(values = c("Significant" = "#d73027",
                                            "Not significant" = "#2c7fb8")) +
      ggplot2::labs(
        title = "Interaction Effects Between Components",
        subtitle = "Positive = synergistic, Negative = antagonistic",
        x = "Interaction",
        y = sprintf("Interaction Effect (%s)", x$sm)
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

  } else {
    # Consistency plot
    plot_data <- x$consistency_check
    max_val <- max(abs(c(plot_data$NMA_Effect, plot_data$CNMA_Predicted)))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = NMA_Effect, y = CNMA_Predicted)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "red", size = 1) +
      ggplot2::geom_point(size = 4, color = "#2c7fb8", alpha = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = Treatment), vjust = -1, size = 3) +
      ggplot2::labs(
        title = "Consistency: Standard NMA vs Component NMA",
        subtitle = sprintf("R² = %.3f", x$model_fit$R_squared),
        x = "Standard NMA Effect",
        y = "Component NMA Predicted Effect"
      ) +
      ggplot2::xlim(-max_val, max_val) +
      ggplot2::ylim(-max_val, max_val) +
      ggplot2::theme_minimal()
  }

  print(p)
  invisible(p)
}

#' Print CNMA Comparison Results
#'
#' @param x cnma_comparison object
#' @param ... Additional arguments
#' @export
print.cnma_comparison <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Component NMA Model Comparison\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat("Model Fit Comparison:\n\n")
  print(x$comparison_table, row.names = FALSE, digits = 3)

  cat("\n")
  cat("F-Test for Interaction Terms:\n")
  cat(sprintf("  F(%d, %d) = %.3f, p = %.4f\n",
             x$f_test$df1, x$f_test$df2,
             x$f_test$F_statistic, x$f_test$p_value))

  cat("\n")
  if (x$f_test$p_value < 0.05) {
    cat("✓ Interaction model significantly improves fit (p < 0.05)\n")
    cat("  Components show synergistic/antagonistic effects.\n")
    cat("  Recommended: Use interaction model\n")
  } else {
    cat("○ No significant improvement with interaction model (p ≥ 0.05)\n")
    cat("  Additive effects may be sufficient.\n")
    cat("  Recommended: Use simpler additive model\n")
  }

  cat(sprintf("\nPreferred model: %s\n", toupper(x$preferred_model)))
  cat("\n")

  invisible(x)
}
