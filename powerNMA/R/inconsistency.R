#' Inconsistency Assessment for Network Meta-Analysis
#'
#' Comprehensive tools for detecting and quantifying inconsistency in NMA.
#' Based on Dias et al. (2010) and Valkenhoef et al. (2016).
#'
#' @name inconsistency
NULL

#' Node-Splitting Analysis
#'
#' Performs node-splitting to assess inconsistency by comparing direct and
#' indirect evidence for each comparison. Based on Dias et al. (2010) and
#' automated approach by Valkenhoef et al. (2016).
#'
#' @param nma_result A netmeta object from network meta-analysis
#' @param data Original pairwise data frame
#' @param comparisons Optional vector of comparisons to split (format: "A:B").
#'   If NULL, splits all comparisons with both direct and indirect evidence.
#' @param p_threshold P-value threshold for flagging inconsistency (default: 0.10)
#' @return node_split object with direct/indirect comparisons and inconsistency tests
#' @export
#' @references
#' Dias S, et al. (2010). Checking consistency in mixed treatment comparison
#' meta-analysis. Statistics in Medicine, 29(7-8):932-944.
#'
#' Valkenhoef G, et al. (2016). Automated generation of node-splitting models
#' for assessment of inconsistency in network meta-analysis.
#' Research Synthesis Methods, 7(1):80-93.
#'
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab,
#'                         data = data, reference.group = "A")
#' node_split_results <- node_splitting(nma, data)
#' print(node_split_results)
#' }
node_splitting <- function(nma_result, data, comparisons = NULL, p_threshold = 0.10) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Identify comparisons with both direct and indirect evidence
  direct_comps <- identify_direct_comparisons(data)

  if (is.null(comparisons)) {
    # Auto-select all comparisons in potentially inconsistent loops
    comparisons <- direct_comps
    message(sprintf("Auto-selected %d comparisons with direct evidence for node-splitting",
                   length(comparisons)))
  }

  results <- list()

  for (comp in comparisons) {
    parts <- strsplit(comp, ":")[[1]]
    if (length(parts) != 2) next

    treat1 <- parts[1]
    treat2 <- parts[2]

    # Run node-split for this comparison
    split_result <- tryCatch({
      perform_single_node_split(nma_result, data, treat1, treat2)
    }, error = function(e) {
      warning(sprintf("Node-splitting failed for %s: %s", comp, e$message))
      NULL
    })

    if (!is.null(split_result)) {
      results[[comp]] <- split_result
    }
  }

  # Summarize results
  summary_table <- do.call(rbind, lapply(names(results), function(comp) {
    r <- results[[comp]]
    data.frame(
      comparison = comp,
      direct_effect = r$direct_TE,
      direct_se = r$direct_seTE,
      indirect_effect = r$indirect_TE,
      indirect_se = r$indirect_seTE,
      difference = r$difference,
      difference_se = r$difference_se,
      p_value = r$p_value,
      inconsistent = r$p_value < p_threshold,
      stringsAsFactors = FALSE
    )
  }))

  output <- list(
    results = results,
    summary = summary_table,
    n_tests = nrow(summary_table),
    n_inconsistent = sum(summary_table$p_value < p_threshold),
    p_threshold = p_threshold
  )

  class(output) <- c("node_split", "list")
  output
}

#' Perform Single Node-Split
#'
#' Internal function to split a single comparison
#'
#' @param nma_result netmeta object
#' @param data Pairwise data
#' @param treat1 First treatment
#' @param treat2 Second treatment
#' @return List with direct/indirect estimates and inconsistency test
#' @keywords internal
perform_single_node_split <- function(nma_result, data, treat1, treat2) {

  # Separate direct evidence for this comparison
  direct_data <- data[
    (data$treat1 == treat1 & data$treat2 == treat2) |
    (data$treat1 == treat2 & data$treat2 == treat1),
  ]

  # Remove direct evidence to get indirect
  indirect_data <- data[
    !((data$treat1 == treat1 & data$treat2 == treat2) |
      (data$treat1 == treat2 & data$treat2 == treat1)),
  ]

  if (nrow(direct_data) == 0) {
    stop("No direct evidence for this comparison")
  }

  if (nrow(indirect_data) < 2) {
    stop("Insufficient indirect evidence")
  }

  # Meta-analysis of direct evidence only
  direct_ma <- tryCatch({
    if (nrow(direct_data) == 1) {
      # Single study - use point estimate
      list(
        TE = direct_data$TE[1],
        seTE = direct_data$seTE[1]
      )
    } else {
      # Multiple studies - meta-analyze
      meta::metagen(
        TE = direct_data$TE,
        seTE = direct_data$seTE,
        studlab = direct_data$studlab,
        sm = nma_result$sm
      )
    }
  }, error = function(e) {
    stop(sprintf("Direct evidence meta-analysis failed: %s", e$message))
  })

  # Network meta-analysis excluding direct evidence
  indirect_nma <- tryCatch({
    netmeta::netmeta(
      TE = TE, seTE = seTE,
      treat1 = treat1, treat2 = treat2,
      studlab = studlab,
      data = indirect_data,
      sm = nma_result$sm,
      reference.group = nma_result$reference.group
    )
  }, error = function(e) {
    stop(sprintf("Indirect NMA failed: %s", e$message))
  })

  # Extract direct estimate
  if (is.list(direct_ma) && !inherits(direct_ma, "meta")) {
    direct_TE <- direct_ma$TE
    direct_seTE <- direct_ma$seTE
  } else {
    direct_TE <- direct_ma$TE.random
    direct_seTE <- direct_ma$seTE.random
  }

  # Extract indirect estimate for this comparison
  trts <- rownames(indirect_nma$TE.random)
  if (!(treat1 %in% trts && treat2 %in% trts)) {
    stop("Treatments not in indirect network")
  }

  indirect_TE <- indirect_nma$TE.random[treat2, treat1]
  indirect_seTE <- indirect_nma$seTE.random[treat2, treat1]

  # Test for inconsistency: difference between direct and indirect
  difference <- direct_TE - indirect_TE
  difference_se <- sqrt(direct_seTE^2 + indirect_seTE^2)
  z_score <- difference / difference_se
  p_value <- 2 * (1 - pnorm(abs(z_score)))

  list(
    treat1 = treat1,
    treat2 = treat2,
    direct_TE = direct_TE,
    direct_seTE = direct_seTE,
    indirect_TE = indirect_TE,
    indirect_seTE = indirect_seTE,
    difference = difference,
    difference_se = difference_se,
    z_score = z_score,
    p_value = p_value
  )
}

#' Design-by-Treatment Inconsistency Model
#'
#' Assesses inconsistency using the design-by-treatment interaction approach.
#' Decomposes total heterogeneity into within-design (heterogeneity) and
#' between-design (inconsistency) components.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data with studlab, treat1, treat2
#' @return List with heterogeneity and inconsistency statistics
#' @export
#' @references
#' Higgins JPT, et al. (2012). Consistency and inconsistency in network
#' meta-analysis: concepts and models for multi-arm studies.
#' Research Synthesis Methods, 3(2):98-110.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' design_inc <- design_inconsistency(nma, data)
#' print(design_inc)
#' }
design_inconsistency <- function(nma_result, data) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Use netmeta's decomposition if available
  if (!is.null(nma_result$Q.decomp)) {
    decomp <- nma_result$Q.decomp

    result <- list(
      Q_total = decomp$Q,
      Q_heterogeneity = decomp$Q.heterogeneity,
      Q_inconsistency = decomp$Q.inconsistency,
      df_total = decomp$df.Q,
      df_heterogeneity = decomp$df.Q.heterogeneity,
      df_inconsistency = decomp$df.Q.inconsistency,
      p_heterogeneity = pchisq(decomp$Q.heterogeneity,
                               decomp$df.Q.heterogeneity,
                               lower.tail = FALSE),
      p_inconsistency = pchisq(decomp$Q.inconsistency,
                               decomp$df.Q.inconsistency,
                               lower.tail = FALSE),
      interpretation = interpret_design_inconsistency(
        pchisq(decomp$Q.inconsistency, decomp$df.Q.inconsistency, lower.tail = FALSE)
      )
    )
  } else {
    # Manual calculation if decomposition not available
    result <- calculate_design_inconsistency_manual(nma_result, data)
  }

  class(result) <- c("design_inconsistency", "list")
  result
}

#' Calculate Design Inconsistency Manually
#'
#' @param nma_result netmeta object
#' @param data Pairwise data
#' @return List with inconsistency statistics
#' @keywords internal
calculate_design_inconsistency_manual <- function(nma_result, data) {

  # Get Q statistic from NMA
  Q_total <- nma_result$Q
  df_total <- nma_result$df.Q

  # Calculate within-design heterogeneity by running separate meta-analyses
  # for each direct comparison
  direct_comparisons <- unique(data[, c("treat1", "treat2")])

  Q_het <- 0
  df_het <- 0

  for (i in seq_len(nrow(direct_comparisons))) {
    t1 <- direct_comparisons$treat1[i]
    t2 <- direct_comparisons$treat2[i]

    comp_data <- data[
      (data$treat1 == t1 & data$treat2 == t2) |
      (data$treat1 == t2 & data$treat2 == t1),
    ]

    if (nrow(comp_data) > 1) {
      ma <- meta::metagen(
        TE = comp_data$TE,
        seTE = comp_data$seTE,
        studlab = comp_data$studlab
      )
      Q_het <- Q_het + ma$Q
      df_het <- df_het + ma$df.Q
    }
  }

  # Inconsistency = Total - Heterogeneity
  Q_inc <- Q_total - Q_het
  df_inc <- df_total - df_het

  p_het <- pchisq(Q_het, df_het, lower.tail = FALSE)
  p_inc <- if (df_inc > 0) pchisq(Q_inc, df_inc, lower.tail = FALSE) else NA

  list(
    Q_total = Q_total,
    Q_heterogeneity = Q_het,
    Q_inconsistency = Q_inc,
    df_total = df_total,
    df_heterogeneity = df_het,
    df_inconsistency = df_inc,
    p_heterogeneity = p_het,
    p_inconsistency = p_inc,
    interpretation = interpret_design_inconsistency(p_inc)
  )
}

#' Loop Inconsistency (Bucher's Method)
#'
#' Identifies and tests for inconsistency in closed loops of the network.
#' Uses Bucher's indirect comparison method.
#'
#' @param nma_result A netmeta object
#' @param data Original pairwise data
#' @return Data frame with loop inconsistency results
#' @export
#' @references
#' Bucher HC, et al. (1997). The results of direct and indirect treatment
#' comparisons in meta-analysis of randomized controlled trials.
#' Journal of Clinical Epidemiology, 50(6):683-691.
#'
#' @examples
#' \dontrun{
#' nma <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, data = data)
#' loops <- loop_inconsistency(nma, data)
#' print(loops)
#' }
loop_inconsistency <- function(nma_result, data) {

  if (!inherits(nma_result, "netmeta")) {
    stop("nma_result must be a netmeta object")
  }

  # Identify all closed triangular loops in the network
  loops <- identify_triangular_loops(data)

  if (length(loops) == 0) {
    message("No triangular loops found in the network")
    return(NULL)
  }

  results <- lapply(loops, function(loop) {
    test_loop_inconsistency(nma_result, data, loop)
  })

  # Combine into data frame
  loop_table <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r)) return(NULL)
    data.frame(
      loop = paste(r$treatments, collapse = "-"),
      IF = r$IF,  # Inconsistency Factor
      se_IF = r$se_IF,
      lower_CI = r$lower_CI,
      upper_CI = r$upper_CI,
      z_score = r$z_score,
      p_value = r$p_value,
      stringsAsFactors = FALSE
    )
  }))

  loop_table$inconsistent <- loop_table$p_value < 0.10

  class(loop_table) <- c("loop_inconsistency", "data.frame")
  loop_table
}

#' Identify Direct Comparisons
#'
#' @param data Pairwise data
#' @return Vector of comparison names
#' @keywords internal
identify_direct_comparisons <- function(data) {
  comps <- unique(paste(
    pmin(as.character(data$treat1), as.character(data$treat2)),
    pmax(as.character(data$treat1), as.character(data$treat2)),
    sep = ":"
  ))
  comps
}

#' Identify Triangular Loops
#'
#' @param data Pairwise data
#' @return List of triangular loops (each a 3-element vector)
#' @keywords internal
identify_triangular_loops <- function(data) {

  # Build adjacency matrix
  trts <- sort(unique(c(as.character(data$treat1), as.character(data$treat2))))
  n_trts <- length(trts)

  adj_matrix <- matrix(FALSE, n_trts, n_trts)
  rownames(adj_matrix) <- colnames(adj_matrix) <- trts

  for (i in seq_len(nrow(data))) {
    t1 <- as.character(data$treat1[i])
    t2 <- as.character(data$treat2[i])
    adj_matrix[t1, t2] <- TRUE
    adj_matrix[t2, t1] <- TRUE
  }

  # Find all triangles
  loops <- list()

  for (i in 1:(n_trts - 2)) {
    for (j in (i + 1):(n_trts - 1)) {
      for (k in (j + 1):n_trts) {
        if (adj_matrix[i, j] && adj_matrix[j, k] && adj_matrix[k, i]) {
          loops[[length(loops) + 1]] <- trts[c(i, j, k)]
        }
      }
    }
  }

  loops
}

#' Test Loop Inconsistency
#'
#' @param nma_result netmeta object
#' @param data Pairwise data
#' @param loop Vector of 3 treatment names
#' @return List with inconsistency factor and statistics
#' @keywords internal
test_loop_inconsistency <- function(nma_result, data, loop) {

  if (length(loop) != 3) return(NULL)

  tryCatch({
    # Get direct estimates for the three sides of the triangle
    AB_direct <- get_direct_estimate(data, loop[1], loop[2])
    BC_direct <- get_direct_estimate(data, loop[2], loop[3])
    AC_direct <- get_direct_estimate(data, loop[1], loop[3])

    # Calculate inconsistency factor: IF = AB + BC - AC
    IF <- AB_direct$TE + BC_direct$TE - AC_direct$TE
    se_IF <- sqrt(AB_direct$seTE^2 + BC_direct$seTE^2 + AC_direct$seTE^2)

    # Test statistic
    z_score <- IF / se_IF
    p_value <- 2 * (1 - pnorm(abs(z_score)))

    # 95% CI
    lower_CI <- IF - 1.96 * se_IF
    upper_CI <- IF + 1.96 * se_IF

    list(
      treatments = loop,
      IF = IF,
      se_IF = se_IF,
      lower_CI = lower_CI,
      upper_CI = upper_CI,
      z_score = z_score,
      p_value = p_value
    )
  }, error = function(e) {
    warning(sprintf("Loop inconsistency test failed for %s: %s",
                   paste(loop, collapse = "-"), e$message))
    NULL
  })
}

#' Get Direct Estimate for Comparison
#'
#' @param data Pairwise data
#' @param treat1 First treatment
#' @param treat2 Second treatment
#' @return List with TE and seTE
#' @keywords internal
get_direct_estimate <- function(data, treat1, treat2) {

  comp_data <- data[
    (data$treat1 == treat1 & data$treat2 == treat2) |
    (data$treat1 == treat2 & data$treat2 == treat1),
  ]

  if (nrow(comp_data) == 0) {
    stop(sprintf("No direct evidence for %s vs %s", treat1, treat2))
  }

  # Ensure correct direction
  comp_data$TE_adj <- ifelse(
    comp_data$treat1 == treat1,
    comp_data$TE,
    -comp_data$TE
  )

  if (nrow(comp_data) == 1) {
    list(TE = comp_data$TE_adj[1], seTE = comp_data$seTE[1])
  } else {
    ma <- meta::metagen(
      TE = comp_data$TE_adj,
      seTE = comp_data$seTE,
      studlab = comp_data$studlab
    )
    list(TE = ma$TE.random, seTE = ma$seTE.random)
  }
}

#' Interpret Design Inconsistency
#'
#' @param p_value P-value from inconsistency test
#' @return Character interpretation
#' @keywords internal
interpret_design_inconsistency <- function(p_value) {
  if (is.na(p_value)) {
    "Cannot assess (insufficient data)"
  } else if (p_value >= 0.10) {
    "No evidence of inconsistency (p >= 0.10)"
  } else if (p_value >= 0.05) {
    "Weak evidence of inconsistency (0.05 <= p < 0.10)"
  } else if (p_value >= 0.01) {
    "Moderate evidence of inconsistency (0.01 <= p < 0.05)"
  } else {
    "Strong evidence of inconsistency (p < 0.01)"
  }
}

#' Print Node-Splitting Results
#'
#' @param x node_split object
#' @param ... Additional arguments
#' @export
print.node_split <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Node-Splitting Analysis for Inconsistency\n")
  cat("══════════════════════════════════════════════════\n\n")

  cat(sprintf("Comparisons tested: %d\n", x$n_tests))
  cat(sprintf("Inconsistent (p < %.2f): %d (%.1f%%)\n\n",
             x$p_threshold,
             x$n_inconsistent,
             100 * x$n_inconsistent / x$n_tests))

  if (nrow(x$summary) > 0) {
    cat("Results:\n")
    print(x$summary, row.names = FALSE)
  }

  if (x$n_inconsistent > 0) {
    cat("\n⚠ WARNING: Inconsistency detected in", x$n_inconsistent, "comparison(s)\n")
    cat("Consider investigating sources of inconsistency:\n")
    cat("  • Clinical/methodological heterogeneity\n")
    cat("  • Violations of transitivity assumption\n")
    cat("  • Effect modification\n")
  } else {
    cat("\n✓ No significant inconsistency detected\n")
  }

  invisible(x)
}

#' Print Design Inconsistency Results
#'
#' @param x design_inconsistency object
#' @param ... Additional arguments
#' @export
print.design_inconsistency <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Design-by-Treatment Inconsistency Model\n")
  cat("══════════════════════════════════════════════════\n\n")

  cat(sprintf("Total Q statistic: %.2f (df = %d, p = %.4f)\n",
             x$Q_total, x$df_total,
             pchisq(x$Q_total, x$df_total, lower.tail = FALSE)))

  cat("\nDecomposition:\n")
  cat(sprintf("  Within-design (heterogeneity): Q = %.2f (df = %d, p = %.4f)\n",
             x$Q_heterogeneity, x$df_heterogeneity, x$p_heterogeneity))
  cat(sprintf("  Between-design (inconsistency): Q = %.2f (df = %d, p = %.4f)\n",
             x$Q_inconsistency, x$df_inconsistency, x$p_inconsistency))

  cat("\nInterpretation:", x$interpretation, "\n")

  invisible(x)
}

#' Print Loop Inconsistency Results
#'
#' @param x loop_inconsistency object
#' @param ... Additional arguments
#' @export
print.loop_inconsistency <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════\n")
  cat("  Loop Inconsistency Analysis (Bucher's Method)\n")
  cat("══════════════════════════════════════════════════\n\n")

  if (is.null(x) || nrow(x) == 0) {
    cat("No closed loops found in the network\n")
    return(invisible(x))
  }

  cat(sprintf("Closed loops tested: %d\n", nrow(x)))
  cat(sprintf("Inconsistent loops (p < 0.10): %d\n\n", sum(x$inconsistent)))

  print(as.data.frame(x), row.names = FALSE)

  if (any(x$inconsistent)) {
    cat("\n⚠ WARNING: Inconsistency detected in some loops\n")
  } else {
    cat("\n✓ No significant loop inconsistency detected\n")
  }

  invisible(x)
}
