# ============================================================================
# PERFORMANCE BENCHMARKING: Large Network Meta-Analyses
# ============================================================================
#
# Purpose: Benchmark powerNMA automatic pathways on large networks
# Focus: Computational efficiency for networks with >20 treatments
# Methods tested: auto_standard_nma, auto_experimental_nma
#
# Benchmarks:
# - Small network: 5 studies, 4 treatments
# - Medium network: 20 studies, 10 treatments
# - Large network: 50 studies, 20 treatments
# - Very large network: 100 studies, 30 treatments
#
# Metrics:
# - Total runtime
# - Memory usage
# - Convergence success rate
# - Scalability (runtime vs network size)
#
# ============================================================================

library(powerNMA)
library(netmeta)
library(dplyr)

cat("=============================================================\n")
cat("PERFORMANCE BENCHMARKING: Large Network Meta-Analyses\n")
cat("=============================================================\n\n")

# ============================================================================
# Benchmark Configuration
# ============================================================================

BENCHMARK_CONFIGS <- list(
  small = list(
    n_studies = 5,
    n_treatments = 4,
    label = "Small (5 studies, 4 treatments)"
  ),
  medium = list(
    n_studies = 20,
    n_treatments = 10,
    label = "Medium (20 studies, 10 treatments)"
  ),
  large = list(
    n_studies = 50,
    n_treatments = 20,
    label = "Large (50 studies, 20 treatments)"
  ),
  very_large = list(
    n_studies = 100,
    n_treatments = 30,
    label = "Very Large (100 studies, 30 treatments)"
  )
)

# Performance targets (seconds)
PERFORMANCE_TARGETS <- list(
  small = 5,
  medium = 30,
  large = 120,
  very_large = 300
)

# ============================================================================
# Helper Functions
# ============================================================================

#' Generate synthetic network for benchmarking
#' @param n_studies Number of studies
#' @param n_treatments Number of treatments
#' @param tau Between-study heterogeneity
simulate_network <- function(n_studies, n_treatments, tau = 0.1) {

  set.seed(20251101)

  # True treatment effects
  true_effects <- c(0, rnorm(n_treatments - 1, mean = 0.3, sd = 0.2))

  # Study effects
  study_effects <- rnorm(n_studies, mean = 0, sd = tau)

  # Create star network (all treatments compared to treatment 1)
  comparisons <- expand.grid(
    study = 1:n_studies,
    treat2 = 2:n_treatments
  )
  comparisons$treat1 <- 1

  # Add some loops for connectivity
  if (n_treatments >= 4 && n_studies >= 10) {
    # Add comparisons between treatments 2 and 3 in some studies
    loop_studies <- sample(1:n_studies, min(5, floor(n_studies/2)))
    for (s in loop_studies) {
      comparisons <- rbind(comparisons, data.frame(
        study = s,
        treat1 = 2,
        treat2 = 3
      ))
    }
  }

  # Calculate effects
  comparisons$TE <- apply(comparisons, 1, function(row) {
    s <- as.numeric(row["study"])
    t1 <- as.numeric(row["treat1"])
    t2 <- as.numeric(row["treat2"])

    true_diff <- true_effects[t2] - true_effects[t1]
    true_diff + study_effects[s] + rnorm(1, 0, 0.05)
  })

  # Sample sizes
  comparisons$n1 <- sample(30:100, nrow(comparisons), replace = TRUE)
  comparisons$n2 <- sample(30:100, nrow(comparisons), replace = TRUE)

  # Standard errors
  comparisons$seTE <- sqrt(1/comparisons$n1 + 1/comparisons$n2)

  # Format
  comparisons <- comparisons %>%
    mutate(
      studlab = paste0("Study", study),
      treat1 = paste0("T", treat1),
      treat2 = paste0("T", treat2)
    ) %>%
    select(studlab, treat1, treat2, TE, seTE)

  return(comparisons)
}


#' Benchmark a single configuration
#' @param config Benchmark configuration
benchmark_configuration <- function(config_name, config) {

  cat("\\n--- Benchmarking:", config$label, "---\\n")

  # Generate network
  cat("  Generating network...")
  data <- simulate_network(config$n_studies, config$n_treatments)
  cat(" Done (", nrow(data), " comparisons)\\n", sep = "")

  # Benchmark auto_standard_nma
  cat("  Running auto_standard_nma...")

  start_time <- Sys.time()
  start_mem <- gc()[2, 2]  # Used memory before

  result_std <- tryCatch({
    auto_standard_nma(
      data = data,
      verbose = FALSE
    )
  }, error = function(e) {
    list(status = "failed", error = e$message)
  })

  end_time <- Sys.time()
  end_mem <- gc()[2, 2]  # Used memory after

  runtime_std <- as.numeric(difftime(end_time, start_time, units = "secs"))
  mem_used_std <- end_mem - start_mem

  cat(" Done (", round(runtime_std, 2), " sec)\\n", sep = "")

  # Check success
  success_std <- !is.null(result_std$primary_analysis) &&
                 result_std$primary_analysis$status == "completed"

  # Benchmark auto_experimental_nma (optional, slower)
  cat("  Running auto_experimental_nma...")

  start_time_exp <- Sys.time()

  result_exp <- tryCatch({
    auto_experimental_nma(
      data = data,
      research_question = "decision_making",
      verbose = FALSE
    )
  }, error = function(e) {
    list(status = "failed", error = e$message)
  })

  end_time_exp <- Sys.time()

  runtime_exp <- as.numeric(difftime(end_time_exp, start_time_exp, units = "secs"))

  cat(" Done (", round(runtime_exp, 2), " sec)\\n", sep = "")

  success_exp <- !is.null(result_exp$experimental_analyses)

  # Return benchmarking results
  return(list(
    config = config_name,
    n_studies = config$n_studies,
    n_treatments = config$n_treatments,
    n_comparisons = nrow(data),
    runtime_std = runtime_std,
    runtime_exp = runtime_exp,
    memory_mb = mem_used_std,
    success_std = success_std,
    success_exp = success_exp
  ))
}


# ============================================================================
# Run All Benchmarks
# ============================================================================

cat("Starting performance benchmarking...\\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

benchmark_results <- list()

for (config_name in names(BENCHMARK_CONFIGS)) {
  config <- BENCHMARK_CONFIGS[[config_name]]

  result <- benchmark_configuration(config_name, config)
  benchmark_results[[config_name]] <- result

  # Check against target
  target <- PERFORMANCE_TARGETS[[config_name]]
  within_target <- result$runtime_std <= target

  cat("  Target runtime: <=", target, "sec | Actual:", round(result$runtime_std, 2), "sec |",
      ifelse(within_target, "PASS", "FAIL"), "\\n")
}


# ============================================================================
# Analyze Results
# ============================================================================

cat("\\n\\n=============================================================\n")
cat("PERFORMANCE BENCHMARK RESULTS\n")
cat("=============================================================\n\n")

# Convert to data frame
results_df <- do.call(rbind, lapply(benchmark_results, function(x) {
  data.frame(
    config = x$config,
    n_studies = x$n_studies,
    n_treatments = x$n_treatments,
    n_comparisons = x$n_comparisons,
    runtime_std = x$runtime_std,
    runtime_exp = x$runtime_exp,
    memory_mb = x$memory_mb,
    success_std = x$success_std,
    success_exp = x$success_exp,
    stringsAsFactors = FALSE
  )
}))

# Print results table
cat("Runtime Summary (seconds):\\n")
cat(sprintf("%-15s %10s %10s %15s %15s\\n",
            "Network Size", "Studies", "Treatments", "Auto Std", "Auto Exp"))
cat(strrep("-", 75), "\\n")

for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  cat(sprintf("%-15s %10d %10d %15.2f %15.2f\\n",
              row$config,
              row$n_studies,
              row$n_treatments,
              row$runtime_std,
              row$runtime_exp))
}

cat("\\n")

# Success rates
cat("Success Rates:\\n")
cat("  auto_standard_nma:", sum(results_df$success_std), "/", nrow(results_df),
    sprintf("(%.1f%%)\\n", 100 * mean(results_df$success_std)))
cat("  auto_experimental_nma:", sum(results_df$success_exp), "/", nrow(results_df),
    sprintf("(%.1f%%)\\n\\n", 100 * mean(results_df$success_exp)))

# Performance targets
cat("Performance Target Achievement:\\n")
for (config_name in names(BENCHMARK_CONFIGS)) {
  target <- PERFORMANCE_TARGETS[[config_name]]
  actual <- results_df[results_df$config == config_name, "runtime_std"]
  within <- actual <= target

  cat(sprintf("  %s: Target <= %d sec | Actual %.2f sec | %s\\n",
              BENCHMARK_CONFIGS[[config_name]]$label,
              target,
              actual,
              ifelse(within, "PASS", "FAIL")))
}

cat("\\n")

# Scalability analysis
cat("Scalability Analysis:\\n")

# Fit linear model: runtime ~ n_studies + n_treatments
if (nrow(results_df) >= 3) {
  lm_model <- lm(runtime_std ~ n_studies + n_treatments, data = results_df)

  cat("  Linear model: runtime = b0 + b1*n_studies + b2*n_treatments\\n")
  cat(sprintf("    Intercept: %.3f\\n", coef(lm_model)[1]))
  cat(sprintf("    Studies coefficient: %.3f sec/study\\n", coef(lm_model)[2]))
  cat(sprintf("    Treatments coefficient: %.3f sec/treatment\\n", coef(lm_model)[3]))
  cat(sprintf("    R-squared: %.3f\\n\\n", summary(lm_model)$r.squared))

  # Predict runtime for even larger networks
  cat("  Predicted runtimes for larger networks:\\n")
  predictions <- data.frame(
    n_studies = c(200, 500, 1000),
    n_treatments = c(50, 100, 200)
  )
  predictions$predicted_runtime <- predict(lm_model, newdata = predictions)

  for (i in 1:nrow(predictions)) {
    cat(sprintf("    %d studies, %d treatments: %.1f sec (%.1f min)\\n",
                predictions$n_studies[i],
                predictions$n_treatments[i],
                predictions$predicted_runtime[i],
                predictions$predicted_runtime[i] / 60))
  }
}

cat("\\n")

# Memory usage
cat("Memory Usage:\\n")
cat(sprintf("  Average: %.1f MB\\n", mean(results_df$memory_mb, na.rm = TRUE)))
cat(sprintf("  Max: %.1f MB\\n\\n", max(results_df$memory_mb, na.rm = TRUE)))

# Overall assessment
all_targets_met <- all(sapply(names(PERFORMANCE_TARGETS), function(config_name) {
  actual <- results_df[results_df$config == config_name, "runtime_std"]
  target <- PERFORMANCE_TARGETS[[config_name]]
  actual <= target
}))

all_successful <- all(results_df$success_std)

cat("=============================================================\n")
if (all_targets_met && all_successful) {
  cat("PERFORMANCE BENCHMARK: PASS\\n")
  cat("=============================================================\n\n")
  cat("All performance targets met.\\n")
  cat("powerNMA automatic pathways scale efficiently to large networks.\\n")
  cat("Suitable for networks with 100+ studies and 30+ treatments.\\n")
} else {
  cat("PERFORMANCE BENCHMARK: REVIEW NEEDED\\n")
  cat("=============================================================\n\n")
  if (!all_successful) {
    cat("Some analyses failed. Review error messages above.\\n")
  }
  if (!all_targets_met) {
    cat("Some performance targets not met. Consider optimization.\\n")
  }
}

cat("\\n")

# Save results
results_file <- "performance_benchmark_results.csv"
write.csv(results_df, results_file, row.names = FALSE)
cat("Results saved to:", results_file, "\\n")

cat("\\n=============================================================\n")
cat("BENCHMARK COMPLETE\n")
cat("=============================================================\n")
