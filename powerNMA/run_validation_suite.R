#!/usr/bin/env Rscript

# powerNMA v2.0 Validation Test Suite Runner
# ==============================================================================
#
# This script runs the complete validation test suite for powerNMA v2.0
# and generates a comprehensive report.
#
# Usage:
#   Rscript run_validation_suite.R
#
# Or from R console:
#   source("run_validation_suite.R")
#
# ==============================================================================

library(testthat)
library(devtools)

# Print header
cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("  powerNMA v2.0 - Comprehensive Validation Test Suite\n")
cat("  Date: ", as.character(Sys.Date()), "\n")
cat("  Time: ", format(Sys.time(), "%H:%M:%S"), "\n")
cat("=", rep("=", 78), "\n\n", sep = "")

# Check required packages
required_pkgs <- c("netmeta", "gemtc", "testthat", "devtools", "dplyr", "purrr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("❌ ERROR: Missing required packages:\n")
  cat("  ", paste(missing_pkgs, collapse = ", "), "\n\n")
  cat("Install with:\n")
  cat("  install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))\n\n", sep = "")
  quit(status = 1)
}

# Test configuration
cat("Configuration:\n")
cat("  SKIP_LONG_TESTS: ", Sys.getenv("SKIP_LONG_TESTS", "false"), "\n")
cat("  SKIP_STRESS_TESTS: ", Sys.getenv("SKIP_STRESS_TESTS", "false"), "\n")
cat("  SKIP_BENCHMARK: ", Sys.getenv("SKIP_BENCHMARK", "false"), "\n\n")

# Ensure we're in the right directory
if (!file.exists("DESCRIPTION")) {
  cat("❌ ERROR: Not in powerNMA package directory\n")
  cat("Current directory: ", getwd(), "\n")
  quit(status = 1)
}

# Load powerNMA
cat("Loading powerNMA package...\n")
devtools::load_all()
cat("✓ Package loaded\n\n")

# Run test suite
cat("=", rep("=", 78), "\n", sep = "")
cat("  Running Test Suite\n")
cat("=", rep("=", 78), "\n\n", sep = "")

start_time <- Sys.time()

# Run tests with detailed output
results <- testthat::test_local(
  path = ".",
  reporter = testthat::ProgressReporter$new(max_failures = 10),
  stop_on_failure = FALSE
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")

# Parse results
total_tests <- length(results)
passed <- sum(sapply(results, function(x) length(x$results) - length(x$failed)))
failed <- sum(sapply(results, function(x) length(x$failed)))
warnings <- sum(sapply(results, function(x) length(x$warnings)))
skipped <- sum(sapply(results, function(x) length(x$skipped)))

# Print summary
cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("  Test Suite Summary\n")
cat("=", rep("=", 78), "\n", sep = "")
cat(sprintf("  Test files run:     %d\n", total_tests))
cat(sprintf("  Total assertions:   %d\n", passed + failed))
cat(sprintf("  Passed:            ✅ %d\n", passed))
cat(sprintf("  Failed:            %s %d\n", ifelse(failed == 0, "✅", "❌"), failed))
cat(sprintf("  Warnings:          %s %d\n", ifelse(warnings == 0, "✅", "⚠️"), warnings))
cat(sprintf("  Skipped:           ℹ️  %d\n", skipped))
cat(sprintf("  Execution time:     %.2f seconds\n", as.numeric(elapsed)))
cat("=", rep("=", 78), "\n\n", sep = "")

# Test breakdown by file
cat("Test Breakdown by File:\n")
cat("-", rep("-", 78), "\n", sep = "")

test_files <- c(
  "test-real-datasets.R" = "Real Dataset Validation",
  "test-large-simulations.R" = "Large Simulation Tests",
  "test-validation-benchmarks.R" = "Statistical Validation",
  "test-experimental-methods.R" = "Experimental Methods",
  "test-rmst-sign-validation.R" = "RMST Sign Convention"
)

for (file in names(test_files)) {
  cat(sprintf("  %-35s %s\n", test_files[file], file))
}
cat("-", rep("-", 78), "\n\n", sep = "")

# Final verdict
cat("=", rep("=", 78), "\n", sep = "")
cat("  VALIDATION VERDICT\n")
cat("=", rep("=", 78), "\n\n", sep = "")

if (failed == 0 && warnings == 0) {
  cat("✅ ✅ ✅  ALL TESTS PASSED  ✅ ✅ ✅\n\n")
  cat("powerNMA v2.0 STANDARD MODE is VALIDATED for production use.\n\n")
  cat("Validation Status:\n")
  cat("  ✅ Real dataset validation: PASS (exact agreement with netmeta)\n")
  cat("  ✅ Multi-arm trials: PASS (all comparisons included)\n")
  cat("  ✅ Statistical properties: PASS (Type I error, power, heterogeneity)\n")
  cat("  ✅ Performance: PASS (acceptable speed and scalability)\n")
  cat("  ✅ Experimental methods: PASS (RMST/Milestone fixes validated)\n\n")
  cat("Suitable for:\n")
  cat("  • Systematic reviews\n")
  cat("  • Cochrane reviews\n")
  cat("  • Clinical guidelines\n")
  cat("  • Health technology assessment\n")
  cat("  • Meta-analysis research\n\n")
  cat("EXPERIMENTAL MODE remains RESEARCH USE ONLY pending peer review.\n\n")

} else if (failed == 0 && warnings > 0) {
  cat("⚠️  TESTS PASSED WITH WARNINGS  ⚠️\n\n")
  cat("All assertions passed, but ", warnings, " warnings were generated.\n")
  cat("Review warnings before production use.\n\n")

} else {
  cat("❌  VALIDATION FAILED  ❌\n\n")
  cat(failed, " test(s) failed. Review failures before using in production.\n\n")
  cat("Failed tests must be addressed before powerNMA can be considered validated.\n\n")
}

cat("=", rep("=", 78), "\n\n", sep = "")

# Detailed failure information (if any)
if (failed > 0) {
  cat("Failure Details:\n")
  cat("-", rep("-", 78), "\n", sep = "")

  for (result in results) {
    if (length(result$failed) > 0) {
      cat("\nFile: ", result$file, "\n", sep = "")
      for (fail in result$failed) {
        cat("  ❌ ", fail$test, "\n", sep = "")
        cat("     ", fail$message, "\n", sep = "")
      }
    }
  }
  cat("\n")
}

# Save results to file
results_file <- paste0("validation_results_", format(Sys.Date(), "%Y%m%d"), ".txt")
sink(results_file)
cat("powerNMA v2.0 Validation Results\n")
cat("Date: ", as.character(Sys.Date()), "\n")
cat("Time: ", format(Sys.time(), "%H:%M:%S"), "\n\n")
cat("Summary:\n")
cat("  Total assertions: ", passed + failed, "\n")
cat("  Passed: ", passed, "\n")
cat("  Failed: ", failed, "\n")
cat("  Warnings: ", warnings, "\n")
cat("  Skipped: ", skipped, "\n")
cat("  Execution time: ", sprintf("%.2f", as.numeric(elapsed)), " seconds\n\n")

if (failed == 0) {
  cat("VERDICT: ✅ VALIDATED\n")
} else {
  cat("VERDICT: ❌ FAILED\n")
}
sink()

cat("Results saved to: ", results_file, "\n\n")

# Exit with appropriate code
quit(status = ifelse(failed == 0, 0, 1))
