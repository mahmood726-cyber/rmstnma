# ============================================================================
# Progress Indicators and Performance Utilities
# ============================================================================

#' Simple Progress Bar for powerNMA Operations
#'
#' Creates a text-based progress indicator for long-running operations.
#'
#' @param total Total number of iterations
#' @param label Optional label for the progress bar
#' @return Progress bar object with update() method
#' @keywords internal
#' @examples
#' \dontrun{
#' pb <- create_progress_bar(100, "Processing studies")
#' for (i in 1:100) {
#'   # Do work
#'   pb$update(i)
#' }
#' pb$finish()
#' }
create_progress_bar <- function(total, label = "Progress") {
  start_time <- Sys.time()
  last_print <- 0

  pb <- list(
    total = total,
    label = label,
    start_time = start_time,

    update = function(current) {
      pct <- round(100 * current / total)

      # Only print every 10% to avoid spam
      if (pct >= last_print + 10 || current == total) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        eta <- if (current > 0) elapsed * (total - current) / current else 0

        cat(sprintf("\r%s: %d%% [%d/%d] (ETA: %.0fs)",
                   label, pct, current, total, eta))

        last_print <<- pct

        if (current == total) {
          cat(" ✓\n")
        }
      }
      flush.console()
    },

    finish = function() {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(sprintf("\r%s: Complete! (%.2fs)              \n", label, elapsed))
      flush.console()
    }
  )

  class(pb) <- "powernma_progress"
  pb
}

#' Time a Code Block
#'
#' Measures execution time and prints result.
#'
#' @param expr Expression to time
#' @param label Optional label for the timing
#' @return Result of expr (invisibly)
#' @export
#' @examples
#' \dontrun{
#' result <- time_operation({
#'   data <- simulate_nma_data(n_studies = 100)
#'   run_powernma(data, data_type = "pairwise")
#' }, label = "Large NMA")
#' }
time_operation <- function(expr, label = "Operation") {
  start_time <- Sys.time()
  result <- force(expr)
  elapsed <- difftime(Sys.time(), start_time, units = "secs")

  cat(sprintf("⏱  %s completed in %.2f seconds\n", label, as.numeric(elapsed)))

  invisible(result)
}

#' Benchmark Multiple NMA Configurations
#'
#' Compares performance of different NMA configurations.
#'
#' @param data NMA data
#' @param configs List of configuration objects
#' @param config_names Names for each configuration
#' @return Data frame with timing results
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_nma_data(n_studies = 30)
#'
#' configs <- list(
#'   setup_powernma(use_bayesian = FALSE, run_sensitivity = FALSE),
#'   setup_powernma(use_bayesian = TRUE, run_sensitivity = FALSE),
#'   setup_powernma(use_bayesian = TRUE, run_sensitivity = TRUE)
#' )
#'
#' results <- benchmark_configs(
#'   data,
#'   configs,
#'   c("Basic", "Bayesian", "Full")
#' )
#' print(results)
#' }
benchmark_configs <- function(data, configs, config_names = NULL) {
  if (is.null(config_names)) {
    config_names <- paste0("Config_", seq_along(configs))
  }

  results <- data.frame(
    Configuration = config_names,
    Time_Seconds = numeric(length(configs)),
    Success = logical(length(configs)),
    stringsAsFactors = FALSE
  )

  cat("\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("  Benchmarking NMA Configurations\n")
  cat("═══════════════════════════════════════════════════\n\n")

  for (i in seq_along(configs)) {
    cat(sprintf("Testing: %s...\n", config_names[i]))

    start_time <- Sys.time()
    result <- tryCatch({
      run_powernma(
        data = data,
        data_type = "pairwise",
        mode = "standard",
        config = configs[[i]]
      )
      TRUE
    }, error = function(e) {
      cat(sprintf("  Error: %s\n", e$message))
      FALSE
    })

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    results$Time_Seconds[i] <- elapsed
    results$Success[i] <- result

    if (result) {
      cat(sprintf("  ✓ Completed in %.2f seconds\n\n", elapsed))
    } else {
      cat("  ✗ Failed\n\n")
    }
  }

  # Summary
  cat("═══════════════════════════════════════════════════\n")
  cat("  Benchmark Results\n")
  cat("═══════════════════════════════════════════════════\n\n")

  print(results, row.names = FALSE)

  if (all(results$Success)) {
    fastest <- which.min(results$Time_Seconds)
    slowest <- which.max(results$Time_Seconds)
    speedup <- results$Time_Seconds[slowest] / results$Time_Seconds[fastest]

    cat(sprintf("\n🏆 Fastest: %s (%.2fs)\n",
                results$Configuration[fastest],
                results$Time_Seconds[fastest]))
    cat(sprintf("🐌 Slowest: %s (%.2fs)\n",
                results$Configuration[slowest],
                results$Time_Seconds[slowest]))
    cat(sprintf("⚡ Speedup: %.1fx\n\n", speedup))
  }

  invisible(results)
}

#' Memory Usage Reporter
#'
#' Reports current memory usage and provides recommendations.
#'
#' @export
#' @examples
#' report_memory_usage()
report_memory_usage <- function() {
  if (.Platform$OS.type == "windows") {
    mem_info <- memory.size()
    mem_max <- memory.limit()
    cat(sprintf("💾 Memory: %.0f MB used (limit: %.0f MB)\n",
                mem_info, mem_max))
  } else {
    # Try to get memory info on Unix-like systems
    mem_info <- gc()
    used_mb <- sum(mem_info[, 2])
    cat(sprintf("💾 Memory: %.0f MB currently allocated\n", used_mb))
  }

  # Run garbage collection
  gc_result <- gc(verbose = FALSE)
  cat("🧹 Garbage collection completed\n")

  invisible(gc_result)
}

#' Estimate Memory Requirements
#'
#' Estimates memory needed for a given network size.
#'
#' @param n_studies Number of studies
#' @param n_treatments Number of treatments
#' @param use_bayesian Whether Bayesian analysis will be used
#' @return Estimated memory in MB
#' @export
#' @examples
#' estimate_memory_needs(n_studies = 100, n_treatments = 10)
estimate_memory_needs <- function(n_studies, n_treatments, use_bayesian = FALSE) {
  # Rough estimates based on typical NMA memory usage
  base_mb <- 50  # Base R + package loading
  data_mb <- (n_studies * n_treatments) / 1000  # Data storage
  analysis_mb <- (n_treatments^2) / 10  # Network computations

  if (use_bayesian) {
    mcmc_mb <- n_treatments * 5  # MCMC chains
  } else {
    mcmc_mb <- 0
  }

  total_mb <- base_mb + data_mb + analysis_mb + mcmc_mb

  cat(sprintf("📊 Estimated memory requirements:\n"))
  cat(sprintf("   Base:           %6.1f MB\n", base_mb))
  cat(sprintf("   Data:           %6.1f MB\n", data_mb))
  cat(sprintf("   Analysis:       %6.1f MB\n", analysis_mb))
  if (use_bayesian) {
    cat(sprintf("   MCMC:           %6.1f MB\n", mcmc_mb))
  }
  cat(sprintf("   ───────────────────────────\n"))
  cat(sprintf("   Total:          %6.1f MB\n\n", total_mb))

  if (total_mb > 1000) {
    cat("⚠️  Warning: Analysis may require >1GB of memory\n")
    cat("   Consider:\n")
    cat("   • Using smaller subsets\n")
    cat("   • Disabling Bayesian analysis\n")
    cat("   • Running on a machine with more RAM\n\n")
  } else if (total_mb > 500) {
    cat("💡 Moderate memory requirements (>500MB)\n")
    cat("   Should be fine on most modern computers\n\n")
  } else {
    cat("✅ Low memory requirements (<500MB)\n")
    cat("   Should run smoothly\n\n")
  }

  invisible(total_mb)
}
