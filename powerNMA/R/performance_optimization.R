#' Performance Optimization Utilities
#'
#' @description
#' Production-grade performance optimization utilities including:
#' - Intelligent memoization for expensive computations
#' - Parallel processing infrastructure
#' - Progress tracking for long operations
#' - Memory-efficient data handling
#'
#' These utilities implement performance best practices identified
#' in comprehensive codebase review.
#'
#' @author powerNMA Development Team
#' @name performance_optimization
NULL

# Global cache environment for memoization
.powernma_memo_cache <- new.env(parent = emptyenv())

#' Clear Memoization Cache
#'
#' @description
#' Clears the global memoization cache to free memory.
#'
#' @param pattern Optional pattern to match function names (regex)
#'
#' @return Invisible number of items removed
#' @export
#'
#' @examples
#' clear_memo_cache()  # Clear all
#' clear_memo_cache("league")  # Clear league table cache
clear_memo_cache <- function(pattern = NULL) {
  if (is.null(pattern)) {
    # Clear all
    n_items <- length(ls(.powernma_memo_cache))
    rm(list = ls(.powernma_memo_cache), envir = .powernma_memo_cache)
  } else {
    # Clear matching pattern
    items <- ls(.powernma_memo_cache)
    matching <- grep(pattern, items, value = TRUE)
    n_items <- length(matching)
    if (n_items > 0) {
      rm(list = matching, envir = .powernma_memo_cache)
    }
  }

  message(sprintf("Cleared %d cached item(s)", n_items))
  invisible(n_items)
}

#' Get Cache Statistics
#'
#' @description
#' Returns statistics about the memoization cache.
#'
#' @return List with cache statistics
#' @export
#'
#' @examples
#' get_cache_stats()
get_cache_stats <- function() {
  items <- ls(.powernma_memo_cache)
  n_items <- length(items)

  if (n_items == 0) {
    return(list(
      n_items = 0,
      total_size_mb = 0,
      items = character(0)
    ))
  }

  # Calculate total size
  sizes <- sapply(items, function(item) {
    obj <- get(item, envir = .powernma_memo_cache)
    object.size(obj)
  })

  total_size_bytes <- sum(sizes)
  total_size_mb <- total_size_bytes / (1024^2)

  list(
    n_items = n_items,
    total_size_mb = round(total_size_mb, 2),
    items = items,
    sizes_mb = round(sizes / (1024^2), 2)
  )
}

#' Memoize Function
#'
#' @description
#' Creates a memoized version of a function that caches results
#' based on input arguments. Useful for expensive computations.
#'
#' @param func Function to memoize
#' @param cache_key_func Optional function to generate cache keys
#'   from arguments. Default uses deparse of all arguments.
#' @param max_cache_size Maximum number of cached results (default: 100)
#' @param expire_seconds Cache expiration time in seconds (default: NULL, never expires)
#'
#' @return Memoized function
#' @export
#'
#' @examples
#' \dontrun{
#' # Memoize an expensive function
#' expensive_calc <- function(x, y) {
#'   Sys.sleep(2)  # Simulate expensive computation
#'   x + y
#' }
#' fast_calc <- memoize(expensive_calc)
#'
#' system.time(fast_calc(5, 3))  # Slow first time
#' system.time(fast_calc(5, 3))  # Fast second time (cached)
#' }
memoize <- function(func, cache_key_func = NULL, max_cache_size = 100,
                    expire_seconds = NULL) {
  func_name <- deparse(substitute(func))

  # Create cache for this function
  cache_name <- paste0("memo_", func_name, "_", as.integer(Sys.time()))

  function(...) {
    # Generate cache key
    if (is.null(cache_key_func)) {
      args_list <- list(...)
      cache_key <- digest::digest(args_list)
    } else {
      cache_key <- cache_key_func(...)
    }

    full_key <- paste0(cache_name, "_", cache_key)

    # Check if result is cached
    if (exists(full_key, envir = .powernma_memo_cache)) {
      cached_item <- get(full_key, envir = .powernma_memo_cache)

      # Check expiration
      if (!is.null(expire_seconds)) {
        time_elapsed <- as.numeric(Sys.time() - cached_item$timestamp)
        if (time_elapsed > expire_seconds) {
          # Expired, remove from cache
          rm(list = full_key, envir = .powernma_memo_cache)
        } else {
          # Valid cache hit
          return(cached_item$result)
        }
      } else {
        # No expiration, return cached result
        return(cached_item$result)
      }
    }

    # Compute result
    result <- func(...)

    # Store in cache
    cached_item <- list(
      result = result,
      timestamp = Sys.time()
    )
    assign(full_key, cached_item, envir = .powernma_memo_cache)

    # Enforce max cache size (remove oldest if exceeded)
    cache_items <- ls(.powernma_memo_cache, pattern = paste0("^", cache_name))
    if (length(cache_items) > max_cache_size) {
      # Get timestamps
      timestamps <- sapply(cache_items, function(item) {
        get(item, envir = .powernma_memo_cache)$timestamp
      })
      # Remove oldest
      oldest <- names(sort(timestamps))[1]
      rm(list = oldest, envir = .powernma_memo_cache)
    }

    result
  }
}

#' Parallel Apply with Progress
#'
#' @description
#' Applies a function in parallel with optional progress tracking.
#' Automatically detects available cores and uses sensible defaults.
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param n_cores Number of cores to use (default: auto-detect - 1)
#' @param show_progress Show progress bar (default: TRUE)
#' @param progress_label Label for progress bar
#' @param ... Additional arguments passed to FUN
#'
#' @return List of results
#' @export
#'
#' @examples
#' \dontrun{
#' # Parallel simulation
#' results <- parallel_apply(1:100, function(i) {
#'   simulate_and_analyze(n_studies = 50)
#' }, show_progress = TRUE)
#' }
parallel_apply <- function(X, FUN, n_cores = NULL, show_progress = TRUE,
                           progress_label = "Processing", ...) {
  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  n_items <- length(X)

  # Use sequential processing if too few items or only 1 core
  if (n_items < n_cores || n_cores == 1) {
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = n_items, style = 3)
      on.exit(close(pb))

      results <- lapply(seq_along(X), function(i) {
        result <- FUN(X[[i]], ...)
        setTxtProgressBar(pb, i)
        result
      })
    } else {
      results <- lapply(X, FUN, ...)
    }
    return(results)
  }

  # Check if future package is available
  if (requireNamespace("future.apply", quietly = TRUE)) {
    # Use future.apply for parallel processing
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    future::plan(future::multisession, workers = n_cores)

    if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
      # Use progressr for progress tracking
      progressr::with_progress({
        p <- progressr::progressor(along = X)
        results <- future.apply::future_lapply(X, function(x) {
          result <- FUN(x, ...)
          p()
          result
        }, future.seed = TRUE)
      })
    } else {
      results <- future.apply::future_lapply(X, FUN, future.seed = TRUE, ...)
    }
  } else if (requireNamespace("parallel", quietly = TRUE)) {
    # Fall back to parallel package
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export required objects
    parallel::clusterExport(cl, varlist = ls(envir = parent.frame()),
                           envir = parent.frame())

    if (show_progress) {
      # Simple progress tracking
      message(sprintf("Processing %d items on %d cores...", n_items, n_cores))
      results <- parallel::parLapply(cl, X, FUN, ...)
      message("Done!")
    } else {
      results <- parallel::parLapply(cl, X, FUN, ...)
    }
  } else {
    # No parallel backend available, use sequential
    warning("No parallel backend available (future.apply or parallel). Using sequential processing.")
    results <- lapply(X, FUN, ...)
  }

  results
}

#' Batch Process with Chunking
#'
#' @description
#' Processes large datasets in chunks to manage memory usage.
#'
#' @param X Vector or list to process
#' @param FUN Function to apply to each chunk
#' @param chunk_size Size of each chunk (default: 1000)
#' @param combine_func Function to combine chunk results (default: c)
#' @param show_progress Show progress (default: TRUE)
#' @param ... Additional arguments passed to FUN
#'
#' @return Combined results from all chunks
#' @export
#'
#' @examples
#' \dontrun{
#' # Process large dataset in chunks
#' large_data <- 1:100000
#' results <- batch_process(large_data, function(chunk) {
#'   mean(chunk)
#' }, chunk_size = 10000)
#' }
batch_process <- function(X, FUN, chunk_size = 1000,
                         combine_func = c, show_progress = TRUE, ...) {
  n_items <- length(X)
  n_chunks <- ceiling(n_items / chunk_size)

  if (show_progress) {
    message(sprintf("Processing %d items in %d chunks of size %d",
                    n_items, n_chunks, chunk_size))
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
    on.exit(close(pb))
  }

  chunk_results <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_items)
    chunk <- X[start_idx:end_idx]

    chunk_results[[i]] <- FUN(chunk, ...)

    if (show_progress) {
      setTxtProgressBar(pb, i)
    }

    # Force garbage collection after each chunk
    if (i %% 10 == 0) {
      gc(verbose = FALSE)
    }
  }

  # Combine results
  do.call(combine_func, chunk_results)
}

#' Optimize Network Structure
#'
#' @description
#' Optimizes network data structure for faster computations.
#' Converts character columns to factors, removes unnecessary columns, etc.
#'
#' @param data Network data frame
#' @param keep_cols Columns to keep (NULL = keep all)
#'
#' @return Optimized data frame
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data_opt <- optimize_network_data(data)
#' }
optimize_network_data <- function(data, keep_cols = NULL) {
  if (!is.null(keep_cols)) {
    data <- data[, keep_cols, drop = FALSE]
  }

  # Convert character columns to factors for memory efficiency
  char_cols <- sapply(data, is.character)
  for (col in names(data)[char_cols]) {
    if (length(unique(data[[col]])) < nrow(data) / 2) {
      # Only convert if it saves memory (< 50% unique values)
      data[[col]] <- factor(data[[col]])
    }
  }

  data
}

#' Measure Execution Time
#'
#' @description
#' Measures execution time of an expression with formatted output.
#'
#' @param expr Expression to time
#' @param label Label for timing output
#' @param silent Suppress output (default: FALSE)
#'
#' @return List with result and timing information
#' @export
#'
#' @examples
#' \dontrun{
#' result <- measure_time({
#'   expensive_operation()
#' }, label = "Expensive operation")
#' }
measure_time <- function(expr, label = "Execution", silent = FALSE) {
  start_time <- Sys.time()

  result <- expr

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (!silent) {
    if (elapsed < 1) {
      time_str <- sprintf("%.0f ms", elapsed * 1000)
    } else if (elapsed < 60) {
      time_str <- sprintf("%.2f sec", elapsed)
    } else {
      time_str <- sprintf("%.2f min", elapsed / 60)
    }

    message(sprintf("%s completed in %s", label, time_str))
  }

  list(
    result = result,
    elapsed_seconds = elapsed,
    start_time = start_time,
    end_time = end_time
  )
}

#' Profile Function Performance
#'
#' @description
#' Profiles a function to identify performance bottlenecks.
#'
#' @param func Function to profile
#' @param ... Arguments to pass to func
#' @param n_runs Number of runs for benchmarking (default: 10)
#'
#' @return List with profiling results
#' @export
#'
#' @examples
#' \dontrun{
#' profile_results <- profile_function(my_function, data = data, n_runs = 100)
#' print(profile_results)
#' }
profile_function <- function(func, ..., n_runs = 10) {
  if (requireNamespace("microbenchmark", quietly = TRUE)) {
    # Use microbenchmark for accurate timing
    benchmark <- microbenchmark::microbenchmark(
      func(...),
      times = n_runs
    )

    list(
      mean_ms = mean(benchmark$time) / 1e6,
      median_ms = median(benchmark$time) / 1e6,
      min_ms = min(benchmark$time) / 1e6,
      max_ms = max(benchmark$time) / 1e6,
      benchmark = benchmark
    )
  } else {
    # Fall back to system.time
    times <- replicate(n_runs, {
      timing <- system.time(func(...))
      timing["elapsed"]
    })

    list(
      mean_ms = mean(times) * 1000,
      median_ms = median(times) * 1000,
      min_ms = min(times) * 1000,
      max_ms = max(times) * 1000,
      times = times
    )
  }
}

#' Memory-Efficient Data Frame Creation
#'
#' @description
#' Creates data frames with appropriate column types to minimize memory usage.
#'
#' @param ... Name-value pairs for columns
#' @param stringsAsFactors Convert strings to factors (default: TRUE for small cardinality)
#'
#' @return Data frame
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' df <- efficient_data_frame(
#'   studlab = rep(paste0("Study", 1:10), each = 3),
#'   treatment = rep(c("A", "B", "C"), 10)
#' )
#' }
efficient_data_frame <- function(..., stringsAsFactors = NULL) {
  df <- data.frame(..., stringsAsFactors = FALSE)

  # Auto-detect if we should use factors
  if (is.null(stringsAsFactors)) {
    char_cols <- sapply(df, is.character)
    for (col in names(df)[char_cols]) {
      # Use factor if < 50% unique values
      if (length(unique(df[[col]])) < nrow(df) / 2) {
        df[[col]] <- factor(df[[col]])
      }
    }
  } else if (stringsAsFactors) {
    char_cols <- sapply(df, is.character)
    for (col in names(df)[char_cols]) {
      df[[col]] <- factor(df[[col]])
    }
  }

  df
}

#' Get System Performance Info
#'
#' @description
#' Returns information about system resources and R session.
#'
#' @return List with system information
#' @export
#'
#' @examples
#' get_system_info()
get_system_info <- function() {
  info <- list(
    r_version = R.version.string,
    platform = R.version$platform,
    n_cores = parallel::detectCores(),
    memory_limit_mb = if (.Platform$OS.type == "windows") {
      memory.limit()
    } else {
      NA
    },
    session_memory_mb = round(sum(gc()[, 2])),
    parallel_backend = if (requireNamespace("future", quietly = TRUE)) {
      "future"
    } else if (requireNamespace("parallel", quietly = TRUE)) {
      "parallel"
    } else {
      "none"
    }
  )

  info
}

#' Estimate Memory Usage
#'
#' @description
#' Estimates memory usage of an object or expression.
#'
#' @param x Object or expression to measure
#'
#' @return Memory size in MB
#' @export
#'
#' @examples
#' data <- matrix(rnorm(10000), ncol = 100)
#' estimate_memory(data)
estimate_memory <- function(x) {
  size_bytes <- object.size(x)
  size_mb <- as.numeric(size_bytes) / (1024^2)
  round(size_mb, 2)
}
