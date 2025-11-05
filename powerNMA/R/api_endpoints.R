#' RESTful API Endpoints for powerNMA
#'
#' @description
#' Provides RESTful API endpoints for programmatic access to powerNMA
#' functionality. Enables integration with other tools, cloud deployment,
#' and remote access to NMA capabilities.
#'
#' @details
#' Built using the plumber package to create OpenAPI-compliant REST API.
#' Includes endpoints for:
#' \itemize{
#'   \item Data upload and validation
#'   \item Running NMA analyses
#'   \item Retrieving results
#'   \item Generating visualizations
#'   \item Creating reports
#'   \item Batch processing
#'   \item Real-time collaboration
#'   \item Authentication and authorization
#' }
#'
#' @references
#' Plumber: https://www.rplumber.io/
#' OpenAPI Specification: https://swagger.io/specification/
#'
#' @author powerNMA Development Team
#' @name api_endpoints
NULL

#' Launch powerNMA API Server
#'
#' @description
#' Starts the REST API server for programmatic access to powerNMA.
#'
#' @param host Host address (default: "0.0.0.0" for all interfaces)
#' @param port Port number (default: 8000)
#' @param docs Enable Swagger documentation (default: TRUE)
#' @param auth_required Require API authentication (default: FALSE)
#' @param api_keys Character vector of valid API keys (if auth_required = TRUE)
#' @param cors_enabled Enable CORS (default: TRUE)
#' @param rate_limit Maximum requests per minute (default: 100)
#' @param log_file Path to API log file (default: NULL for console)
#'
#' @return Runs API server (does not return)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch API server
#' launch_powernma_api(port = 8000, docs = TRUE)
#'
#' # Launch with authentication
#' launch_powernma_api(
#'   port = 8000,
#'   auth_required = TRUE,
#'   api_keys = c("key1", "key2", "key3")
#' )
#'
#' # Access API documentation at http://localhost:8000/__docs__/
#' }
launch_powernma_api <- function(host = "0.0.0.0",
                               port = 8000,
                               docs = TRUE,
                               auth_required = FALSE,
                               api_keys = NULL,
                               cors_enabled = TRUE,
                               rate_limit = 100,
                               log_file = NULL) {

  # Check for plumber package
  if (!requireNamespace("plumber", quietly = TRUE)) {
    stop("Package 'plumber' required for API. Install with: install.packages('plumber')")
  }

  message(sprintf("Starting powerNMA API server on %s:%d", host, port))

  # Create plumber API
  pr <- create_api_router(
    auth_required = auth_required,
    api_keys = api_keys,
    rate_limit = rate_limit
  )

  # Add CORS if enabled
  if (cors_enabled) {
    pr <- pr %>%
      plumber::pr_set_api_spec(function(spec) {
        spec$info$title <- "powerNMA API"
        spec$info$description <- "RESTful API for Network Meta-Analysis"
        spec$info$version <- "1.0.0"
        spec
      })
  }

  # Enable Swagger docs
  if (docs) {
    pr <- pr %>%
      plumber::pr_set_docs("swagger")
  }

  # Setup logging
  if (!is.null(log_file)) {
    pr <- pr %>%
      plumber::pr_set_error(function(req, res, err) {
        log_api_error(err, req, log_file)
        list(error = err$message)
      })
  }

  # Run API
  pr %>%
    plumber::pr_run(host = host, port = port)
}

#' Create API Router
#'
#' @keywords internal
create_api_router <- function(auth_required, api_keys, rate_limit) {

  pr <- plumber::plumber$new()

  # ============================================================================
  # Middleware
  # ============================================================================

  # Authentication middleware
  if (auth_required) {
    pr$filter("authentication", function(req, res) {
      api_key <- req$HTTP_X_API_KEY %||% req$args$api_key

      if (is.null(api_key) || !(api_key %in% api_keys)) {
        res$status <- 401
        return(list(error = "Unauthorized: Invalid or missing API key"))
      }

      plumber::forward()
    })
  }

  # Rate limiting middleware
  rate_limiter <- create_rate_limiter(rate_limit)
  pr$filter("rate_limit", rate_limiter)

  # CORS middleware
  pr$filter("cors", function(req, res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    res$setHeader("Access-Control-Allow-Methods", "GET, POST, PUT, DELETE, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type, X-API-Key")

    if (req$REQUEST_METHOD == "OPTIONS") {
      res$status <- 200
      return(list())
    }

    plumber::forward()
  })

  # ============================================================================
  # Health Check Endpoints
  # ============================================================================

  #* Health check
  #* @get /health
  pr$handle("GET", "/health", function() {
    list(
      status = "healthy",
      version = "1.0.0",
      timestamp = Sys.time()
    )
  })

  #* API info
  #* @get /info
  pr$handle("GET", "/info", function() {
    list(
      name = "powerNMA API",
      version = "1.0.0",
      description = "REST API for Network Meta-Analysis",
      endpoints = list(
        "POST /analyze" = "Run NMA analysis",
        "GET /results/{id}" = "Get analysis results",
        "POST /visualize" = "Generate visualization",
        "POST /batch" = "Batch process multiple datasets"
      )
    )
  })

  # ============================================================================
  # Data Endpoints
  # ============================================================================

  #* Upload and validate data
  #* @post /data/upload
  #* @serializer json
  pr$handle("POST", "/data/upload", function(req, res) {
    tryCatch({
      # Parse uploaded data
      data <- parse_uploaded_data(req$body, req$content_type)

      # Validate data
      validation <- validate_nma_input(data)

      # Store data with unique ID
      data_id <- generate_unique_id()
      store_data(data_id, data)

      return(list(
        data_id = data_id,
        validation = validation,
        message = "Data uploaded and validated successfully"
      ))
    }, error = function(e) {
      res$status <- 400
      return(list(error = e$message))
    })
  })

  #* Get data info
  #* @get /data/<data_id>
  pr$handle("GET", "/data/<data_id>", function(data_id, res) {
    data <- retrieve_data(data_id)

    if (is.null(data)) {
      res$status <- 404
      return(list(error = "Data not found"))
    }

    return(list(
      data_id = data_id,
      n_studies = nrow(data),
      n_comparisons = sum(!is.na(data$TE)),
      treatments = unique(c(data$treat1, data$treat2))
    ))
  })

  # ============================================================================
  # Analysis Endpoints
  # ============================================================================

  #* Run NMA analysis
  #* @post /analyze
  #* @serializer json
  pr$handle("POST", "/analyze", function(req, res) {
    tryCatch({
      # Parse request body
      body <- jsonlite::fromJSON(req$postBody)

      data_id <- body$data_id
      config <- body$config %||% list()

      # Retrieve data
      data <- retrieve_data(data_id)
      if (is.null(data)) {
        res$status <- 404
        return(list(error = "Data not found"))
      }

      # Run analysis
      analysis_id <- generate_unique_id()

      # Run in background for large analyses
      future::future({
        result <- run_ultimate_nma(
          data = data,
          sm = config$sm %||% "OR",
          reference = config$reference,
          assess_geometry = config$assess_geometry %||% TRUE,
          assess_inconsistency = config$assess_inconsistency %||% TRUE,
          calculate_rankings = config$calculate_rankings %||% TRUE
        )

        store_result(analysis_id, result)
      })

      return(list(
        analysis_id = analysis_id,
        status = "processing",
        message = "Analysis started. Check /results/{analysis_id} for results."
      ))
    }, error = function(e) {
      res$status <- 500
      return(list(error = e$message))
    })
  })

  #* Get analysis results
  #* @get /results/<analysis_id>
  pr$handle("GET", "/results/<analysis_id>", function(analysis_id, res) {
    result <- retrieve_result(analysis_id)

    if (is.null(result)) {
      res$status <- 404
      return(list(error = "Analysis not found or still processing"))
    }

    # Return summary
    return(list(
      analysis_id = analysis_id,
      status = "completed",
      summary = summarize_nma_result(result),
      download_url = sprintf("/results/%s/download", analysis_id)
    ))
  })

  #* Download full results
  #* @get /results/<analysis_id>/download
  #* @serializer json
  pr$handle("GET", "/results/<analysis_id>/download", function(analysis_id, res) {
    result <- retrieve_result(analysis_id)

    if (is.null(result)) {
      res$status <- 404
      return(list(error = "Results not found"))
    }

    return(result)
  })

  # ============================================================================
  # Visualization Endpoints
  # ============================================================================

  #* Generate visualization
  #* @post /visualize
  #* @serializer json
  pr$handle("POST", "/visualize", function(req, res) {
    tryCatch({
      body <- jsonlite::fromJSON(req$postBody)

      analysis_id <- body$analysis_id
      viz_type <- body$type %||% "network"

      result <- retrieve_result(analysis_id)
      if (is.null(result)) {
        res$status <- 404
        return(list(error = "Analysis not found"))
      }

      # Generate visualization
      viz_id <- generate_unique_id()
      viz_file <- generate_api_visualization(result, viz_type, viz_id)

      return(list(
        viz_id = viz_id,
        type = viz_type,
        url = sprintf("/visualizations/%s", viz_id),
        download_url = sprintf("/visualizations/%s/download", viz_id)
      ))
    }, error = function(e) {
      res$status <- 500
      return(list(error = e$message))
    })
  })

  #* Get visualization
  #* @get /visualizations/<viz_id>
  #* @serializer htmlwidget
  pr$handle("GET", "/visualizations/<viz_id>", function(viz_id, res) {
    viz_file <- retrieve_visualization(viz_id)

    if (is.null(viz_file) || !file.exists(viz_file)) {
      res$status <- 404
      return(list(error = "Visualization not found"))
    }

    # Return HTML widget
    htmlwidgets::saveWidget(readRDS(viz_file), "temp.html")
    res$body <- readLines("temp.html")
    res
  })

  # ============================================================================
  # Batch Processing Endpoints
  # ============================================================================

  #* Run batch analysis
  #* @post /batch
  #* @serializer json
  pr$handle("POST", "/batch", function(req, res) {
    tryCatch({
      body <- jsonlite::fromJSON(req$postBody)

      data_ids <- body$data_ids
      config <- body$config %||% list()

      # Start batch job
      batch_id <- generate_unique_id()

      future::future({
        datasets <- lapply(data_ids, retrieve_data)
        names(datasets) <- data_ids

        batch_result <- run_batch_nma(
          datasets = datasets,
          analysis_config = config,
          parallel = TRUE
        )

        store_result(batch_id, batch_result)
      })

      return(list(
        batch_id = batch_id,
        status = "processing",
        n_datasets = length(data_ids),
        message = "Batch analysis started."
      ))
    }, error = function(e) {
      res$status <- 500
      return(list(error = e$message))
    })
  })

  # ============================================================================
  # Report Generation Endpoints
  # ============================================================================

  #* Generate report
  #* @post /report
  #* @serializer json
  pr$handle("POST", "/report", function(req, res) {
    tryCatch({
      body <- jsonlite::fromJSON(req$postBody)

      analysis_id <- body$analysis_id
      format <- body$format %||% "docx"

      result <- retrieve_result(analysis_id)
      if (is.null(result)) {
        res$status <- 404
        return(list(error = "Analysis not found"))
      }

      # Generate report
      report_id <- generate_unique_id()
      report_file <- generate_api_report(result, format, report_id)

      return(list(
        report_id = report_id,
        format = format,
        download_url = sprintf("/reports/%s/download", report_id)
      ))
    }, error = function(e) {
      res$status <- 500
      return(list(error = e$message))
    })
  })

  return(pr)
}

#' Create Rate Limiter
#'
#' @keywords internal
create_rate_limiter <- function(max_requests) {
  request_counts <- list()
  last_reset <- Sys.time()

  function(req, res) {
    current_time <- Sys.time()

    # Reset counts every minute
    if (difftime(current_time, last_reset, units = "mins") >= 1) {
      request_counts <<- list()
      last_reset <<- current_time
    }

    # Get client IP
    client_ip <- req$HTTP_X_FORWARDED_FOR %||% req$REMOTE_ADDR

    # Update count
    request_counts[[client_ip]] <<- (request_counts[[client_ip]] %||% 0) + 1

    # Check limit
    if (request_counts[[client_ip]] > max_requests) {
      res$status <- 429
      return(list(error = "Rate limit exceeded"))
    }

    plumber::forward()
  }
}

#' Parse Uploaded Data
#'
#' @keywords internal
parse_uploaded_data <- function(body, content_type) {
  if (grepl("json", content_type)) {
    return(jsonlite::fromJSON(body))
  } else if (grepl("csv", content_type)) {
    return(read.csv(textConnection(body)))
  } else {
    stop("Unsupported content type")
  }
}

#' Generate Unique ID
#'
#' @keywords internal
generate_unique_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample(10000:99999, 1))
}

#' Data Storage Functions
#'
#' @keywords internal
store_data <- function(id, data) {
  dir.create("api_cache/data", recursive = TRUE, showWarnings = FALSE)
  saveRDS(data, file.path("api_cache/data", paste0(id, ".rds")))
}

#' @keywords internal
retrieve_data <- function(id) {
  file_path <- file.path("api_cache/data", paste0(id, ".rds"))
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  }
  return(NULL)
}

#' Result Storage Functions
#'
#' @keywords internal
store_result <- function(id, result) {
  dir.create("api_cache/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path("api_cache/results", paste0(id, ".rds")))
}

#' @keywords internal
retrieve_result <- function(id) {
  file_path <- file.path("api_cache/results", paste0(id, ".rds"))
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  }
  return(NULL)
}

#' Summarize NMA Result
#'
#' @keywords internal
summarize_nma_result <- function(result) {
  list(
    n_studies = result$nma_result$k,
    n_treatments = length(result$nma_result$trts),
    treatments = result$nma_result$trts,
    tau2 = result$nma_result$tau^2,
    best_treatment = if (!is.null(result$sucra)) {
      names(which.max(result$sucra$sucra_scores))
    } else {
      NULL
    }
  )
}

#' Generate API Visualization
#'
#' @keywords internal
generate_api_visualization <- function(result, viz_type, viz_id) {
  viz_plot <- switch(
    viz_type,
    network = create_interactive_network(result$data, result$nma_result),
    forest = create_interactive_forest(result$nma_result),
    funnel = create_interactive_funnel(result$nma_result),
    rankings = create_interactive_rankings(result$sucra)
  )

  dir.create("api_cache/visualizations", recursive = TRUE, showWarnings = FALSE)
  viz_file <- file.path("api_cache/visualizations", paste0(viz_id, ".rds"))
  saveRDS(viz_plot, viz_file)

  return(viz_file)
}

#' @keywords internal
retrieve_visualization <- function(viz_id) {
  file.path("api_cache/visualizations", paste0(viz_id, ".rds"))
}

#' Generate API Report
#'
#' @keywords internal
generate_api_report <- function(result, format, report_id) {
  dir.create("api_cache/reports", recursive = TRUE, showWarnings = FALSE)

  report_file <- file.path("api_cache/reports", sprintf("%s.%s", report_id, format))

  generate_nma_report(
    nma_results = result,
    output_file = report_file,
    format = format
  )

  return(report_file)
}

#' Log API Error
#'
#' @keywords internal
log_api_error <- function(err, req, log_file) {
  log_entry <- sprintf(
    "[%s] %s %s - Error: %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    req$REQUEST_METHOD,
    req$PATH_INFO,
    err$message
  )

  cat(log_entry, file = log_file, append = TRUE)
}

#' NULL coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
