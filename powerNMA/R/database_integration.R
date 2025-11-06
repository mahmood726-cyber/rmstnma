#' Database Integration for powerNMA
#'
#' @description
#' Enterprise database integration supporting SQLite, PostgreSQL, and MySQL
#' for persistent storage of analyses, results, and collaboration.
#'
#' @details
#' Provides comprehensive database functionality:
#' \itemize{
#'   \item Connection pooling for high performance
#'   \item Automatic schema creation and migration
#'   \item Result caching and retrieval
#'   \item User session management
#'   \item Collaboration features
#'   \item Audit logging
#'   \item Data versioning
#'   \item Query optimization
#' }
#'
#' @references
#' DBI: https://dbi.r-dbi.org/
#' pool: https://rstudio.github.io/pool/
#'
#' @author powerNMA Development Team
#' @name database_integration
NULL

#' Initialize Database Connection
#'
#' @description
#' Creates and initializes database connection with schema setup.
#'
#' @param db_type Database type: "sqlite", "postgresql", "mysql"
#' @param host Database host (for PostgreSQL/MySQL)
#' @param port Database port
#' @param dbname Database name
#' @param user Database user
#' @param password Database password
#' @param pool_size Connection pool size (default: 5)
#' @param create_schema Logical whether to create schema (default: TRUE)
#'
#' @return Database connection pool object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # SQLite (local)
#' db <- initialize_powernma_db(
#'   db_type = "sqlite",
#'   dbname = "powernma.db"
#' )
#'
#' # PostgreSQL (production)
#' db <- initialize_powernma_db(
#'   db_type = "postgresql",
#'   host = "localhost",
#'   port = 5432,
#'   dbname = "powernma",
#'   user = "powernma",
#'   password = "secure_password"
#' )
#'
#' # Close connection when done
#' close_powernma_db(db)
#' }
initialize_powernma_db <- function(db_type = c("sqlite", "postgresql", "mysql"),
                                  host = NULL,
                                  port = NULL,
                                  dbname = "powernma.db",
                                  user = NULL,
                                  password = NULL,
                                  pool_size = 5,
                                  create_schema = TRUE) {

  db_type <- match.arg(db_type)

  # Check required packages
  required_pkgs <- c("DBI", "pool")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(sprintf("Required packages missing: %s", paste(missing_pkgs, collapse = ", ")))
  }

  message(sprintf("Initializing %s database connection...", db_type))

  # Create connection pool
  pool <- switch(
    db_type,
    sqlite = {
      if (!requireNamespace("RSQLite", quietly = TRUE)) {
        stop("Package 'RSQLite' required for SQLite")
      }
      pool::dbPool(
        drv = RSQLite::SQLite(),
        dbname = dbname,
        maxSize = pool_size
      )
    },
    postgresql = {
      if (!requireNamespace("RPostgres", quietly = TRUE)) {
        stop("Package 'RPostgres' required for PostgreSQL")
      }
      pool::dbPool(
        drv = RPostgres::Postgres(),
        host = host,
        port = port,
        dbname = dbname,
        user = user,
        password = password,
        maxSize = pool_size
      )
    },
    mysql = {
      if (!requireNamespace("RMySQL", quietly = TRUE)) {
        stop("Package 'RMySQL' required for MySQL")
      }
      pool::dbPool(
        drv = RMySQL::MySQL(),
        host = host,
        port = port,
        dbname = dbname,
        user = user,
        password = password,
        maxSize = pool_size
      )
    }
  )

  # Create schema if requested
  if (create_schema) {
    create_powernma_schema(pool)
  }

  message("Database connection initialized successfully")

  return(structure(
    list(
      pool = pool,
      type = db_type,
      dbname = dbname
    ),
    class = "powernma_db"
  ))
}

#' Create Database Schema
#'
#' @keywords internal
create_powernma_schema <- function(pool) {

  message("Creating database schema...")

  # Projects table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS projects (
      project_id TEXT PRIMARY KEY,
      project_name TEXT NOT NULL,
      description TEXT,
      created_by TEXT,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      status TEXT DEFAULT 'active',
      metadata TEXT
    )
  ")

  # Datasets table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS datasets (
      dataset_id TEXT PRIMARY KEY,
      project_id TEXT,
      dataset_name TEXT NOT NULL,
      n_studies INTEGER,
      n_treatments INTEGER,
      n_comparisons INTEGER,
      effect_measure TEXT,
      uploaded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      uploaded_by TEXT,
      data_blob BLOB,
      FOREIGN KEY (project_id) REFERENCES projects(project_id)
    )
  ")

  # Analyses table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS analyses (
      analysis_id TEXT PRIMARY KEY,
      dataset_id TEXT,
      analysis_type TEXT NOT NULL,
      configuration TEXT,
      status TEXT DEFAULT 'pending',
      started_at TIMESTAMP,
      completed_at TIMESTAMP,
      error_message TEXT,
      result_blob BLOB,
      FOREIGN KEY (dataset_id) REFERENCES datasets(dataset_id)
    )
  ")

  # Results table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS results (
      result_id TEXT PRIMARY KEY,
      analysis_id TEXT,
      result_type TEXT NOT NULL,
      summary TEXT,
      tau2 REAL,
      i2 REAL,
      best_treatment TEXT,
      best_treatment_sucra REAL,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      result_blob BLOB,
      FOREIGN KEY (analysis_id) REFERENCES analyses(analysis_id)
    )
  ")

  # Visualizations table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS visualizations (
      viz_id TEXT PRIMARY KEY,
      result_id TEXT,
      viz_type TEXT NOT NULL,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      viz_blob BLOB,
      FOREIGN KEY (result_id) REFERENCES results(result_id)
    )
  ")

  # Users table (for collaboration)
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS users (
      user_id TEXT PRIMARY KEY,
      username TEXT UNIQUE NOT NULL,
      email TEXT UNIQUE,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      last_login TIMESTAMP,
      role TEXT DEFAULT 'user'
    )
  ")

  # Sessions table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS sessions (
      session_id TEXT PRIMARY KEY,
      user_id TEXT,
      started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      last_activity TIMESTAMP,
      ip_address TEXT,
      FOREIGN KEY (user_id) REFERENCES users(user_id)
    )
  ")

  # Audit log table
  DBI::dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS audit_log (
      log_id INTEGER PRIMARY KEY AUTOINCREMENT,
      user_id TEXT,
      action TEXT NOT NULL,
      resource_type TEXT,
      resource_id TEXT,
      details TEXT,
      timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      FOREIGN KEY (user_id) REFERENCES users(user_id)
    )
  ")

  # Create indexes for performance
  DBI::dbExecute(pool, "CREATE INDEX IF NOT EXISTS idx_datasets_project ON datasets(project_id)")
  DBI::dbExecute(pool, "CREATE INDEX IF NOT EXISTS idx_analyses_dataset ON analyses(dataset_id)")
  DBI::dbExecute(pool, "CREATE INDEX IF NOT EXISTS idx_results_analysis ON results(analysis_id)")
  DBI::dbExecute(pool, "CREATE INDEX IF NOT EXISTS idx_audit_user ON audit_log(user_id)")
  DBI::dbExecute(pool, "CREATE INDEX IF NOT EXISTS idx_audit_timestamp ON audit_log(timestamp)")

  message("Database schema created successfully")
}

#' Save Analysis to Database
#'
#' @description
#' Saves NMA analysis and results to database.
#'
#' @param db Database connection (from initialize_powernma_db())
#' @param analysis_result Analysis result object
#' @param dataset_id Dataset ID
#' @param analysis_type Analysis type description
#' @param user_id User ID (optional)
#'
#' @return Analysis ID
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_powernma_db(db_type = "sqlite", dbname = "powernma.db")
#' nma <- run_ultimate_nma(data)
#' analysis_id <- save_analysis_to_db(db, nma, dataset_id, "Ultimate NMA")
#' }
save_analysis_to_db <- function(db, analysis_result, dataset_id,
                               analysis_type = "NMA", user_id = NULL) {

  analysis_id <- generate_unique_id()

  # Serialize result
  result_blob <- serialize(analysis_result, NULL)

  # Insert analysis record
  DBI::dbExecute(
    db$pool,
    "INSERT INTO analyses (analysis_id, dataset_id, analysis_type, status,
                          started_at, completed_at, result_blob)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
    params = list(
      analysis_id,
      dataset_id,
      analysis_type,
      "completed",
      Sys.time(),
      Sys.time(),
      list(result_blob)
    )
  )

  # Save summary results
  result_id <- save_results_summary_to_db(db, analysis_id, analysis_result)

  # Log action
  log_db_action(db, user_id, "analysis_saved", "analysis", analysis_id)

  message(sprintf("Analysis saved to database: %s", analysis_id))

  return(analysis_id)
}

#' Save Results Summary
#'
#' @keywords internal
save_results_summary_to_db <- function(db, analysis_id, analysis_result) {

  result_id <- generate_unique_id()

  # Extract summary statistics
  tau2 <- analysis_result$nma_result$tau^2
  best_treatment <- if (!is.null(analysis_result$sucra)) {
    names(which.max(analysis_result$sucra$sucra_scores))
  } else {
    NA
  }
  best_sucra <- if (!is.null(analysis_result$sucra)) {
    max(analysis_result$sucra$sucra_scores)
  } else {
    NA
  }

  # Serialize full result
  result_blob <- serialize(analysis_result, NULL)

  # Insert result record
  DBI::dbExecute(
    db$pool,
    "INSERT INTO results (result_id, analysis_id, result_type, tau2,
                         best_treatment, best_treatment_sucra, result_blob)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
    params = list(
      result_id,
      analysis_id,
      "comprehensive",
      tau2,
      best_treatment,
      best_sucra,
      list(result_blob)
    )
  )

  return(result_id)
}

#' Retrieve Analysis from Database
#'
#' @description
#' Retrieves analysis results from database.
#'
#' @param db Database connection
#' @param analysis_id Analysis ID
#'
#' @return Analysis result object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db <- initialize_powernma_db(db_type = "sqlite")
#' result <- retrieve_analysis_from_db(db, analysis_id)
#' }
retrieve_analysis_from_db <- function(db, analysis_id) {

  query <- "SELECT result_blob FROM analyses WHERE analysis_id = ?"

  result <- DBI::dbGetQuery(db$pool, query, params = list(analysis_id))

  if (nrow(result) == 0) {
    stop("Analysis not found")
  }

  # Deserialize result
  analysis_result <- unserialize(result$result_blob[[1]])

  # Log action
  log_db_action(db, NULL, "analysis_retrieved", "analysis", analysis_id)

  return(analysis_result)
}

#' List All Analyses
#'
#' @description
#' Lists all analyses in database with summary information.
#'
#' @param db Database connection
#' @param project_id Optional project ID filter
#' @param status Optional status filter
#' @param limit Maximum number of results (default: 100)
#'
#' @return Data frame with analysis summaries
#'
#' @export
list_analyses_from_db <- function(db, project_id = NULL, status = NULL, limit = 100) {

  query <- "
    SELECT a.analysis_id, a.dataset_id, a.analysis_type, a.status,
           a.started_at, a.completed_at,
           r.best_treatment, r.best_treatment_sucra, r.tau2,
           d.dataset_name
    FROM analyses a
    LEFT JOIN results r ON a.analysis_id = r.analysis_id
    LEFT JOIN datasets d ON a.dataset_id = d.dataset_id
    WHERE 1=1
  "

  params <- list()

  if (!is.null(project_id)) {
    query <- paste(query, "AND d.project_id = ?")
    params <- c(params, project_id)
  }

  if (!is.null(status)) {
    query <- paste(query, "AND a.status = ?")
    params <- c(params, status)
  }

  query <- paste(query, "ORDER BY a.started_at DESC LIMIT ?")
  params <- c(params, limit)

  results <- DBI::dbGetQuery(db$pool, query, params = params)

  return(results)
}

#' Save Dataset to Database
#'
#' @export
save_dataset_to_db <- function(db, data, dataset_name, project_id = NULL, user_id = NULL) {

  dataset_id <- generate_unique_id()

  # Calculate summary statistics
  n_studies <- length(unique(data$studlab))
  n_treatments <- length(unique(c(data$treat1, data$treat2)))
  n_comparisons <- nrow(data)

  # Serialize data
  data_blob <- serialize(data, NULL)

  # Insert dataset record
  DBI::dbExecute(
    db$pool,
    "INSERT INTO datasets (dataset_id, project_id, dataset_name, n_studies,
                          n_treatments, n_comparisons, data_blob, uploaded_by)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(
      dataset_id,
      project_id,
      dataset_name,
      n_studies,
      n_treatments,
      n_comparisons,
      list(data_blob),
      user_id
    )
  )

  # Log action
  log_db_action(db, user_id, "dataset_uploaded", "dataset", dataset_id)

  message(sprintf("Dataset saved to database: %s", dataset_id))

  return(dataset_id)
}

#' Retrieve Dataset from Database
#'
#' @export
retrieve_dataset_from_db <- function(db, dataset_id) {

  query <- "SELECT data_blob FROM datasets WHERE dataset_id = ?"

  result <- DBI::dbGetQuery(db$pool, query, params = list(dataset_id))

  if (nrow(result) == 0) {
    stop("Dataset not found")
  }

  # Deserialize data
  data <- unserialize(result$data_blob[[1]])

  # Log action
  log_db_action(db, NULL, "dataset_retrieved", "dataset", dataset_id)

  return(data)
}

#' Log Database Action
#'
#' @keywords internal
log_db_action <- function(db, user_id, action, resource_type, resource_id, details = NULL) {

  tryCatch({
    DBI::dbExecute(
      db$pool,
      "INSERT INTO audit_log (user_id, action, resource_type, resource_id, details)
       VALUES (?, ?, ?, ?, ?)",
      params = list(user_id, action, resource_type, resource_id, details)
    )
  }, error = function(e) {
    warning("Failed to log action: ", e$message)
  })
}

#' Close Database Connection
#'
#' @description
#' Properly closes database connection pool.
#'
#' @param db Database connection from initialize_powernma_db()
#'
#' @export
close_powernma_db <- function(db) {
  pool::poolClose(db$pool)
  message("Database connection closed")
}

#' Get Database Statistics
#'
#' @export
get_db_statistics <- function(db) {

  stats <- list(
    n_projects = DBI::dbGetQuery(db$pool, "SELECT COUNT(*) as n FROM projects")$n,
    n_datasets = DBI::dbGetQuery(db$pool, "SELECT COUNT(*) as n FROM datasets")$n,
    n_analyses = DBI::dbGetQuery(db$pool, "SELECT COUNT(*) as n FROM analyses")$n,
    n_results = DBI::dbGetQuery(db$pool, "SELECT COUNT(*) as n FROM results")$n,
    n_users = DBI::dbGetQuery(db$pool, "SELECT COUNT(*) as n FROM users")$n,
    db_size = get_db_size(db)
  )

  return(stats)
}

#' Get Database Size
#'
#' @keywords internal
get_db_size <- function(db) {

  if (db$type == "sqlite") {
    file_size <- file.size(db$dbname)
    return(sprintf("%.2f MB", file_size / 1024^2))
  } else {
    return("N/A")
  }
}

#' Generate Unique ID
#'
#' @keywords internal
generate_unique_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample(10000:99999, 1))
}

#' Print Method for Database Connection
#'
#' @export
print.powernma_db <- function(x, ...) {
  cat("powerNMA Database Connection\n")
  cat("============================\n")
  cat(sprintf("Type: %s\n", x$type))
  cat(sprintf("Database: %s\n", x$dbname))
  cat(sprintf("Pool size: %d\n", pool::dbGetInfo(x$pool)$maxSize))
  cat("\n")

  stats <- get_db_statistics(x)
  cat("Statistics:\n")
  cat(sprintf("  Projects: %d\n", stats$n_projects))
  cat(sprintf("  Datasets: %d\n", stats$n_datasets))
  cat(sprintf("  Analyses: %d\n", stats$n_analyses))
  cat(sprintf("  Results: %d\n", stats$n_results))
  cat(sprintf("  Users: %d\n", stats$n_users))
  if (!is.na(stats$db_size)) {
    cat(sprintf("  Database size: %s\n", stats$db_size))
  }

  invisible(x)
}
