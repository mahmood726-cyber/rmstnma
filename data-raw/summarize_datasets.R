# Summarize and validate raw datasets without saving .rda

library(readr)

src_dir <- file.path("data-raw", "raw")
if (!dir.exists(src_dir)) {
  stop("No directory ", src_dir)
}

files <- list.files(src_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE, ignore.case = TRUE)
out <- data.frame(
  file = character(),
  ok = logical(),
  n_rows = integer(),
  n_studies = integer(),
  n_treatments = integer(),
  time_min = numeric(),
  time_max = numeric(),
  issues = character(),
  stringsAsFactors = FALSE
)

for (f in files) {
  ext <- tools::file_ext(f)
  df <- switch(tolower(ext),
               csv = suppressMessages(readr::read_csv(f, show_col_types = FALSE)),
               tsv = suppressMessages(readr::read_tsv(f, show_col_types = FALSE)))
  df <- standardize_nma_columns(df)
  v <- validate_nma_dataset(df, strict = TRUE, quiet = TRUE)
  smry <- v$summary
  out <- rbind(out, data.frame(
    file = basename(f),
    ok = v$ok,
    n_rows = if (!is.null(smry)) smry$n_rows else NA_integer_,
    n_studies = if (!is.null(smry)) smry$n_studies else NA_integer_,
    n_treatments = if (!is.null(smry)) smry$n_treatments else NA_integer_,
    time_min = if (!is.null(smry)) smry$time_range[1] else NA_real_,
    time_max = if (!is.null(smry)) smry$time_range[2] else NA_real_,
    issues = if (length(v$issues) > 0) paste(v$issues, collapse = "; ") else ""
  ))
}

out_path <- file.path("data-raw", "validation_report.csv")
utils::write.csv(out, out_path, row.names = FALSE)
message("Wrote ", out_path)

