# Build standardized NMA datasets from CSV/TSV files
#
# Place raw files in data-raw/raw/ (one file per dataset). Supported formats:
# - .csv (comma separated)
# - .tsv (tab separated)
#
# This script will attempt to standardize column names, validate the content,
# and save each dataset into data/ as an .rda object with the same base name as
# the source file. It will stop on validation errors.

library(readr)
library(dplyr)

src_dir <- file.path("data-raw", "raw")
out_dir <- "data"
if (!dir.exists(src_dir)) dir.create(src_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Scanning ", src_dir)
files <- list.files(src_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0L) {
  message("No raw files found. Place CSV/TSV files in ", src_dir)
}

for (f in files) {
  message("\nProcessing ", basename(f))
  ext <- tools::file_ext(f)
  df <- switch(tolower(ext),
               csv = suppressMessages(readr::read_csv(f, show_col_types = FALSE)),
               tsv = suppressMessages(readr::read_tsv(f, show_col_types = FALSE)),
               stop("Unsupported extension: ", ext))

  df <- standardize_nma_columns(df)
  v <- validate_nma_dataset(df, strict = FALSE, quiet = TRUE)
  if (!v$ok) {
    msg <- paste(v$issues, collapse = "\n - ")
    stop("Validation failed for ", basename(f), ":\n - ", msg)
  }

  # Save to data/ with object name = base filename
  obj_name <- make.names(sub("\\.[^.]+$", "", basename(f)))
  assign(obj_name, df)
  save(list = obj_name, file = file.path(out_dir, paste0(obj_name, ".rda")), compress = "bzip2")
  message("Saved: ", file.path(out_dir, paste0(obj_name, ".rda")))
}

