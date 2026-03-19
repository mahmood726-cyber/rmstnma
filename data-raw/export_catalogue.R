in_path <- 'data-raw/embedded_validation.csv'
out_dir <- file.path('inst','extdata')
out_path <- file.path(out_dir, 'catalogue.csv')

if (!file.exists(in_path)) stop('Validation file not found: ', in_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
x <- utils::read.csv(in_path, stringsAsFactors = FALSE)
utils::write.csv(x, out_path, row.names = FALSE)
message('Exported: ', out_path)
