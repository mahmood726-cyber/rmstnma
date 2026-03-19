cat_file <- 'inst/extdata/catalogue.csv'
reg_file <- 'inst/registry/registry.csv'
out_file <- 'inst/metadata/datasets.md'

dir.create('inst/metadata', showWarnings = FALSE, recursive = TRUE)

if (!file.exists(cat_file)) stop('Catalogue not found: ', cat_file)
x <- utils::read.csv(cat_file, stringsAsFactors = FALSE)
if (file.exists(reg_file)) r <- utils::read.csv(reg_file, stringsAsFactors = FALSE) else r <- data.frame()

# Merge registry info if present
if (nrow(r)) {
  keep <- c('id','source_url','license','citation')
  r2 <- r[, intersect(keep, names(r)), drop = FALSE]
  x2 <- merge(x, r2, by = 'id', all.x = TRUE)
} else {
  x2 <- x
}

header <- c(
  '# Datasets Overview',
  '',
  'This page summarizes the curated datasets shipped with the package, including their',
  'structure, counts, and available provenance where known.',
  '',
  '| id | arm | contrast | n_rows_arm | n_treatments | n_studies | moderators | source_url | license | citation |',
  '|---|:---:|:---:|---:|---:|---:|---|---|---|---|'
)

lines <- header
for (i in seq_len(nrow(x2))) {
  row <- x2[i, ]
  arm <- ifelse(isTRUE(row$has_arm), 'Y', '')
  con <- ifelse(isTRUE(row$has_contrast), 'Y', '')
  mods <- if (!is.na(row$covariates) && nzchar(row$covariates)) gsub('\n', ' ', row$covariates) else ''
  src <- if ('source_url' %in% names(row)) as.character(row$source_url) else ''
  lic <- if ('license' %in% names(row)) as.character(row$license) else ''
  cit <- if ('citation' %in% names(row)) as.character(row$citation) else ''
  lines <- c(lines, sprintf('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |',
                            row$id, arm, con, row$n_rows_arm, row$n_treatments, row$n_studies, mods, src, lic, cit))
}

writeLines(lines, out_file)
message('Wrote ', out_file)
