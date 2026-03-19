options(repos = 'https://cloud.r-project.org')
ensure <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss)
}
ensure(c('devtools'))

message('Running devtools::document() ...')
devtools::document()

message('Running devtools::check() ...')
res <- devtools::check(document = FALSE, cran = FALSE)
print(res)

if (length(res$errors)) stop('R CMD check reported errors')
if (length(res$warnings)) warning('R CMD check reported warnings')
message('devtools::check completed.')
