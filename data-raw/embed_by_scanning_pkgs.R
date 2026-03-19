pkgs <- c('netmeta','gemtc','multinma','MBNMAdose','pcnetmeta','bnma','BUGSnet','nmaINLA')

sanitize_id <- function(pkg, item) {
  base <- paste0(pkg, '__', item)
  gsub('[^A-Za-z0-9_]+','_', base)
}

ensure_src <- function(pkg, item) {
  e <- new.env(parent = emptyenv())
  ok <- try(suppressWarnings(utils::data(list=item, package=pkg, envir=e)), silent=TRUE)
  if (!inherits(ok, 'try-error') && exists(item, envir=e, inherits=FALSE)) return(get(item, envir=e, inherits=FALSE))
  item2 <- sub("\\s*\\(.*$", "", item)
  ok <- try(suppressWarnings(utils::data(list=item2, package=pkg, envir=e)), silent=TRUE)
  if (!inherits(ok, 'try-error') && exists(item2, envir=e, inherits=FALSE)) return(get(item2, envir=e, inherits=FALSE))
  NULL
}

source('R/normalize_more.R', local = FALSE)
source('R/upstream.R', local = FALSE)
as_nma <- as_nma_tables
write_local <- write_embedded_local

dir.create(file.path('inst','extdata'), recursive = TRUE, showWarnings = FALSE)

embedded <- character(); failed <- character()
cat('Scanning packages: ', paste(pkgs, collapse=', '), "\n", sep='')
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) next
  ddf <- try(utils::data(package = pkg)$results, silent = TRUE)
  if (inherits(ddf, 'try-error') || is.null(ddf) || !NROW(ddf)) next
  items <- unique(as.character(as.data.frame(ddf, stringsAsFactors = FALSE)[['Item']]))
  cat('Package ', pkg, ': ', length(items), ' items\n', sep='')
  for (it in items) {
    id <- sanitize_id(pkg, it)
    # Skip if already embedded
    if (file.exists(file.path('inst','extdata', paste0(id, '_studies.csv'))) ||
        file.exists(file.path('inst','extdata', paste0(id, '_contrast.csv')))) next
    obj <- ensure_src(pkg, it)
    if (is.null(obj)) next
    tables <- NULL
    if (is.data.frame(obj)) {
      # Try pairwise -> arm, then contrast
      tables <- tryCatch({ out <- to_arm_from_pairwise_df(obj); if (!is.null(out)) list(studies = out$arm, nodes = out$nodes) else NULL }, error=function(e) NULL)
      if (is.null(tables)) tables <- tryCatch({ out <- to_contrast_from_effects_df(obj); if (!is.null(out)) list(contrast = out$contrast, nodes = out$nodes) else NULL }, error=function(e) NULL)
    }
    if (is.null(tables)) tables <- tryCatch(as_nma(obj), error = function(e) NULL)
    if (is.null(tables)) { failed <- c(failed, id); next }
    ok <- tryCatch(write_local(id, tables), error = function(e) FALSE)
    if (isTRUE(ok)) embedded <- c(embedded, id) else failed <- c(failed, id)
  }
}
cat('Embedded by scanning: ', length(embedded), ' datasets\n', sep='')
