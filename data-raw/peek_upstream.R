pkgs <- c('multinma','netmeta','BUGSnet','MBNMAdose','pcnetmeta','nmaINLA','bnma','gemtc')
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) next
  cat('--- PACKAGE:', pkg, '---\n')
  dd <- try(utils::data(package = pkg)$results, silent = TRUE)
  if (inherits(dd, 'try-error') || is.null(dd)) next
  items <- unique(as.data.frame(dd, stringsAsFactors = FALSE)[['Item']])
  for (nm in head(items, 8)) {
    e <- new.env(parent = emptyenv()); suppressWarnings(utils::data(list = nm, package = pkg, envir = e))
    if (!exists(nm, envir = e, inherits = FALSE)) next
    obj <- get(nm, envir = e, inherits = FALSE)
    cat(pkg, '::', nm, '\n', sep='')
    str(obj)
    cat('\n')
  }
}
