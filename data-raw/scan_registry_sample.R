reg <- read.csv('inst/registry/registry.csv', stringsAsFactors = FALSE)
pkgs <- c('netmeta','gemtc','multinma','MBNMAdose','pcnetmeta','bnma')
rows <- subset(reg, package %in% pkgs)[, c('id','package','object')]
out <- data.frame(id=character(), package=character(), object=character(), class=character(), is_df=logical(), nrows=integer(), ncols=integer(), stringsAsFactors=FALSE)
for (i in seq_len(nrow(rows))) {
  pkg <- rows$package[i]; obj <- rows$object[i]
  e <- new.env(parent = emptyenv())
  ok <- try(suppressWarnings(utils::data(list=obj, package=pkg, envir=e)), silent=TRUE)
  if (inherits(ok, 'try-error') || !exists(obj, envir=e, inherits=FALSE)) {
    # try sanitized name
    obj2 <- sub("\\s*\\(.*$", "", obj)
    ok <- try(suppressWarnings(utils::data(list=obj2, package=pkg, envir=e)), silent=TRUE)
    nm <- if (exists(obj2, envir=e, inherits=FALSE)) obj2 else obj
  } else nm <- obj
  if (exists(nm, envir=e, inherits=FALSE)) {
    x <- get(nm, envir=e, inherits=FALSE)
    cl <- paste(class(x), collapse=';')
    isdf <- is.data.frame(x)
    nr <- if (isdf) nrow(x) else NA_integer_
    nc <- if (isdf) ncol(x) else NA_integer_
    out <- rbind(out, data.frame(id=rows$id[i], package=pkg, object=nm, class=cl, is_df=isdf, nrows=nr, ncols=nc, stringsAsFactors=FALSE))
  } else {
    out <- rbind(out, data.frame(id=rows$id[i], package=pkg, object=obj, class='NOT_FOUND', is_df=NA, nrows=NA, ncols=NA, stringsAsFactors=FALSE))
  }
}
write.csv(out, file='data-raw/scan_registry_sample.csv', row.names=FALSE)
cat('Wrote data-raw/scan_registry_sample.csv with ', nrow(out), ' rows.\n', sep='')
