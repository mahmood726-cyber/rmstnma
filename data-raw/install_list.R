pkgs <- c('netmeta','gemtc','multinma','MBNMAdose','pcnetmeta','bnma','httr','readxl')
ok <- character(); fail <- character()
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat('Installing ', p, ' ...\n', sep='')
    tryCatch(install.packages(p, repos='https://cloud.r-project.org'), error = function(e) fail <<- c(fail, p))
  }
  if (requireNamespace(p, quietly = TRUE)) ok <- c(ok, p) else fail <- c(fail, p)
}
cat('Installed/available: ', paste(ok, collapse=', '), '\n', sep='')
if (length(fail)) cat('Unavailable (skipped): ', paste(unique(fail), collapse=', '), '\n', sep='')
