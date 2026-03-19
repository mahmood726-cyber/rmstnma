repos <- 'https://cloud.r-project.org'
reg <- read.csv('inst/registry/registry.csv', stringsAsFactors = FALSE)
ext <- if (dir.exists('inst/extdata')) 'inst/extdata' else 'inst/extdata'
have <- unique(sub('(_nodes|_studies|_contrast)?\\.csv$','', list.files(ext, pattern='\\.csv$', full.names=FALSE)))
todo <- subset(reg, !(id %in% have) & !is.na(package) & !is.na(object))
pkgs <- unique(todo$package)
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  install.packages(miss, repos = repos)
  cat('Installed:', paste(miss, collapse = ', '), '\n')
} else {
  cat('No installs needed.\n')
}

