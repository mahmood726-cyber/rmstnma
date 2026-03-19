# Hand-pick some upstream datasets and embed them with robust pairwise mapping

source('data-raw/recipes/02_embed_known_pkgs.R')

embed_one <- function(pkg, item, id = NULL){
  if (!requireNamespace(pkg, quietly = TRUE)) return(FALSE)
  e <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(list=item, package=pkg, envir=e))
  if (!exists(item, envir=e, inherits=FALSE)) return(FALSE)
  x <- get(item, envir=e, inherits=FALSE)
  if (!is.data.frame(x)) return(FALSE)
  out <- to_arm_from_pairwise_df(x)
  if (is.null(out)) return(FALSE)
  if (is.null(id)) id <- paste0(pkg, '__', item)
  write_tables(id, arm = out$arm, nodes = out$nodes)
  TRUE
}

dir.create('inst/extdata', recursive = TRUE, showWarnings = FALSE)

targets <- list(
  list(pkg='netmeta', item='smokingcessation'),
  list(pkg='netmeta', item='dietaryfat'),
  list(pkg='netmeta', item='Franchini2012'),
  list(pkg='netmeta', item='Woods2010'),
  list(pkg='netmeta', item='parkinson')
)

ok <- 0L
for (t in targets) {
  id <- paste0(t$pkg, '__', t$item)
  if (file.exists(file.path('inst','extdata', paste0(id, '_studies.csv')))) next
  res <- embed_one(t$pkg, t$item, id)
  if (isTRUE(res)) ok <- ok + 1L
}
cat('Handpicked embedded: ', ok, '\n', sep='')

