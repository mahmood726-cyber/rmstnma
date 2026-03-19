args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: Rscript data-raw/peek_pairwise70.R <path-to-.rda>')
p <- args[[1]]
e <- new.env(parent = emptyenv())
load(p, envir = e)
nms <- ls(e)
cat('Objects in', p, ':', paste(nms, collapse = ', '), '\n')
for (nm in nms) {
  x <- get(nm, envir = e)
  cat('---', nm, 'class =', paste(class(x), collapse = '/'), '\n')
  if (is.data.frame(x)) {
    cat('Columns:', paste(names(x), collapse = ', '), '\n')
    print(utils::head(x, 3))
  } else {
    str(x)
  }
}
