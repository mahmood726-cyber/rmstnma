suppressPackageStartupMessages({
  library(jsonlite)
  have_httr <- requireNamespace('httr', quietly = TRUE)
  have_readxl <- requireNamespace('readxl', quietly = TRUE)
})

source_if <- function(f) if (file.exists(f)) source(f, local = FALSE)
source_if('R/normalize_more.R')
source_if('R/upstream.R')
source_if('R/registry_utils.R')
source_if('data-raw/remote_recipes/apply_recipes.R')

here_ext <- function() {
  p <- file.path('inst','extdata'); dir.create(p, recursive = TRUE, showWarnings = FALSE); p
}

registry_path <- function() {
  p <- system.file('registry','registry.csv', package = 'nmadatasets', mustWork = FALSE)
  if (!nzchar(p)) p <- file.path('inst','registry','registry.csv')
  p
}

regfile <- registry_path()
if (!file.exists(regfile)) stop('No registry found at ', regfile)
reg <- utils::read.csv(regfile, stringsAsFactors = FALSE)

is_remote <- function(row) {
  sr <- tolower(as.character(row$source))
  su <- as.character(row$source_url)
  (nzchar(sr) && sr == 'remote' && nzchar(su))
}

rows <- reg[apply(reg, 1, function(r){ is_remote(as.list(r)) }), , drop = FALSE]
if (!nrow(rows)) {
  message('No remote-only entries found in registry.')
  quit(save='no')
}

zenodo_files <- function(url) {
  rid <- NA_character_
  m <- regexec('zenodo\\.org/(record|records)/([0-9]+)', url)
  mr <- regmatches(url, m)[[1]]
  if (length(mr) >= 3) rid <- mr[3]
  if (is.na(rid)) return(character())
  api <- paste0('https://zenodo.org/api/records/', rid)
  js <- tryCatch(jsonlite::fromJSON(api, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(js) || is.null(js$files)) return(character())
  files <- js$files
  urls <- files$links$download
  names(urls) <- files$filename
  keep <- grepl('\\.(csv|tsv|xlsx|xls|rds|rdata)$', tolower(names(urls)))
  urls[keep]
}

parse_github_repo <- function(url) {
  m <- regexec('github\\.com/([^/]+)/([^/]+)', url)
  mr <- regmatches(url, m)[[1]]
  if (length(mr) >= 3) list(owner = mr[2], repo = sub('\\.git$', '', mr[3])) else NULL
}

gh_get <- function(url, token = Sys.getenv('GITHUB_TOKEN', ''), accept = 'application/vnd.github+json') {
  if (!have_httr) return(NULL)
  headers <- httr::add_headers(Accept = accept)
  if (nzchar(token)) headers <- httr::add_headers(Authorization = paste('Bearer', token), Accept = accept)
  res <- tryCatch(httr::GET(url, headers, httr::user_agent('nmadatasets-remote/0.1')),
                  error = function(e) NULL)
  if (is.null(res) || httr::http_error(res)) return(NULL)
  txt <- httr::content(res, as = 'text', encoding = 'UTF-8')
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

github_files <- function(url) {
  info <- parse_github_repo(url)
  if (is.null(info)) return(character())
  api_repo <- paste0('https://api.github.com/repos/', info$owner, '/', info$repo)
  js <- gh_get(api_repo)
  if (is.null(js)) return(character())
  branch <- if (!is.null(js$default_branch)) js$default_branch else 'main'
  api_tree <- paste0('https://api.github.com/repos/', info$owner, '/', info$repo, '/git/trees/', branch, '?recursive=1')
  tree <- gh_get(api_tree)
  if (is.null(tree) || is.null(tree$tree)) return(character())
  paths <- vapply(tree$tree$path, as.character, '')
  keep <- grepl('\\.(csv|tsv|xlsx|xls|rds|rdata)$', tolower(paths))
  paths <- paths[keep]
  if (!length(paths)) return(character())
  ord <- order(!(grepl('(^|/)data|(^|/)inst/|(^|/)extdata', tolower(paths))))
  paths <- paths[ord]
  raw_base <- paste0('https://raw.githubusercontent.com/', info$owner, '/', info$repo, '/', branch, '/')
  urls <- paste0(raw_base, paths)
  names(urls) <- basename(paths)
  urls
}

download_one <- function(id, url) {
  ext <- tolower(tools::file_ext(url))
  dest <- tempfile(paste0('nmadl_', id, '_'), fileext = paste0('.', ext))
  ok <- tryCatch(utils::download.file(url, destfile = dest, mode = 'wb', quiet = TRUE), error = function(e) 1L)
  if (is.integer(ok) && ok != 0L) return(NULL)
  dest
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c('csv','tsv')) return(tryCatch(utils::read.csv(path, check.names = FALSE), error=function(e) NULL))
  if (ext %in% c('xlsx','xls')) {
    if (!have_readxl) return(NULL)
    return(tryCatch(readxl::read_excel(path), error=function(e) NULL))
  }
  if (ext == 'rds') return(tryCatch(readRDS(path), error=function(e) NULL))
  if (ext == 'rdata') {
    e <- new.env(parent = emptyenv());
    tryCatch({ load(path, envir = e); objs <- ls(envir = e); for (nm in objs) { x <- get(nm, envir = e); if (is.data.frame(x)) return(x) } ; NULL }, error=function(e) NULL)
  } else NULL
}

standardize_and_write <- function(id, path) {
  obj <- read_any(path)
  if (is.null(obj)) return(FALSE)
  # Try per-repo recipes first
  if (exists('apply_remote_recipe', inherits = TRUE)) {
    rec <- tryCatch(apply_remote_recipe(id, obj), error=function(e) NULL)
    if (!is.null(rec)) obj <- rec
  }
  # try data.frame first, then list via heuristics
  tabs <- tryCatch(as_nma_tables(obj), error = function(e) NULL)
  if (is.null(tabs) && is.list(obj)) tabs <- tryCatch(as_nma_tables(obj), error=function(e) NULL)
  if (is.null(tabs)) return(FALSE)
  write_embedded_local(id, tabs, dir = here_ext())
  TRUE
}

embedded <- character(); failed <- character()

for (i in seq_len(nrow(rows))) {
  id <- rows$id[i]
  url <- rows$source_url[i]
  if (file.exists(file.path(here_ext(), paste0(id, '_studies.csv'))) ||
      file.exists(file.path(here_ext(), paste0(id, '_contrast.csv')))) next
  urls <- character()
  if (grepl('zenodo\\.org', url, ignore.case = TRUE)) {
    urls <- zenodo_files(url)
  } else if (grepl('github\\.com', url, ignore.case = TRUE)) {
    urls <- github_files(url)
  } else if (grepl('\\.(csv|tsv|xlsx|xls|rds|rdata)$', tolower(url))) {
    urls <- url
  }
  if (!length(urls)) next
  ok <- FALSE
  for (u in urls) {
    tmp <- download_one(id, u)
    if (is.null(tmp)) next
    if (standardize_and_write(id, tmp)) { ok <- TRUE; break }
  }
  if (ok) embedded <- c(embedded, id) else failed <- c(failed, id)
}

message('Remote embedded: ', length(embedded))
if (length(failed)) utils::write.csv(data.frame(id = failed), file = file.path('data-raw','download_remote_failed.csv'), row.names = FALSE)
