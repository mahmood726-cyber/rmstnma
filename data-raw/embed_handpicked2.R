source('R/normalize_more.R', local = FALSE)
source('R/upstream.R', local = FALSE)

ids <- c(
  'bnma__thrombolytic','bnma__blocker','bnma__parkinsons','bnma__parkinsons_contrast','bnma__smoking','bnma__statins',
  'multinma__social_anxiety','multinma__hta_psoriasis','multinma__plaque_psoriasis_agd','multinma__plaque_psoriasis_ipd',
  'BUGSnet__diabetes','BUGSnet__diabetes.sim','BUGSnet__parkinsons',
  'nmaINLA__Diabetesdat','nmaINLA__Dietaryfatdat','nmaINLA__Parkinsondat'
)

ok <- 0L; fail <- character()
for (id in ids) {
  if (file.exists(file.path('inst','extdata', paste0(id, '_studies.csv'))) ||
      file.exists(file.path('inst','extdata', paste0(id, '_contrast.csv')))) next
  x <- tryCatch(get_upstream_raw(id), error=function(e) NULL)
  if (is.null(x)) { fail <- c(fail, id); next }
  tabs <- tryCatch(as_nma_tables(x), error=function(e) NULL)
  if (is.null(tabs)) { fail <- c(fail, id); next }
  write_embedded_local(id, tabs, dir = file.path('inst','extdata'))
  ok <- ok + 1L
}
cat('Handpicked2 embedded: ', ok, '; failed: ', length(fail), '\n', sep='')
if (length(fail)) write.csv(data.frame(id=fail), file='data-raw/handpicked2_failed.csv', row.names=FALSE)

