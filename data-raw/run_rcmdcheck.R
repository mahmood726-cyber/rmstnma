options(repos = 'https://cloud.r-project.org')
if (!requireNamespace('rcmdcheck', quietly = TRUE)) install.packages('rcmdcheck')
res <- rcmdcheck::rcmdcheck(args = c('--no-manual','--ignore-vignettes'), error_on = 'warning')
print(res)
