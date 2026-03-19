options(repos = 'https://cloud.r-project.org')

# Register additional candidates (CRAN/GitHub/Zenodo)
if (file.exists('data-raw/find_more.R')) source('data-raw/find_more.R')

# Download from remote (GitHub/Zenodo) and embed
source('data-raw/download_remote.R')

# Validate all embedded datasets
source('data-raw/validate_embedded.R')

cat('Done. Wrote data-raw/embedded_validation.csv\n')
