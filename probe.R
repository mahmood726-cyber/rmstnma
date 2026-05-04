# Baseline probe for rmstnma (Restricted Mean Survival Time NMA package).
#
# Loads the canonical example_network and emits structural / first-pass
# numerical signals from rmst_network(). All signals are deterministic
# given a fixed dataset.
#
# Run: Rscript probe.R

suppressPackageStartupMessages({
  if (!"rmstnma" %in% loadedNamespaces()) {
    pkgload::load_all(rprojroot::find_package_root_file(), quiet = TRUE)
  }
})

data(example_network)
net <- rmst_network(
  data = example_network,
  study = "study_id",
  trt = "treatment",
  time = "time",
  surv = "survival"
)

out <- list(
  n_studies = net$n_studies,
  n_treatments = net$n_treatments,
  n_rows = nrow(example_network),
  time_min = round(min(example_network$time), 6),
  time_max = round(max(example_network$time), 6),
  surv_min = round(min(example_network$survival), 6),
  surv_max = round(max(example_network$survival), 6),
  n_treatment_levels = length(unique(example_network$treatment))
)

cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 6))
cat("\n")
