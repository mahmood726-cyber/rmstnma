## Meta–Meta Analysis Example (fixed-effect illustration)
##
## This script aggregates all standardized contrast datasets shipped in
## inst/extdata (those ending with _contrast.csv). For each dataset, it
## performs an inverse-variance fixed-effect meta-analysis to obtain a dataset-
## level pooled effect (TE, seTE). It then performs a second-stage fixed-effect
## meta-analysis across datasets to obtain an overall pooled estimate.
##
## Tuning/filters to improve stability:
## - Coerce TE and seTE to numeric (defensive against "NA" strings)
## - Per-dataset, drop rows with non-finite values and extreme seTE
##   (default: remove seTE above 95th percentile within-dataset, and absolute seTE > 5)
## - Drop dataset-level results with pooled seTE > 5
##
## Note: This is an illustrative example. You may replace with a random-effects
## method (e.g., DerSimonian–Laird) if desired; or stratify by outcome/measure.

suppressWarnings(suppressMessages({
  library(utils)
}))

ext <- file.path("inst", "extdata")
files <- list.files(ext, pattern = "_contrast\\.csv$", full.names = TRUE)

summ <- data.frame(
  id = character(), n_effects = integer(), TE = numeric(), seTE = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {
  id <- sub("_contrast\\.csv$", "", basename(f))
  df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) next
  needed <- c("TE", "seTE")
  if (!all(needed %in% names(df))) next
  # Coerce to numeric in case columns were read as character
  suppressWarnings({
    df$TE <- as.numeric(df$TE)
    df$seTE <- as.numeric(df$seTE)
  })
  # Drop non-finite and non-positive seTE
  dd <- subset(df, is.finite(TE) & is.finite(seTE) & seTE > 0)
  if (!nrow(dd)) next
  # Drop extreme seTE within dataset (95th percentile) and absolute cap (5)
  q95 <- stats::quantile(dd$seTE, probs = 0.95, na.rm = TRUE, names = FALSE)
  cap <- max(5, q95)
  dd <- dd[dd$seTE <= cap, , drop = FALSE]
  if (!nrow(dd)) next
  w <- 1 / (dd$seTE^2)
  TEp <- sum(w * dd$TE) / sum(w)
  seTEp <- sqrt(1 / sum(w))
  # Skip dataset-level if pooled seTE is still extreme
  if (is.finite(seTEp) && seTEp <= 5) {
    summ <- rbind(summ, data.frame(id = id, n_effects = nrow(dd), TE = TEp, seTE = seTEp))
  }
}

if (!nrow(summ)) {
  warning("No contrast datasets with usable TE/seTE found after filtering.")
} else {
  # Second-stage pooled estimate across datasets
  w2 <- 1 / (summ$seTE^2)
  overall_TE <- sum(w2 * summ$TE) / sum(w2)
  overall_seTE <- sqrt(1 / sum(w2))
  overall <- data.frame(id = "OVERALL", n_effects = sum(summ$n_effects), TE = overall_TE, seTE = overall_seTE)
  out <- rbind(summ, overall)
  dir.create(ext, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, file.path(ext, "meta_meta_results.csv"), row.names = FALSE)
  message("Wrote ", file.path(ext, "meta_meta_results.csv"))
}
