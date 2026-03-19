
test_that("rmst_network creates valid object", {

  # Use example data
  data(example_network)

  # Create network
  net <- rmst_network(
    data = example_network,
    study = "study_id",
    trt = "treatment",
    time = "time",
    surv = "survival"
  )

  # Test structure
  expect_s3_class(net, "rmst_network")
  expect_equal(net$n_studies, 3)
  expect_equal(net$n_treatments, 3)
})

test_that("diagnose_tau works", {

  data(example_network)

  net <- rmst_network(
    data = example_network,
    study = "study_id",
    trt = "treatment",
    time = "time",
    surv = "survival"
  )

  diag <- diagnose_tau(net)

  expect_s3_class(diag, "tau_diagnostics")
  expect_true(length(diag$suggested_tau) > 0)
})

test_that("validate_nma_dataset accepts example data", {
  data(example_network)
  v <- validate_nma_dataset(example_network, strict = TRUE, quiet = TRUE)
  expect_true(v$ok)
  expect_true(v$summary$n_studies >= 1)
})

test_that("standardize_nma_columns maps common synonyms", {
  df <- data.frame(
    study = c("A", "A"),
    trt = c("X", "Y"),
    t = c(0, 1),
    S = c(1, 0.9),
    stringsAsFactors = FALSE
  )
  out <- standardize_nma_columns(df)
  expect_true(all(c("study_id", "treatment", "time", "survival") %in% names(out)))
  v <- validate_nma_dataset(out, strict = FALSE, quiet = TRUE)
  expect_true(v$ok)
})
