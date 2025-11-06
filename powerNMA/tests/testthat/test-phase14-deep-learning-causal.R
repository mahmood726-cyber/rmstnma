# Tests for Phase 14: Deep Learning, Causal Inference & Advanced Methods

context("Phase 14: Graph Neural Networks and Causal Inference")

# Test 1: Graph Neural Networks
test_that("GNN training works", {
  skip_on_cran()
  
  data <- simulate_nma_data(n_studies = 20, n_treatments = 6)
  
  # Create node features
  treatments <- unique(c(data$treat1, data$treat2))
  node_features <- matrix(rnorm(length(treatments) * 10), length(treatments), 10)
  rownames(node_features) <- treatments
  
  gnn_model <- train_gnn_nma(
    nma_data = data,
    node_features = node_features,
    architecture = "GCN",
    n_layers = 2,
    hidden_dim = 32,
    n_epochs = 10,
    task = "link_prediction"
  )
  
  expect_s3_class(gnn_model, "gnn_nma")
  expect_equal(gnn_model$architecture, "GCN")
})

# Test 2: GNN Embeddings
test_that("Treatment embeddings extraction works", {
  data <- simulate_nma_data(n_studies = 15, n_treatments = 5)
  treatments <- unique(c(data$treat1, data$treat2))
  node_features <- matrix(rnorm(length(treatments) * 8), length(treatments), 8)
  rownames(node_features) <- treatments
  
  gnn_model <- train_gnn_nma(
    nma_data = data,
    node_features = node_features,
    architecture = "GAT",
    n_layers = 2,
    n_epochs = 5
  )
  
  embeddings <- extract_treatment_embeddings(gnn_model)
  
  expect_true(is.matrix(embeddings))
  expect_equal(nrow(embeddings), length(treatments))
})

# Test 3: G-Formula
test_that("G-formula causal inference works", {
  data <- simulate_nma_data(n_studies = 20)
  data$age <- rnorm(nrow(data), 60, 10)
  data$sex <- sample(c("M", "F"), nrow(data), replace = TRUE)
  
  gformula_result <- run_gformula_nma(
    nma_data = data,
    treatment_var = "treat1",
    outcome_var = "TE",
    confounders = c("age", "sex"),
    method = "parametric"
  )
  
  expect_s3_class(gformula_result, "gformula_nma")
  expect_true(!is.null(gformula_result$causal_effects))
})

# Test 4: Quantile NMA
test_that("Quantile NMA works", {
  skip_if_not_installed("quantreg")
  
  data <- simulate_nma_data(n_studies = 25)
  
  quantile_result <- run_quantile_nma(
    nma_data = data,
    quantiles = c(0.25, 0.5, 0.75),
    method = "quantile_regression"
  )
  
  expect_s3_class(quantile_result, "quantile_nma")
  expect_equal(length(quantile_result$quantiles), 3)
})

# Test 5: Competing Risks NMA
test_that("Competing risks NMA works", {
  data <- simulate_nma_data(n_studies = 15)
  data$time <- rexp(nrow(data), 0.1)
  data$event <- sample(c("death", "progression", "censored"), nrow(data), replace = TRUE)
  
  cr_result <- run_competing_risks_nma(
    nma_data = data,
    event_types = c("death", "progression"),
    time_var = "time",
    event_var = "event",
    method = "cause_specific"
  )
  
  expect_s3_class(cr_result, "competing_risks_nma")
  expect_equal(length(cr_result$event_results), 2)
})

# Test 6: Crossover Trial Synthesis
test_that("Crossover trial synthesis works", {
  skip_if_not_installed("netmeta")
  
  data <- simulate_nma_data(n_studies = 10)
  
  crossover_result <- synthesize_crossover_trials(
    crossover_data = data,
    within_patient_correlation = 0.6,
    adjust_carryover = TRUE
  )
  
  expect_s3_class(crossover_result, "crossover_nma")
  expect_equal(crossover_result$within_patient_correlation, 0.6)
})

# Test 7: Finite Mixture NMA
test_that("Finite mixture NMA works", {
  skip_if_not_installed("flexmix")
  skip_if_not_installed("netmeta")
  
  data <- simulate_nma_data(n_studies = 30)
  
  mixture_result <- run_finite_mixture_nma(
    nma_data = data,
    n_classes = 2,
    method = "latent_class"
  )
  
  expect_s3_class(mixture_result, "finite_mixture_nma")
  expect_equal(mixture_result$n_classes, 2)
})

# Test 8: GNN Clustering
test_that("GNN-based clustering works", {
  data <- simulate_nma_data(n_studies = 20, n_treatments = 8)
  treatments <- unique(c(data$treat1, data$treat2))
  node_features <- matrix(rnorm(length(treatments) * 10), length(treatments), 10)
  rownames(node_features) <- treatments
  
  gnn_model <- train_gnn_nma(
    nma_data = data,
    node_features = node_features,
    architecture = "GCN",
    n_epochs = 10
  )
  
  clustering <- cluster_treatments_gnn(
    gnn_model = gnn_model,
    n_clusters = 3,
    method = "kmeans"
  )
  
  expect_s3_class(clustering, "gnn_clustering")
  expect_equal(clustering$n_clusters, 3)
})

message("Phase 14 Deep Learning & Causal Inference Tests Complete")
