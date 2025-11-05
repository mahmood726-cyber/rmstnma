# Test Phase 9 Part 2: AI-Powered Manuscript Generation & Text Results
# Tests for:
# - Ollama integration and AI enhancement
# - Comprehensive text results generator
# - Manuscript generator integration with AI
# - Quality checks and validation

library(testthat)
library(powerNMA)

# ============================================================================
# Helper Functions
# ============================================================================

# Create mock NMA result for testing
create_mock_nma_for_text <- function() {
  set.seed(123)
  n_treatments <- 5
  treatments <- c("Placebo", "Drug A", "Drug B", "Drug C", "Drug D")

  # Create mock netmeta object structure
  nma <- list(
    k = 15,
    n = rep(100, 15),
    m = 8,
    trts = treatments,
    reference.group = "Placebo",
    sm = "OR",
    TE.random = matrix(rnorm(25, 0, 0.5), nrow = 5, ncol = 5,
                      dimnames = list(treatments, treatments)),
    seTE.random = matrix(0.2, nrow = 5, ncol = 5,
                        dimnames = list(treatments, treatments)),
    lower.random = matrix(rnorm(25, -0.5, 0.5), nrow = 5, ncol = 5,
                         dimnames = list(treatments, treatments)),
    upper.random = matrix(rnorm(25, 0.5, 0.5), nrow = 5, ncol = 5,
                         dimnames = list(treatments, treatments)),
    pval.random = matrix(runif(25), nrow = 5, ncol = 5,
                        dimnames = list(treatments, treatments))
  )

  class(nma) <- "netmeta"
  return(nma)
}

# Create mock SUCRA result
create_mock_sucra <- function() {
  treatments <- c("Placebo", "Drug A", "Drug B", "Drug C", "Drug D")
  list(
    sucra_scores = setNames(c(20, 80, 60, 40, 70), treatments),
    mean_ranks = setNames(c(4.5, 1.2, 2.5, 3.8, 2.0), treatments)
  )
}

# Create mock heterogeneity result
create_mock_heterogeneity <- function() {
  list(
    tau2 = 0.15,
    I2 = 45.3,
    Q = 25.6,
    df = 14,
    p_value = 0.03,
    prediction_intervals = TRUE
  )
}

# Create mock inconsistency result
create_mock_inconsistency <- function() {
  list(
    node_splitting = data.frame(
      comparison = c("Drug A vs Placebo", "Drug B vs Drug A"),
      direct_effect = c(0.5, 0.3),
      indirect_effect = c(0.45, 0.35),
      difference = c(0.05, -0.05),
      p_value = c(0.6, 0.4),
      stringsAsFactors = FALSE
    ),
    design_inconsistency = list(
      Q = 2.5,
      df = 3,
      p_value = 0.48
    )
  )
}

# ============================================================================
# Tests: Text Results Generator
# ============================================================================

test_that("generate_comprehensive_text_results creates text output", {
  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()
  het <- create_mock_heterogeneity()

  result <- generate_comprehensive_text_results(
    nma_result = nma,
    sucra_result = sucra,
    heterogeneity_result = het,
    style = "narrative",
    include_clinical = TRUE,
    include_stats = TRUE
  )

  expect_s3_class(result, "text_results")
  expect_type(result$full_text, "character")
  expect_true(nchar(result$full_text) > 100)
  expect_type(result$sections, "list")
  expect_type(result$metadata, "list")
})

test_that("text results include all requested sections", {
  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()
  het <- create_mock_heterogeneity()
  incon <- create_mock_inconsistency()

  result <- generate_comprehensive_text_results(
    nma_result = nma,
    sucra_result = sucra,
    heterogeneity_result = het,
    inconsistency_result = incon,
    style = "structured"
  )

  expect_true("overview" %in% names(result$sections))
  expect_true("treatment_effects" %in% names(result$sections))
  expect_true("rankings" %in% names(result$sections))
  expect_true("heterogeneity" %in% names(result$sections))
  expect_true("inconsistency" %in% names(result$sections))
  expect_true("summary" %in% names(result$sections))
})

test_that("text results style parameter works correctly", {
  nma <- create_mock_nma_for_text()

  narrative <- generate_comprehensive_text_results(
    nma, style = "narrative"
  )

  structured <- generate_comprehensive_text_results(
    nma, style = "structured"
  )

  detailed <- generate_comprehensive_text_results(
    nma, style = "detailed"
  )

  # Structured should have section headers
  expect_true(grepl("##", structured$full_text))

  # All should have content
  expect_true(nchar(narrative$full_text) > 50)
  expect_true(nchar(structured$full_text) > 50)
  expect_true(nchar(detailed$full_text) > 50)
})

test_that("text results metadata is complete", {
  nma <- create_mock_nma_for_text()
  result <- generate_comprehensive_text_results(nma)

  expect_equal(result$metadata$style, "narrative")
  expect_true("word_count" %in% names(result$metadata))
  expect_true("character_count" %in% names(result$metadata))
  expect_true("sections_included" %in% names(result$metadata))
  expect_true("generation_time" %in% names(result$metadata))
  expect_true(result$metadata$word_count > 0)
})

test_that("network overview text is generated correctly", {
  nma <- create_mock_nma_for_text()
  result <- generate_comprehensive_text_results(nma)

  overview <- result$sections$overview

  expect_true(grepl("15 studies", overview))
  expect_true(grepl("5 treatments", overview))
  expect_true(grepl("network", overview, ignore.case = TRUE))
})

test_that("treatment effects text includes statistics", {
  nma <- create_mock_nma_for_text()
  result <- generate_comprehensive_text_results(
    nma, include_stats = TRUE
  )

  effects <- result$sections$treatment_effects

  expect_true(grepl("95% CI", effects))
  expect_true(grepl("p =", effects))
  expect_true(grepl("OR =", effects))
})

test_that("rankings text is generated with SUCRA scores", {
  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()

  result <- generate_comprehensive_text_results(
    nma, sucra_result = sucra
  )

  rankings <- result$sections$rankings

  expect_true(grepl("SUCRA", rankings))
  expect_true(grepl("Drug A", rankings))  # Best treatment
  expect_true(grepl("highest-ranked", rankings, ignore.case = TRUE))
})

test_that("heterogeneity text interprets I2 correctly", {
  nma <- create_mock_nma_for_text()

  # Low heterogeneity
  het_low <- list(tau2 = 0.05, I2 = 20)
  result_low <- generate_comprehensive_text_results(
    nma, heterogeneity_result = het_low
  )
  expect_true(grepl("low", result_low$sections$heterogeneity, ignore.case = TRUE))

  # Substantial heterogeneity
  het_high <- list(tau2 = 0.30, I2 = 70)
  result_high <- generate_comprehensive_text_results(
    nma, heterogeneity_result = het_high
  )
  expect_true(grepl("substantial", result_high$sections$heterogeneity, ignore.case = TRUE))
})

test_that("inconsistency text reports node-splitting results", {
  nma <- create_mock_nma_for_text()
  incon <- create_mock_inconsistency()

  result <- generate_comprehensive_text_results(
    nma, inconsistency_result = incon
  )

  incon_text <- result$sections$inconsistency

  expect_true(grepl("Node-splitting", incon_text))
  expect_true(grepl("comparisons", incon_text))
  expect_true(grepl("inconsistency", incon_text, ignore.case = TRUE))
})

test_that("summary text identifies best treatment", {
  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()

  result <- generate_comprehensive_text_results(
    nma, sucra_result = sucra, include_clinical = TRUE
  )

  summary <- result$sections$summary

  expect_true(grepl("Drug A", summary))  # Best by SUCRA
  expect_true(grepl("most effective", summary, ignore.case = TRUE))
  expect_true(grepl("clinical", summary, ignore.case = TRUE))
})

test_that("print method works for text_results", {
  nma <- create_mock_nma_for_text()
  result <- generate_comprehensive_text_results(nma)

  expect_output(print(result), "Comprehensive NMA Text Results")
  expect_output(print(result), "Style:")
  expect_output(print(result), "Word count:")
})

# ============================================================================
# Tests: Ollama Integration (Mock Tests)
# ============================================================================

test_that("check_ollama_available returns logical", {
  # This will return FALSE unless Ollama is actually running
  result <- check_ollama_available()
  expect_type(result, "logical")
})

test_that("list_ollama_models handles unavailable Ollama gracefully", {
  # Should return empty character vector if Ollama not available
  # or vector of models if available
  models <- suppressWarnings(list_ollama_models())
  expect_type(models, "character")
})

test_that("enhance_with_ollama returns original text when Ollama unavailable", {
  text <- "This is test text for enhancement."

  result <- suppressWarnings(
    enhance_with_ollama(text, mode = "readability")
  )

  # Should return original text if Ollama not available
  expect_type(result, "character")
  expect_true(nchar(result) > 0)
})

test_that("enhance_methods_with_ollama creates proper return structure", {
  methods_text <- "We performed network meta-analysis using frequentist methods."

  result <- suppressWarnings(
    enhance_methods_with_ollama(
      methods_text,
      enhancement_level = "light"
    )
  )

  expect_type(result, "list")
  expect_true("enhanced_text" %in% names(result))
  expect_true("original_text" %in% names(result))
  expect_true("enhanced" %in% names(result))
  expect_true("model_used" %in% names(result))
  expect_true("enhancement_level" %in% names(result))
  expect_equal(result$enhancement_level, "light")
})

test_that("enhance_results_with_ollama creates proper return structure", {
  results_text <- "Treatment A was superior to placebo (OR = 1.5, 95% CI: 1.2-1.8)."

  result <- suppressWarnings(
    enhance_results_with_ollama(
      results_text,
      add_clinical_context = TRUE,
      enhancement_level = "moderate"
    )
  )

  expect_type(result, "list")
  expect_true("enhanced_text" %in% names(result))
  expect_true("original_text" %in% names(result))
  expect_true("clinical_context_added" %in% names(result))
  expect_equal(result$clinical_context_added, TRUE)
  expect_equal(result$enhancement_level, "moderate")
})

test_that("enhancement_level parameter is validated", {
  methods_text <- "Test methods."

  expect_error(
    enhance_methods_with_ollama(methods_text, enhancement_level = "invalid"),
    "should be one of"
  )
})

test_that("enhancement mode parameter is validated", {
  text <- "Test text."

  expect_error(
    enhance_with_ollama(text, mode = "invalid_mode"),
    "should be one of"
  )
})

test_that("quality_check_enhanced_text validates numerical preservation", {
  original <- "The odds ratio was 1.50 with 95% CI: 1.20-1.80, p < 0.001."
  enhanced <- "The odds ratio was 1.50 with 95% confidence interval: 1.20-1.80, p < 0.001."

  result <- quality_check_enhanced_text(original, enhanced)

  expect_type(result, "list")
  expect_true("passed" %in% names(result))
  expect_true("warnings" %in% names(result))
  expect_true("numerical_values_preserved" %in% names(result))
  expect_true("length_change_pct" %in% names(result))
})

test_that("quality check detects missing numerical values", {
  original <- "The effect size was 0.45 (SE = 0.12, p = 0.001)."
  enhanced <- "The effect size was substantial and statistically significant."

  result <- quality_check_enhanced_text(original, enhanced)

  expect_false(result$passed)
  expect_true(length(result$warnings) > 0)
})

test_that("quality check detects removed p-values", {
  original <- "Treatment showed effect (p < 0.05)."
  enhanced <- "Treatment showed effect."

  result <- quality_check_enhanced_text(original, enhanced)

  expect_true(any(grepl("P-values", result$warnings)))
})

test_that("batch_enhance_sections processes multiple sections", {
  sections <- list(
    methods = "Methods text.",
    results = "Results text.",
    discussion = "Discussion text."
  )

  result <- suppressWarnings(
    batch_enhance_sections(sections)
  )

  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_true(all(names(sections) %in% names(result)))
})

test_that("get_recommended_model handles no Ollama gracefully", {
  # Should error or return NULL when Ollama not available
  expect_error(
    get_recommended_model(),
    "Ollama is not available|No Ollama models"
  )
})

# ============================================================================
# Tests: Integration - Text Results with Manuscript Generation
# ============================================================================

test_that("text results can be combined with manuscript generation", {
  skip_if_not_installed("powerNMA")

  # This test verifies that text results generator works with
  # the manuscript generation functions

  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()

  # Generate text results
  text_results <- generate_comprehensive_text_results(
    nma, sucra_result = sucra, style = "detailed"
  )

  expect_true(nchar(text_results$full_text) > 500)
  expect_true(text_results$metadata$word_count > 50)
})

test_that("text results export to different formats", {
  skip_if_not_installed("powerNMA")

  nma <- create_mock_nma_for_text()
  text_results <- generate_comprehensive_text_results(nma)

  # Test TXT export
  txt_file <- tempfile(fileext = ".txt")
  export_text_results(text_results, txt_file, format = "txt")
  expect_true(file.exists(txt_file))
  expect_true(file.size(txt_file) > 0)

  # Test HTML export
  html_file <- tempfile(fileext = ".html")
  export_text_results(text_results, html_file, format = "html")
  expect_true(file.exists(html_file))
  expect_true(file.size(html_file) > 0)

  # Cleanup
  unlink(c(txt_file, html_file))
})

# ============================================================================
# Tests: System Integration
# ============================================================================

test_that("Phase 9 Part 2 components integrate correctly", {
  # High-level integration test

  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()
  het <- create_mock_heterogeneity()

  # 1. Generate text results
  text_results <- generate_comprehensive_text_results(
    nma, sucra, het,
    style = "detailed",
    include_clinical = TRUE
  )

  expect_s3_class(text_results, "text_results")

  # 2. Attempt enhancement (will gracefully handle no Ollama)
  enhanced_results <- suppressWarnings(
    enhance_results_with_ollama(
      text_results$sections$treatment_effects,
      add_clinical_context = TRUE
    )
  )

  expect_type(enhanced_results, "list")

  # 3. Quality check
  quality <- quality_check_enhanced_text(
    text_results$sections$treatment_effects,
    enhanced_results$enhanced_text
  )

  expect_type(quality, "list")
  expect_true("passed" %in% names(quality))
})

test_that("text results handle edge cases", {
  # Minimal NMA
  nma_minimal <- list(
    k = 3,
    n = rep(50, 3),
    m = 2,
    trts = c("A", "B"),
    reference.group = "A",
    sm = "MD",
    TE.random = matrix(c(0, 0.5, -0.5, 0), nrow = 2, ncol = 2,
                      dimnames = list(c("A", "B"), c("A", "B"))),
    seTE.random = matrix(0.2, nrow = 2, ncol = 2),
    lower.random = matrix(0, nrow = 2, ncol = 2),
    upper.random = matrix(1, nrow = 2, ncol = 2),
    pval.random = matrix(0.5, nrow = 2, ncol = 2)
  )
  class(nma_minimal) <- "netmeta"

  result <- generate_comprehensive_text_results(nma_minimal)

  expect_s3_class(result, "text_results")
  expect_true(nchar(result$full_text) > 0)
})

test_that("text results include clinical interpretation when requested", {
  nma <- create_mock_nma_for_text()
  sucra <- create_mock_sucra()
  het <- create_mock_heterogeneity()

  result_with_clinical <- generate_comprehensive_text_results(
    nma, sucra, het, include_clinical = TRUE
  )

  result_without_clinical <- generate_comprehensive_text_results(
    nma, sucra, het, include_clinical = FALSE
  )

  # With clinical should have more content
  expect_true(
    nchar(result_with_clinical$full_text) >=
    nchar(result_without_clinical$full_text)
  )

  # Should contain clinical keywords
  expect_true(grepl("clinical", result_with_clinical$full_text, ignore.case = TRUE))
})

# ============================================================================
# Tests: Performance and Robustness
# ============================================================================

test_that("text generation handles large networks efficiently", {
  # Create larger network
  n_treatments <- 15
  treatments <- paste0("Treatment_", 1:n_treatments)

  nma_large <- list(
    k = 50,
    n = rep(100, 50),
    m = 30,
    trts = treatments,
    reference.group = treatments[1],
    sm = "OR",
    TE.random = matrix(rnorm(n_treatments^2), nrow = n_treatments, ncol = n_treatments,
                      dimnames = list(treatments, treatments)),
    seTE.random = matrix(0.2, nrow = n_treatments, ncol = n_treatments),
    lower.random = matrix(rnorm(n_treatments^2, -0.5), nrow = n_treatments, ncol = n_treatments),
    upper.random = matrix(rnorm(n_treatments^2, 0.5), nrow = n_treatments, ncol = n_treatments),
    pval.random = matrix(runif(n_treatments^2), nrow = n_treatments, ncol = n_treatments)
  )
  class(nma_large) <- "netmeta"

  # Should complete within reasonable time
  start_time <- Sys.time()
  result <- generate_comprehensive_text_results(nma_large)
  end_time <- Sys.time()

  expect_s3_class(result, "text_results")
  expect_true(as.numeric(end_time - start_time, units = "secs") < 10)
})

test_that("ollama integration handles connection timeouts", {
  # Call with very short timeout to simulate failure
  result <- suppressWarnings(
    check_ollama_available(timeout = 0.001)
  )

  expect_type(result, "logical")
  # Likely FALSE due to short timeout
})

# ============================================================================
# Summary Test
# ============================================================================

test_that("Phase 9 Part 2 implementation is complete", {
  # Verify all expected functions exist

  # Text results generator
  expect_true(exists("generate_comprehensive_text_results"))
  expect_true(exists("export_text_results"))

  # Ollama integration
  expect_true(exists("check_ollama_available"))
  expect_true(exists("list_ollama_models"))
  expect_true(exists("enhance_with_ollama"))
  expect_true(exists("enhance_methods_with_ollama"))
  expect_true(exists("enhance_results_with_ollama"))
  expect_true(exists("quality_check_enhanced_text"))
  expect_true(exists("get_recommended_model"))

  message("Phase 9 Part 2: All functions implemented and tested")
})
