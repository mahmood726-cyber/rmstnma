#' Ollama Integration for AI-Powered Text Enhancement
#'
#' @description
#' This module provides integration with local Ollama LLM instances for
#' AI-powered enhancement of generated manuscript sections. It can improve
#' readability, add clinical context, polish language, and ensure quality.
#'
#' @details
#' The module requires a local Ollama instance running. Install Ollama from
#' https://ollama.ai/ and ensure the desired model is pulled (e.g., llama2,
#' mistral, or specialized medical models).
#'
#' @references
#' Ollama: https://ollama.ai/
#'
#' @author powerNMA Development Team
#' @name ollama_integration
NULL

#' Check Ollama Availability
#'
#' @description
#' Checks if Ollama is installed and running on the local system.
#'
#' @param host Character string specifying Ollama host (default: "http://localhost:11434")
#' @param timeout Numeric timeout in seconds for connection check
#'
#' @return Logical indicating if Ollama is available
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (check_ollama_available()) {
#'   message("Ollama is ready!")
#' }
#' }
check_ollama_available <- function(host = "http://localhost:11434", timeout = 5) {

  tryCatch({
    # Try to connect to Ollama API
    response <- httr::GET(
      paste0(host, "/api/tags"),
      httr::timeout(timeout)
    )

    if (httr::status_code(response) == 200) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }, error = function(e) {
    return(FALSE)
  })
}

#' List Available Ollama Models
#'
#' @description
#' Lists all models available in the local Ollama instance.
#'
#' @param host Character string specifying Ollama host
#'
#' @return Character vector of available model names
#'
#' @export
#'
#' @examples
#' \dontrun{
#' models <- list_ollama_models()
#' print(models)
#' }
list_ollama_models <- function(host = "http://localhost:11434") {

  if (!check_ollama_available(host)) {
    warning("Ollama is not available. Please ensure it is installed and running.")
    return(character(0))
  }

  tryCatch({
    response <- httr::GET(paste0(host, "/api/tags"))
    content <- httr::content(response, as = "parsed")

    if (!is.null(content$models)) {
      models <- sapply(content$models, function(m) m$name)
      return(models)
    } else {
      return(character(0))
    }
  }, error = function(e) {
    warning("Failed to retrieve models: ", e$message)
    return(character(0))
  })
}

#' Call Ollama API
#'
#' @description
#' Low-level function to make API calls to Ollama.
#'
#' @param prompt Character string with the prompt
#' @param model Character string specifying model name
#' @param host Character string specifying Ollama host
#' @param temperature Numeric temperature for sampling (0-1)
#' @param max_tokens Maximum number of tokens to generate
#' @param system_prompt Optional system prompt
#'
#' @return Character string with model response
#'
#' @keywords internal
call_ollama_api <- function(prompt, model = "llama2",
                           host = "http://localhost:11434",
                           temperature = 0.7, max_tokens = 2000,
                           system_prompt = NULL) {

  if (!check_ollama_available(host)) {
    stop("Ollama is not available. Please ensure it is installed and running.")
  }

  # Prepare request body
  body <- list(
    model = model,
    prompt = prompt,
    stream = FALSE,
    options = list(
      temperature = temperature,
      num_predict = max_tokens
    )
  )

  if (!is.null(system_prompt)) {
    body$system <- system_prompt
  }

  # Make request
  tryCatch({
    response <- httr::POST(
      paste0(host, "/api/generate"),
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      httr::content_type_json(),
      httr::timeout(120)
    )

    if (httr::status_code(response) == 200) {
      content <- httr::content(response, as = "parsed")
      return(content$response)
    } else {
      stop("Ollama API returned error: ", httr::status_code(response))
    }
  }, error = function(e) {
    stop("Failed to call Ollama API: ", e$message)
  })
}

#' Enhance Text with Ollama
#'
#' @description
#' Main function to enhance any text using Ollama AI. Can improve readability,
#' add context, or polish language based on specified mode.
#'
#' @param text Character string with text to enhance
#' @param mode Character string specifying enhancement mode:
#'   \itemize{
#'     \item "readability" - Improve clarity and flow
#'     \item "clinical" - Add clinical context and interpretation
#'     \item "polish" - Polish language and grammar
#'     \item "expand" - Expand with additional details
#'     \item "concise" - Make more concise
#'     \item "academic" - Enhance academic tone
#'   }
#' @param model Character string specifying Ollama model to use
#' @param host Character string specifying Ollama host
#' @param temperature Numeric temperature for sampling
#' @param preserve_structure Logical whether to preserve original structure
#'
#' @return Character string with enhanced text
#'
#' @export
#'
#' @examples
#' \dontrun{
#' text <- "We found treatment A was better than B (OR=1.5, p<0.05)."
#' enhanced <- enhance_with_ollama(text, mode = "clinical")
#' cat(enhanced)
#' }
enhance_with_ollama <- function(text,
                               mode = c("readability", "clinical", "polish",
                                       "expand", "concise", "academic"),
                               model = "llama2",
                               host = "http://localhost:11434",
                               temperature = 0.7,
                               preserve_structure = TRUE) {

  mode <- match.arg(mode)

  if (!check_ollama_available(host)) {
    warning("Ollama not available. Returning original text.")
    return(text)
  }

  # Build system prompt based on mode
  system_prompt <- build_enhancement_system_prompt(mode, preserve_structure)

  # Build user prompt
  user_prompt <- sprintf(
    "Please enhance the following text according to the instructions:\n\n%s",
    text
  )

  # Call Ollama
  enhanced_text <- call_ollama_api(
    prompt = user_prompt,
    model = model,
    host = host,
    temperature = temperature,
    max_tokens = nchar(text) * 2,  # Allow expansion
    system_prompt = system_prompt
  )

  # Clean up response
  enhanced_text <- clean_ollama_response(enhanced_text)

  return(enhanced_text)
}

#' Build Enhancement System Prompt
#'
#' @description
#' Creates appropriate system prompt based on enhancement mode.
#'
#' @param mode Enhancement mode
#' @param preserve_structure Whether to preserve structure
#'
#' @return Character string with system prompt
#'
#' @keywords internal
build_enhancement_system_prompt <- function(mode, preserve_structure) {

  base_prompt <- "You are an expert medical writer specializing in systematic reviews and meta-analyses. "

  structure_instruction <- if (preserve_structure) {
    "Preserve the overall structure and organization of the text. "
  } else {
    "You may reorganize the text if needed. "
  }

  mode_instructions <- switch(
    mode,
    readability = paste0(
      "Improve the clarity and readability of the text. ",
      "Use clear, straightforward language while maintaining scientific accuracy. ",
      "Improve flow and coherence between sentences. ",
      "Keep all numerical results and statistical values exactly as provided."
    ),
    clinical = paste0(
      "Add clinical context and interpretation to the statistical results. ",
      "Explain what the findings mean for clinical practice. ",
      "Maintain all statistical values exactly as provided. ",
      "Add relevant clinical implications where appropriate."
    ),
    polish = paste0(
      "Polish the language and grammar. ",
      "Improve sentence structure and word choice. ",
      "Ensure consistency in terminology and style. ",
      "Maintain all numerical values exactly as provided."
    ),
    expand = paste0(
      "Expand the text with additional relevant details. ",
      "Add context and explanation where helpful. ",
      "Maintain scientific rigor and accuracy. ",
      "Keep all original numerical values and add contextual information."
    ),
    concise = paste0(
      "Make the text more concise while preserving all key information. ",
      "Remove redundancy and wordiness. ",
      "Maintain all critical statistical results. ",
      "Improve clarity through brevity."
    ),
    academic = paste0(
      "Enhance the academic tone and style. ",
      "Use appropriate academic language. ",
      "Ensure proper citation style and terminology. ",
      "Maintain scientific rigor and precision."
    )
  )

  return(paste0(base_prompt, structure_instruction, mode_instructions))
}

#' Clean Ollama Response
#'
#' @description
#' Cleans and normalizes Ollama API response.
#'
#' @param text Character string with raw response
#'
#' @return Cleaned character string
#'
#' @keywords internal
clean_ollama_response <- function(text) {

  # Remove potential markdown formatting
  text <- gsub("```.*?```", "", text, perl = TRUE)

  # Remove excessive whitespace
  text <- gsub("\\s+", " ", text)
  text <- gsub("\\n\\s*\\n\\s*\\n+", "\n\n", text)

  # Trim
  text <- trimws(text)

  return(text)
}

#' Enhance Methods Section with Ollama
#'
#' @description
#' Specialized function to enhance a generated Methods section using Ollama AI.
#'
#' @param methods_text Character string with Methods section text
#' @param enhancement_level Character string: "light", "moderate", "heavy"
#' @param model Character string specifying Ollama model
#' @param host Character string specifying Ollama host
#'
#' @return List with enhanced text and metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' methods <- generate_methods_section(nma_result)
#' enhanced_methods <- enhance_methods_with_ollama(methods$text)
#' }
enhance_methods_with_ollama <- function(methods_text,
                                       enhancement_level = c("light", "moderate", "heavy"),
                                       model = "llama2",
                                       host = "http://localhost:11434") {

  enhancement_level <- match.arg(enhancement_level)

  if (!check_ollama_available(host)) {
    warning("Ollama not available. Returning original text.")
    return(list(
      enhanced_text = methods_text,
      original_text = methods_text,
      enhanced = FALSE,
      model_used = NA,
      enhancement_level = enhancement_level
    ))
  }

  # Apply enhancements based on level
  enhanced <- methods_text

  if (enhancement_level %in% c("light", "moderate", "heavy")) {
    # Improve readability
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "readability",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  if (enhancement_level %in% c("moderate", "heavy")) {
    # Polish language
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "polish",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  if (enhancement_level == "heavy") {
    # Enhance academic tone
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "academic",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  return(list(
    enhanced_text = enhanced,
    original_text = methods_text,
    enhanced = TRUE,
    model_used = model,
    enhancement_level = enhancement_level,
    character_increase = nchar(enhanced) - nchar(methods_text),
    character_increase_pct = 100 * (nchar(enhanced) - nchar(methods_text)) / nchar(methods_text)
  ))
}

#' Enhance Results Section with Ollama
#'
#' @description
#' Specialized function to enhance a generated Results section using Ollama AI.
#'
#' @param results_text Character string with Results section text
#' @param add_clinical_context Logical whether to add clinical interpretation
#' @param enhancement_level Character string: "light", "moderate", "heavy"
#' @param model Character string specifying Ollama model
#' @param host Character string specifying Ollama host
#'
#' @return List with enhanced text and metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- generate_results_section(nma_result, sucra_result)
#' enhanced_results <- enhance_results_with_ollama(results$text,
#'                                                  add_clinical_context = TRUE)
#' }
enhance_results_with_ollama <- function(results_text,
                                       add_clinical_context = TRUE,
                                       enhancement_level = c("light", "moderate", "heavy"),
                                       model = "llama2",
                                       host = "http://localhost:11434") {

  enhancement_level <- match.arg(enhancement_level)

  if (!check_ollama_available(host)) {
    warning("Ollama not available. Returning original text.")
    return(list(
      enhanced_text = results_text,
      original_text = results_text,
      enhanced = FALSE,
      model_used = NA,
      enhancement_level = enhancement_level,
      clinical_context_added = FALSE
    ))
  }

  # Apply enhancements
  enhanced <- results_text

  # Add clinical context if requested
  if (add_clinical_context) {
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "clinical",
      model = model,
      host = host,
      temperature = 0.6
    )
  }

  # Apply level-based enhancements
  if (enhancement_level %in% c("light", "moderate", "heavy")) {
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "readability",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  if (enhancement_level %in% c("moderate", "heavy")) {
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "polish",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  if (enhancement_level == "heavy") {
    enhanced <- enhance_with_ollama(
      enhanced,
      mode = "academic",
      model = model,
      host = host,
      temperature = 0.5
    )
  }

  return(list(
    enhanced_text = enhanced,
    original_text = results_text,
    enhanced = TRUE,
    model_used = model,
    enhancement_level = enhancement_level,
    clinical_context_added = add_clinical_context,
    character_increase = nchar(enhanced) - nchar(results_text),
    character_increase_pct = 100 * (nchar(enhanced) - nchar(results_text)) / nchar(results_text)
  ))
}

#' Batch Enhance Multiple Sections
#'
#' @description
#' Enhances multiple manuscript sections in batch.
#'
#' @param sections Named list of text sections to enhance
#' @param modes Named list of enhancement modes for each section (optional)
#' @param model Character string specifying Ollama model
#' @param host Character string specifying Ollama host
#' @param progress Logical whether to show progress
#'
#' @return Named list of enhanced sections with metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sections <- list(
#'   methods = "Methods text...",
#'   results = "Results text..."
#' )
#' enhanced_sections <- batch_enhance_sections(sections)
#' }
batch_enhance_sections <- function(sections,
                                  modes = NULL,
                                  model = "llama2",
                                  host = "http://localhost:11434",
                                  progress = TRUE) {

  if (!check_ollama_available(host)) {
    warning("Ollama not available. Returning original sections.")
    return(lapply(sections, function(s) {
      list(enhanced_text = s, original_text = s, enhanced = FALSE)
    }))
  }

  # Default modes if not specified
  if (is.null(modes)) {
    modes <- setNames(rep("readability", length(sections)), names(sections))
  }

  # Enhance each section
  enhanced_sections <- list()

  for (i in seq_along(sections)) {
    section_name <- names(sections)[i]

    if (progress) {
      message(sprintf("Enhancing section %d of %d: %s",
                     i, length(sections), section_name))
    }

    enhanced_sections[[section_name]] <- enhance_with_ollama(
      sections[[i]],
      mode = modes[[section_name]],
      model = model,
      host = host
    )
  }

  return(enhanced_sections)
}

#' Quality Check Enhanced Text
#'
#' @description
#' Performs quality checks on AI-enhanced text to ensure it maintains
#' statistical accuracy and scientific rigor.
#'
#' @param original Character string with original text
#' @param enhanced Character string with enhanced text
#'
#' @return List with quality metrics and warnings
#'
#' @export
#'
#' @examples
#' \dontrun{
#' original <- "OR = 1.50, 95% CI: 1.20-1.80, p < 0.001"
#' enhanced <- enhance_with_ollama(original)
#' quality <- quality_check_enhanced_text(original, enhanced)
#' }
quality_check_enhanced_text <- function(original, enhanced) {

  warnings <- character(0)

  # Check if numerical values are preserved
  original_numbers <- gregexpr("[0-9]+\\.?[0-9]*", original)[[1]]
  enhanced_numbers <- gregexpr("[0-9]+\\.?[0-9]*", enhanced)[[1]]

  if (length(original_numbers) != length(enhanced_numbers)) {
    warnings <- c(warnings, "Warning: Number of numerical values changed")
  }

  # Check p-value reporting
  original_pvals <- grepl("p\\s*[<>=]\\s*[0-9.]+", original, ignore.case = TRUE)
  enhanced_pvals <- grepl("p\\s*[<>=]\\s*[0-9.]+", enhanced, ignore.case = TRUE)

  if (original_pvals && !enhanced_pvals) {
    warnings <- c(warnings, "Warning: P-values may have been removed or altered")
  }

  # Check confidence intervals
  original_ci <- grepl("confidence interval|\\bCI\\b|95%", original, ignore.case = TRUE)
  enhanced_ci <- grepl("confidence interval|\\bCI\\b|95%", enhanced, ignore.case = TRUE)

  if (original_ci && !enhanced_ci) {
    warnings <- c(warnings, "Warning: Confidence intervals may have been removed")
  }

  # Check length change
  length_change_pct <- 100 * (nchar(enhanced) - nchar(original)) / nchar(original)

  if (abs(length_change_pct) > 100) {
    warnings <- c(warnings, sprintf("Warning: Large length change (%.1f%%)", length_change_pct))
  }

  return(list(
    passed = length(warnings) == 0,
    warnings = warnings,
    original_length = nchar(original),
    enhanced_length = nchar(enhanced),
    length_change_pct = length_change_pct,
    numerical_values_preserved = length(original_numbers) == length(enhanced_numbers)
  ))
}

#' Get Recommended Ollama Model
#'
#' @description
#' Recommends the best available Ollama model for medical text enhancement.
#'
#' @param host Character string specifying Ollama host
#' @param prefer_medical Logical whether to prefer medical models
#'
#' @return Character string with recommended model name
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model <- get_recommended_model()
#' print(paste("Using model:", model))
#' }
get_recommended_model <- function(host = "http://localhost:11434",
                                 prefer_medical = TRUE) {

  if (!check_ollama_available(host)) {
    stop("Ollama is not available")
  }

  available_models <- list_ollama_models(host)

  if (length(available_models) == 0) {
    stop("No Ollama models found. Please pull a model first.")
  }

  # Preference order (medical models first if preferred)
  if (prefer_medical) {
    preference_order <- c(
      "meditron", "medllama", "biomistral",  # Medical-specific
      "mixtral", "llama3", "llama2",          # General large models
      "mistral", "gemma"                      # Smaller models
    )
  } else {
    preference_order <- c(
      "mixtral", "llama3", "llama2",
      "mistral", "gemma",
      "meditron", "medllama", "biomistral"
    )
  }

  # Find first available model from preference list
  for (preferred in preference_order) {
    matching <- available_models[grepl(preferred, available_models, ignore.case = TRUE)]
    if (length(matching) > 0) {
      return(matching[1])
    }
  }

  # If no preferred model found, return first available
  return(available_models[1])
}
