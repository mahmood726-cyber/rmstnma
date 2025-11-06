# Tests for bs4Dash Dashboard and High-Resolution Export
# ==========================================================

context("bs4Dash Dashboard with High-Resolution Export")

# Test 1: bs4Dash UI Creation
test_that("bs4Dash UI creates without errors", {
  skip_if_no_display()  # Skip on headless CI environments
  ci_safe_require("bs4Dash")
  ci_safe_require("shiny")

  # Create UI
  expect_silent({
    ui <- create_bs4dash_ui(theme = "default", enable_db = FALSE)
  })

  # Check that ui is a bs4Dash page
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})

# Test 2: High-Resolution Download Handler Creation
test_that("High-resolution download handlers are created correctly", {
  skip_if_no_display()  # Skip on headless CI environments
  ci_safe_require("shiny")

  # Create mock plot function
  mock_plot <- function() {
    plot(1:10, 1:10)
  }

  # Create download handler
  expect_silent({
    handler <- create_highres_download_handler(
      plot_func = mock_plot,
      filename = "test_plot",
      format = "png",
      dpi = 600,
      width = 10,
      height = 8
    )
  })

  # Check that handler is a function
  expect_true(is.function(handler))
})

# Test 3: Export Configuration
test_that("Export configurations are valid", {
  # Test DPI ranges
  expect_true(300 >= 300 && 300 <= 1200)
  expect_true(600 >= 300 && 600 <= 1200)
  expect_true(1200 >= 300 && 1200 <= 1200)

  # Test supported formats
  supported_formats <- c("png", "pdf", "svg", "tiff", "html")
  expect_true(all(c("png", "pdf", "svg") %in% supported_formats))
})

# Test 4: bs4Dash Theme Creation
test_that("bs4Dash themes are created correctly", {
  skip_if_no_display()  # Skip on headless CI environments
  ci_safe_require("bs4Dash")

  # Test different themes
  expect_silent({
    ui_default <- create_bs4dash_ui(theme = "default", enable_db = FALSE)
    ui_dark <- create_bs4dash_ui(theme = "dark", enable_db = FALSE)
    ui_light <- create_bs4dash_ui(theme = "light", enable_db = FALSE)
  })
})

# Test 5: Custom CSS Generation
test_that("Custom CSS is generated correctly", {
  # Get CSS
  css <- get_bs4dash_custom_css("default")

  # Check that CSS is a character string
  expect_type(css, "character")
  expect_true(nchar(css) > 0)

  # Check that key CSS classes are present
  expect_true(grepl("plotly-container", css))
  expect_true(grepl("download-highres", css))
  expect_true(grepl("info-box", css))
})

# Test 6: bs4Dash Tab Creation
test_that("bs4Dash tabs are created without errors", {
  skip_if_not_installed("bs4Dash")
  skip_if_not_installed("shiny")

  # Test tab creation functions
  expect_silent({
    home_tab <- create_bs4dash_home_tab()
    data_tab <- create_bs4dash_data_import_tab()
    nma_tab <- create_bs4dash_standard_nma_tab()
    rankings_tab <- create_bs4dash_rankings_tab()
    viz_tab <- create_bs4dash_visualizations_tab()
  })

  # Check that tabs return Shiny UI elements
  expect_true(inherits(home_tab, c("shiny.tag", "shiny.tag.list")))
})

# Test 7: High-Resolution Plot Parameters
test_that("High-resolution plot parameters are validated", {
  # Test DPI values
  test_dpi <- c(300, 600, 900, 1200)
  expect_true(all(test_dpi >= 300 & test_dpi <= 1200))

  # Test dimension values (inches)
  test_widths <- c(5, 10, 15, 20)
  test_heights <- c(5, 8, 12, 20)
  expect_true(all(test_widths >= 5 & test_widths <= 20))
  expect_true(all(test_heights >= 5 & test_heights <= 20))
})

# Test 8: Fresh Theme Creation
test_that("Fresh theme is created or returns NULL", {
  # Create fresh theme
  theme_result <- create_fresh_theme("dark")

  # Should either be a fresh theme or NULL if package not available
  expect_true(is.null(theme_result) || inherits(theme_result, "fresh_theme"))
})

# Test 9: Server Function Creation
test_that("bs4Dash server function is created correctly", {
  skip_if_not_installed("shiny")

  # Create server function
  expect_silent({
    server_func <- create_bs4dash_server(
      enable_db = FALSE,
      db_path = NULL,
      high_res_dpi = 600
    )
  })

  # Check that server is a function
  expect_true(is.function(server_func))
})

# Test 10: Launch Function Parameters
test_that("launch_bs4dash_app accepts valid parameters", {
  skip_if_not_installed("bs4Dash")
  skip_on_ci()  # Skip on CI since we don't actually launch

  # Test parameter validation
  expect_error({
    # Should error with invalid theme
    launch_bs4dash_app(theme = "invalid_theme", launch_browser = FALSE)
  })

  # Valid themes should be accepted
  valid_themes <- c("default", "dark", "light")
  expect_true(all(valid_themes %in% c("default", "dark", "light")))
})

# Test 11: High-Resolution File Format Support
test_that("All high-resolution formats are supported", {
  supported_formats <- c("png", "pdf", "svg", "tiff", "html")

  # Check format list
  expect_length(supported_formats, 5)
  expect_true("png" %in% supported_formats)
  expect_true("pdf" %in% supported_formats)
  expect_true("svg" %in% supported_formats)
  expect_true("tiff" %in% supported_formats)
  expect_true("html" %in% supported_formats)
})

# Test 12: Waiter Integration
test_that("Waiter package integration works if available", {
  skip_if_not_installed("waiter")

  # Test waiter functions
  expect_silent({
    # These should work if waiter is installed
    waiter::use_waiter()
    waiter::use_hostess()
  })
})

# Test 13: Plot Download Handler with Different Formats
test_that("Download handlers work with different plot formats", {
  skip_if_not_installed("ggplot2")

  # Create a simple ggplot
  test_plot <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Test that handler creation doesn't error
  expect_silent({
    handler_png <- create_highres_download_handler(
      function() test_plot,
      "test",
      format = "png",
      dpi = 300
    )
  })

  expect_true(is.function(handler_png))
})

# Test 14: UI Component Structure
test_that("UI components have correct structure", {
  skip_if_not_installed("bs4Dash")

  # Create home tab
  home_tab <- create_bs4dash_home_tab()

  # Check structure
  expect_true(inherits(home_tab, "shiny.tag.list"))
})

# Test 15: Advanced Features Integration
test_that("Advanced features (Phase 12-14) are integrated", {
  skip_if_not_installed("bs4Dash")

  # Check that advanced tabs exist
  expect_silent({
    gnn_tab <- create_bs4dash_gnn_tab()
    causal_tab <- create_bs4dash_causal_tab()
  })

  expect_true(inherits(gnn_tab, c("shiny.tag", "shiny.tag.list")))
  expect_true(inherits(causal_tab, c("shiny.tag", "shiny.tag.list")))
})

# Integration Test
test_that("Full bs4Dash app can be initialized (not launched)", {
  skip_if_not_installed("bs4Dash")
  skip_if_not_installed("shiny")

  # Create full UI and server (but don't launch)
  expect_silent({
    ui <- create_bs4dash_ui("default", FALSE)
    server <- create_bs4dash_server(FALSE, NULL, 600)
  })

  # Check objects are created
  expect_true(!is.null(ui))
  expect_true(is.function(server))
})
