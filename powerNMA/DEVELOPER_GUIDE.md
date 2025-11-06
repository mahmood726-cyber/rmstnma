# powerNMA Developer Guide

## Package Architecture

### Core Modules

```
powerNMA/
├── R/
│   ├── utils.R                  # Core utilities and helpers
│   ├── data_functions.R         # Data generation and validation
│   ├── nma_core.R              # Core NMA functions
│   ├── powernma.R              # Main analysis runner
│   └── powernma-package.R      # Package documentation
├── inst/
│   ├── extdata/                # Example datasets
│   └── shiny/                  # Shiny dashboard (future)
├── tests/
│   └── testthat/               # Unit tests
├── man/                        # Documentation (auto-generated)
└── vignettes/                  # Long-form documentation
```

## Development Workflow

### Setup

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "testthat", "usethis"))

# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

### Adding New Features

1. **Write the function** in appropriate `R/*.R` file
2. **Add roxygen2 documentation** with `#'` comments
3. **Export if public** with `@export` tag
4. **Write tests** in `tests/testthat/test-*.R`
5. **Update documentation** with `devtools::document()`
6. **Run checks** with `devtools::check()`

### Documentation Standards

```r
#' Function Title
#'
#' Detailed description of what the function does.
#'
#' @param param1 Description of param1
#' @param param2 Description of param2
#' @return Description of return value
#' @export
#' @examples
#' \dontrun{
#' # Example code
#' result <- my_function(param1 = "value")
#' }
my_function <- function(param1, param2 = NULL) {
  # Function body
}
```

## Code Style

### Naming Conventions

- **Functions**: `snake_case` (e.g., `run_powernma`)
- **Variables**: `snake_case` (e.g., `ref_treatment`)
- **Classes**: `snake_case` with `_` (e.g., `powernma_result`)
- **Constants**: `UPPER_CASE` (e.g., `DEFAULT_SEED`)

### Code Organization

```r
# 1. Function documentation
#' Function title
#' @param x Description
#' @export

# 2. Input validation
function_name <- function(x, y = NULL) {
  if (!is.numeric(x)) stop("x must be numeric")

  # 3. Main logic
  result <- process_data(x, y)

  # 4. Return with appropriate class
  class(result) <- c("custom_class", "list")
  result
}
```

### Error Handling

```r
# Use informative error messages
if (condition) {
  stop("Descriptive error message with context")
}

# Use warnings for non-fatal issues
if (potential_issue) {
  warning("Warning message explaining the issue")
}

# Use messages for user information
msg("Processing %d studies...", n_studies)
```

## Testing

### Test Structure

```r
test_that("descriptive test name", {
  # Setup
  data <- generate_test_data()

  # Execute
  result <- my_function(data)

  # Assert
  expect_equal(result$value, expected_value)
  expect_true(is.numeric(result$metric))
  expect_s3_class(result, "expected_class")
})
```

### Test Coverage Goals

- **Core functions**: 100% coverage
- **Utility functions**: 90%+ coverage
- **Edge cases**: Explicitly tested
- **Error conditions**: Explicitly tested

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific file
devtools::test_file("tests/testthat/test-basic.R")

# Check coverage
covr::package_coverage()
```

## Adding Dependencies

### Required Packages (Imports)

```r
# In DESCRIPTION:
Imports:
    newpackage (>= 1.0.0)

# In code:
newpackage::function_name()  # Explicit namespace
```

### Optional Packages (Suggests)

```r
# In DESCRIPTION:
Suggests:
    optionalpackage

# In code:
if (has_pkg("optionalpackage")) {
  optionalpackage::advanced_feature()
} else {
  message("Install 'optionalpackage' for advanced features")
}
```

## Version Control

### Commit Messages

Follow conventional commits:

```
type(scope): subject

body (optional)

footer (optional)
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `test`: Tests
- `refactor`: Code refactoring
- `perf`: Performance improvement
- `chore`: Maintenance

Examples:
```
feat(nma): add transportability weighting
fix(rmst): correct sign convention in treatment effects
docs(readme): update installation instructions
test(milestone): add zero-event edge case tests
```

## Release Checklist

### Pre-Release

- [ ] All tests pass (`devtools::test()`)
- [ ] No R CMD check errors/warnings (`devtools::check()`)
- [ ] Documentation is up to date (`devtools::document()`)
- [ ] NEWS.md updated with changes
- [ ] Version number bumped in DESCRIPTION
- [ ] Examples run successfully
- [ ] Vignettes build without errors

### Release Process

1. Update version in DESCRIPTION (semantic versioning)
2. Update NEWS.md with release notes
3. Run full check: `devtools::check()`
4. Build package: `devtools::build()`
5. Tag release in git: `git tag v1.0.0`
6. Push to GitHub: `git push --tags`

## Performance Optimization

### Profiling

```r
# Profile code
prof <- profvis::profvis({
  results <- run_powernma(large_dataset)
})
print(prof)
```

### Common Optimizations

1. **Vectorize operations** instead of loops
2. **Pre-allocate vectors/lists** before filling
3. **Use appropriate data structures** (data.table for large data)
4. **Cache expensive computations** with memoization
5. **Parallelize independent operations** where beneficial

## Troubleshooting

### Common Issues

**Issue**: Package won't load
```r
# Solution: Clean and rebuild
devtools::clean_dll()
devtools::load_all()
```

**Issue**: Documentation not updating
```r
# Solution: Force regeneration
devtools::document(roclets = c('rd', 'collate', 'namespace'))
```

**Issue**: Tests failing intermittently
```r
# Solution: Set seed, check for timing issues
set.seed(42)
```

## Contributing Guidelines

### Pull Request Process

1. Fork the repository
2. Create feature branch: `git checkout -b feature/my-feature`
3. Write code following style guide
4. Add tests for new functionality
5. Update documentation
6. Run checks: `devtools::check()`
7. Commit with conventional commit messages
8. Push and create pull request

### Code Review Criteria

- [ ] Follows code style guidelines
- [ ] Includes appropriate tests
- [ ] Documentation is complete
- [ ] No decrease in test coverage
- [ ] Passes R CMD check
- [ ] Commit messages are clear

## Useful Resources

### R Package Development

- [R Packages (2e)](https://r-pkgs.org/) by Hadley Wickham & Jennifer Bryan
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/R-exts.html)
- [testthat documentation](https://testthat.r-lib.org/)
- [roxygen2 documentation](https://roxygen2.r-lib.org/)

### Statistical Methods

- [netmeta package](https://cran.r-project.org/package=netmeta)
- [gemtc package](https://cran.r-project.org/package=gemtc)
- [Cochrane Handbook](https://training.cochrane.org/handbook)

### Development Tools

```r
# Essential development packages
install.packages(c(
  "devtools",      # Package development
  "roxygen2",      # Documentation
  "testthat",      # Testing
  "usethis",       # Package setup
  "pkgdown",       # Website generation
  "covr",          # Test coverage
  "lintr",         # Code linting
  "goodpractice",  # Best practices
  "profvis",       # Profiling
  "bench"          # Benchmarking
))
```

## Contact

- **Maintainer**: research@example.com
- **Issues**: https://github.com/your-org/powerNMA/issues
- **Discussions**: https://github.com/your-org/powerNMA/discussions
