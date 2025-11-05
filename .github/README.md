# CI/CD Configuration for powerNMA

This directory contains GitHub Actions workflows and configuration for automated testing, building, and deployment of the powerNMA R package.

## 🚀 Workflows

### R-CMD-check.yaml
**Purpose**: Comprehensive package checking across multiple R versions and operating systems

**Triggers**:
- Push to `main`, `master`, `develop`, or `claude/**` branches
- Pull requests to `main`, `master`, `develop`
- Daily scheduled runs at 3 AM UTC

**Test Matrix**:
- ✅ macOS (R release)
- ✅ Windows (R release)
- ✅ Ubuntu (R devel, release, oldrel-1)

**What it does**:
1. Runs `R CMD check` on the package
2. Executes all testthat tests
3. Checks for warnings and errors
4. Uploads check results if failures occur

**Status Badge**:
```markdown
[![R-CMD-check](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/R-CMD-check.yaml)
```

---

### test-coverage.yaml
**Purpose**: Code coverage analysis and reporting

**Triggers**:
- Push to `main`, `master`, `develop`
- Pull requests to `main`, `master`, `develop`

**What it does**:
1. Runs test suite with coverage tracking
2. Uploads coverage report to Codecov
3. Comments on PRs with coverage changes

**Status Badge**:
```markdown
[![codecov](https://codecov.io/gh/mahmood726-cyber/rmstnma/branch/main/graph/badge.svg)](https://codecov.io/gh/mahmood726-cyber/rmstnma)
```

---

### pkgdown.yaml
**Purpose**: Automated documentation website generation

**Triggers**:
- Push to `main`, `master`
- Pull requests to `main`, `master`
- Releases
- Manual dispatch

**What it does**:
1. Builds pkgdown website from package documentation
2. Deploys to GitHub Pages (`gh-pages` branch)
3. Makes documentation available at: `https://mahmood726-cyber.github.io/rmstnma/`

**How to view**:
Visit the documentation site after enabling GitHub Pages in repository settings.

---

### pr-commands.yaml
**Purpose**: Helpful commands for pull request maintainers

**Triggers**:
- Comments on pull requests by MEMBERS or OWNERS

**Available Commands**:
- `/document` - Regenerate documentation with `roxygen2::roxygenise()`
- `/style` - Auto-format code with `styler::style_pkg()`

**Usage**:
Comment `/document` or `/style` on any PR to trigger the respective action.

---

### lint.yaml
**Purpose**: Code style and quality checking

**Triggers**:
- Push to `main`, `master`, `develop`
- Pull requests to `main`, `master`, `develop`

**What it does**:
1. Runs `lintr::lint_package()` on all R code
2. Checks for style violations
3. Reports issues in workflow logs

**Customization**:
Create `.lintr` file in `powerNMA/` directory to customize linting rules.

---

### release.yaml
**Purpose**: Automated release creation and package distribution

**Triggers**:
- Push of version tags (e.g., `v1.0.0`, `v2.1.3`)

**What it does**:
1. Builds source package (`.tar.gz`)
2. Creates GitHub Release with changelog
3. Uploads package as release asset
4. Provides installation instructions

**Usage**:
```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```

---

### benchmark.yaml
**Purpose**: Performance monitoring and regression detection

**Triggers**:
- Weekly schedule (Sundays at 2 AM UTC)
- Push to `main`/`master` affecting R code
- Manual dispatch

**What it does**:
1. Runs performance benchmarks
2. Stores results as artifacts (90-day retention)
3. Comments benchmark summary on commits

**Custom Benchmarks**:
Add custom benchmark scripts to `powerNMA/tests/benchmark/run_benchmarks.R`

---

## ⚙️ Configuration Files

### dependabot.yml
**Purpose**: Automatic dependency updates

**What it monitors**:
- GitHub Actions versions (weekly updates)
- Workflow dependencies

**How it works**:
- Creates automated PRs for dependency updates
- Labels PRs appropriately
- Limits to 10 open PRs at a time

---

## 📊 Status Badges

Add these badges to your main README.md:

```markdown
# powerNMA

[![R-CMD-check](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/mahmood726-cyber/rmstnma/branch/main/graph/badge.svg)](https://codecov.io/gh/mahmood726-cyber/rmstnma)
[![lint](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/lint.yaml/badge.svg)](https://github.com/mahmood726-cyber/rmstnma/actions/workflows/lint.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
```

---

## 🛠️ Setup Instructions

### 1. Enable GitHub Actions
GitHub Actions should be enabled by default. Check:
- Repository → Settings → Actions → Allow all actions

### 2. Enable GitHub Pages (for pkgdown)
- Repository → Settings → Pages
- Source: Deploy from a branch
- Branch: `gh-pages` / `root`

### 3. Add Codecov Integration (Optional)
1. Sign up at https://codecov.io
2. Add repository
3. No additional secrets needed (uses GITHUB_TOKEN)

### 4. Branch Protection (Recommended)
Set up branch protection rules for `main`/`master`:
- Require status checks to pass before merging
  - R-CMD-check
  - test-coverage
  - lint
- Require pull request reviews
- Require branches to be up to date before merging

---

## 🔍 Troubleshooting

### Workflow fails with "Package check failed"
1. Check the workflow logs for specific errors
2. Run `devtools::check("powerNMA")` locally to reproduce
3. Fix errors in the package code
4. Push fixes to trigger re-run

### Coverage report not uploading
1. Ensure repository is public or Codecov is configured for private repos
2. Check that `GITHUB_TOKEN` has appropriate permissions
3. Verify Codecov integration is active

### pkgdown site not deploying
1. Check that GitHub Pages is enabled
2. Ensure `gh-pages` branch exists
3. Verify workflow has `contents: write` permission
4. Check workflow logs for deployment errors

### Benchmark workflow fails
1. Create `powerNMA/tests/benchmark/` directory if missing
2. Ensure benchmark scripts don't require interactive input
3. Check that all dependencies are installed

---

## 📝 Best Practices

### For Contributors
1. **Always run checks locally before pushing**:
   ```r
   devtools::check("powerNMA")
   devtools::test("powerNMA")
   ```

2. **Write tests for new features**:
   - Add tests in `powerNMA/tests/testthat/`
   - Aim for >80% code coverage

3. **Document all functions**:
   - Use roxygen2 comments
   - Include examples
   - Run `/document` command on PRs

4. **Follow code style**:
   - Use tidyverse style guide
   - Run `/style` command on PRs
   - Fix lintr warnings

### For Maintainers
1. **Review CI results before merging PRs**
2. **Create releases with semantic versioning** (MAJOR.MINOR.PATCH)
3. **Monitor benchmark results** for performance regressions
4. **Keep dependencies updated** (review Dependabot PRs)
5. **Maintain changelog** in `NEWS.md`

---

## 🔄 Continuous Integration Flow

```
┌─────────────────┐
│  Code Push/PR   │
└────────┬────────┘
         │
         ├──────────────────────────────────────┐
         │                                      │
         ▼                                      ▼
┌─────────────────┐                   ┌──────────────┐
│  R-CMD-check    │                   │    lint      │
│  (Multi-OS/R)   │                   │  (Style)     │
└────────┬────────┘                   └──────┬───────┘
         │                                    │
         ▼                                    ▼
┌─────────────────┐                   ┌──────────────┐
│ test-coverage   │                   │   PASS?      │
│   (Codecov)     │                   │              │
└────────┬────────┘                   └──────┬───────┘
         │                                    │
         │                    NO              │  YES
         ├────────────────────────────────────┤
         │                                    │
         ▼                                    ▼
┌─────────────────┐                   ┌──────────────┐
│  Fix Issues &   │                   │  Merge PR    │
│  Re-push        │                   │              │
└─────────────────┘                   └──────┬───────┘
                                              │
                                              ▼
                                      ┌──────────────┐
                                      │  pkgdown     │
                                      │  (Docs)      │
                                      └──────────────┘
```

---

## 📚 Additional Resources

- [GitHub Actions for R](https://github.com/r-lib/actions)
- [usethis::use_github_action()](https://usethis.r-lib.org/reference/github_actions.html)
- [R Packages Book - CI/CD](https://r-pkgs.org/software-development-practices.html#sec-sw-dev-practices-ci)
- [pkgdown documentation](https://pkgdown.r-lib.org/)
- [Codecov for R](https://docs.codecov.com/docs/r)

---

**Last Updated**: 2025-11-05
**Maintained by**: powerNMA Development Team
