## Test environments
* local: Windows 11, R 4.5.1
* win-builder: R-devel (2025-09-17 r88852)
* win-builder: R-release (4.4.2)
* win-builder: R-oldrelease (4.3.3)

## R CMD check results
0 errors | 0 warnings | 1 note

### NOTE:
* New submission
* Possibly misspelled words in DESCRIPTION:
  - Kaplan (11:39) - This is correct: Kaplan-Meier is a standard statistical term
  - RMST (9:20) - This is correct: Restricted Mean Survival Time is the core concept of the package

## Additional information
* Package uses optional dependency 'cmdstanr' from Additional_repositories (https://mc-stan.org/r-packages/)
* SystemRequirements: CmdStan (>= 2.33.0) is optional - package functions check for its availability
* All examples requiring Stan/cmdstanr are wrapped in \dontrun{} to prevent CRAN check timeouts
* The package has been tested without cmdstanr installed and handles its absence gracefully

## Reverse dependencies
This is a new release, so there are no reverse dependencies.

## Notes for CRAN reviewers
* The terms "Kaplan-Meier" and "RMST" flagged as possibly misspelled are standard statistical terminology
* cmdstanr is correctly specified in Additional_repositories as per CRAN policy
* The package provides meaningful functionality even without the optional Stan dependencies
