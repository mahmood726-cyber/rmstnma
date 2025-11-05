# ============================================================================
# Package Initialization
# ============================================================================

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("powerNMA")

  packageStartupMessage(
    "\n",
    "═══════════════════════════════════════════════════════════\n",
    "  powerNMA v", version, " - Network Meta-Analysis Suite\n",
    "═══════════════════════════════════════════════════════════\n",
    "\n",
    "📊 Quick Start:\n",
    "  • data <- simulate_nma_data(n_studies = 20)\n",
    "  • results <- run_powernma(data, data_type = 'pairwise')\n",
    "\n",
    "📖 Documentation:\n",
    "  • ?run_powernma     - Main analysis function\n",
    "  • ?diagnose_nma_data - Data quality diagnostics\n",
    "  • browseVignettes('powerNMA') - Tutorials\n",
    "\n",
    "⚠️  Mode Guide:\n",
    "  • mode='standard' - Validated methods (default, use for publications)\n",
    "  • mode='experimental' - Novel methods (research use only)\n",
    "\n",
    "💡 Pro Tips:\n",
    "  • Always run diagnose_nma_data() before analysis\n",
    "  • Use validate_nma_input() to check data format\n",
    "  • See ?setup_powernma for configuration options\n",
    "\n",
    "══════════════════════════════════════════════════════════\n"
  )
}

.onLoad <- function(libname, pkgname) {
  # Register S3 methods if needed
  invisible()
}
