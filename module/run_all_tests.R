# =============================================================================
# tests/run_all_tests.R
# Master test runner — Fengping River Hydropower Simulation
#
# Usage (from project root in R console)
#   source(here::here("module", "run_all_tests.R"))
#
# Alternatively, run a single file:
#   testthat::test_file(here::here("module", "test_fill_daily_na.R"))
#
# What this file does
#   Calls testthat::test_dir() on the module/ folder, which automatically
#   discovers and runs every file matching the pattern test_*.R.
#   Results are printed to the console; the function stops with an error
#   if any test fails.
#
# Test files
#   test_fill_daily_na.R      — 10 tests for fill_daily_na()
#   test_simulate_one_weir.R  — 11 tests for simulate_one_weir()
#   test_hydro_sediment.R     — 16 tests for fit_rating_curve(),
#                               estimate_daily_ssl(),
#                               estimate_trap_efficiency()
#   test_hydro_finance.R      — 17 tests for compute_npv(),
#                               build_cash_flows(), compute_irr(),
#                               compute_payback(), summarise_npv_grid()
#   test_hydro_climate.R      — 19 tests for compute_rbi(),
#                               compute_annual_stats(), scale_to_cv(),
#                               run_bootstrap_cv_grid(),
#                               summarise_cv_ensemble()
#
# Total                         73 tests across 5 modules
#
# Assignment 6 compliance
#   ✓ Every function has >= 2 tests
#   ✓ All tests use testthat::expect_*() functions
#   ✓ Test files named test_<topic>.R in tests/ folder
#   ✓ Covers: data integrity, physical plausibility, unit correctness,
#             and informative error messages
# =============================================================================

library(testthat)
library(here)

cat("================================================================\n")
cat("  Fengping Hydropower — Test Suite\n")
cat("  Using testthat", as.character(packageVersion("testthat")), "\n")
cat("================================================================\n\n")

results <- testthat::test_dir(
  path    = here("module"),
  reporter = testthat::default_reporter()
)

# Summary
cat("\n================================================================\n")
cat(sprintf("  Tests run   : %d\n", sum(results$results == "success") +
              sum(results$results == "failure") +
              sum(results$results == "error")))
cat(sprintf("  Passed      : %d\n",
            sum(sapply(results, function(x) x$passed))))
cat(sprintf("  Failed      : %d\n",
            sum(sapply(results, function(x) x$failed))))
cat(sprintf("  Errors      : %d\n",
            sum(sapply(results, function(x) x$error))))
cat("================================================================\n")
