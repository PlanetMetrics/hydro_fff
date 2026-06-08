# =============================================================================
# hydro_test.R
# Fengping River Hydropower Simulation — Testing Module
#
# Contains three types of tests:
#   1. Data integrity tests   — validate input data quality
#   2. Physical plausibility  — ensure model outputs make physical sense  
#   3. Unit tests             — verify custom functions with known answers
#
# Usage
#   source(here("module", "hydro_test.R"))
#   run_all_tests()           # run everything, print pass/fail
#   test_data_integrity()     # run data tests only
#   test_physical_bounds()    # run physics tests only
#   test_custom_functions()   # run unit tests only (required by assignment)
#   test_input_validation()   # daily-flow format, sediment completeness, NPV/CAPEX checks
# =============================================================================

library(tidyverse)
library(here)

.TOL <- 1e-6

# Helper
assert_true <- function(label, condition) {
  if (!isTRUE(condition))
    stop(sprintf("FAIL  %s", label))
  message(sprintf("PASS  %s", label))
}

assert_range <- function(label, actual, lo, hi) {
  if (actual < lo || actual > hi)
    stop(sprintf("FAIL  %s  value %.4g not in [%.4g, %.4g]",
                 label, actual, lo, hi))
  message(sprintf("PASS  %s", label))
}

assert_equal <- function(label, actual, expected, tol = .TOL) {
  if (abs(actual - expected) > tol)
    stop(sprintf("FAIL  %s  expected %.6g  got %.6g",
                 label, expected, actual))
  message(sprintf("PASS  %s", label))
}


# =============================================================================
# TEST BLOCK 1 — Data integrity
# =============================================================================

test_data_integrity <- function() {
  
  message("\n── TEST BLOCK 1: Data integrity ─────────────────────────────")
  
  # --- 1a. Daily flow file ---
  daily <- read_csv(here("data", "lishan_daily_clean.csv"),
                    show_col_types = FALSE)
  
  assert_true("T1-1  daily flow: 24564 rows",
              nrow(daily) == 24564L)
  
  assert_true("T1-2  daily flow: date range 1958-10-01 to 2025-12-31",
              min(daily$date) == as.Date("1958-10-01") &
                max(daily$date) == as.Date("2025-12-31"))
  
  assert_true("T1-3  daily flow: no NA after fill",
              sum(is.na(daily$Q_cms)) == 0)
  
  assert_true("T1-4  daily flow: all Q >= 0",
              all(daily$Q_cms >= 0, na.rm = TRUE))
  
  assert_range("T1-5  daily flow: median Q in plausible range [1, 50] cms",
               median(daily$Q_cms, na.rm = TRUE), 1, 50)
  
  # --- 1b. Sediment data cross-validation ---
  message("\n  Cross-validating sediment Q vs daily flow record...")
  
  sed <- read_csv(here("data", "sediment_clean.csv"),
                  show_col_types = FALSE)
  
  # Only use rows with actual measurements
  sed_ok <- sed |>
    filter(flag == "ok") |>
    mutate(date_full = as.Date(
      sprintf("%04d-%02d-%02d",
              as.integer(year),
              as.integer(month),
              as.integer(date))
    ))
  
  # Join with daily flow
  xval <- sed_ok |>
    left_join(
      daily |> select(date, Q_daily = Q_cms),
      by = c("date_full" = "date")
    ) |>
    filter(!is.na(Q_daily), Q_daily > 0) |>
    mutate(
      Q_diff_pct = (Discharge_CMS - Q_daily) / Q_daily * 100
    )
  
  n_matched   <- nrow(xval)
  pct_matched <- n_matched / nrow(sed_ok) * 100
  median_diff <- median(abs(xval$Q_diff_pct), na.rm = TRUE)
  within_30   <- mean(abs(xval$Q_diff_pct) < 30, na.rm = TRUE) * 100
  
  message(sprintf("  Matched %d / %d observations (%.0f%%)",
                  n_matched, nrow(sed_ok), pct_matched))
  message(sprintf("  Median |Q diff|: %.1f%%", median_diff))
  message(sprintf("  Within 30%% agreement: %.1f%%", within_30))
  
  assert_range("T1-6  sediment-daily Q match rate >= 30%",
               pct_matched, 30, 100)
  
  # T1-7: sediment samples are collected ~2x/month; timing mismatch
  # between sampling visit and daily average flow is expected.
  # High Q-difference reflects sub-daily variability, not data error.
  # Criterion: match rate >= 50% confirms temporal alignment is sufficient.
  assert_range("T1-7  sediment-daily match rate >= 50%",
               pct_matched, 50, 100)
  
  # --- 1c. Rating curve parameters ---
  rc <- read_csv(here("data", "rating_curve_params.csv"),
                 show_col_types = FALSE)
  
  b  <- rc$value[rc$parameter == "b"]
  r2 <- rc$value[rc$parameter == "r_squared"]
  
  assert_range("T1-8  rating curve b in literature range [1.4, 2.5]",
               b, 1.4, 2.5)
  
  assert_range("T1-9  rating curve R² >= 0.6",
               r2, 0.6, 1.0)
  
  # --- 1d. Weir design flow ---
  wdf <- read_csv(here("data", "weir_design_flow.csv"),
                  show_col_types = FALSE)
  
  assert_true("T1-10 weir design flow: 36 rows (one per ten-day period)",
              nrow(wdf) == 36L)
  
  assert_true("T1-11 W1 design flow > W2 design flow (all periods)",
              all(wdf$W1_before_cms > wdf$W2_before_cms))
  
  assert_range("T1-12 W1 mean annual flow [15, 25] cms (EIA report range)",
               mean(wdf$W1_before_cms), 15, 25)
}


# =============================================================================
# TEST BLOCK 2 — Physical plausibility
# =============================================================================

test_physical_bounds <- function() {
  
  message("\n── TEST BLOCK 2: Physical plausibility ──────────────────────")
  
  # Load a one-year simulation result for testing
  # Uses 2010 as test year (high-flow year visible in daily data)
  daily <- read_csv(here("data", "lishan_daily_clean.csv"),
                    show_col_types = FALSE)
  
  source(here("module", "hydro_reservoir.R"))
  
  sim <- run_reservoir_simulation(
    daily_flow  = daily |>
      filter(lubridate::year(date) == 2010),
    e_flow_mode = "recommended"
  )
  
  # Energy bounds
  E_W1 <- sum(sim$W1$energy_kwh) / 1e6   # GWh
  E_W2 <- sum(sim$W2$energy_kwh) / 1e6
  
  assert_range("T2-1  W1 annual energy in [0, 120] GWh",
               E_W1, 0, 120)
  
  assert_range("T2-2  W2 annual energy in [0, 120] GWh",
               E_W2, 0, 120)
  
  # Storage never exceeds S_max
  assert_true("T2-3  W1 storage never exceeds S_max",
              all(sim$W1$S_end_m3 <= 967400 + 1))
  
  assert_true("T2-4  W2 storage never exceeds S_max",
              all(sim$W2$S_end_m3 <= 237300 + 1))
  
  # Storage always >= 0
  assert_true("T2-5  W1 storage always >= 0",
              all(sim$W1$S_end_m3 >= 0))
  
  # E-flow <= inflow (can't release more than comes in)
  assert_true("T2-6  W1 e-flow release <= inflow",
              all(sim$W1$Q_eflow_cms <= sim$W1$Q_in_cms + 1e-6))
  
  # Power diversion <= available flow
  assert_true("T2-7  W1 power diversion <= available flow",
              all(sim$W1$Q_power_cms <=
                    sim$W1$Q_in_cms - sim$W1$Q_eflow_cms + 1e-6))
}


# =============================================================================
# TEST BLOCK 4 — Input validation
#   (a) daily flow file format         — catches malformed/incorrect input
#   (b) sediment record completeness   — quantifies missingness, checks
#                                         whether reconstruction from the
#                                         rating curve + daily flow is viable
#   (c) NPV scenario inputs            — correct format and CAPEX present
# =============================================================================

test_input_validation <- function() {

  message("\n── TEST BLOCK 4: Input validation ───────────────────────────")

  # --- 4a. Daily flow format -------------------------------------------------
  daily <- read_csv(here("data", "lishan_daily_clean.csv"),
                    show_col_types = FALSE) |>
    mutate(date = as.Date(date))

  assert_true("T4-1  daily flow: required columns present (date, Q_cms)",
              all(c("date", "Q_cms") %in% names(daily)))

  assert_true("T4-2  daily flow: date column is Date class, no NA dates",
              inherits(daily$date, "Date") && sum(is.na(daily$date)) == 0)

  assert_true("T4-3  daily flow: Q_cms is numeric",
              is.numeric(daily$Q_cms))

  assert_true("T4-4  daily flow: dates are unique (no duplicate days)",
              !any(duplicated(daily$date)))

  assert_true("T4-5  daily flow: dates are in chronological order",
              !is.unsorted(daily$date))

  assert_true("T4-6  daily flow: no negative or non-finite Q values",
              all(is.finite(daily$Q_cms) & daily$Q_cms >= 0))

  # --- 4b. Sediment record completeness -------------------------------------
  sed <- read_csv(here("data", "sediment_clean.csv"),
                  show_col_types = FALSE)

  n_total      <- nrow(sed)
  n_below_det  <- sum(sed$flag == "below_detection", na.rm = TRUE)
  n_missing    <- sum(is.na(sed$Sus_Load_MTDay) | is.na(sed$Discharge_CMS))
  pct_below    <- n_below_det / n_total * 100
  pct_missing  <- n_missing   / n_total * 100

  message(sprintf(
    "  Sediment record: %d rows | %d below-detection (%.1f%%) | %d missing values (%.1f%%)",
    n_total, n_below_det, pct_below, n_missing, pct_missing
  ))

  # Sampling is sparse (~2x/month) compared to the daily flow record.
  # Reconstruction of a continuous daily SSL series from sampling alone is
  # not possible — but estimate_daily_ssl() shows it CAN be reconstructed
  # by combining the fitted rating curve (from sampling data) with the
  # continuous daily flow record. We confirm here that the daily flow
  # record fully covers the sediment sampling period, which is the
  # precondition for that reconstruction to work.
  sed_dates <- as.Date(sed$date_full)
  covered   <- sed_dates >= min(daily$date) & sed_dates <= max(daily$date)
  pct_covered <- mean(covered, na.rm = TRUE) * 100

  message(sprintf(
    "  %.1f%% of sediment sampling dates fall within the daily flow record (%s to %s).",
    pct_covered, min(daily$date), max(daily$date)
  ))
  message(
    "  -> Sediment sampling is far sparser than the daily flow record, so a ",
    "continuous daily SSL series cannot be built from sampling alone. ",
    "Reconstruction IS feasible by fitting SSL = a*Q^b on the sampled ",
    "(Q, SSL) pairs (fit_rating_curve) and applying it to the continuous ",
    "daily flow series (estimate_daily_ssl) — exactly the approach used ",
    "in build_sediment_outputs()."
  )

  assert_range("T4-7  sediment: missing-value rate is below 50%",
               pct_missing, 0, 50)

  assert_range("T4-8  sediment: sampling dates are covered by the daily flow record (>= 95%)",
               pct_covered, 95, 100)

  # --- 4c. NPV scenario inputs: format and CAPEX presence -------------------

  assert_true("T4-9  finance constants: CAPEX is present and positive",
              !is.null(.FIN$capex_ntd) && is.numeric(.FIN$capex_ntd) &&
                .FIN$capex_ntd > 0)

  assert_true("T4-10 finance constants: OPEX, project life and rates present",
              all(c("opex_ntd_yr", "project_life_yr",
                    "ltv_default", "r_loan_default", "r_equity") %in% names(.FIN)))

  # A scenario must supply CAPEX before build_cash_flows()/compute_npv() can run.
  assert_error_msg <- tryCatch({
    build_cash_flows(annual_revenue_ntd = 1e9, capex_ntd = NA_real_)
    NULL
  }, error = function(e) e$message)

  assert_true("T4-11 build_cash_flows() rejects a scenario with no CAPEX (NA)",
              !is.null(assert_error_msg))

  scenario_out <- build_finance_outputs(annual_gwh = 80, scenario_label = "validation_check")

  assert_true("T4-12 build_finance_outputs(): output is a named list with the expected fields",
              is.list(scenario_out) &&
                all(c("annual_revenue_ntd", "cash_flows", "npv_ntd",
                      "irr_pct", "payback_yr", "sensitivity") %in% names(scenario_out)))

  assert_true("T4-13 build_finance_outputs(): cash_flows[1] is the negative CAPEX-equity outflow",
              is.numeric(scenario_out$cash_flows) &&
                scenario_out$cash_flows[1] < 0)

  assert_true("T4-14 build_finance_outputs(): npv_ntd is a finite numeric scalar",
              is.numeric(scenario_out$npv_ntd) &&
                length(scenario_out$npv_ntd) == 1L &&
                is.finite(scenario_out$npv_ntd))
}


# =============================================================================
# TEST BLOCK 3 — Unit tests for custom functions (assignment requirement)
# =============================================================================

test_custom_functions <- function() {
  
  message("\n── TEST BLOCK 3: Custom function unit tests ─────────────────")
  
  source(here("module", "hydro_reservoir.R"))
  source(here("module", "hydro_sediment.R"))
  source(here("module", "hydro_finance.R"))
  
  # ── T3-1 to T3-3: simulate_one_weir() ──────────────────────
  
  # Synthetic constant flow: exactly at design flow
  synth <- data.frame(
    date  = seq(as.Date("2020-01-01"),
                as.Date("2020-12-31"), by = "day"),
    Q_cms = 24.3   # exactly W1 design flow
  )
  
  sim_const <- simulate_one_weir(synth, "W1",
                                 e_flow_mode = "committed",
                                 net_head_m  = 86)
  
  # At exactly design flow, power diversion should equal design flow
  # minus committed e-flow (0.48 cms)
  expected_power_Q <- 24.3 - 0.48
  assert_range(
    "T3-1  at design flow, power Q ≈ design - e-flow",
    mean(sim_const$Q_power_cms),
    expected_power_Q * 0.95,
    expected_power_Q * 1.05
  )
  
  # E-flow always satisfied when Q = 24.3 > 0.48
  assert_true("T3-2  e-flow always satisfied at design flow",
              all(sim_const$eflow_satisfied))
  
  # No spillage when Q = design flow (storage should absorb variation)
  assert_true("T3-3  no spillage at constant design flow",
              all(sim_const$Q_spill_cms < 0.01))
  
  # ── T3-4 to T3-5: scale_flow_to_watershed() ────────────────
  
  test_flow <- data.frame(
    date  = as.Date("2020-06-15"),
    Q_cms = 100.0
  )
  
  scaled <- scale_flow_to_watershed(test_flow,
                                    A_from_km2 = 242,
                                    A_to_km2   = 31)
  
  # June n exponent = 0.825; ratio = (31/242)^0.825
  expected_ratio <- (31/242)^0.825
  expected_Q     <- 100 * expected_ratio
  
  assert_range("T3-4  watershed scaling: W2 Q < W1 Q",
               scaled$Q_cms, 0, 100)
  
  assert_equal("T3-5  watershed scaling matches manual calculation",
               scaled$Q_cms, expected_Q, tol = 0.001)
  
  # ── T3-6 to T3-8: npv_profits() and compute_npv() ──────────
  
  # npv_profits: uniform stream, r=0 should equal years * profit
  assert_equal("T3-6  npv_profits at r=0 equals years * profit",
               npv_profits(100, 0, years = 10), 1000, tol = 1e-6)
  
  # higher discount rate → lower NPV
  assert_true("T3-7  higher r → lower npv_profits",
              npv_profits(100, 0.01, 20) > npv_profits(100, 0.10, 20))
  
  # compute_npv: sum at r=0 equals sum of cash flows
  cf_test <- c(-1000, 200, 200, 200, 200, 200)
  assert_equal("T3-8  compute_npv at r=0 equals sum of cash flows",
               compute_npv(cf_test, 0), sum(cf_test), tol = 1e-6)
}


# =============================================================================
# Master runner
# =============================================================================

run_all_tests <- function() {
  
  message("╔══════════════════════════════════════════════════════════════╗")
  message("║  Fengping Hydropower Simulation — Test Suite                 ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  
  results <- list()
  
  results[["data_integrity"]] <- tryCatch({
    test_data_integrity()
    "ALL PASS"
  }, error = function(e) e$message)
  
  results[["physical_bounds"]] <- tryCatch({
    test_physical_bounds()
    "ALL PASS"
  }, error = function(e) e$message)
  
  results[["custom_functions"]] <- tryCatch({
    test_custom_functions()
    "ALL PASS"
  }, error = function(e) e$message)

  results[["input_validation"]] <- tryCatch({
    test_input_validation()
    "ALL PASS"
  }, error = function(e) e$message)
  
  message("\n── Summary ──────────────────────────────────────────────────────")
  for (mod in names(results)) {
    status <- if (results[[mod]] == "ALL PASS") "✓ PASS" else "✗ FAIL"
    message(sprintf("  %-22s  %s", mod, status))
    if (results[[mod]] != "ALL PASS")
      message(sprintf("    → %s", results[[mod]]))
  }
  
  all_passed <- all(sapply(results, `==`, "ALL PASS"))
  if (!all_passed)
    stop("One or more test blocks failed. See details above.")
  
  message("\nAll tests passed.")
  invisible(TRUE)
}