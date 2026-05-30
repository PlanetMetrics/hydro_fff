# =============================================================================
# tests/test_hydro_climate.R
# testthat unit tests for hydro_climate.R
#
# Functions covered
#   compute_rbi()               — Richards-Baker Flashiness Index
#   compute_annual_stats()      — year-by-year summary statistics
#   scale_to_cv()               — Layer 2 CV scaling
#   run_bootstrap_cv_grid()     — two-layer bootstrap over CV grid
#   summarise_cv_ensemble()     — percentile summary of grid output
#
# Coverage
#   compute_rbi
#     T-01  Constant flow gives RBI = 0
#     T-02  Highly variable flow gives RBI > 0
#     T-03  RBI is always in [0, 1]
#     T-04  Vector of length 1 triggers an error
#     T-05  All-NA vector returns NA (not an error)
#
#   compute_annual_stats
#     T-06  Returns one row per calendar year
#     T-07  cv_Q = sd_Q / mu_Q for every row
#     T-08  Missing date column triggers an error
#
#   scale_to_cv
#     T-09  Achieved CV is within 5% of cv_target
#     T-10  Annual mean is preserved within 1% after scaling
#     T-11  No values below q_floor after scaling
#     T-12  cv_target <= 0 triggers an error
#
#   run_bootstrap_cv_grid
#     T-13  Returns length(cv_grid) * n_sim rows
#     T-14  cv_achieved is always positive
#     T-15  annual_Q_m3 is consistent with mu_yr_cms
#     T-16  All required columns are present
#     T-17  n_sim < 10 triggers an error
#
#   summarise_cv_ensemble
#     T-18  Returns one row per cv_level
#     T-19  mu_P50 is within mu_P10 and mu_P90
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_climate.R"))

# ---------------------------------------------------------------------------
# Shared synthetic daily flow (3 full years, known statistics)
# ---------------------------------------------------------------------------

set.seed(2025L)
n_days <- 3 * 365
synth_daily <- data.frame(
  date  = seq(as.Date("2000-01-01"), by = "day", length.out = n_days),
  Q_cms = pmax(rnorm(n_days, mean = 15, sd = 8), 0.5)
)

# ---------------------------------------------------------------------------
# compute_rbi — T-01 to T-05
# ---------------------------------------------------------------------------

test_that("T-01: constant flow vector gives RBI = 0", {

  Q_const <- rep(10, 100)

  expect_equal(
    compute_rbi(Q_const), 0,
    tolerance = 1e-9,
    label = "RBI must be 0 for constant flow (no variability)"
  )
})

test_that("T-02: alternating high-low flow gives RBI > 0", {

  Q_alt <- rep(c(1, 100), 50)   # extreme alternation

  expect_gt(
    compute_rbi(Q_alt), 0,
    label = "RBI must be positive for variable flow"
  )
})

test_that("T-03: RBI is always in [0, 1]", {

  set.seed(99L)
  for (i in 1:10) {
    Q_rand <- pmax(rnorm(200, mean = 10, sd = 5), 0.01)
    rbi    <- compute_rbi(Q_rand)
    expect_true(
      rbi >= 0 && rbi <= 1,
      label = paste("RBI =", round(rbi, 4), "is outside [0, 1] on iteration", i)
    )
  }
})

test_that("T-04: vector of length 1 triggers an error", {

  expect_error(
    compute_rbi(10),
    regexp = NULL,
    label  = "Should error when Q_vec has fewer than 2 elements"
  )
})

test_that("T-05: all-NA vector returns NA without error", {

  Q_na <- rep(NA_real_, 50)

  expect_equal(
    compute_rbi(Q_na), NA_real_,
    label = "All-NA input should return NA (not throw an error)"
  )
})

# ---------------------------------------------------------------------------
# compute_annual_stats — T-06 to T-08
# ---------------------------------------------------------------------------

ann <- compute_annual_stats(synth_daily)

test_that("T-06: compute_annual_stats() returns one row per calendar year", {

  n_years_expected <- length(unique(lubridate::year(synth_daily$date)))

  expect_equal(
    nrow(ann), n_years_expected,
    label = "One row per year expected in annual_stats output"
  )
})

test_that("T-07: cv_Q equals sd_Q / mu_Q for every row", {

  cv_manual <- ann$sd_Q / ann$mu_Q

  expect_equal(
    ann$cv_Q, cv_manual,
    tolerance = 1e-9,
    label = "cv_Q must equal sd_Q / mu_Q by definition"
  )
})

test_that("T-08: missing 'date' column triggers an error", {

  bad_df <- data.frame(Q_cms = c(10, 20, 30))

  expect_error(
    compute_annual_stats(bad_df),
    regexp = NULL,
    label  = "Should error when 'date' column is absent"
  )
})

# ---------------------------------------------------------------------------
# scale_to_cv — T-09 to T-12
# ---------------------------------------------------------------------------

set.seed(77L)
Q_test   <- pmax(rnorm(365, mean = 12, sd = 4), 0.1)
cv_target <- 2.5
Q_scaled  <- scale_to_cv(Q_test, cv_target = cv_target)

test_that("T-09: achieved CV is within 5% of cv_target", {

  cv_out <- sd(Q_scaled) / mean(Q_scaled)

  expect_equal(
    cv_out, cv_target,
    tolerance = cv_target * 0.05,
    label = paste("Achieved CV =", round(cv_out, 3),
                  "must be within 5% of target", cv_target)
  )
})

test_that("T-10: annual mean is preserved within 1% after CV scaling", {

  mu_in  <- mean(Q_test)
  mu_out <- mean(Q_scaled)
  rel_err <- abs(mu_out - mu_in) / mu_in

  expect_lt(
    rel_err, 0.01,
    label = paste("Relative mean error =", round(rel_err * 100, 3),
                  "% must be < 1% (mean-preserving property)")
  )
})

test_that("T-11: no values fall below q_floor after scaling", {

  q_floor  <- 0.01
  Q_scaled2 <- scale_to_cv(Q_test, cv_target = 3.0, q_floor = q_floor)

  expect_true(
    all(Q_scaled2 >= q_floor),
    label = paste("All scaled values must be >= q_floor =", q_floor)
  )
})

test_that("T-12: cv_target <= 0 triggers an error", {

  expect_error(
    scale_to_cv(Q_test, cv_target = 0),
    regexp = NULL,
    label  = "cv_target = 0 is not meaningful and must error"
  )
})

# ---------------------------------------------------------------------------
# run_bootstrap_cv_grid — T-13 to T-17
# ---------------------------------------------------------------------------

grid_small <- run_bootstrap_cv_grid(
  daily_flow = synth_daily,
  cv_grid    = c(1.5, 2.5),
  n_sim      = 20L,
  seed       = 42L
)

test_that("T-13: run_bootstrap_cv_grid() returns length(cv_grid) * n_sim rows", {

  expect_equal(
    nrow(grid_small), 2L * 20L,
    label = "Total rows must equal length(cv_grid) * n_sim"
  )
})

test_that("T-14: cv_achieved is always strictly positive", {

  expect_true(
    all(grid_small$cv_achieved > 0),
    label = "cv_achieved must be positive for all rows"
  )
})

test_that("T-15: annual_Q_m3 is consistent with mu_yr_cms (within 0.1%)", {

  expected_vol <- grid_small$mu_yr_cms * 365 * 86400
  rel_err      <- abs(grid_small$annual_Q_m3 - expected_vol) / expected_vol

  expect_true(
    all(rel_err < 0.001),
    label = "annual_Q_m3 must equal mu_yr_cms * 365 * 86400 within 0.1%"
  )
})

test_that("T-16: all required output columns are present", {

  required <- c("cv_level", "sim_id", "year_src",
                "mu_yr_cms", "cv_achieved", "rbi", "annual_Q_m3")

  expect_true(
    all(required %in% names(grid_small)),
    label = paste("Missing columns:",
                  paste(setdiff(required, names(grid_small)), collapse = ", "))
  )
})

test_that("T-17: n_sim < 10 triggers an error", {

  expect_error(
    run_bootstrap_cv_grid(synth_daily, cv_grid = c(2.0), n_sim = 5L),
    regexp = NULL,
    label  = "n_sim = 5 is below minimum of 10 and must error"
  )
})

# ---------------------------------------------------------------------------
# summarise_cv_ensemble — T-18 to T-19
# ---------------------------------------------------------------------------

summ <- summarise_cv_ensemble(grid_small)

test_that("T-18: summarise_cv_ensemble() returns one row per cv_level", {

  expect_equal(
    nrow(summ), length(unique(grid_small$cv_level)),
    label = "One summary row per cv_level expected"
  )
})

test_that("T-19: mu_P50 is between mu_P10 and mu_P90 for all rows", {

  expect_true(
    all(summ$mu_P10 <= summ$mu_P50 & summ$mu_P50 <= summ$mu_P90),
    label = "mu_P50 must lie between mu_P10 and mu_P90 (ordered percentiles)"
  )
})
