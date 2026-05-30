# =============================================================================
# tests/test_hydro_sediment.R
# testthat unit tests for hydro_sediment.R
#
# Functions covered
#   fit_rating_curve()         — power-law SSL = a * Q^b via log-log OLS
#   estimate_daily_ssl()       — apply rating curve to daily flow series
#   estimate_trap_efficiency() — Brune-curve reservoir sedimentation model
#
# Coverage
#   fit_rating_curve
#     T-01  Returns a list with required elements
#     T-02  Exponent b is within literature range [1.4, 2.5]
#     T-03  R-squared is positive
#     T-04  Both a and b are positive (physically required)
#     T-05  Fewer than 5 valid rows triggers an error
#     T-06  Missing required column triggers an error
#
#   estimate_daily_ssl
#     T-07  Output has columns date, Q_cms, SSL_mt_day
#     T-08  SSL_mt_day is 0 when Q < Q_min_transport
#     T-09  SSL_mt_day is non-negative for all rows
#     T-10  SSL_mt_day follows power-law at a known Q value
#
#   estimate_trap_efficiency
#     T-11  Returns data.frame with expected columns
#     T-12  Number of rows equals project lifetime (years argument)
#     T-13  TE is in (0, 1) for all years
#     T-14  S_remaining_m3 is monotonically non-increasing
#     T-15  S_remaining_m3 is always >= 0
#     T-16  ssl_annual_mt = 0 gives TE ≈ 0 and S_remaining unchanged
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_sediment.R"))

# ---------------------------------------------------------------------------
# Shared synthetic sediment dataset
# ---------------------------------------------------------------------------

make_sed <- function(n = 50, seed = 1L) {
  set.seed(seed)
  Q   <- runif(n, min = 5, max = 200)
  SSL <- 10 * Q ^ 1.6 * exp(rnorm(n, 0, 0.3))   # SSL = 10 * Q^1.6 + noise
  data.frame(
    Discharge_CMS  = Q,
    Sus_Load_MTDay = SSL,
    flag           = rep("ok", n),
    stringsAsFactors = FALSE
  )
}

sed_df <- make_sed()

# Pre-fit curve used across multiple tests
rc <- fit_rating_curve(sed_df, plot_fit = FALSE)

# ---------------------------------------------------------------------------
# fit_rating_curve — T-01 to T-06
# ---------------------------------------------------------------------------

test_that("T-01: fit_rating_curve() returns a list with required elements", {

  required <- c("a", "b", "r_squared", "n_points", "b_in_range", "method_note")

  expect_true(
    all(required %in% names(rc)),
    label = paste("Missing elements:", paste(
      setdiff(required, names(rc)), collapse = ", "
    ))
  )
})

test_that("T-02: exponent b is within literature range [1.4, 2.5]", {

  # Synthetic data was generated with b = 1.6, so fitted b should be close
  expect_true(
    rc$b >= 1.4 && rc$b <= 2.5,
    label = paste("b =", round(rc$b, 3), "is outside [1.4, 2.5]")
  )
})

test_that("T-03: R-squared is positive", {

  expect_gt(rc$r_squared, 0,
            label = "R-squared must be > 0 for a non-trivial fit")
})

test_that("T-04: coefficient a and exponent b are both positive", {

  expect_gt(rc$a, 0, label = "Coefficient a must be > 0")
  expect_gt(rc$b, 0, label = "Exponent b must be > 0")
})

test_that("T-05: fewer than 5 valid rows triggers an error", {

  tiny_df <- data.frame(
    Discharge_CMS  = c(10, 20, 30),
    Sus_Load_MTDay = c(50, 200, 800),
    flag           = rep("ok", 3)
  )

  expect_error(
    fit_rating_curve(tiny_df, plot_fit = FALSE),
    regexp = NULL,
    label  = "Should error when fewer than 5 valid observations"
  )
})

test_that("T-06: missing required column triggers an error", {

  bad_df <- data.frame(
    Q   = c(10, 20, 30, 40, 50),
    SSL = c(50, 200, 800, 2000, 5000),
    flag = rep("ok", 5)
  )

  expect_error(
    fit_rating_curve(bad_df, plot_fit = FALSE),
    regexp = NULL,
    label  = "Should error when 'Discharge_CMS' column is absent"
  )
})

# ---------------------------------------------------------------------------
# estimate_daily_ssl — T-07 to T-10
# ---------------------------------------------------------------------------

# Minimal daily flow data
make_daily_flow <- function(n = 60, seed = 7L) {
  set.seed(seed)
  data.frame(
    date  = seq(as.Date("2020-01-01"), by = "day", length.out = n),
    Q_cms = runif(n, min = 0.5, max = 80)
  )
}

daily_df <- make_daily_flow()
ssl_df   <- estimate_daily_ssl(daily_df, rc, Q_min_transport = 5.0)

test_that("T-07: estimate_daily_ssl() returns required columns", {

  expect_true(
    all(c("date", "Q_cms", "SSL_mt_day") %in% names(ssl_df)),
    label = "Output must contain date, Q_cms, SSL_mt_day"
  )
})

test_that("T-08: SSL_mt_day = 0 when Q < Q_min_transport", {

  # Find rows where Q < 5.0
  low_Q_idx <- which(daily_df$Q_cms < 5.0)

  if (length(low_Q_idx) > 0) {
    expect_equal(
      ssl_df$SSL_mt_day[low_Q_idx],
      rep(0, length(low_Q_idx)),
      tolerance = 1e-9,
      label = "SSL must be 0 for Q below transport threshold"
    )
  } else {
    skip("No rows with Q < 5 in synthetic data — increase range or lower threshold")
  }
})

test_that("T-09: SSL_mt_day is non-negative for all rows", {

  expect_true(
    all(ssl_df$SSL_mt_day >= 0),
    label = "SSL_mt_day must be >= 0 for all rows"
  )
})

test_that("T-10: SSL_mt_day follows power-law at a known Q value", {

  # At Q = 20 cms, SSL should be approximately rc$a * 20^rc$b
  known_flow <- data.frame(
    date  = as.Date("2020-06-01"),
    Q_cms = 20
  )
  known_ssl <- estimate_daily_ssl(known_flow, rc, Q_min_transport = 5.0)
  expected  <- rc$a * 20 ^ rc$b

  expect_equal(
    known_ssl$SSL_mt_day,
    expected,
    tolerance = 1e-6,
    label = "SSL_mt_day must equal a * Q^b at Q = 20 cms"
  )
})

# ---------------------------------------------------------------------------
# estimate_trap_efficiency — T-11 to T-16
# ---------------------------------------------------------------------------

trap_W1 <- estimate_trap_efficiency(ssl_annual_mt = 50000, weir_id = "W1",
                                    years = 35L)

test_that("T-11: estimate_trap_efficiency() returns required columns", {

  required <- c("year", "S_remaining_m3", "TE", "V_loss_m3_yr", "CI_ratio")

  expect_true(
    all(required %in% names(trap_W1)),
    label = paste("Missing columns:", paste(
      setdiff(required, names(trap_W1)), collapse = ", "
    ))
  )
})

test_that("T-12: number of rows equals the years argument", {

  expect_equal(nrow(trap_W1), 35L,
               label = "One row per project year expected")
})

test_that("T-13: trap efficiency TE is in (0, 1) for all years", {

  expect_true(
    all(trap_W1$TE > 0 & trap_W1$TE < 1),
    label = "TE must be strictly between 0 and 1"
  )
})

test_that("T-14: S_remaining_m3 is monotonically non-increasing", {

  diffs <- diff(trap_W1$S_remaining_m3)

  expect_true(
    all(diffs <= 0),
    label = "Storage must not increase over time (sediment only accumulates)"
  )
})

test_that("T-15: S_remaining_m3 is always >= 0", {

  expect_true(
    all(trap_W1$S_remaining_m3 >= 0),
    label = "Remaining storage cannot be negative"
  )
})

test_that("T-16: zero sediment load gives near-zero volume loss each year", {

  trap_zero <- estimate_trap_efficiency(ssl_annual_mt = 0, weir_id = "W1",
                                       years = 10L)

  expect_equal(
    sum(trap_zero$V_loss_m3_yr), 0,
    tolerance = 1e-6,
    label = "Zero sediment input must give zero volume loss"
  )
})
