# =============================================================================
# tests/test_estimate_daily_ssl.R
# testthat unit tests for estimate_daily_ssl()
#
# Function contract
#   Input  : daily_flow  data.frame(date = Date, Q_cms = numeric)
#            curve       list with elements $a, $b (from fit_rating_curve())
#            Q_min_transport  numeric, cms below which SSL = 0
#   Output : data.frame(date, Q_cms, SSL_mt_day)
#            SSL_mt_day = a * Q_cms^b  when Q_cms >= Q_min_transport, else 0
#
# This is the function that lets us "match" sparse sediment sampling data
# (used to fit a*Q^b via fit_rating_curve) to the continuous daily flow
# record, producing a daily suspended-sediment-load estimate.
#
# Coverage
#   T-01  Output has the required columns: date, Q_cms, SSL_mt_day
#   T-02  SSL_mt_day = 0 when Q_cms < Q_min_transport
#   T-03  SSL_mt_day is non-negative for all rows
#   T-04  SSL_mt_day follows the power law a * Q^b at a known Q value
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_sediment.R"))

# ---------------------------------------------------------------------------
# Shared synthetic data
# ---------------------------------------------------------------------------

# Synthetic sediment sampling data: SSL = 10 * Q^1.6 + noise
make_sed <- function(n = 50, seed = 1L) {
  set.seed(seed)
  Q   <- runif(n, min = 5, max = 200)
  SSL <- 10 * Q ^ 1.6 * exp(rnorm(n, 0, 0.3))
  data.frame(
    Discharge_CMS  = Q,
    Sus_Load_MTDay = SSL,
    flag           = rep("ok", n),
    stringsAsFactors = FALSE
  )
}

# Synthetic daily flow record
make_daily_flow <- function(n = 60, seed = 7L) {
  set.seed(seed)
  data.frame(
    date  = seq(as.Date("2020-01-01"), by = "day", length.out = n),
    Q_cms = runif(n, min = 0.5, max = 80)
  )
}

rc       <- fit_rating_curve(make_sed(), plot_fit = FALSE)
daily_df <- make_daily_flow()
ssl_df   <- estimate_daily_ssl(daily_df, rc, Q_min_transport = 5.0)

# ---------------------------------------------------------------------------
# T-01: required output columns
# ---------------------------------------------------------------------------

test_that("T-01: estimate_daily_ssl() returns required columns", {

  expect_true(
    all(c("date", "Q_cms", "SSL_mt_day") %in% names(ssl_df)),
    label = "Output must contain date, Q_cms, SSL_mt_day"
  )
})

# ---------------------------------------------------------------------------
# T-02: SSL = 0 below the transport threshold
# ---------------------------------------------------------------------------

test_that("T-02: SSL_mt_day = 0 when Q_cms < Q_min_transport", {

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

# ---------------------------------------------------------------------------
# T-03: SSL is never negative
# ---------------------------------------------------------------------------

test_that("T-03: SSL_mt_day is non-negative for all rows", {

  expect_true(
    all(ssl_df$SSL_mt_day >= 0),
    label = "SSL_mt_day must be >= 0 for all rows"
  )
})

# ---------------------------------------------------------------------------
# T-04: SSL follows the fitted power law a * Q^b
# ---------------------------------------------------------------------------

test_that("T-04: SSL_mt_day follows power-law at a known Q value", {

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
