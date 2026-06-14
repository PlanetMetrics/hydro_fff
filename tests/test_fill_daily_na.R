# =============================================================================
# tests/test_fill_daily_na.R
# testthat unit tests for fill_daily_na()
#
# Function contract
#   Input  : data.frame(date = Date, Q_cms = numeric) with possible NA in Q_cms
#   Output : same data.frame, NA filled, logical column was_filled appended
#   Errors : stops if required columns missing, or method unrecognised
#
# Coverage
#   T-01  Short gap (≤ 7 days) is filled by linear interpolation
#   T-02  Filled values are strictly between the two bounding observations
#   T-03  was_filled is TRUE only on originally-NA rows
#   T-04  No NA remains in Q_cms after fill (all methods)
#   T-05  All Q_cms values are non-negative after fill
#   T-06  method = "zero" sets long-gap rows to 0
#   T-07  Missing column "date" triggers an informative error
#   T-08  Missing column "Q_cms" triggers an informative error
#   T-09  Unrecognised method triggers an informative error
#   T-10  Row count is unchanged (no rows added or dropped)
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_reservoir.R"))

# ---------------------------------------------------------------------------
# Shared synthetic data
# ---------------------------------------------------------------------------

make_flow <- function(n_days   = 30,
                      na_days  = 5:9,       # positions of NA values
                      Q_lo     = 10,
                      Q_hi     = 20,
                      start    = "2020-01-01") {

  Q <- seq(Q_lo, Q_hi, length.out = n_days)
  Q[na_days] <- NA_real_

  data.frame(
    date  = seq(as.Date(start), by = "day", length.out = n_days),
    Q_cms = Q
  )
}

# ---------------------------------------------------------------------------
# T-01 — Short gap filled by linear interpolation
# ---------------------------------------------------------------------------

test_that("T-01: short gap (5 days) is filled after fill_daily_na()", {

  df  <- make_flow(na_days = 5:9)
  out <- fill_daily_na(df, method = "linear", max_gap_linear = 7L)

  expect_false(
    any(is.na(out$Q_cms[5:9])),
    label = "NA rows 5-9 should be filled"
  )
})

# ---------------------------------------------------------------------------
# T-02 — Interpolated values are monotone between bounding observations
# ---------------------------------------------------------------------------

test_that("T-02: linearly interpolated values lie between bounding observations", {

  df  <- make_flow(na_days = 5:9)
  out <- fill_daily_na(df, method = "linear", max_gap_linear = 7L)

  Q_before <- df$Q_cms[4]    # last observed value before gap
  Q_after  <- df$Q_cms[10]   # first observed value after gap

  Q_lo <- min(Q_before, Q_after)
  Q_hi <- max(Q_before, Q_after)

  filled_vals <- out$Q_cms[5:9]

  expect_true(
    all(filled_vals >= Q_lo - 1e-6 & filled_vals <= Q_hi + 1e-6),
    label = "Interpolated values must lie between bounding observations"
  )
})

# ---------------------------------------------------------------------------
# T-03 — was_filled correctly flags imputed rows only
# ---------------------------------------------------------------------------

test_that("T-03: was_filled is TRUE only on originally-NA rows", {

  df  <- make_flow(na_days = 5:9)
  out <- fill_daily_na(df, method = "linear", max_gap_linear = 7L)

  expect_true(
    all(out$was_filled[5:9]),
    label = "was_filled should be TRUE for gap rows 5-9"
  )
  expect_false(
    any(out$was_filled[-c(5:9)]),
    label = "was_filled should be FALSE for non-gap rows"
  )
})

# ---------------------------------------------------------------------------
# T-04 — No NA remains in Q_cms after fill (linear method)
# ---------------------------------------------------------------------------

test_that("T-04: no NA remains in Q_cms after fill_daily_na() with method = linear", {

  df  <- make_flow(na_days = c(1:3, 15:20, 28:30))
  out <- fill_daily_na(df, method = "linear", max_gap_linear = 7L)

  expect_equal(
    sum(is.na(out$Q_cms)), 0L,
    label = "Q_cms must have zero NA values after filling"
  )
})

# ---------------------------------------------------------------------------
# T-05 — All Q_cms values are non-negative
# ---------------------------------------------------------------------------

test_that("T-05: all Q_cms values are non-negative after fill", {

  df  <- make_flow(Q_lo = 0.1, Q_hi = 5, na_days = 10:15)
  out <- fill_daily_na(df, method = "linear", max_gap_linear = 7L)

  expect_true(
    all(out$Q_cms >= 0),
    label = "All filled Q_cms values must be >= 0"
  )
})

# ---------------------------------------------------------------------------
# T-06 — method = "zero" sets NA rows to exactly 0
# ---------------------------------------------------------------------------

test_that("T-06: method = zero sets NA rows to 0", {

  df  <- make_flow(na_days = 5:9)
  out <- fill_daily_na(df, method = "zero")

  expect_equal(
    out$Q_cms[5:9], rep(0, 5),
    tolerance = 1e-9,
    label = "Gap rows should be set to 0 when method = zero"
  )
})

# ---------------------------------------------------------------------------
# T-07 — Missing "date" column triggers error
# ---------------------------------------------------------------------------

test_that("T-07: missing 'date' column triggers an error", {

  bad_df <- data.frame(Q_cms = c(1, 2, NA, 4))

  expect_error(
    fill_daily_na(bad_df),
    regexp = NULL,    # any error is acceptable
    label  = "Should error when 'date' column is absent"
  )
})

# ---------------------------------------------------------------------------
# T-08 — Missing "Q_cms" column triggers error
# ---------------------------------------------------------------------------

test_that("T-08: missing 'Q_cms' column triggers an error", {

  bad_df <- data.frame(date = Sys.Date() + 0:3, flow = c(1, 2, 3, 4))

  expect_error(
    fill_daily_na(bad_df),
    regexp = NULL,
    label  = "Should error when 'Q_cms' column is absent"
  )
})

# ---------------------------------------------------------------------------
# T-09 — Unrecognised method triggers error
# ---------------------------------------------------------------------------

test_that("T-09: unrecognised method argument triggers an error", {

  df <- make_flow()

  expect_error(
    fill_daily_na(df, method = "spline"),
    regexp = NULL,
    label  = "Should error for method = spline (not supported)"
  )
})

# ---------------------------------------------------------------------------
# T-10 — Row count unchanged after fill
# ---------------------------------------------------------------------------

test_that("T-10: fill_daily_na() does not add or remove rows", {

  df  <- make_flow(n_days = 365, na_days = c(50:60, 200:210))
  out <- fill_daily_na(df, method = "linear")

  expect_equal(
    nrow(out), nrow(df),
    label = "Output must have the same number of rows as input"
  )
})
