# =============================================================================
# tests/test_simulate_one_weir.R
# testthat unit tests for simulate_one_weir()
#
# Function contract
#   Input  : daily_flow (date, Q_cms), weir_id, e_flow_mode, S_init, net_head_m
#   Output : data.frame — one row per day with energy, flow components, storage
#   Errors : stops if required columns missing or weir_id unrecognised
#
# Coverage
#   T-01  Output has correct column names
#   T-02  Output has the same number of rows as input
#   T-03  Storage S_end_m3 never exceeds S_max
#   T-04  Storage S_end_m3 is always >= 0
#   T-05  Q_power_cms <= Q_in_cms - Q_eflow_cms (power cannot exceed available)
#   T-06  Q_eflow_cms <= Q_in_cms (cannot release more than what arrives)
#   T-07  energy_kwh >= 0 always
#   T-08  At constant design flow, e-flow is always satisfied (committed mode)
#   T-09  At zero inflow, energy_kwh = 0 (no generation without water)
#   T-10  Missing "Q_cms" column triggers an error
#   T-11  Unrecognised weir_id triggers an error
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_reservoir.R"))

# ---------------------------------------------------------------------------
# Shared synthetic flow series
# ---------------------------------------------------------------------------

make_synth <- function(Q_val  = 24.3,
                       n_days = 365,
                       start  = "2020-01-01") {
  data.frame(
    date  = seq(as.Date(start), by = "day", length.out = n_days),
    Q_cms = Q_val
  )
}

# Run baseline simulation once and reuse
sim_W1 <- simulate_one_weir(
  make_synth(Q_val = 24.3),
  weir_id    = "W1",
  e_flow_mode = "committed",
  net_head_m  = 86
)

# ---------------------------------------------------------------------------
# T-01 — Required output columns are present
# ---------------------------------------------------------------------------

test_that("T-01: simulate_one_weir() returns all required columns", {

  required_cols <- c(
    "date", "Q_in_cms", "Q_eflow_cms", "Q_power_cms", "Q_spill_cms",
    "S_start_m3", "S_end_m3", "energy_kwh",
    "eflow_satisfied", "power_satisfied"
  )

  expect_true(
    all(required_cols %in% names(sim_W1)),
    label = paste("Missing columns:", paste(
      setdiff(required_cols, names(sim_W1)), collapse = ", "
    ))
  )
})

# ---------------------------------------------------------------------------
# T-02 — Row count matches input
# ---------------------------------------------------------------------------

test_that("T-02: output has the same number of rows as input", {

  df <- make_synth(n_days = 100)
  out <- simulate_one_weir(df, "W1", e_flow_mode = "committed", net_head_m = 86)

  expect_equal(nrow(out), 100L,
               label = "Row count must equal input length")
})

# ---------------------------------------------------------------------------
# T-03 — Storage never exceeds S_max
# ---------------------------------------------------------------------------

test_that("T-03: S_end_m3 never exceeds S_max (W1 = 967400 m3)", {

  S_max_W1 <- 967400

  expect_true(
    all(sim_W1$S_end_m3 <= S_max_W1 + 1e-3),
    label = "S_end_m3 must not exceed S_max for W1"
  )
})

# ---------------------------------------------------------------------------
# T-04 — Storage never goes negative
# ---------------------------------------------------------------------------

test_that("T-04: S_end_m3 is always >= 0", {

  expect_true(
    all(sim_W1$S_end_m3 >= 0),
    label = "Storage cannot be negative"
  )
})

# ---------------------------------------------------------------------------
# T-05 — Power diversion never exceeds available flow
# ---------------------------------------------------------------------------

test_that("T-05: Q_power_cms <= Q_in_cms - Q_eflow_cms (water balance)", {

  Q_avail <- sim_W1$Q_in_cms - sim_W1$Q_eflow_cms

  expect_true(
    all(sim_W1$Q_power_cms <= Q_avail + 1e-6),
    label = "Power diversion must not exceed available flow after e-flow"
  )
})

# ---------------------------------------------------------------------------
# T-06 — E-flow release never exceeds inflow
# ---------------------------------------------------------------------------

test_that("T-06: Q_eflow_cms <= Q_in_cms (cannot release more than arrives)", {

  expect_true(
    all(sim_W1$Q_eflow_cms <= sim_W1$Q_in_cms + 1e-6),
    label = "E-flow release must not exceed daily inflow"
  )
})

# ---------------------------------------------------------------------------
# T-07 — Energy is always non-negative
# ---------------------------------------------------------------------------

test_that("T-07: energy_kwh >= 0 for all days", {

  expect_true(
    all(sim_W1$energy_kwh >= 0),
    label = "Daily energy cannot be negative"
  )
})

# ---------------------------------------------------------------------------
# T-08 — At constant design flow, e-flow always satisfied (committed mode)
# ---------------------------------------------------------------------------

test_that("T-08: e-flow always satisfied when Q_in = design flow (committed)", {

  # W1 design flow = 24.3 cms >> committed e-flow = 0.48 cms
  expect_true(
    all(sim_W1$eflow_satisfied),
    label = "eflow_satisfied must be TRUE every day when Q = 24.3 cms"
  )
})

# ---------------------------------------------------------------------------
# T-09 — At zero inflow, energy = 0
# ---------------------------------------------------------------------------

test_that("T-09: energy_kwh = 0 when Q_in = 0 (no water, no generation)", {

  zero_flow <- make_synth(Q_val = 0, n_days = 30)
  out <- simulate_one_weir(zero_flow, "W1",
                           e_flow_mode = "committed",
                           net_head_m  = 86,
                           S_init      = 0)

  expect_equal(
    sum(out$energy_kwh), 0,
    tolerance = 1e-6,
    label = "Total energy must be 0 when inflow and initial storage are both 0"
  )
})

# ---------------------------------------------------------------------------
# T-10 — Missing Q_cms column triggers error
# ---------------------------------------------------------------------------

test_that("T-10: missing 'Q_cms' column triggers an error", {

  bad_df <- data.frame(date = Sys.Date() + 0:9, flow = rep(10, 10))

  expect_error(
    simulate_one_weir(bad_df, "W1"),
    regexp = NULL,
    label  = "Should error when Q_cms column is absent"
  )
})

# ---------------------------------------------------------------------------
# T-11 — Unrecognised weir_id triggers error
# ---------------------------------------------------------------------------

test_that("T-11: unrecognised weir_id triggers an error", {

  df <- make_synth()

  expect_error(
    simulate_one_weir(df, weir_id = "W3"),
    regexp = NULL,
    label  = "Should error for weir_id = W3 (not defined)"
  )
})
