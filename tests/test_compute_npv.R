# =============================================================================
# tests/test_compute_npv.R
# testthat unit tests for compute_npv()
#
# Function contract
#   Input  : cash_flows     numeric vector, cash_flows[1] = year-0 CAPEX
#                            outflow (negative); remaining elements are
#                            annual net cash flows (years 1..T)
#            discount_rate  numeric, equity discount rate (default .FIN$r_equity)
#   Output : numeric scalar, NPV = sum(cash_flows[t] / (1 + r)^t)
#   Errors : stops if cash_flows has fewer than 2 elements
#
# This is the core NPV calculation used by every finance scenario
# (build_cash_flows() -> compute_npv()), and is where CAPEX enters the
# calculation as the negative year-0 cash flow.
#
# Coverage
#   T-01  At r = 0, NPV equals the simple sum of all cash flows
#   T-02  A higher discount rate gives a lower NPV (all else equal)
#   T-03  All-negative cash flows (CAPEX with no revenue) give a negative NPV
#   T-04  Fewer than 2 cash-flow values triggers an error
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_finance.R"))

# ---------------------------------------------------------------------------
# T-01: NPV at r = 0 equals the sum of cash flows
# ---------------------------------------------------------------------------

test_that("T-01: at r = 0, NPV equals sum of all cash flows", {

  # -1000 CAPEX at year 0, then 300/yr for 4 years
  cf <- c(-1000, 300, 300, 300, 300)

  expect_equal(
    compute_npv(cf, discount_rate = 0),
    sum(cf),
    tolerance = 1e-6,
    label = "NPV at r=0 must equal simple sum of cash flows"
  )
})

# ---------------------------------------------------------------------------
# T-02: higher discount rate gives lower NPV
# ---------------------------------------------------------------------------

test_that("T-02: higher discount rate gives lower NPV", {

  cf <- c(-500, 200, 200, 200)

  npv_low  <- compute_npv(cf, discount_rate = 0.01)
  npv_high <- compute_npv(cf, discount_rate = 0.15)

  expect_gt(npv_low, npv_high,
            label = "NPV at r=1% must exceed NPV at r=15%")
})

# ---------------------------------------------------------------------------
# T-03: all-negative cash flows give a negative NPV
# ---------------------------------------------------------------------------

test_that("T-03: all-negative cash flows give a negative NPV", {

  # CAPEX outflow with no offsetting revenue (e.g. project never generates)
  cf <- c(-500, -100, -100, -100)

  expect_lt(
    compute_npv(cf, discount_rate = 0.05), 0,
    label = "NPV must be negative when all cash flows are negative"
  )
})

# ---------------------------------------------------------------------------
# T-04: fewer than 2 cash-flow values triggers an error
# ---------------------------------------------------------------------------

test_that("T-04: fewer than 2 cash flow values triggers an error", {

  expect_error(
    compute_npv(c(-1000), discount_rate = 0.05),
    regexp = NULL,
    label  = "Should error with only one cash flow element"
  )
})
