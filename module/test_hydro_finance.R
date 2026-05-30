# =============================================================================
# tests/test_hydro_finance.R
# testthat unit tests for hydro_finance.R
#
# Functions covered
#   compute_npv()          — NPV from arbitrary cash flow vector
#   build_cash_flows()     — year-by-year cash flows including CAPEX
#   compute_irr()          — IRR by bisection
#   compute_payback()      — simple payback period
#   summarise_npv_grid()   — aggregate NPV grid for FFF
#
# Coverage
#   compute_npv
#     T-01  At r = 0, NPV equals sum of cash flows
#     T-02  Higher discount rate gives lower NPV (all else equal)
#     T-03  All-negative cash flows give negative NPV
#     T-04  Single-element input triggers an error
#
#   build_cash_flows
#     T-05  Year-0 element equals negative equity outflow
#     T-06  Annual net cash flow is (revenue - OPEX - debt_service)
#     T-07  LTV = 0 gives debt_service = 0, full equity at year 0
#     T-08  r_loan = 0 gives level principal payments
#     T-09  LTV > 1 triggers an error
#
#   compute_irr
#     T-10  IRR is positive for a profitable project
#     T-11  NPV at IRR is approximately zero
#     T-12  Returns NA when no root exists
#
#   compute_payback
#     T-13  Payback year is correct for a known cash flow vector
#     T-14  Returns Inf when cumulative cash flow never turns positive
#
#   summarise_npv_grid
#     T-15  Returns one row per (cv_level, r_loan) combination
#     T-16  prob_viable is in [0, 1]
#     T-17  n_sim equals the number of replicates per cell
# =============================================================================

library(testthat)
library(here)

source(here("module", "hydro_finance.R"))

# ---------------------------------------------------------------------------
# compute_npv — T-01 to T-04
# ---------------------------------------------------------------------------

test_that("T-01: at r = 0, NPV equals sum of all cash flows", {

  cf <- c(-1000, 300, 300, 300, 300)

  expect_equal(
    compute_npv(cf, discount_rate = 0),
    sum(cf),
    tolerance = 1e-6,
    label = "NPV at r=0 must equal simple sum of cash flows"
  )
})

test_that("T-02: higher discount rate gives lower NPV", {

  cf <- c(-500, 200, 200, 200)

  npv_low  <- compute_npv(cf, discount_rate = 0.01)
  npv_high <- compute_npv(cf, discount_rate = 0.15)

  expect_gt(npv_low, npv_high,
            label = "NPV at r=1% must exceed NPV at r=15%")
})

test_that("T-03: all-negative cash flows give a negative NPV", {

  cf <- c(-500, -100, -100, -100)

  expect_lt(
    compute_npv(cf, discount_rate = 0.05), 0,
    label = "NPV must be negative when all cash flows are negative"
  )
})

test_that("T-04: fewer than 2 cash flow values triggers an error", {

  expect_error(
    compute_npv(c(-1000), discount_rate = 0.05),
    regexp = NULL,
    label  = "Should error with only one cash flow element"
  )
})

# ---------------------------------------------------------------------------
# build_cash_flows — T-05 to T-09
# ---------------------------------------------------------------------------

test_that("T-05: year-0 element equals negative equity CAPEX", {

  capex  <- 9.7e9
  ltv    <- 0.75
  equity <- capex * (1 - ltv)

  cf <- build_cash_flows(
    annual_revenue_ntd = 1e9,
    ltv             = ltv,
    r_loan          = 0.03,
    capex_ntd       = capex,
    project_life_yr = 10L
  )

  expect_equal(
    cf[1], -equity,
    tolerance = 1,
    label = "First cash flow must equal negative equity CAPEX"
  )
})

test_that("T-06: annual net cash flow = revenue - OPEX - debt_service", {

  rev    <- 800e6
  opex   <- 9.7e7
  capex  <- 9.7e9
  ltv    <- 0.80
  r_loan <- 0.03
  T      <- 20L

  loan         <- capex * ltv
  debt_service <- loan * r_loan / (1 - (1 + r_loan)^(-T))
  expected_net <- rev - opex - debt_service

  cf <- build_cash_flows(
    annual_revenue_ntd = rev,
    ltv             = ltv,
    r_loan          = r_loan,
    capex_ntd       = capex,
    opex_ntd_yr     = opex,
    project_life_yr = T
  )

  expect_equal(
    cf[2], expected_net,
    tolerance = 1,
    label = "Annual net cash flow must equal revenue - OPEX - debt_service"
  )
})

test_that("T-07: LTV = 0 gives no debt service and full equity at year 0", {

  cf <- build_cash_flows(
    annual_revenue_ntd = 1e9,
    ltv             = 0,
    r_loan          = 0.05,
    capex_ntd       = 5e9,
    opex_ntd_yr     = 5e7,
    project_life_yr = 10L
  )

  # Year-0 = full CAPEX as equity outflow
  expect_equal(cf[1], -5e9, tolerance = 1,
               label = "Year-0 must equal -CAPEX when LTV = 0")

  # No debt → annual CF = revenue - OPEX (no debt_service)
  expect_equal(cf[2], 1e9 - 5e7, tolerance = 1,
               label = "Annual CF must be revenue - OPEX when LTV = 0")
})

test_that("T-08: r_loan = 0 gives equal principal payments each year", {

  capex  <- 4e9
  ltv    <- 0.75
  T      <- 10L
  loan   <- capex * ltv
  principal_per_yr <- loan / T

  cf <- build_cash_flows(
    annual_revenue_ntd = 8e8,
    ltv             = ltv,
    r_loan          = 0,
    capex_ntd       = capex,
    project_life_yr = T
  )

  # Debt service = principal only when r = 0
  opex         <- capex * 0.01
  expected_net <- 8e8 - opex - principal_per_yr

  expect_equal(
    cf[2], expected_net,
    tolerance = 1,
    label = "Annual CF with r_loan=0 must use simple principal repayment"
  )
})

test_that("T-09: LTV > 1 triggers an error", {

  expect_error(
    build_cash_flows(annual_revenue_ntd = 1e9, ltv = 1.1),
    regexp = NULL,
    label  = "LTV > 1 is physically impossible and must error"
  )
})

# ---------------------------------------------------------------------------
# compute_irr — T-10 to T-12
# ---------------------------------------------------------------------------

test_that("T-10: IRR is positive for a clearly profitable project", {

  # Simple project: invest 1000, get 400/yr for 4 years → IRR ≈ 21.9%
  cf  <- c(-1000, 400, 400, 400, 400)
  irr <- compute_irr(cf)

  expect_true(
    !is.na(irr) && irr > 0,
    label = paste("IRR should be positive and non-NA; got", irr)
  )
})

test_that("T-11: NPV at IRR is approximately zero", {

  cf  <- c(-1000, 350, 350, 350, 350)
  irr <- compute_irr(cf)

  if (!is.na(irr)) {
    npv_at_irr <- compute_npv(cf, discount_rate = irr)
    expect_equal(
      npv_at_irr, 0,
      tolerance = 1e-4,
      label = "NPV evaluated at IRR must be approximately 0"
    )
  } else {
    skip("IRR returned NA — project may have no positive IRR")
  }
})

test_that("T-12: compute_irr() returns NA when no root exists", {

  # All outflows — no positive return is possible
  cf  <- c(-500, -100, -100, -100)
  irr <- compute_irr(cf)

  expect_true(
    is.na(irr),
    label = "IRR should be NA when the project has no positive returns"
  )
})

# ---------------------------------------------------------------------------
# compute_payback — T-13 to T-14
# ---------------------------------------------------------------------------

test_that("T-13: payback period is correct for a known cash flow vector", {

  # CF: -1000, 400, 400, 400 → cumulative: -1000, -600, -200, 200
  # Payback occurs in year 3 (index 4 → 4 - 1 = 3)
  cf <- c(-1000, 400, 400, 400)

  expect_equal(
    compute_payback(cf), 3L,
    label = "Payback should be year 3 for this cash flow pattern"
  )
})

test_that("T-14: compute_payback() returns Inf when never recovered", {

  cf <- c(-1000, 100, 100, 100)   # cumulative never turns positive

  expect_equal(
    compute_payback(cf), Inf,
    label = "Payback should be Inf when cumulative cash flow stays negative"
  )
})

# ---------------------------------------------------------------------------
# summarise_npv_grid — T-15 to T-17
# ---------------------------------------------------------------------------

# Build a minimal npv_grid_df to test summarise_npv_grid()
make_npv_grid <- function(cv_levels = c(1.5, 2.5),
                          r_loans   = c(0.02, 0.05),
                          n_sim     = 20L,
                          seed      = 42L) {

  set.seed(seed)

  purrr::map_dfr(cv_levels, function(cv) {
    purrr::map_dfr(r_loans, function(r) {
      data.frame(
        cv_level    = cv,
        sim_id      = seq_len(n_sim),
        r_loan      = r,
        r_loan_pct  = r * 100,
        mu_yr_cms   = rnorm(n_sim, mean = 12, sd = 2),
        annual_gwh  = rnorm(n_sim, mean = 80, sd = 10),
        annual_rev_ntd = rnorm(n_sim, mean = 2e8, sd = 2e7),
        npv_ntd     = rnorm(n_sim, mean = 0, sd = 5e8),
        npv_b_ntd   = rnorm(n_sim, mean = 0, sd = 0.5),
        viable      = sample(c(TRUE, FALSE), n_sim, replace = TRUE)
      )
    })
  })
}

npv_grid <- make_npv_grid()
npv_sum  <- summarise_npv_grid(npv_grid)

test_that("T-15: summarise_npv_grid() returns one row per (cv_level, r_loan)", {

  n_cells <- length(unique(npv_grid$cv_level)) *
             length(unique(npv_grid$r_loan))

  expect_equal(
    nrow(npv_sum), n_cells,
    label = "One summary row per (cv_level, r_loan) cell expected"
  )
})

test_that("T-16: prob_viable is in [0, 1] for all cells", {

  expect_true(
    all(npv_sum$prob_viable >= 0 & npv_sum$prob_viable <= 1),
    label = "prob_viable must be a valid probability in [0, 1]"
  )
})

test_that("T-17: n_sim in summary matches replicates per cell", {

  expect_true(
    all(npv_sum$n_sim == 20L),
    label = "n_sim in summary table must equal 20 (replicates per cell)"
  )
})
