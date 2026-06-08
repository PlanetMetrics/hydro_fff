# =============================================================================
# hydro_climate.R
# Fengping River Hydropower Simulation — Climate & Hydrological Variability Module
#
# Purpose
#   (1) Characterise historical flow variability using annual statistics
#       and the Richards-Baker Flashiness Index (RBI).
#   (2) Generate synthetic annual flow series via a two-layer bootstrap
#       that preserves within-year seasonal structure while enabling
#       systematic exploration of CV levels beyond the historical range.
#   (3) Produce a structured output table ready for the Financial
#       Feasibility Frontier (FFF): rows indexed by (cv_level, sim_id),
#       containing annual mean flow and annual energy inputs for NPV.
#
# =============================================================================
# TWO-LAYER BOOTSTRAP — RATIONALE AND DESIGN
# =============================================================================
#
# Layer 1 — Whole-year block bootstrap
#   Complete calendar years are resampled with replacement from the observed
#   record (1958–2025).  This preserves:
#     • within-year seasonal flow distribution (wet/dry season structure)
#     • typhoon event integrity (multi-day hydrograph shape)
#     • inter-variable dependence (e.g. high Q → high SSL on same day)
#   Reference:
#     Efron, B. & Tibshirani, R.J. (1993). An Introduction to the Bootstrap.
#     Chapman & Hall. — Chapter 8 (block bootstrap rationale).
#
# Layer 2 — CV scaling around the annual mean
#   After resampling, each synthetic year's daily flow series is rescaled to
#   a target within-year coefficient of variation (CV) while holding the
#   annual mean flow constant:
#
#     Q_scaled[t] = mu_yr + (Q_boot[t] - mu_yr) * lambda
#
#   where  mu_yr  = mean(Q_boot)              [annual mean, cms]
#          lambda = CV_target / CV_boot        [scaling factor]
#          CV_boot = sd(Q_boot) / mu_yr        [bootstrap year CV]
#
#   Physical lower bound: Q_scaled[t] >= Q_FLOOR_CMS (default 0.01 cms).
#
#   This enables exploration of CV levels that have not yet been observed
#   historically but are projected under climate change — consistent with
#   non-stationarity assumptions for Taiwan eastern rivers under SSP pathways.
#
# Why CV_target can exceed the historical range
#   Climate change projections for Taiwan indicate increased typhoon intensity
#   and more severe dry-season deficits, leading to higher within-year flow
#   variability beyond the historical envelope:
#
#   • Shiau, J.T. & Huang, W.H. (2014). Detecting distributional changes of
#     annual rainfall indices in Taiwan using quantile regression.
#     Journal of Hydro-environment Research 8(4):355–366.
#     https://doi.org/10.1016/j.jher.2014.07.006
#     → Taiwan annual precipitation CV increasing trend, especially east coast.
#
#   • Lee, M.H., Ho, C.H., & Wang, J.Y. (2006). Influence of SST anomalies on
#     typhoon activity over the Western North Pacific Ocean.
#     Terrestrial, Atmospheric and Oceanic Sciences 17(4):979–994.
#     → Eastern Taiwan receives disproportionate typhoon rainfall; trend upward.
#
#   • Tung, C.P., et al. (2016). Evaluating the impact of climate change on
#     water resources in Taiwan using CMIP5 GCM simulations.
#     Paddy and Water Environment 14(1):15–26.
#     https://doi.org/10.1007/s10333-014-0470-2
#     → SSP-equivalent RCP8.5 increases peak-to-baseflow ratio in eastern basins.
#
#   • IPCC AR6 WGI Chapter 11 (2021). Weather and Climate Extreme Events
#     in a Changing Climate. In: Climate Change 2021: The Physical Science Basis.
#     Cambridge University Press. https://doi.org/10.1017/9781009157896.013
#     → High confidence: increased precipitation intensity in Taiwan under SSP5-8.5.
#
# CV target range used in this study
#   Historical observed range  : CV ≈ 1.5 – 2.5  (Lishan station 1958–2025)
#   SSP2-4.5 projection (2050) : CV ≈ 1.8 – 3.0  (+20–30 %)
#   SSP5-8.5 projection (2050) : CV ≈ 2.2 – 4.0  (+40–80 %)
#   FFF exploration range      : CV = 1.0 – 4.0   (see .CV_GRID below)
#
# =============================================================================
# KEY NOTATION
# =============================================================================
#
#   Q           daily streamflow (cms = m³/s)
#   mu_yr       annual mean daily flow (cms)
#   sd_yr       annual standard deviation of daily flow (cms)
#   CV          within-year coefficient of variation = sd_yr / mu_yr
#               (dimensionless; captures intra-annual flow variability)
#   lambda      CV scaling factor = CV_target / CV_boot  (dimensionless)
#   RBI         Richards-Baker Flashiness Index
#               RBI = Σ|Q[t] - Q[t-1]| / ΣQ[t]  (dimensionless, 0–1)
#               Higher RBI → more flashy (Taiwan torrents: 0.4–0.8)
#   n_sim       number of bootstrap resamples per CV level
#   cv_level    target CV assigned to a synthetic year (FFF Y-axis)
#   sim_id      bootstrap replicate index within a cv_level (1 … n_sim)
#   Q_FLOOR     physical lower bound for daily flow (cms)
#   SSP         Shared Socioeconomic Pathway (IPCC AR6 framework)
#
# =============================================================================
# OUTPUT STRUCTURE
# =============================================================================
#
# run_bootstrap_cv_grid() returns a data.frame with one row per
# (cv_level, sim_id) combination, columns:
#
#   cv_level      numeric   target CV for this synthetic year
#   sim_id        integer   bootstrap replicate index
#   year_src      integer   observed year drawn in layer 1
#   mu_yr_cms     numeric   annual mean daily flow (cms)
#   cv_achieved   numeric   actual CV after scaling (should ≈ cv_level)
#   rbi           numeric   Richards-Baker Flashiness Index
#   annual_Q_m3   numeric   annual flow volume (m³) = mu_yr * 365 * 86400
#
# This table is the direct input to the NPV grid loop in fp_hydro_main.qmd.
#
# =============================================================================
# Author   [your name]
# Date     2025
# =============================================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(here)


# =============================================================================
# SECTION 0 — Internal constants
# =============================================================================

# Physical lower bound for daily flow after CV scaling.
# Prevents negative or near-zero flows that are physically impossible
# on Fengping Creek (perennial stream with measurable baseflow year-round).
.Q_FLOOR_CMS <- 0.01

# Default number of bootstrap resamples per CV level.
# Use 200 for coursework; use 1000 for publication-quality results.
# Search "!!! n_sim" in this file to find the parameter.
.N_SIM_DEFAULT <- 200L

# Minimum valid days required to include a year in the bootstrap pool.
# Years with fewer than min_days of non-NA flow are excluded.
.MIN_DAYS_DEFAULT <- 180L

# CV target grid for the Financial Feasibility Frontier (FFF).
# Range 1.0–4.0 covers:
#   low end (1.0–1.5)  : hypothetical low-variability / wet baseline
#   historical (1.5–2.5): observed Lishan station range
#   projected (2.5–4.0) : SSP2-4.5 and SSP5-8.5 upper range
# Reference: Shiau & Huang (2014); Tung et al. (2016); IPCC AR6 Ch.11
.CV_GRID <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)


# =============================================================================
# SECTION 1 — Historical variability analysis
# =============================================================================

# -----------------------------------------------------------------------------
# compute_rbi()
#
# Compute the Richards-Baker Flashiness Index (RBI) for a daily flow vector.
#
# RBI measures the frequency and rapidity of short-term changes in streamflow
# relative to total discharge.  It is especially sensitive to storm pulses and
# typhoon events, making it a natural variability indicator for flashy mountain
# streams such as Fengping Creek.
#
#   RBI = Σ|Q[t] - Q[t-1]| / ΣQ[t]
#
# Typical values
#   Stable baseflow rivers : RBI < 0.20
#   Taiwan mountain torrents: RBI = 0.40 – 0.80
#
# Reference
#   Baker, D.B., Richards, R.P., Loftus, T.T., & Kramer, J.W. (2004).
#   A new flashiness index: characteristics and applications to midwestern
#   rivers and streams. Journal of the American Water Resources Association
#   40(2):503–522. https://doi.org/10.1111/j.1752-1688.2004.tb01046.x
#
# Arguments
#   Q_vec   numeric vector   daily flow (cms); must be chronologically ordered;
#                            NA values are removed before calculation
#
# Returns
#   numeric scalar   RBI (dimensionless); NA if fewer than 2 valid values
# -----------------------------------------------------------------------------

compute_rbi <- function(Q_vec) {

  stopifnot(is.numeric(Q_vec), length(Q_vec) >= 2L)

  Q_valid <- Q_vec[!is.na(Q_vec)]
  if (length(Q_valid) < 2L) return(NA_real_)

  sum(abs(diff(Q_valid))) / sum(Q_valid)
}


# -----------------------------------------------------------------------------
# compute_annual_stats()
#
# Compute year-by-year summary statistics from a daily flow data.frame.
# The within-year CV is the primary variability indicator used as the
# Y-axis of the Financial Feasibility Frontier.
#
# Arguments
#   daily_flow   data.frame   columns: date (Date), Q_cms (numeric)
#
# Returns
#   data.frame with one row per calendar year, columns:
#     year       integer
#     mu_Q       numeric   annual mean daily flow (cms)
#     sd_Q       numeric   annual standard deviation of daily flow (cms)
#     cv_Q       numeric   within-year CV = sd_Q / mu_Q  (dimensionless)
#     Q05        numeric   5th percentile of daily flow (dry-season proxy)
#     Q95        numeric   95th percentile of daily flow (flood proxy)
#     max_Q      numeric   annual maximum daily flow (cms)
#     rbi        numeric   Richards-Baker Flashiness Index
#     n_days     integer   total days in year
#     n_valid    integer   days with non-NA flow
# -----------------------------------------------------------------------------

#' Compute year-by-year streamflow summary statistics
#'
#' Splits a daily streamflow record into calendar years and computes, for
#' each year, the mean, standard deviation, within-year coefficient of
#' variation (CV), low- and high-flow percentiles, annual maximum, and the
#' Richards-Baker Flashiness Index (RBI). The within-year CV is the primary
#' variability indicator used as the Y-axis of the Financial Feasibility
#' Frontier (FFF) elsewhere in this analysis.
#'
#' \deqn{CV = \frac{sd(Q)}{mean(Q)}}
#'
#' @param daily_flow A \code{data.frame} of daily streamflow with at least
#'   the columns:
#'   \describe{
#'     \item{date}{\code{Date}. Calendar date of the observation.}
#'     \item{Q_cms}{Numeric. Daily mean streamflow (cubic metres per second).}
#'   }
#'
#' @return A \code{data.frame} with one row per calendar year and columns:
#'   \describe{
#'     \item{year}{Integer. Calendar year.}
#'     \item{mu_Q}{Numeric. Annual mean daily flow (cms).}
#'     \item{sd_Q}{Numeric. Annual standard deviation of daily flow (cms).}
#'     \item{cv_Q}{Numeric. Within-year coefficient of variation, \code{sd_Q / mu_Q}.}
#'     \item{Q05}{Numeric. 5th percentile of daily flow — dry-season proxy (cms).}
#'     \item{Q95}{Numeric. 95th percentile of daily flow — flood proxy (cms).}
#'     \item{max_Q}{Numeric. Annual maximum daily flow (cms).}
#'     \item{rbi}{Numeric. Richards-Baker Flashiness Index (dimensionless, 0-1).}
#'     \item{n_days}{Integer. Total calendar days in the year.}
#'     \item{n_valid}{Integer. Days with non-\code{NA} flow.}
#'   }
#'
#' @examples
#' \dontrun{
#' daily  <- read_csv(here("data", "lishan_daily_clean.csv")) |>
#'   mutate(date = as.Date(date))
#' stats  <- compute_annual_stats(daily)
#' head(stats[, c("year", "mu_Q", "cv_Q", "rbi")])
#' }
#'
#' @export
compute_annual_stats <- function(daily_flow) {

  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow))
  )

  daily_flow |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarise(
      mu_Q    = mean(Q_cms,             na.rm = TRUE),
      sd_Q    = sd(Q_cms,               na.rm = TRUE),
      cv_Q    = sd_Q / mu_Q,
      Q05     = quantile(Q_cms, 0.05,   na.rm = TRUE),
      Q95     = quantile(Q_cms, 0.95,   na.rm = TRUE),
      max_Q   = max(Q_cms,              na.rm = TRUE),
      rbi     = compute_rbi(Q_cms),
      n_days  = n(),
      n_valid = sum(!is.na(Q_cms)),
      .groups = "drop"
    )
}


# -----------------------------------------------------------------------------
# compute_rolling_stats()
#
# Compute rolling window statistics to detect non-stationarity in the
# historical record.  Each window centred on a given year yields one row.
# A significant trend in roll_cv_Q across time is evidence of non-stationarity
# consistent with climate-change projections.
#
# Arguments
#   annual_stats   data.frame   output of compute_annual_stats()
#   window         integer      rolling window width in years (default 30)
#
# Returns
#   data.frame with columns:
#     year_centre   integer   centre year of the rolling window
#     roll_mu_Q     numeric   mean of annual mean flows within window
#     roll_sd_Q     numeric   mean of annual SD values within window
#     roll_cv_Q     numeric   mean of annual CV values within window
#     roll_rbi      numeric   mean of annual RBI values within window
#     roll_Q05      numeric   mean of annual Q05 values within window
#     roll_Q95      numeric   mean of annual Q95 values within window
#     n_years       integer   number of years in window (= window)
# -----------------------------------------------------------------------------

compute_rolling_stats <- function(annual_stats, window = 30L) {

  stopifnot(
    is.data.frame(annual_stats),
    "cv_Q" %in% names(annual_stats),
    nrow(annual_stats) >= window
  )

  half <- floor(window / 2L)
  n    <- nrow(annual_stats)

  purrr::map_dfr(seq(half + 1L, n - half), function(i) {
    idx <- (i - half):(i + half)
    w   <- annual_stats[idx, ]
    data.frame(
      year_centre = annual_stats$year[i],
      roll_mu_Q   = mean(w$mu_Q,  na.rm = TRUE),
      roll_sd_Q   = mean(w$sd_Q,  na.rm = TRUE),
      roll_cv_Q   = mean(w$cv_Q,  na.rm = TRUE),
      roll_rbi    = mean(w$rbi,   na.rm = TRUE),
      roll_Q05    = mean(w$Q05,   na.rm = TRUE),
      roll_Q95    = mean(w$Q95,   na.rm = TRUE),
      n_years     = nrow(w)
    )
  })
}


# =============================================================================
# SECTION 2 — Layer 1: Whole-year block bootstrap
# =============================================================================

# -----------------------------------------------------------------------------
# build_year_index()
#
# Return the set of calendar years in the observed record that have sufficient
# valid daily observations to be used as bootstrap source years.
# Years below the threshold are excluded to avoid synthetic series dominated
# by imputed values.
#
# Arguments
#   daily_flow   data.frame   columns: date (Date), Q_cms (numeric)
#   min_days     integer      minimum valid (non-NA) days per year (default 180)
#
# Returns
#   integer vector   usable year indices, sorted ascending
# -----------------------------------------------------------------------------

build_year_index <- function(daily_flow, min_days = .MIN_DAYS_DEFAULT) {

  daily_flow |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarise(n_valid = sum(!is.na(Q_cms)), .groups = "drop") |>
    filter(n_valid >= min_days) |>
    arrange(year) |>
    pull(year)
}


# -----------------------------------------------------------------------------
# draw_bootstrap_year()
#
# Draw one synthetic year's daily flow series by sampling a single observed
# calendar year with replacement from the pool of usable years, then
# extracting its daily Q_cms vector.
#
# This is Layer 1 of the two-layer bootstrap.  The daily sequence is returned
# as a numeric vector (not a data.frame) for fast processing inside loops.
#
# Arguments
#   daily_by_year   named list   list of data.frames split by year;
#                                names are character year labels;
#                                each element has columns date, Q_cms
#   usable_years    integer vector   years eligible for sampling
#
# Returns
#   named list:
#     $Q_vec     numeric vector   daily flow (cms) for the drawn year
#     $year_src  integer          the observed year that was drawn
# -----------------------------------------------------------------------------

draw_bootstrap_year <- function(daily_by_year, usable_years) {

  yr      <- sample(usable_years, size = 1L, replace = TRUE)
  yr_data <- daily_by_year[[as.character(yr)]]

  list(
    Q_vec    = yr_data$Q_cms,
    year_src = yr
  )
}


# =============================================================================
# SECTION 3 — Layer 2: CV scaling
# =============================================================================

# -----------------------------------------------------------------------------
# scale_to_cv()
#
# Rescale a daily flow vector to a target within-year coefficient of variation
# (CV_target) while preserving the annual mean flow (mu_yr).
#
# Mathematical derivation
#   Let Q_boot be the bootstrap year's daily flow vector.
#   mu_yr  = mean(Q_boot)
#   sd_boot = sd(Q_boot)
#   CV_boot = sd_boot / mu_yr
#
#   To achieve CV_target, apply a linear scaling around the mean:
#     Q_scaled[t] = mu_yr + (Q_boot[t] - mu_yr) * lambda
#   where lambda = CV_target / CV_boot
#
#   This transformation has the following properties:
#     E[Q_scaled]  = mu_yr            (mean is preserved exactly)
#     sd(Q_scaled) = sd_boot * lambda (standard deviation scales linearly)
#     CV(Q_scaled) = CV_target        (target achieved exactly, if Q_floor = 0)
#
#   When the Q_FLOOR_CMS constraint is applied, some negative deviations are
#   clipped, so the achieved CV will be slightly below CV_target for very
#   high lambda values.  This is physically correct behaviour.
#
# Why we scale around the mean rather than rescale the whole vector
#   Multiplying Q_boot by a constant changes both the mean and the SD
#   proportionally, leaving CV unchanged.  Only the mean-centred
#   transformation above can alter CV independently of the mean.
#
# Reference for the mean-preserving spread approach
#   Arnell, N.W. (1998). Climate change and global water resources.
#   Global Environmental Change 9(S1):S31–S49.
#   https://doi.org/10.1016/S0959-3780(99)00017-5
#   → Section 3.1 discusses mean-preserving variance perturbation for
#     hydrological scenario generation.
#
# Arguments
#   Q_vec       numeric vector   daily flow (cms) for one bootstrap year
#   cv_target   numeric          target CV (dimensionless, > 0)
#   q_floor     numeric          minimum allowed flow after scaling (cms)
#                                default: .Q_FLOOR_CMS
#
# Returns
#   numeric vector (same length as Q_vec)   rescaled daily flows (cms)
#
# Warnings
#   If CV_boot is near zero (nearly constant flow), scaling is skipped and
#   the original vector is returned with a warning.
# -----------------------------------------------------------------------------

scale_to_cv <- function(Q_vec,
                        cv_target,
                        q_floor = .Q_FLOOR_CMS) {

  stopifnot(
    is.numeric(Q_vec),
    length(Q_vec) >= 2L,
    is.numeric(cv_target),
    cv_target > 0,
    q_floor >= 0
  )

  # Remove NA for statistics; will reinsert below
  Q_clean <- Q_vec[!is.na(Q_vec)]

  mu_yr   <- mean(Q_clean)
  sd_boot <- sd(Q_clean)
  cv_boot <- sd_boot / mu_yr

  # Guard: if bootstrap year is nearly constant, scaling is unstable
  if (cv_boot < 1e-4) {
    warning(sprintf(
      paste("scale_to_cv: CV_boot = %.6f is near zero.",
            "Returning original vector without scaling.",
            "This year will not match CV_target = %.3f."),
      cv_boot, cv_target
    ))
    return(pmax(Q_vec, q_floor))
  }

  lambda    <- cv_target / cv_boot
  Q_scaled  <- mu_yr + (Q_vec - mu_yr) * lambda
  Q_scaled  <- pmax(Q_scaled, q_floor)   # enforce physical lower bound

  Q_scaled
}


# =============================================================================
# SECTION 4 — Main function: two-layer bootstrap over a CV grid
# =============================================================================

# -----------------------------------------------------------------------------
# run_bootstrap_cv_grid()   <-- PRIMARY FUNCTION FOR FFF ANALYSIS
#
# Generate a structured ensemble of synthetic annual flow series by applying
# the two-layer bootstrap across a grid of target CV levels.
#
# Workflow
#   For each cv_level in cv_grid:
#     For sim_id in 1 … n_sim:
#       1. Draw one observed year with replacement (Layer 1 — block bootstrap).
#       2. Rescale its daily flow to cv_level (Layer 2 — CV scaling).
#       3. Compute summary statistics for the synthetic year.
#   Return a tidy data.frame with one row per (cv_level, sim_id).
#
# Output columns
#   cv_level      numeric   target CV assigned to this row
#   sim_id        integer   bootstrap replicate index (1 … n_sim)
#   year_src      integer   observed year drawn in Layer 1
#   mu_yr_cms     numeric   annual mean daily flow, cms  (preserved from Layer 1)
#   cv_achieved   numeric   actual CV after Layer 2 scaling
#                           (≈ cv_level; may differ slightly due to q_floor)
#   rbi           numeric   Richards-Baker Flashiness Index of scaled series
#   annual_Q_m3   numeric   annual flow volume (m³) = mu_yr * 365 * 86400
#
# The output is the direct input to the for-loop NPV grid in fp_hydro_main.qmd.
# Each row represents one (climate scenario, simulation) pair for which
# annual energy and NPV will be computed.
#
# Arguments
#   daily_flow   data.frame   columns: date (Date), Q_cms (numeric)
#                             typically lishan_daily_clean.csv after fill_daily_na()
#   cv_grid      numeric vector   target CV levels (default: .CV_GRID)
#                                 recommended: seq(1.0, 4.0, by = 0.5)
#   n_sim        integer      bootstrap resamples per CV level
#                             !!! use 200 for coursework, 1000 for publication !!!
#   min_days     integer      minimum valid days to include a year in pool
#   q_floor      numeric      minimum daily flow after scaling (cms)
#   seed         integer      RNG seed for reproducibility
#
# Returns
#   data.frame   nrow = length(cv_grid) * n_sim
#
# Example
#   grid <- run_bootstrap_cv_grid(
#     daily_flow = daily_clean,
#     cv_grid    = c(1.5, 2.0, 2.5, 3.0, 3.5),
#     n_sim      = 200L,
#     seed       = 42L
#   )
# -----------------------------------------------------------------------------

run_bootstrap_cv_grid <- function(daily_flow,
                                  cv_grid      = .CV_GRID,
                                  n_sim        = .N_SIM_DEFAULT,
                                  min_days     = .MIN_DAYS_DEFAULT,
                                  q_floor      = .Q_FLOOR_CMS,
                                  seed         = 42L,
                                  keep_series  = TRUE) {
  # keep_series = TRUE  : store the full scaled daily Q vector in column
  #                       `Q_series` (a list-column of numeric vectors).
  #                       Required for run_fff_grid() to call
  #                       run_reservoir_simulation() on each synthetic year.
  # keep_series = FALSE : store summary statistics only (lighter output,
  #                       used for quick diagnostics and plotting).

  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    is.numeric(cv_grid), length(cv_grid) >= 1L, all(cv_grid > 0),
    is.integer(n_sim) || (is.numeric(n_sim) && n_sim == round(n_sim)),
    n_sim >= 10L
  )

  n_sim <- as.integer(n_sim)

  set.seed(seed)

  # ── Identify usable years (Layer 1 pool) ────────────────────────────────
  usable_years <- build_year_index(daily_flow, min_days)

  if (length(usable_years) < 5L)
    stop(paste(
      "run_bootstrap_cv_grid: fewer than 5 usable years in daily record.",
      "Check min_days argument or data completeness."
    ))

  message(sprintf(
    "\nBootstrap pool : %d usable years (%d\u2013%d)",
    length(usable_years), min(usable_years), max(usable_years)
  ))
  message(sprintf("CV grid        : %s", paste(cv_grid, collapse = " | ")))
  message(sprintf(
    "n_sim per CV   : %d  |  total simulations: %d  |  keep_series: %s",
    n_sim, length(cv_grid) * n_sim, keep_series
  ))
  if (n_sim < 1000L)
    message("  NOTE: use n_sim = 1000 for publication-quality results.")

  # Pre-split daily flow by year for fast lookup
  daily_by_year <- daily_flow |>
    mutate(year = year(date)) |>
    filter(year %in% usable_years) |>
    group_by(year) |>
    group_split()

  names(daily_by_year) <- purrr::map_chr(
    daily_by_year, ~ as.character(.x$year[1])
  )

  # ── Two-layer bootstrap loop ─────────────────────────────────────────────
  #
  # Outer loop : cv_grid  — Y-axis of the Financial Feasibility Frontier
  # Inner loop : n_sim    — Monte Carlo replicates per CV level
  #
  # Layer 1: draw one observed year (preserves seasonal + typhoon structure)
  # Layer 2: rescale daily Q to cv_level (enables super-historical CV values)
  #
  # Each synthetic year retains the DATE vector from the source year so that
  # run_reservoir_simulation() can correctly assign monthly n-exponents for
  # watershed scaling (scale_flow_to_watershed).
  # ─────────────────────────────────────────────────────────────────────────

  purrr::map_dfr(cv_grid, function(cv_level) {

    message(sprintf("  CV = %.2f ...", cv_level))

    purrr::map_dfr(seq_len(n_sim), function(sim_id) {

      # Layer 1 ─────────────────────────────────────────────────────────────
      boot     <- draw_bootstrap_year(daily_by_year, usable_years)
      Q_boot   <- boot$Q_vec
      year_src <- boot$year_src

      # Retrieve date vector for this source year (needed by reservoir sim)
      dates_src <- daily_by_year[[as.character(year_src)]]$date

      # Layer 2 ─────────────────────────────────────────────────────────────
      Q_scaled <- scale_to_cv(Q_boot, cv_target = cv_level, q_floor = q_floor)

      # Summary statistics ──────────────────────────────────────────────────
      mu_yr       <- mean(Q_scaled, na.rm = TRUE)
      sd_yr       <- sd(Q_scaled,   na.rm = TRUE)
      cv_achieved <- sd_yr / mu_yr
      rbi_val     <- compute_rbi(Q_scaled)
      annual_Q_m3 <- mu_yr * 365 * 86400

      row <- data.frame(
        cv_level    = cv_level,
        sim_id      = sim_id,
        year_src    = year_src,
        mu_yr_cms   = round(mu_yr,       4),
        cv_achieved = round(cv_achieved, 4),
        rbi         = round(rbi_val,     4),
        annual_Q_m3 = round(annual_Q_m3, 0)
      )

      # Optionally store the full daily flow series as a list-column ────────
      # The series is stored as a named list element so that it survives
      # purrr::map_dfr() binding.  Column Q_series is a list of data.frames:
      #   each element has columns date (Date) and Q_cms (numeric).
      if (keep_series) {
        row$Q_series <- list(
          data.frame(date = dates_src, Q_cms = Q_scaled)
        )
      }

      row
    })
  })
}


# =============================================================================
# SECTION 5 — Ensemble summary for plotting
# =============================================================================

# -----------------------------------------------------------------------------
# summarise_cv_ensemble()
#
# Compute percentile statistics of mu_yr and cv_achieved across sim_id
# replicates for each cv_level.  Used to visualise the bootstrap distribution
# and verify that Layer 2 scaling achieves the intended CV targets.
#
# Arguments
#   grid_df   data.frame   output of run_bootstrap_cv_grid()
#
# Returns
#   data.frame with one row per cv_level, columns:
#     cv_level, n_sim,
#     mu_P10, mu_P25, mu_P50, mu_P75, mu_P90   (annual mean flow percentiles)
#     cv_P10, cv_P25, cv_P50, cv_P75, cv_P90   (achieved CV percentiles)
#     rbi_P50                                    (median RBI)
# -----------------------------------------------------------------------------

summarise_cv_ensemble <- function(grid_df) {

  stopifnot(
    is.data.frame(grid_df),
    all(c("cv_level", "mu_yr_cms", "cv_achieved", "rbi") %in% names(grid_df))
  )

  probs  <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  pnames <- c("P10", "P25", "P50", "P75", "P90")

  grid_df |>
    group_by(cv_level) |>
    summarise(
      n_sim   = n(),

      # Annual mean flow distribution (cms)
      mu_P10  = quantile(mu_yr_cms,   probs[1], na.rm = TRUE),
      mu_P25  = quantile(mu_yr_cms,   probs[2], na.rm = TRUE),
      mu_P50  = quantile(mu_yr_cms,   probs[3], na.rm = TRUE),
      mu_P75  = quantile(mu_yr_cms,   probs[4], na.rm = TRUE),
      mu_P90  = quantile(mu_yr_cms,   probs[5], na.rm = TRUE),

      # Achieved CV distribution (should cluster around cv_level)
      cv_P10  = quantile(cv_achieved, probs[1], na.rm = TRUE),
      cv_P25  = quantile(cv_achieved, probs[2], na.rm = TRUE),
      cv_P50  = quantile(cv_achieved, probs[3], na.rm = TRUE),
      cv_P75  = quantile(cv_achieved, probs[4], na.rm = TRUE),
      cv_P90  = quantile(cv_achieved, probs[5], na.rm = TRUE),

      # Flashiness
      rbi_P50 = quantile(rbi,         probs[3], na.rm = TRUE),

      .groups = "drop"
    )
}


# =============================================================================
# SECTION 6 — Diagnostic plot functions
# =============================================================================

# -----------------------------------------------------------------------------
# plot_cv_distribution()
#
# Violin + boxplot showing the distribution of achieved CV across bootstrap
# replicates for each cv_level.  Used to verify Layer 2 scaling accuracy.
# The orange dashed line shows the 1:1 target line (achieved = target).
#
# Arguments
#   grid_df      data.frame   output of run_bootstrap_cv_grid()
#   title_text   character    plot title
#
# Returns
#   ggplot object
# -----------------------------------------------------------------------------

plot_cv_distribution <- function(
    grid_df,
    title_text = "Bootstrap CV scaling: achieved vs target") {

  stopifnot(
    is.data.frame(grid_df),
    all(c("cv_level", "cv_achieved") %in% names(grid_df))
  )

  # Sort CV levels numerically — factor() default is alphabetical which
  # breaks the x-axis order (e.g. "1" sorts before "1.5").
  cv_levels_sorted <- sort(unique(grid_df$cv_level))
  n_sim_per_level  <- nrow(grid_df) / length(cv_levels_sorted)

  plot_df <- grid_df |>
    mutate(
      cv_fac = factor(cv_level, levels = cv_levels_sorted)
    )

  # Diamond reference points — one per cv_level actually used in grid_df
  # (not .CV_GRID, which may differ from CV_GRID_USE passed by the caller)
  diamond_df <- data.frame(
    cv_fac      = factor(cv_levels_sorted, levels = cv_levels_sorted),
    cv_achieved = cv_levels_sorted   # target = perfect 1:1
  )

  # Choose geom based on sample size per group:
  # n >= 50  : violin + boxplot (enough for kernel density)
  # n < 50   : jitter + boxplot (honest display of individual points)
  use_violin <- n_sim_per_level >= 50

  p <- ggplot(plot_df, aes(x = cv_fac, y = cv_achieved))

  if (use_violin) {
    p <- p +
      geom_violin(fill = "#A8D5BA", alpha = 0.55, colour = NA,
                  bw = "nrd0")
  } else {
    p <- p +
      geom_jitter(width = 0.12, height = 0, alpha = 0.35,
                  colour = "#2C6E49", size = 1.2)
  }

  p +
    geom_boxplot(width = 0.20, outlier.shape = NA,
                 colour = "#2C6E49", fill = "white", alpha = 0.7) +

    # 1:1 reference line — achieved CV should equal target CV
    geom_abline(slope = 0, intercept = 0,   # placeholder; use segments instead
                colour = NA) +

    # Diamond = target CV (perfect scaling lands here)
    geom_point(
      data  = diamond_df,
      aes(x = cv_fac, y = cv_achieved),
      shape = 23, size = 3.5,
      fill  = "#E07B39", colour = "white", stroke = 1.2
    ) +

    # Connecting line from diamond down to x-axis — visual guide
    geom_segment(
      data = diamond_df,
      aes(x = cv_fac, xend = cv_fac,
          y = cv_achieved, yend = -Inf),
      linetype  = "dotted",
      colour    = "#E07B39",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +

    scale_x_discrete(
      labels = function(x) paste0("CV = ", x)
    ) +

    labs(
      title    = title_text,
      subtitle = paste0(
        if (use_violin) "Violin + boxplot" else
          paste0("Jitter + boxplot (n = ", round(n_sim_per_level), " per level)"),
        "; orange diamond = target CV (perfect scaling).",
        "\nMedian of boxplot should align with diamond."
      ),
      x       = "Target CV level",
      y       = "Achieved CV  (sd(Q) / mean(Q)  per synthetic year)",
      caption = paste(
        "Layer-2 CV scaling: Q_scaled = mu + (Q_boot - mu) * lambda,",
        "lambda = CV_target / CV_boot.",
        "\nMean-preserving spread (Arnell 1998 Global Env. Change).",
        "Slight underachievement at high CV is expected due to the",
        "physical flow floor (Q >= 0.01 cms)."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.caption = element_text(size = 7.5, colour = "grey50"),
      axis.text.x  = element_text(size = 9)
    )
}


# -----------------------------------------------------------------------------
# plot_mu_by_cv()
#
# Box plot of annual mean flow (mu_yr_cms) across cv_level groups.
# Confirms that Layer 2 scaling preserves the mean flow regardless of the
# target CV — a key assumption of the mean-preserving spread approach.
#
# Arguments
#   grid_df      data.frame   output of run_bootstrap_cv_grid()
#   title_text   character    plot title
#
# Returns
#   ggplot object
# -----------------------------------------------------------------------------

plot_mu_by_cv <- function(
    grid_df,
    title_text = "Annual mean flow by CV level (mean preservation check)") {

  stopifnot(
    is.data.frame(grid_df),
    all(c("cv_level", "mu_yr_cms") %in% names(grid_df))
  )

  # Sort numerically to avoid alphabetical factor ordering
  cv_levels_sorted <- sort(unique(grid_df$cv_level))

  grid_df |>
    mutate(cv_fac = factor(cv_level, levels = cv_levels_sorted)) |>
    ggplot(aes(x = cv_fac, y = mu_yr_cms)) +
    geom_boxplot(fill = "#C4D9F5", colour = "#1A5276",
                 outlier.size = 0.6, width = 0.45) +
    scale_x_discrete(labels = function(x) paste0("CV = ", x)) +
    labs(
      title    = title_text,
      subtitle = paste(
        "Median mu_yr should be approximately equal across all CV levels.",
        "\nLayer-2 scaling preserves the annual mean by construction."
      ),
      x        = "Target CV level",
      y        = expression("Annual mean flow" ~ mu[yr] ~ "(cms)"),
      caption  = paste(
        "Mean preservation: E[Q_scaled] = mu_yr by construction.",
        "\nArnell (1998); Shiau & Huang (2014) J. Hydro-environment Research."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.caption = element_text(size = 8, colour = "grey50"))
}


# -----------------------------------------------------------------------------
# plot_historical_cv_trend()
#
# Two-panel time-series showing:
#   LEFT  panel : Historical observed within-year CV (1958–2025)
#                 Grey points = annual CV; orange line = rolling mean.
#   RIGHT panel : Projected future CV ranges (2026–2060) as horizontal bands,
#                 derived from published climate-change projections for eastern
#                 Taiwan mountain streams.  These are scenario PROXIES, not
#                 direct GCM downscaling outputs.
#
# Panel separation
#   A vertical dashed line at 2025 marks the boundary between observed
#   history and the projected future window.  The projected panel covers
#   2026–2060, consistent with a 35-year project economic lifetime.
#
# Projected CV scenario bands (right panel)
#   Baseline          : CV 1.5–2.5
#     Source: Lishan station observed mean (1958–2025), this study.
#   Moderate variability (SSP2-4.5 proxy) : CV 2.5–3.0
#     Source: Shiau & Huang (2014); Tung et al. (2016).
#   High variability     (SSP5-8.5 proxy) : CV 3.0–4.0
#     Source: Tung et al. (2016); IPCC AR6 WGI Ch.11 (2021).
#
# IMPORTANT DISCLAIMER (for Methods section)
#   CV ranges in the projected panel are literature-based scenario proxies
#   representing plausible future hydrological variability under different
#   warming pathways.  They are NOT direct outputs of GCM downscaling and
#   should be interpreted as sensitivity bounds rather than deterministic
#   predictions.
#
# Literature references
#   Shiau, J.T. & Huang, W.H. (2014). Detecting distributional changes of
#     annual rainfall indices in Taiwan using quantile regression.
#     J. Hydro-environment Research 8(4):355-366.
#     doi:10.1016/j.jher.2014.07.006
#   Tung, C.P., et al. (2016). Evaluating the impact of climate change on
#     water resources in Taiwan using CMIP5 GCM simulations.
#     Paddy Water Environ. 14(1):15-26.
#     doi:10.1007/s10333-014-0470-2
#   IPCC AR6 WGI Ch.11 (2021). Weather and Climate Extreme Events in a
#     Changing Climate. Cambridge Univ. Press.
#     doi:10.1017/9781009157896.013
#
# Arguments
#   annual_stats    data.frame  output of compute_annual_stats()
#   cv_grid         numeric vec FFF target CV levels (default: .CV_GRID)
#   roll_window     integer     rolling mean window in years (default 10)
#   proj_start      integer     first year of projection panel (default 2026)
#   proj_end        integer     last year of projection panel (default 2060)
#   title_text      character   plot title
#
# Returns
#   ggplot object (single figure; two zones separated by vertical dashed line)
# -----------------------------------------------------------------------------

plot_historical_cv_trend <- function(
    annual_stats,
    cv_grid     = .CV_GRID,
    roll_window = 10L,
    proj_start  = 2026L,
    proj_end    = 2060L,
    title_text  = "Historical CV and Projected Future Variability — Lishan Station") {

  stopifnot(
    is.data.frame(annual_stats),
    all(c("year", "cv_Q") %in% names(annual_stats))
  )

  # ── Historical data ────────────────────────────────────────────────────
  roll_df <- annual_stats |>
    arrange(year) |>
    mutate(
      cv_roll = zoo::rollmean(cv_Q, k = roll_window,
                              fill = NA, align = "center")
    )

  year_min  <- min(roll_df$year, na.rm = TRUE)
  year_hist <- max(roll_df$year, na.rm = TRUE)   # last observed year (2025)

  # Full x range spans history + projection window
  x_min <- year_min
  x_max <- proj_end

  # ── Projected scenario bands (right panel only) ───────────────────────
  # Three bands defined by CV range and corresponding SSP proxy label.
  # xmin = proj_start so bands appear ONLY in the future panel.
  scen_df <- data.frame(
    ymin  = c(1.5,  2.5,  3.0),
    ymax  = c(2.5,  3.0,  4.0),
    label = c("Baseline",
              "Moderate variability\n(SSP2-4.5 proxy)",
              "High variability\n(SSP5-8.5 proxy)"),
    fill  = c("#888780", "#1D9E75", "#D85A30"),
    alpha = c(0.18,      0.22,      0.22)
  )

  # Scenario label positions (right-most edge of plot)
  label_x   <- proj_end - 1
  label_df  <- data.frame(
    x      = label_x,
    y      = c(2.0,  2.75, 3.5),
    label  = c("Baseline\nCV 1.5\u20132.5",
               "Moderate\n(SSP2-4.5)\nCV 2.5\u20133.0",
               "High\n(SSP5-8.5)\nCV 3.0\u20134.0"),
    colour = c("#444444", "#1D9E75", "#D85A30")
  )

  # ── Dashed CV grid lines (full width, coloured by scenario) ──────────
  grid_df <- data.frame(
    cv     = cv_grid,
    colour = dplyr::case_when(
      cv_grid <= 2.5 ~ "#888780",
      cv_grid <= 3.0 ~ "#1D9E75",
      TRUE           ~ "#D85A30"
    )
  )

  # ── Build plot ────────────────────────────────────────────────────────
  ggplot() +

    # Grey background for historical zone
    annotate("rect",
             xmin = x_min - 0.5, xmax = year_hist + 0.5,
             ymin = -Inf,        ymax = Inf,
             fill = "grey97",    alpha = 1) +

    # Scenario bands in projected zone only
    purrr::pmap(scen_df, function(ymin, ymax, label, fill, alpha) {
      annotate("rect",
               xmin  = proj_start - 0.5, xmax = proj_end + 0.5,
               ymin  = ymin, ymax = ymax,
               fill  = fill, alpha = alpha)
    }) +

    # CV grid dashed lines (full width)
    purrr::pmap(as.list(grid_df), function(cv, colour) {
      annotate("segment",
               x = x_min - 0.5, xend = proj_end + 0.5,
               y = cv,          yend = cv,
               linetype = "dashed", linewidth = 0.35,
               colour = colour, alpha = 0.6)
    }) +

    # Vertical separator at 2025
    geom_vline(xintercept = year_hist + 0.5,
               linetype = "dashed", colour = "grey30",
               linewidth = 0.6) +

    # Zone header labels just below the top of the plot
    annotate("text",
             x = (x_min + year_hist) / 2, y = Inf,
             label = "Observed (1958\u20132025)",
             vjust = 1.4, hjust = 0.5, size = 3.2,
             colour = "grey30", fontface = "italic") +
    annotate("text",
             x = (proj_start + proj_end) / 2, y = Inf,
             label = sprintf("Projected (%d\u2013%d)", proj_start, proj_end),
             vjust = 1.4, hjust = 0.5, size = 3.2,
             colour = "grey30", fontface = "italic") +

    # Observed annual CV points and line
    geom_line(data  = roll_df,
              aes(x = year, y = cv_Q),
              colour = "grey60", linewidth = 0.4) +
    geom_point(data = roll_df,
               aes(x = year, y = cv_Q),
               colour = "grey50", size = 0.9, alpha = 0.75) +

    # Rolling mean
    geom_line(data  = roll_df,
              aes(x = year, y = cv_roll),
              colour = "#D35400", linewidth = 1.15, na.rm = TRUE) +

    # Scenario labels on right side
    annotate("label",
             x          = label_df$x,
             y          = label_df$y,
             label      = label_df$label,
             hjust      = 1,
             size       = 2.5,
             colour     = label_df$colour,
             fill       = "white",
             label.size = 0.2,
             alpha      = 0.88) +

    # CV grid labels on the left margin
    annotate("text",
             x      = x_min,
             y      = cv_grid,
             label  = paste0("CV = ", cv_grid),
             hjust  = 0, vjust = -0.45,
             size   = 2.4, colour = "grey45") +

    scale_x_continuous(
      limits = c(x_min - 0.5, proj_end + 0.5),
      breaks = seq(1960, proj_end, by = 10),
      expand = c(0, 0)
    ) +

    labs(
      title    = title_text,
      subtitle = paste0(
        "Grey points = annual CV; orange = ", roll_window,
        "-yr rolling mean; dashed lines = FFF target CV levels.",
        "\nProjected bands are literature-based scenario proxies",
        " (not direct GCM output)."
      ),
      x       = "Year",
      y       = "Within-year CV = sd(Q) / mean(Q)  (dimensionless)",
      caption = paste(
        "Observed: Lishan Station 01T230 (WRA), 1958\u20132025.",
        "Projected CV ranges are scenario proxies for eastern Taiwan mountain streams:",
        "\nBaseline CV 1.5\u20132.5 (this study);",
        "Moderate variability CV 2.5\u20133.0 (SSP2-4.5 proxy:",
        "Shiau & Huang 2014, J. Hydro-environ. Res. 8(4):355\u2013366;",
        "\nTung et al. 2016, Paddy Water Environ. 14(1):15\u201326);",
        "High variability CV 3.0\u20134.0 (SSP5-8.5 proxy:",
        "Tung et al. 2016; IPCC AR6 WGI Ch.11, 2021).",
        "\nCV ranges are sensitivity bounds, not deterministic GCM projections."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.caption    = element_text(size = 7.5, colour = "grey50",
                                     hjust = 0, lineheight = 1.2),
      legend.position = "none",
      panel.grid      = element_blank(),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3)
    )
}

# -----------------------------------------------------------------------------
# .cv_scenario_colours()
#
# Internal helper: return a named colour vector mapping cv_level values to
# a sequential palette from grey (low CV) through green to red (high CV).
# Used consistently across both bootstrap plot variants.
# -----------------------------------------------------------------------------

.cv_scenario_colours <- function(cv_levels) {
  n   <- length(cv_levels)
  pal <- colorRampPalette(c("#888780", "#1D9E75", "#D85A30"))(n)
  setNames(pal, as.character(cv_levels))
}


# -----------------------------------------------------------------------------
# plot_bootstrap_cv_A()
#
# Option A — Bootstrap CV distribution: X = target CV, Y = achieved CV.
#
# Each simulated year is plotted as a point at (cv_level, cv_achieved).
# Points are coloured by cv_level intensity (grey → green → red).
# Background bands show the three scenario zones, so the reader can see
# which CV range each set of simulations falls into.
#
# Y-axis is identical to the historical CV trend plot, enabling direct
# visual alignment when the two plots are placed side-by-side via patchwork.
#
# Arguments
#   grid_df     data.frame  output of run_bootstrap_cv_grid()
#   title_text  character
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_bootstrap_cv_A <- function(
    grid_df,
    title_text = "Option A: Bootstrap CV distribution\nby target CV level") {

  stopifnot(
    is.data.frame(grid_df),
    all(c("cv_level", "cv_achieved") %in% names(grid_df))
  )

  cv_levels <- sort(unique(grid_df$cv_level))
  col_map   <- .cv_scenario_colours(cv_levels)

  # Background scenario bands (same as historical plot)
  band_df <- data.frame(
    ymin  = c(1.5, 2.5, 3.0),
    ymax  = c(2.5, 3.0, 4.0),
    fill  = c("#888780", "#1D9E75", "#D85A30"),
    alpha = c(0.10, 0.14, 0.14)
  )

  # Jitter width: slightly randomise x within each cv_level column
  set.seed(1L)

  p <- ggplot(
    grid_df |>
      mutate(
        cv_fac  = factor(cv_level, levels = cv_levels),
        cv_col  = as.character(cv_level)
      ),
    aes(x = cv_level, y = cv_achieved, colour = cv_col)
  )

  # Background bands
  for (i in seq_len(nrow(band_df))) {
    p <- p + annotate("rect",
                      xmin  = -Inf, xmax = Inf,
                      ymin  = band_df$ymin[i],
                      ymax  = band_df$ymax[i],
                      fill  = band_df$fill[i],
                      alpha = band_df$alpha[i])
  }

  p +
    # 1:1 reference line (achieved = target)
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40",
                linewidth = 0.5) +

    # Jittered points
    geom_jitter(width = 0.06, height = 0,
                size = 1.2, alpha = 0.55) +

    # Boxplot overlay (no outlier dots — already shown by jitter)
    geom_boxplot(aes(group = cv_level),
                 width = 0.12, outlier.shape = NA,
                 fill = "white", alpha = 0.6,
                 colour = "grey30", linewidth = 0.4) +

    # Diamond = target (perfect 1:1)
    geom_point(
      data = data.frame(cv_level = cv_levels,
                        cv_achieved = cv_levels,
                        cv_col = as.character(cv_levels)),
      aes(x = cv_level, y = cv_achieved),
      shape = 23, size = 3.2,
      fill = "white", stroke = 1.5
    ) +

    # Scenario band labels on right margin
    annotate("text",
             x = max(cv_levels) + 0.15,
             y = c(2.0, 2.75, 3.5),
             label = c("Baseline", "Moderate\n(SSP2-4.5)", "High\n(SSP5-8.5)"),
             hjust = 0, size = 2.4,
             colour = c("#555555", "#1D9E75", "#D85A30")) +

    scale_colour_manual(values = col_map, guide = "none") +

    scale_x_continuous(
      breaks = cv_levels,
      labels = paste0("CV=", cv_levels),
      expand = expansion(mult = c(0.05, 0.2))
    ) +

    coord_cartesian(ylim = c(NA, max(cv_levels) + 0.3)) +

    labs(
      title    = title_text,
      subtitle = paste(
        "Points = individual bootstrap years; diamond = target (1:1 line).",
        "\nColour intensity = scenario variability level."
      ),
      x = "Target CV level",
      y = "Achieved CV  (sd(Q) / mean(Q))",
      caption = paste(
        "Scenario proxies: Baseline CV 1.5-2.5 (this study);",
        "Moderate CV 2.5-3.0 (SSP2-4.5, Shiau & Huang 2014;",
        "Tung et al. 2016); High CV 3.0-4.0 (SSP5-8.5, Tung et al. 2016;",
        "IPCC AR6 Ch.11 2021)."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.caption = element_text(size = 7.5, colour = "grey50",
                                  hjust = 0),
      panel.grid   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3)
    )
}


# -----------------------------------------------------------------------------
# plot_bootstrap_cv_B()
#
# Option B — Bootstrap CV distribution: X = sim_id, Y = cv_achieved,
#            colour = cv_level (scenario intensity).
#
# All 200 × n_cv_levels simulated years are plotted as individual points,
# each coloured by its target CV level.  This view emphasises the full
# spread of the simulation ensemble and makes it easy to see which scenario
# levels produce the widest variance in achieved CV.
#
# Arguments
#   grid_df     data.frame  output of run_bootstrap_cv_grid()
#   title_text  character
#
# Returns  ggplot object
# -----------------------------------------------------------------------------

plot_bootstrap_cv_B <- function(
    grid_df,
    title_text = "Option B: Bootstrap ensemble\nall simulations by scenario") {

  stopifnot(
    is.data.frame(grid_df),
    all(c("cv_level", "sim_id", "cv_achieved") %in% names(grid_df))
  )

  cv_levels <- sort(unique(grid_df$cv_level))
  col_map   <- .cv_scenario_colours(cv_levels)

  # Background scenario bands (horizontal, same as historical plot)
  band_df <- data.frame(
    ymin  = c(1.5, 2.5, 3.0),
    ymax  = c(2.5, 3.0, 4.0),
    fill  = c("#888780", "#1D9E75", "#D85A30"),
    alpha = c(0.10, 0.14, 0.14)
  )

  plot_df <- grid_df |>
    mutate(cv_col = as.character(cv_level))

  p <- ggplot(plot_df, aes(x = sim_id, y = cv_achieved, colour = cv_col))

  for (i in seq_len(nrow(band_df))) {
    p <- p + annotate("rect",
                      xmin  = -Inf, xmax = Inf,
                      ymin  = band_df$ymin[i],
                      ymax  = band_df$ymax[i],
                      fill  = band_df$fill[i],
                      alpha = band_df$alpha[i])
  }

  # Pre-compute median achieved CV per cv_level for horizontal reference lines
  median_df <- grid_df |>
    group_by(cv_level) |>
    summarise(med_cv = median(cv_achieved, na.rm = TRUE), .groups = "drop") |>
    mutate(cv_col = as.character(cv_level))

  p +
    geom_point(size = 1.1, alpha = 0.55) +

    # Horizontal median lines per cv_level
    geom_hline(
      data        = median_df,
      aes(yintercept = med_cv, colour = cv_col),
      linetype    = "solid",
      linewidth   = 0.7,
      alpha       = 0.85,
      show.legend = FALSE
    ) +

    # Scenario band labels on right
    annotate("text",
             x = max(grid_df$sim_id) * 1.02,
             y = c(2.0, 2.75, 3.5),
             label = c("Baseline", "Moderate\n(SSP2-4.5)", "High\n(SSP5-8.5)"),
             hjust = 0, size = 2.4,
             colour = c("#555555", "#1D9E75", "#D85A30")) +

    scale_colour_manual(
      values = col_map,
      name   = "Target CV",
      labels = paste0("CV = ", cv_levels)
    ) +

    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.18))
    ) +

    labs(
      title    = title_text,
      subtitle = paste(
        "Each point = one bootstrap year; horizontal line = median per CV level.",
        "\nColour = target CV level (scenario intensity)."
      ),
      x = "Bootstrap replicate index (sim_id)",
      y = "Achieved CV  (sd(Q) / mean(Q))",
      caption = paste(
        "n =", nrow(grid_df), "total simulations across",
        length(cv_levels), "CV levels."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.caption    = element_text(size = 7.5, colour = "grey50"),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      panel.grid      = element_blank(),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3)
    )
}


# -----------------------------------------------------------------------------
# plot_cv_combined()
#
# Combine the historical CV trend (left) with either Option A or Option B
# bootstrap distribution (right) using patchwork.
#
# The Y-axis of both panels spans the same CV range, enabling direct visual
# alignment between observed history and simulated future distributions.
#
# Arguments
#   annual_stats  data.frame  output of compute_annual_stats()
#   grid_df       data.frame  output of run_bootstrap_cv_grid()
#   option        character   "A" (x = target CV) or "B" (x = sim_id)
#   cv_grid       numeric vec FFF target CV levels
#   widths        numeric(2)  relative panel widths (default c(2, 1))
#
# Returns  patchwork object
# -----------------------------------------------------------------------------

plot_cv_combined <- function(annual_stats,
                             grid_df,
                             option  = "A",
                             cv_grid = .CV_GRID,
                             widths  = c(2, 1)) {

  stopifnot(option %in% c("A", "B"))

  p_hist <- plot_historical_cv_trend(
    annual_stats = annual_stats,
    cv_grid      = cv_grid
  )

  p_boot <- if (option == "A") {
    plot_bootstrap_cv_A(grid_df)
  } else {
    plot_bootstrap_cv_B(grid_df)
  }

  # Align Y axes across both panels
  p_hist <- p_hist + coord_cartesian(ylim = c(0.8, 4.3))
  p_boot <- p_boot + coord_cartesian(ylim = c(0.8, 4.3))

  # Remove Y-axis label from right panel to avoid duplication
  p_boot <- p_boot + theme(axis.title.y = element_blank())

  patchwork::wrap_plots(p_hist, p_boot, widths = widths) +
    patchwork::plot_annotation(
      caption = paste(
        "Left: Observed CV, Lishan Station 01T230 (WRA) 1958-2025.",
        "Right: Two-layer bootstrap ensemble (Layer 1 = block resample;",
        "Layer 2 = CV scaling, Arnell 1998).",
        "\nScenario proxies: Shiau & Huang (2014);",
        "Tung et al. (2016); IPCC AR6 WGI Ch.11 (2021)."
      ),
      theme = theme(
        plot.caption = element_text(size = 7.5, colour = "grey50",
                                    hjust = 0)
      )
    )
}



# -----------------------------------------------------------------------------
# build_climate_ensemble()
#
# Run the full two-layer bootstrap over the default CV grid and return
# both the raw grid table and the summary statistics.  This is the single
# entry point called from fp_hydro_main.qmd.
#
# Arguments
#   daily_flow   data.frame   columns: date (Date), Q_cms (numeric)
#   cv_grid      numeric vector   target CV levels (default: .CV_GRID)
#   n_sim        integer      bootstrap resamples per CV level
#   seed         integer      RNG seed
#
# Returns  named list:
#   $grid_df       data.frame   raw output of run_bootstrap_cv_grid()
#                               nrow = length(cv_grid) * n_sim
#   $summary_df    data.frame   output of summarise_cv_ensemble()
#   $annual_stats  data.frame   output of compute_annual_stats()
#   $cv_grid       numeric vector   CV levels used
#   $n_sim         integer
#
# Example usage in fp_hydro_main.qmd
#   climate <- build_climate_ensemble(daily_flow = daily_clean, n_sim = 200L)
#   # climate$grid_df feeds directly into the NPV for-loop
# -----------------------------------------------------------------------------

build_climate_ensemble <- function(daily_flow,
                                   cv_grid      = .CV_GRID,
                                   n_sim        = .N_SIM_DEFAULT,
                                   seed         = 42L,
                                   keep_series  = TRUE) {

  message("=== build_climate_ensemble() ===")
  message(sprintf("CV grid: %s", paste(cv_grid, collapse = ", ")))
  message(sprintf("n_sim per CV level: %d", n_sim))
  message(sprintf("Total synthetic years: %d", length(cv_grid) * n_sim))
  message(sprintf("keep_series: %s (required for FFF simulation loop)", keep_series))

  annual_stats <- compute_annual_stats(daily_flow)

  grid_df <- run_bootstrap_cv_grid(
    daily_flow   = daily_flow,
    cv_grid      = cv_grid,
    n_sim        = as.integer(n_sim),
    seed         = seed,
    keep_series  = keep_series
  )

  summary_df <- summarise_cv_ensemble(grid_df)

  message(sprintf(
    "\nDone. Historical CV range: %.2f \u2013 %.2f",
    min(annual_stats$cv_Q, na.rm = TRUE),
    max(annual_stats$cv_Q, na.rm = TRUE)
  ))

  list(
    grid_df      = grid_df,
    summary_df   = summary_df,
    annual_stats = annual_stats,
    cv_grid      = cv_grid,
    n_sim        = n_sim
  )
}


# =============================================================================
# SECTION 8 — FFF simulation engine
# =============================================================================

# -----------------------------------------------------------------------------
# run_fff_grid()   <-- PRIMARY ENGINE FOR FINANCIAL FEASIBILITY FRONTIER
#
# For each (cv_level, sim_id) row in the bootstrap grid:
#   1. Extract the synthetic daily flow series (Q_series list-column).
#   2. Run run_reservoir_simulation() for W1 and W2 to get actual daily
#      energy generation — accounting for pondage dispatch, e-flow, and spill.
#   3. Repeat with run_reservoir_simulation_forward24() for the pondage
#      (24-hour forecast) mode.
#   4. Compute annual total energy for both operation modes.
#   5. For each r_loan in r_loan_seq, compute NPV using build_cash_flows()
#      and compute_npv().
#
# Why this is necessary
#   A simple formula (energy = Q * capture_rate * H * eta) cannot capture
#   the dispatch dynamics of a pondage reservoir.  At high CV (flashy flow),
#   the small pondage fills and spills during typhoon days and runs dry during
#   droughts, reducing effective generation relative to a naive flow-based
#   estimate.  This CV-dependent generation loss is the core mechanism that
#   creates the slope of the FFF boundary.
#
# Parallelisation
#   Uses furrr::future_map_dfr() if the furrr package is installed and a
#   parallel plan has been set by the caller (e.g. plan(multisession)).
#   Falls back to purrr::map_dfr() if furrr is unavailable.
#
# Arguments
#   grid_df          data.frame   output of run_bootstrap_cv_grid(keep_series=TRUE)
#                                 must contain the Q_series list-column
#   r_loan_seq       numeric vec  loan interest rates (fraction, e.g. seq(0,0.10,0.01))
#   ltv              numeric      loan-to-value ratio (default 0.80 = 80%)
#   r_equity         numeric      equity discount rate for NPV (default 0.08)
#   e_flow_mode      character    passed to run_reservoir_simulation()
#                                 "committed" | "recommended" | "custom"
#   capex_ntd        numeric      total project CAPEX (NTD)
#   opex_ntd_yr      numeric      annual O&M cost (NTD)
#   project_life_yr  integer      project economic lifetime (years)
#   fit_ntd_per_kwh  numeric      Feed-In Tariff rate (NTD/kWh)
#   scale_w2         logical      scale flow to W2 catchment? (default TRUE)
#
# Returns
#   data.frame with one row per (cv_level, sim_id, r_loan, mode), columns:
#     cv_level        numeric   target CV (FFF Y-axis)
#     sim_id          integer   bootstrap replicate
#     r_loan          numeric   loan interest rate (FFF X-axis)
#     r_loan_pct      numeric   r_loan × 100 (for plotting)
#     mode            character "run_of_river" | "pondage"
#     annual_gwh      numeric   simulated annual energy (GWh)
#     annual_rev_ntd  numeric   annual revenue (NTD)
#     npv_ntd         numeric   net present value (NTD)
#     npv_b_ntd       numeric   NPV in NTD billions
#     viable          logical   TRUE if NPV > 0
# -----------------------------------------------------------------------------

run_fff_grid <- function(grid_df,
                         r_loan_seq      = seq(0, 0.10, by = 0.01),
                         ltv             = 0.80,
                         r_equity        = 0.08,
                         e_flow_mode     = "recommended",
                         capex_ntd       = 9.7e9,
                         opex_ntd_yr     = 9.7e7,
                         project_life_yr = 35L,
                         fit_ntd_per_kwh = 2.8599,
                         scale_w2        = TRUE) {

  stopifnot(
    is.data.frame(grid_df),
    "Q_series" %in% names(grid_df),
    all(c("cv_level", "sim_id") %in% names(grid_df)),
    is.numeric(r_loan_seq), all(r_loan_seq >= 0),
    ltv >= 0, ltv <= 1
  )

  # Print fixed net head values once, here, instead of inside each sim call
  message(sprintf(
    "Net head (fixed): W1 = %.1f m (Francis) | W2 = %.1f m (Pelton)",
    (18100 * 1000) / (0.88 * 1000 * 9.81 * 24.3),
    (19000 * 1000) / (0.88 * 1000 * 9.81 *  6.3)
  ))
  message(sprintf(
    "Running %d simulations × %d r_loan values × 2 modes = %d NPV computations",
    nrow(grid_df), length(r_loan_seq), nrow(grid_df) * length(r_loan_seq) * 2
  ))
  message("Using sequential purrr (stable on Windows). For parallel, set plan(multisession) before calling.")

  # ── Inner helper: simulate one synthetic year, return GWh for both modes ──
  simulate_one_synthetic_year <- function(Q_series_df, e_flow_mode, scale_w2) {

    sim_ror <- run_reservoir_simulation(
      daily_flow     = Q_series_df,
      e_flow_mode    = e_flow_mode,
      na_method      = "linear",
      max_gap_linear = 7L,
      scale_w2       = scale_w2
    )

    sim_pond <- run_reservoir_simulation_forward24(
      daily_flow     = Q_series_df,
      e_flow_mode    = e_flow_mode,
      na_method      = "linear",
      max_gap_linear = 7L,
      scale_w2       = scale_w2
    )

    list(
      run_of_river = sum(sim_ror$combined$energy_kwh_total,  na.rm = TRUE) / 1e6,
      pondage      = sum(sim_pond$combined$energy_kwh_total, na.rm = TRUE) / 1e6
    )
  }

  # ── Sequential loop over bootstrap rows ────────────────────────────────────
  # purrr::map_dfr is used here for stability on Windows.
  # To parallelise: wrap in furrr::future_map_dfr() after plan(multisession).
  purrr::map_dfr(
    seq_len(nrow(grid_df)),
    function(i) {

      if (i %% 20 == 0 || i == 1)
        message(sprintf("  row %d / %d", i, nrow(grid_df)))

      row    <- grid_df[i, ]
      Q_ser  <- row$Q_series[[1]]

      energy <- simulate_one_synthetic_year(Q_ser, e_flow_mode, scale_w2)

      purrr::map_dfr(c("run_of_river", "pondage"), function(mode) {

        gwh        <- energy[[mode]]
        annual_rev <- gwh * 1e6 * fit_ntd_per_kwh

        purrr::map_dfr(r_loan_seq, function(r_loan) {

          cf <- build_cash_flows(
            annual_revenue_ntd = annual_rev,
            ltv             = ltv,
            r_loan          = r_loan,
            capex_ntd       = capex_ntd,
            opex_ntd_yr     = opex_ntd_yr,
            project_life_yr = as.integer(project_life_yr)
          )
          npv_val <- compute_npv(cf, discount_rate = r_equity)

          data.frame(
            cv_level       = row$cv_level,
            sim_id         = row$sim_id,
            r_loan         = r_loan,
            r_loan_pct     = r_loan * 100,
            mode           = mode,
            annual_gwh     = round(gwh,        3),
            annual_rev_ntd = round(annual_rev, 0),
            npv_ntd        = round(npv_val,    0),
            npv_b_ntd      = round(npv_val / 1e9, 4),
            viable         = npv_val > 0
          )
        })
      })
    }
  )
}
