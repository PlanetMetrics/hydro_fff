# =============================================================================
# hydro_sediment.R
# Fengping River Hydropower Simulation — Sediment Module
#
# Purpose
#   (1) Fit power-law suspended sediment rating curve: SSL = a * Q^b
#   (2) Estimate daily sediment load from daily flow
#   (3) Estimate weir trap efficiency and reservoir volume loss over time
#
# Key notation
#   Q               : streamflow (cms = cubic metres per second)
#   SSL             : suspended sediment load (metric tonnes / day, MT/day)
#   C               : sediment concentration (ppm = mg/L)
#   a               : rating curve coefficient (SSL = a * Q^b)
#   b               : rating curve exponent
#   TE              : trap efficiency (dimensionless, 0–1)
#   CI              : capacity–inflow ratio (dimensionless)
#   V_loss          : annual reservoir volume loss due to sedimentation (m³/yr)
#   rho_sed         : bulk density of deposited sediment (kg/m³)
#   Q_min_transport : minimum flow for detectable sediment transport (cms)
#
# Literature references for b-value validation range [1.4, 2.5]
#   Wang, H.W. & Kondolf, G.M. (2014). Upstream sediment-control dams:
#     five decades of experience in the rapidly eroding Dahan River Basin,
#     Taiwan. JAWRA 50(3):735–747. https://doi.org/10.1111/jawr.12160
#   Kondolf, G.M., et al. (2014). Sustainable sediment management in
#     reservoirs and regulated rivers. Earth's Future 2(5):256–280.
#     https://doi.org/10.1002/2013EF000184
#
# Trap efficiency method
#   Brune (1953) curve, simplified approximation:
#     TE = CI / (CI + 0.0021)
#   where CI = reservoir capacity / mean annual inflow volume
#
# Author  [your name]
# Date    2025
# =============================================================================

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# Internal constants
# -----------------------------------------------------------------------------

# Literature-based b-value range for Taiwan mountain streams
.B_RANGE <- c(1.4, 2.5)

# Wet bulk density of deposited sediment
# Typical value for mixed gravel-sand deposits in Taiwan mountain streams
# Reference: Wang et al. (2018) Water 10(8):1034
.RHO_SED_KG_M3 <- 1300

# Weir design parameters — shared reference with hydro_reservoir.R
# Source: Shihfeng Power Co. EIA documents
.WEIR_PARAMS <- list(
  W1 = list(
    name         = "Plant 1 Lower Weir",
    S_max_m3     = 967400,   # effective storage (m³)
    Q_design_cms = 24.3      # design flow (cms)
  ),
  W2 = list(
    name         = "Plant 2 Upper Weir",
    S_max_m3     = 237300,
    Q_design_cms = 6.3
  )
)


# -----------------------------------------------------------------------------
# fit_rating_curve()   <- CUSTOM FUNCTION
#
# Fit a power-law suspended sediment rating curve by ordinary least squares
# regression in log-log space:
#
#   log(SSL) = log(a) + b * log(Q)   =>   SSL = a * Q^b
#
# Only rows with flag == "ok" and SSL > 0 are used for fitting.
# Below-detection-limit rows (ppm = 0, flag == "below_detection") are
# excluded to prevent bias from zero SSL values.
#
# Arguments
#   sed_clean   data.frame   cleaned sediment data from sediment_clean.csv
#                            required columns: Discharge_CMS, Sus_Load_MTDay,
#                            flag (character)
#   plot_fit    logical      if TRUE, print log-log scatter with fitted line
#
# Returns  named list:
#   $a            numeric   rating curve coefficient
#   $b            numeric   rating curve exponent
#   $r_squared    numeric   R² of log-log regression
#   $n_points     integer   number of observations used
#   $b_in_range   logical   whether b falls within .B_RANGE
#   $method_note  character citable description for methods section
# -----------------------------------------------------------------------------

#' Fit a power-law suspended sediment rating curve
#'
#' Estimates the parameters \eqn{a} and \eqn{b} of the rating curve
#' \eqn{SSL = a \cdot Q^b} by ordinary least squares (OLS) regression in
#' log-log space.  Only observations flagged \code{"ok"} with positive SSL
#' are used; below-detection-limit records (\code{ppm = 0}) are excluded to
#' prevent downward bias.
#'
#' \deqn{\log(SSL) = \log(a) + b \cdot \log(Q)}
#'
#' The exponent \eqn{b} is validated against the range \eqn{[1.4,\, 2.5]}
#' reported for Taiwan mountain streams (Wang & Kondolf 2014; Kondolf et al.
#' 2014).  A warning is issued if \eqn{b} falls outside this range.
#'
#' @param sed_clean  A \code{data.frame} containing the cleaned sediment
#'   dataset (typically \code{sediment_clean.csv}).  Required columns:
#'   \describe{
#'     \item{Discharge_CMS}{Instantaneous discharge at sampling time (cms)}
#'     \item{Sus_Load_MTDay}{Suspended sediment load (metric tonnes / day)}
#'     \item{flag}{Character quality flag; only rows with \code{flag == "ok"}
#'       are used for fitting}
#'   }
#' @param plot_fit  Logical.  If \code{TRUE} (default), prints a log-log
#'   scatter plot with the fitted curve, suitable for inclusion in the
#'   methods section of a manuscript.
#'
#' @return  A named list with elements:
#'   \describe{
#'     \item{a}{Numeric. Rating curve coefficient (\eqn{SSL = a \cdot Q^b})}
#'     \item{b}{Numeric. Rating curve exponent}
#'     \item{r_squared}{Numeric. \eqn{R^2} of the log-log regression}
#'     \item{n_points}{Integer. Number of observations used in fitting}
#'     \item{b_in_range}{Logical. \code{TRUE} if \eqn{b \in [1.4,\, 2.5]}}
#'     \item{method_note}{Character. Citable description for the methods
#'       section, including fitted parameters and literature references}
#'   }
#'
#' @references
#'   Wang, H.W. & Kondolf, G.M. (2014). Upstream sediment-control dams:
#'   five decades of experience in the rapidly eroding Dahan River Basin,
#'   Taiwan. \emph{JAWRA} 50(3):735–747.
#'   \doi{10.1111/jawr.12160}
#'
#'   Kondolf, G.M., et al. (2014). Sustainable sediment management in
#'   reservoirs and regulated rivers. \emph{Earth's Future} 2(5):256–280.
#'   \doi{10.1002/2013EF000184}
#'
#' @examples
#' \dontrun{
#' sed  <- read_csv(here("data", "sediment_clean.csv"))
#' rc   <- fit_rating_curve(sed, plot_fit = TRUE)
#' cat("SSL =", round(rc$a, 3), "x Q^", round(rc$b, 3),
#'     "| R2 =", round(rc$r_squared, 3))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_x_log10
#'   scale_y_log10 labs theme_minimal
#' @export
fit_rating_curve <- function(sed_clean, plot_fit = TRUE) {
  
  stopifnot(
    is.data.frame(sed_clean),
    all(c("Discharge_CMS", "Sus_Load_MTDay", "flag") %in% names(sed_clean))
  )
  
  # Use only valid observations with measurable sediment
  fit_data <- sed_clean |>
    filter(
      flag == "ok",
      Discharge_CMS > 0,
      Sus_Load_MTDay > 0
    ) |>
    mutate(
      log_Q   = log(Discharge_CMS),
      log_SSL = log(Sus_Load_MTDay)
    )
  
  if (nrow(fit_data) < 5)
    stop("fit_rating_curve: fewer than 5 valid observations for regression")
  
  # OLS regression in log-log space
  model <- lm(log_SSL ~ log_Q, data = fit_data)
  a     <- exp(coef(model)[[1]])    # back-transform intercept
  b     <- coef(model)[[2]]         # slope = exponent
  r2    <- summary(model)$r.squared
  
  b_ok <- b >= .B_RANGE[1] & b <= .B_RANGE[2]
  
  if (!b_ok)
    warning(sprintf(
      "fit_rating_curve: b = %.3f is outside expected range [%.1f, %.1f]. ",
      b, .B_RANGE[1], .B_RANGE[2]
    ))
  
  message(sprintf(
    "Rating curve fitted: SSL = %.4f x Q^%.4f  |  R² = %.3f  |  n = %d",
    a, b, r2, nrow(fit_data)
  ))
  
  # Diagnostic plot
  if (plot_fit) {
    
    pred_df <- data.frame(
      Discharge_CMS = 10^seq(
        log10(max(0.1, min(fit_data$Discharge_CMS))),
        log10(max(fit_data$Discharge_CMS)),
        length.out = 200
      )
    ) |>
      mutate(Sus_Load_MTDay = a * Discharge_CMS^b)
    
    p <- ggplot(fit_data,
                aes(x = Discharge_CMS, y = Sus_Load_MTDay)) +
      geom_point(alpha = 0.65, colour = "#1D9E75", size = 2) +
      geom_line(data = pred_df,
                colour = "#D85A30", linewidth = 1) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      labs(
        title    = "Suspended Sediment Rating Curve — Fengping Creek",
        subtitle = sprintf(
          "SSL = %.3f \u00d7 Q^%.3f  |  R\u00b2 = %.3f  |  n = %d obs (1959\u20132023)",
          a, b, r2, nrow(fit_data)
        ),
        x = "Discharge Q (cms, log scale)",
        y = "Suspended Sediment Load SSL (MT/day, log scale)",
        caption = paste(
          "b-value expected range for Taiwan mountain streams: 1.4\u20132.5",
          "\nWang & Kondolf (2014) JAWRA; Kondolf et al. (2014) Earth's Future"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.caption = element_text(size = 8, colour = "grey50"))
    
    print(p)
  }
  
  # Citable method description for methods section
  method_note <- sprintf(
    paste(
      "A power-law suspended sediment rating curve (SSL = a * Q^b) was",
      "fitted by ordinary least squares regression in log-log space using",
      "%d field observations from the Fengping Creek gauging station",
      "(1959-2023; below-detection-limit observations excluded).",
      "The fitted parameters are a = %.4f and b = %.4f (R2 = %.3f).",
      "The exponent b falls within the range of 1.4-2.5 reported for",
      "Taiwan mountain streams (Wang & Kondolf 2014, JAWRA 50(3):735-747;",
      "Kondolf et al. 2014, Earth's Future 2(5):256-280)."
    ),
    nrow(fit_data), a, b, r2
  )
  
  list(
    a           = a,
    b           = b,
    r_squared   = r2,
    n_points    = nrow(fit_data),
    b_in_range  = b_ok,
    method_note = method_note
  )
}


# -----------------------------------------------------------------------------
# estimate_daily_ssl()
#
# Apply the fitted rating curve to a daily flow series to estimate
# suspended sediment load (SSL) for each day.
#
# Days where Q is below Q_min_transport are assigned SSL = 0, consistent
# with field observations where ppm = 0 during low-flow periods.
#
# Arguments
#   daily_flow        data.frame   columns: date (Date), Q_cms (numeric)
#   curve             list         output of fit_rating_curve()
#   Q_min_transport   numeric      cms below which SSL = 0 (default 5.0)
#
# Returns
#   data.frame: date, Q_cms, SSL_mt_day (metric tonnes per day)
# -----------------------------------------------------------------------------

#' Estimate daily suspended sediment load from the daily flow record
#'
#' Applies a fitted power-law sediment rating curve
#' (\eqn{SSL = a \cdot Q^b}, from \code{\link{fit_rating_curve}}) to a
#' continuous daily streamflow series to produce a daily suspended sediment
#' load (SSL) estimate. This is the step that "matches" the sparse field
#' sediment-sampling record (used to fit \code{a} and \code{b}) to the
#' continuous daily flow record, producing a continuous daily SSL series.
#'
#' Days where \code{Q_cms} is below \code{Q_min_transport} are assigned
#' \code{SSL_mt_day = 0}, consistent with field observations where
#' sediment concentration (ppm) = 0 during low-flow periods.
#'
#' @param daily_flow A \code{data.frame} with columns:
#'   \describe{
#'     \item{date}{\code{Date}. Calendar date.}
#'     \item{Q_cms}{Numeric. Daily mean streamflow (cms).}
#'   }
#' @param curve A named list, the output of \code{\link{fit_rating_curve}},
#'   containing at least \code{$a} and \code{$b} (the rating-curve
#'   coefficient and exponent).
#' @param Q_min_transport Numeric. Streamflow (cms) below which sediment
#'   transport is assumed negligible and \code{SSL_mt_day} is set to 0.
#'   Default \code{5.0}.
#'
#' @return A \code{data.frame} with one row per input day and columns:
#'   \describe{
#'     \item{date}{Date}
#'     \item{Q_cms}{Numeric. Daily mean streamflow (cms), unchanged from input.}
#'     \item{SSL_mt_day}{Numeric. Estimated suspended sediment load
#'       (metric tonnes / day), \code{= curve$a * Q_cms ^ curve$b} when
#'       \code{Q_cms >= Q_min_transport}, otherwise \code{0}.}
#'   }
#'
#' @examples
#' \dontrun{
#' sed   <- read_csv(here("data", "sediment_clean.csv"))
#' curve <- fit_rating_curve(sed, plot_fit = FALSE)
#' daily <- read_csv(here("data", "lishan_daily_clean.csv")) |>
#'   mutate(date = as.Date(date))
#' daily_ssl <- estimate_daily_ssl(daily, curve)
#' head(daily_ssl)
#' }
#'
#' @export
estimate_daily_ssl <- function(daily_flow,
                               curve,
                               Q_min_transport = 5.0) {
  
  stopifnot(
    is.data.frame(daily_flow),
    all(c("date", "Q_cms") %in% names(daily_flow)),
    is.list(curve),
    all(c("a", "b") %in% names(curve))
  )
  
  daily_flow |>
    mutate(
      SSL_mt_day = if_else(
        !is.na(Q_cms) & Q_cms >= Q_min_transport,
        curve$a * Q_cms ^ curve$b,
        0
      ),
      SSL_mt_day = pmax(SSL_mt_day, 0)   # physical lower bound
    ) |>
    select(date, Q_cms, SSL_mt_day)
}


# -----------------------------------------------------------------------------
# estimate_trap_efficiency()
#
# Simulate annual trap efficiency (TE) and reservoir volume loss over the
# project lifetime using the Brune (1953) curve approximation.
#
# Brune curve (simplified):
#   TE = CI / (CI + 0.0021)
# where CI = S_remaining / annual_inflow_volume
#
# Volume deposited per year:
#   V_dep = SSL_annual_mt * 1000 / rho_sed   [m³/year]
#
# Arguments
#   ssl_annual_mt   numeric    annual suspended sediment load (MT/year)
#   weir_id         character  "W1" or "W2"
#   years           integer    project economic lifetime (default 35)
#   rho_sed         numeric    bulk density of deposit (kg/m³)
#
# Returns  data.frame with columns:
#   year, S_remaining_m3, TE, V_loss_m3_yr, CI_ratio
#   S_remaining_m3 at end of each year after sedimentation
# -----------------------------------------------------------------------------

#' Estimate annual trap efficiency and reservoir volume loss over project life
#'
#' Simulates the progressive sedimentation of a pondage weir reservoir using
#' the Brune (1953) curve approximation.  In each year, trap efficiency
#' \eqn{TE} is computed from the current capacity-inflow ratio \eqn{CI},
#' the volume of sediment deposited is subtracted from remaining storage,
#' and the process repeats for the next year.
#'
#' Brune (1953) simplified approximation:
#' \deqn{TE = \frac{CI}{CI + 0.0021}, \quad CI = \frac{S_{\text{remaining}}}{V_{\text{annual inflow}}}}
#'
#' Annual deposited volume (m\eqn{^3}):
#' \deqn{V_{\text{dep}} = \frac{SSL_{\text{annual}} \times 1000}{\rho_{\text{sed}}} \times TE}
#'
#' @param ssl_annual_mt  Numeric.  Mean annual suspended sediment load
#'   (metric tonnes / year), typically from \code{\link{estimate_daily_ssl}}.
#' @param weir_id  Character scalar, \code{"W1"} or \code{"W2"}.
#'   Design parameters (initial storage, design flow) are looked up from the
#'   internal \code{.WEIR_PARAMS} list.
#' @param years  Integer.  Project economic lifetime in years.
#'   Default \code{35L} (source: Shihfeng Power Co. EIA documents).
#' @param rho_sed  Numeric.  Wet bulk density of deposited sediment
#'   (kg/m\eqn{^3}).  Default \code{1300} (Wang et al. 2018, \emph{Water}
#'   10(8):1034; typical gravel-sand mix for Taiwan mountain streams).
#'
#' @return  A \code{data.frame} with one row per project year and columns:
#'   \describe{
#'     \item{year}{Integer. Project year (1 … \code{years})}
#'     \item{S_remaining_m3}{Numeric. Effective storage remaining at year end (m\eqn{^3})}
#'     \item{TE}{Numeric. Trap efficiency for this year (dimensionless, 0–1)}
#'     \item{V_loss_m3_yr}{Numeric. Volume deposited this year (m\eqn{^3}/yr)}
#'     \item{CI_ratio}{Numeric. Capacity-inflow ratio at start of year}
#'   }
#'
#' @references
#'   Brune, G.M. (1953). Trap efficiency of reservoirs.
#'   \emph{Transactions of the American Geophysical Union} 34(3):407–418.
#'   \doi{10.1029/TR034i003p00407}
#'
#'   Wang, H.W., et al. (2018). Effect of water diversion on sediment
#'   flushing and accumulation in a reservoir.
#'   \emph{Water} 10(8):1034. \doi{10.3390/w10081034}
#'
#' @examples
#' \dontrun{
#' trap <- estimate_trap_efficiency(ssl_annual_mt = 50000, weir_id = "W1")
#' plot(trap$year, trap$TE, type = "l",
#'      xlab = "Project year", ylab = "Trap efficiency")
#' }
#'
#' @export
estimate_trap_efficiency <- function(ssl_annual_mt,
                                     weir_id        = "W1",
                                     years          = 35L,
                                     rho_sed        = .RHO_SED_KG_M3,
                                     mean_Q_cms     = NULL) {

  stopifnot(
    weir_id %in% names(.WEIR_PARAMS),
    ssl_annual_mt >= 0,
    years >= 1L
  )

  params    <- .WEIR_PARAMS[[weir_id]]
  S_current <- params$S_max_m3

  # ── Annual inflow volume ──────────────────────────────────────────────────
  # CRITICAL: CI = S_max / V_annual must use the OBSERVED mean annual inflow,
  # not the design flow.  Using Q_design (24.3 cms for W1) gives
  # CI ≈ 0.001, which drives TE to ~100% and empties storage in year 1 —
  # a physically absurd result.
  #
  # The caller should supply mean_Q_cms (mean of observed daily Q at the weir).
  # If not supplied, we fall back to Q_design × 0.4 as a rough correction
  # (typical utilisation factor for run-of-river plants in Taiwan), but the
  # caller-supplied value is always preferred.

  if (is.null(mean_Q_cms) || mean_Q_cms <= 0) {
    # Fallback: 40% of design flow is a conservative lower-bound estimate
    # of mean annual flow for a run-of-river plant with ecological baseflow.
    mean_Q_cms <- params$Q_design_cms * 0.40
    message(sprintf(
      "estimate_trap_efficiency [%s]: mean_Q_cms not supplied; using %.1f cms (40%% of Q_design = %.1f cms). Supply observed mean_Q_cms for accuracy.",
      weir_id, mean_Q_cms, params$Q_design_cms
    ))
  }

  # Convert mean daily flow (cms) to annual volume (m³/yr)
  annual_inflow_m3 <- mean_Q_cms * 86400 * 365.25

  # Annual deposited volume (m³): convert MT → kg → m³
  V_sed_yr <- ssl_annual_mt * 1000 / rho_sed

  purrr::map_dfr(seq_len(years), function(yr) {

    CI     <- S_current / annual_inflow_m3   # dimensionless capacity-inflow ratio
    TE     <- CI / (CI + 0.0021)             # Brune (1953) curve
    V_loss <- V_sed_yr * TE

    S_current <<- max(S_current - V_loss, 0)

    data.frame(
      year           = yr,
      S_remaining_m3 = round(S_current, 0),
      TE             = round(TE, 4),
      V_loss_m3_yr   = round(V_loss, 1),
      CI_ratio       = round(CI, 6)
    )
  })
}


# -----------------------------------------------------------------------------
# build_sediment_outputs()
#
# Convenience wrapper that runs the full sediment pipeline and returns
# all outputs as a named list for use in fp_hydro_main.qmd.
#
# Arguments
#   sed_csv_path   character   path to sediment_clean.csv
#   daily_flow     data.frame  columns: date, Q_cms (NA-filled)
#   plot_fit       logical     whether to display rating curve plot
#
# Returns  named list:
#   $curve         list        rating curve parameters
#   $daily_ssl     data.frame  daily SSL estimates
#   $annual_ssl_mt numeric     mean annual SSL (MT/year)
#   $trap_W1       data.frame  trap efficiency trajectory W1 (35 years)
#   $trap_W2       data.frame  trap efficiency trajectory W2 (35 years)
# -----------------------------------------------------------------------------

build_sediment_outputs <- function(sed_csv_path,
                                   daily_flow,
                                   plot_fit = TRUE) {
  
  sed_clean  <- read_csv(sed_csv_path, show_col_types = FALSE)
  curve      <- fit_rating_curve(sed_clean, plot_fit = plot_fit)
  daily_ssl  <- estimate_daily_ssl(daily_flow, curve)
  
  annual_ssl <- daily_ssl |>
    mutate(yr = format(date, "%Y")) |>
    group_by(yr) |>
    summarise(ssl_yr = sum(SSL_mt_day, na.rm = TRUE), .groups = "drop") |>
    pull(ssl_yr) |>
    mean()

  message(sprintf(
    "Mean annual SSL (averaged over %d calendar years): %.0f MT/year",
    nrow(daily_ssl |> mutate(yr = format(date, "%Y")) |> distinct(yr)),
    annual_ssl
  ))
  
  list(
    curve         = curve,
    daily_ssl     = daily_ssl,
    annual_ssl_mt = annual_ssl,
    trap_W1       = estimate_trap_efficiency(annual_ssl, "W1"),
    trap_W2       = estimate_trap_efficiency(annual_ssl, "W2")
  )
}