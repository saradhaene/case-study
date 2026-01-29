# ============================================================
# JOINT WEEKLY GAM — Exceedance across multiple stations (1970–2020)
# - Reads ECA&D-style daily RR from multiple stations
# - Defines station-specific exceedance using wet-day percentile threshold
# - Aggregates to weekly (ISO week) indicators per station
# - Defines joint exceedance: at least K stations in a week
# - Fits GAM: joint ~ s(time) + s(woy, cyclic)
# - Optional lag-1 model for persistence
# - Produces: summaries, gam.check PNG, time-trend plot, season plot
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(ggplot2)

DATA_DIR <- "C:/Users/Sara/Downloads/Case study"
setwd(DATA_DIR)

# ----------------------------
# 1) Read + clean one station (daily)
# ----------------------------
read_ecad_rr <- function(path, station_name,
                         date_min = as.Date("1970-01-01"),
                         date_max = as.Date("2020-12-31")) {
  read_excel(path, col_names = FALSE) %>%
    rename(raw = 1) %>%
    separate(raw, into = c("STAID", "SOUID", "DATE", "RR", "Q_RR"),
             sep = ",", convert = TRUE) %>%
    mutate(
      station = station_name,
      DATE    = as.Date(as.character(DATE), format = "%Y%m%d"),
      RR      = suppressWarnings(as.numeric(RR)),
      Q_RR    = suppressWarnings(as.numeric(Q_RR)),
      rr_mm   = RR / 10
    ) %>%
    filter(!is.na(DATE), RR != -9999, Q_RR == 0) %>%
    filter(DATE >= date_min, DATE <= date_max) %>%
    arrange(DATE) %>%
    select(DATE, station, rr_mm)
}

# ----------------------------
# 2) Daily exceedance per station
# ----------------------------
make_daily_exceedance <- function(df_station, wet_day_mm = 1, p_extreme = 0.90) {
  thr <- quantile(df_station$rr_mm[df_station$rr_mm > wet_day_mm],
                  probs = p_extreme, type = 8, na.rm = TRUE)

  df_exc <- df_station %>%
    mutate(exc = as.integer(rr_mm > thr)) %>%
    select(DATE, station, exc)

  list(df = df_exc, thr = as.numeric(thr))
}

# ----------------------------
# 3) Weekly aggregation: any exceedance per station
# ----------------------------
to_weekly_any <- function(df_daily_exc) {
  df_daily_exc %>%
    mutate(
      iso_year = isoyear(DATE),
      iso_week = isoweek(DATE)
    ) %>%
    group_by(iso_year, iso_week, station) %>%
    summarise(week_exc = as.integer(any(exc == 1)), .groups = "drop")
}

# ----------------------------
# 4) Build weekly panel (all stations aligned)
# ----------------------------
build_weekly_panel <- function(stations, wet_day_mm = 1, p_extreme = 0.90) {
  thr_tbl <- tibble(station = names(stations), thr_mm = NA_real_)
  weekly_list <- list()

  for (st in names(stations)) {
    df_raw <- read_ecad_rr(stations[[st]], st)
    tmp <- make_daily_exceedance(df_raw, wet_day_mm = wet_day_mm, p_extreme = p_extreme)
    thr_tbl$thr_mm[thr_tbl$station == st] <- tmp$thr

    w <- to_weekly_any(tmp$df) %>%
      filter(station == st) %>%
      select(iso_year, iso_week, week_exc) %>%
      rename(!!paste0("week_exc_", tolower(st)) := week_exc)

    weekly_list[[st]] <- w
  }

  df_week <- Reduce(
    function(x, y) inner_join(x, y, by = c("iso_year", "iso_week")),
    weekly_list
  ) %>%
    arrange(iso_year, iso_week) %>%
    mutate(
      week_start = as.Date(paste0(iso_year, "-", sprintf("%02d", iso_week), "-1"),
                           format = "%G-%V-%u"),
      time = row_number(),
      woy = as.integer(iso_week)
    )

  list(df_week = df_week, thresholds = thr_tbl)
}

# ----------------------------
# 5) Joint outcome: at least K stations exceed in a week
# ----------------------------
add_joint_outcome <- function(df_week, k_joint = 2) {
  exc_cols <- grep("^week_exc_", names(df_week), value = TRUE)

  df_week %>%
    mutate(
      n_exc = rowSums(across(all_of(exc_cols))),
      exc_joint = as.integer(n_exc >= k_joint),
      lag1 = dplyr::lag(exc_joint)
    ) %>%
    filter(!is.na(lag1))
}

# ----------------------------
# 6) Fit GAMs
# ----------------------------
fit_joint_gam <- function(data, k_time = 30, k_woy = 20) {
  gam(
    exc_joint ~ s(time, k = k_time) + s(woy, bs = "cc", k = k_woy),
    data = data,
    family = binomial(link = "logit"),
    method = "REML",
    knots = list(woy = c(0.5, 53.5))
  )
}

fit_joint_gam_lag <- function(data, k_time = 30, k_woy = 20) {
  gam(
    exc_joint ~ s(time, k = k_time) + s(woy, bs = "cc", k = k_woy) + lag1,
    data = data,
    family = binomial(link = "logit"),
    method = "REML",
    knots = list(woy = c(0.5, 53.5))
  )
}

# ----------------------------
# 7) Plot helpers (time trend + seasonality)
# ----------------------------
plot_time_trend <- function(data, model, title, woy_fixed = 26) {
  nd <- data.frame(
    time = seq(min(data$time), max(data$time), length.out = 500),
    woy  = woy_fixed,
    lag1 = 0
  )

  pred <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  nd$p  <- plogis(pred$fit)
  nd$lo <- plogis(pred$fit - 1.96 * pred$se.fit)
  nd$hi <- plogis(pred$fit + 1.96 * pred$se.fit)
  nd$week_start <- data$week_start[1] + (nd$time - 1) * 7

  ggplot(nd, aes(x = week_start, y = p)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = "P(joint exceedance)", title = title) +
    theme_minimal()
}

plot_seasonality <- function(data, model, title, time_fixed = NULL) {
  if (is.null(time_fixed)) {
    time_fixed <- median(data$time)
  }

  nd <- data.frame(
    time = time_fixed,
    woy = 1:53,
    lag1 = 0
  )

  pred <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  nd$p  <- plogis(pred$fit)
  nd$lo <- plogis(pred$fit - 1.96 * pred$se.fit)
  nd$hi <- plogis(pred$fit + 1.96 * pred$se.fit)

  ggplot(nd, aes(x = woy, y = p)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = seq(1, 53, by = 4)) +
    labs(x = "ISO week", y = "P(joint exceedance)", title = title) +
    theme_minimal()
}

# ============================================================
# RUN
# ============================================================

stations <- c(
  Madrid    = "madrid.xlsx",
  Barcelona = "barcelona.xlsx",
  Asturias  = "asturias.xlsx",
  Sevilla   = "sevilla.xlsx",
  Vigo      = "vigo.xlsx",
  Valencia  = "valencia.xlsx",
  Alicante  = "alicante.xlsx"
)

wet_day_mm <- 1
p_extreme  <- 0.90
K_joint    <- 2

panel <- build_weekly_panel(stations, wet_day_mm = wet_day_mm, p_extreme = p_extreme)
print(panel$thresholds)

df_week <- add_joint_outcome(panel$df_week, k_joint = K_joint)

cat("Weeks:", nrow(df_week), "\n")
cat("Joint exceedance rate:", round(mean(df_week$exc_joint), 4), "\n")

# Base model
m_joint <- fit_joint_gam(df_week)
print(summary(m_joint))

png(paste0("joint_gamcheck_K", K_joint, ".png"), width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_joint)
dev.off()

# Optional lag model
m_joint_lag <- fit_joint_gam_lag(df_week)
print(summary(m_joint_lag))

png(paste0("joint_gamcheck_lag_K", K_joint, ".png"), width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_joint_lag)
dev.off()

# Plot time trend + seasonality (base model)
p_time <- plot_time_trend(
  df_week,
  m_joint,
  title = paste0("Joint weekly exceedance (K=", K_joint, ") — time trend")
)
print(p_time)

p_season <- plot_seasonality(
  df_week,
  m_joint,
  title = paste0("Joint weekly exceedance (K=", K_joint, ") — seasonality")
)
print(p_season)
