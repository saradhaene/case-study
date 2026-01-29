# ============================================================
# ALICANTE ONLY — Daily exceedance GAM (1970–2020)
# - Reads ECA&D-style daily RR from alicante.xlsx
# - Keeps only valid observations (Q_RR==0) in 1970-01-01..2020-12-31
# - Defines exceedance using wet-day percentile threshold (recommended)
# - Fits GAM: exc ~ s(time) + s(doy, cyclic)
# - Produces: summaries, gam.check PNG, time-trend plot, season plot
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(ggplot2)

setwd("C:/Users/Sara/Downloads/Case study")

# ----------------------------
# 1) Read + clean Alicante
# ----------------------------
read_ecad_rr <- function(path, station_name = "asturias",
                         date_min = as.Date("1970-01-01"),
                         date_max = as.Date("2020-12-31")) {
  
  read_excel(path, col_names = FALSE) %>%
    rename(raw = 1) %>%
    separate(raw, into = c("STAID","SOUID","DATE","RR","Q_RR"),
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

df <- read_ecad_rr("asturias.xlsx", "Asturias")

cat("Rows:", nrow(df), "\n")
cat("Date range:", format(min(df$DATE)), "to", format(max(df$DATE)), "\n")
print(summary(df$rr_mm))

# ----------------------------
# 2) Define exceedance (recommended: wet-day percentile)
# ----------------------------
wet_day_mm <- 1
p_extreme  <- 0.95   # set to 0.90 if you want more events/power

thr <- quantile(df$rr_mm[df$rr_mm > wet_day_mm],
                probs = p_extreme, type = 8, na.rm = TRUE)

cat("\nThreshold (mm) for exceedance =", round(thr, 2),
    " (wet_day_mm >", wet_day_mm, ", p =", p_extreme, ")\n")

df <- df %>%
  mutate(
    exc  = as.integer(rr_mm > thr),
    time = as.numeric(DATE - min(DATE)),
    doy  = yday(DATE)
  )

cat("\nExceedance rate:", round(mean(df$exc), 4), "\n")
print(table(df$exc))

# ----------------------------
# 3) Fit GAM (binomial-logit; smooth time + cyclic season)
# ----------------------------
m_indiv <- gam(
  exc ~ s(time, k = 30) + s(doy, bs = "cc", k = 20),
  data   = df,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

print(summary(m_indiv))

png("individual_gamcheck.png", width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_indiv)
dev.off()

# ----------------------------
# 4) Plot: time trend (season fixed)
# ----------------------------
nd_time <- data.frame(
  time = seq(min(df$time), max(df$time), length.out = 800),
  doy  = 180
)

pred_t <- predict(m_indiv, newdata = nd_time, type = "link", se.fit = TRUE)
nd_time$p  <- plogis(pred_t$fit)
nd_time$lo <- plogis(pred_t$fit - 1.96 * pred_t$se.fit)
nd_time$hi <- plogis(pred_t$fit + 1.96 * pred_t$se.fit)
nd_time$date <- min(df$DATE) + round(nd_time$time)

ggplot(nd_time, aes(x = date, y = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "P(exceedance)",
    title = paste0("Alicante: exceedance probability over time (doy fixed at 180)\n",
                   "threshold=", round(thr, 1), " mm, p=", p_extreme)
  ) +
  theme_minimal() -> p_time

print(p_time)

# ----------------------------
# 5) Plot: seasonality (time fixed)
# ----------------------------
nd_season <- data.frame(
  doy  = 1:366,
  time = median(df$time)
)

pred_s <- predict(m_indiv, newdata = nd_season, type = "link", se.fit = TRUE)
nd_season$p  <- plogis(pred_s$fit)
nd_season$lo <- plogis(pred_s$fit - 1.96 * pred_s$se.fit)
nd_season$hi <- plogis(pred_s$fit + 1.96 * pred_s$se.fit)

ggplot(nd_season, aes(x = doy, y = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 366, by = 30)) +
  labs(
    x = "Day of year",
    y = "P(exceedance)",
    title = paste0("Seasonal exceedance probability (time fixed)\n",
                   "threshold=", round(thr, 1), " mm, p=", p_extreme)
  ) +
  theme_minimal() -> p_season

print(p_season)

# ----------------------------
# OPTIONAL: compare against a fixed absolute threshold (10 mm)
# (Keep only if you really want it; otherwise delete this block)
# ----------------------------
# THRESHOLD_MM <- 10
# df_abs <- df %>%
#   mutate(exc_abs = as.integer(rr_mm > THRESHOLD_MM))
#
# m_alic_abs <- gam(
#   exc_abs ~ s(time, k = 30) + s(doy, bs = "cc", k = 20),
#   data   = df_abs,
#   family = binomial(link = "logit"),
#   method = "REML",
#   knots  = list(doy = c(0.5, 366.5))
# )
# summary(m_alic_abs)


# lag
# ============================================================
# INDIVIDUAL STATION (daily) — GAM with 1-day persistence (lag)
# - Reads ECA&D-style daily RR from <station>.xlsx
# - Restricts to 1970-01-01..2020-12-31 (Q_RR==0)
# - Defines exceedance: wet-day percentile threshold
# - Fits GAM: exc ~ s(time) + s(doy, cyclic) + lag1
# - Produces: summary, gam.check PNG, time plot, season plot, lag odds ratio
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(ggplot2)

setwd("C:/Users/Sara/Downloads/Case study")

# ----------------------------
# 1) Read + clean one station
# ----------------------------
read_ecad_rr <- function(path, station_name,
                         date_min = as.Date("1970-01-01"),
                         date_max = as.Date("2020-12-31")) {
  
  read_excel(path, col_names = FALSE) %>%
    rename(raw = 1) %>%
    separate(raw, into = c("STAID","SOUID","DATE","RR","Q_RR"),
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
# 2) Choose station here
# ----------------------------
station_name <- "Asturias"
station_file <- "asturias.xlsx"

df <- read_ecad_rr(station_file, station_name)

cat("Rows:", nrow(df), "\n")
cat("Date range:", format(min(df$DATE)), "to", format(max(df$DATE)), "\n")

# ----------------------------
# 3) Define exceedance (wet-day percentile)
# ----------------------------
wet_day_mm <- 1
p_extreme  <- 0.95   # try 0.90 if exceedances are too rare

thr <- quantile(df$rr_mm[df$rr_mm > wet_day_mm],
                probs = p_extreme, type = 8, na.rm = TRUE)

cat("Threshold (mm):", round(thr, 2), "\n")

df <- df %>%
  mutate(
    exc  = as.integer(rr_mm > thr),
    time = as.numeric(DATE - min(DATE)),
    doy  = yday(DATE)
  )

cat("Exceedance rate:", round(mean(df$exc), 4), "\n")
print(table(df$exc))

# ----------------------------
# 4) Add 1-day lag (persistence)
# ----------------------------
df_lag <- df %>%
  arrange(DATE) %>%
  mutate(lag1 = dplyr::lag(exc)) %>%
  filter(!is.na(lag1))

# ----------------------------
# 5) Fit GAM with lag
# ----------------------------
m_lag <- gam(
  exc ~ s(time, k = 30) + s(doy, bs = "cc", k = 20) + lag1,
  data   = df_lag,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

print(summary(m_lag))

# Odds ratio for lag effect
b_lag <- coef(m_lag)["lag1"]
cat("\nLag1 odds ratio (exp(beta)) =", round(exp(b_lag), 3), "\n")

png(paste0("gamcheck_", tolower(station_name), "_lag1.png"),
    width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_lag)
dev.off()

# ----------------------------
# 6) Plot time trend (season fixed, lag fixed)
#     Fix lag1=0 to interpret "baseline" probability
# ----------------------------
nd_time <- data.frame(
  time = seq(min(df_lag$time), max(df_lag$time), length.out = 800),
  doy  = 180,
  lag1 = 0
)

pred_t <- predict(m_lag, newdata = nd_time, type = "link", se.fit = TRUE)
nd_time$p  <- plogis(pred_t$fit)
nd_time$lo <- plogis(pred_t$fit - 1.96 * pred_t$se.fit)
nd_time$hi <- plogis(pred_t$fit + 1.96 * pred_t$se.fit)
nd_time$date <- min(df_lag$DATE) + round(nd_time$time)

p_time <- ggplot(nd_time, aes(x = date, y = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "P(exceedance)",
    title = paste0(station_name, ": time trend (season fixed at doy=180, lag1=0)\n",
                   "thr=", round(thr, 1), " mm (p=", p_extreme, ")")
  ) +
  theme_minimal()

print(p_time)

# ----------------------------
# 7) Plot seasonality (time fixed, lag fixed)
# ----------------------------
nd_season <- data.frame(
  doy  = 1:366,
  time = median(df_lag$time),
  lag1 = 0
)

pred_s <- predict(m_lag, newdata = nd_season, type = "link", se.fit = TRUE)
nd_season$p  <- plogis(pred_s$fit)
nd_season$lo <- plogis(pred_s$fit - 1.96 * pred_s$se.fit)
nd_season$hi <- plogis(pred_s$fit + 1.96 * pred_s$se.fit)

p_season <- ggplot(nd_season, aes(x = doy, y = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 366, by = 30)) +
  labs(
    x = "Day of year",
    y = "P(exceedance)",
    title = paste0(station_name, ": seasonal pattern (time fixed, lag1=0)\n",
                   "thr=", round(thr, 1), " mm (p=", p_extreme, ")")
  ) +
  theme_minimal()

print(p_season)


# nao
station_name <- "sevilla"
station_file <- "sevilla.xlsx"

wet_day_mm <- 1
p_extreme  <- 0.95

df <- read_ecad_rr(
  path = station_file,
  station_name = station_name
)

cat("Station:", station_name, "\n")
cat("Rows:", nrow(df), "\n")
cat("Date range:", min(df$DATE), "to", max(df$DATE), "\n")

thr <- quantile(
  df$rr_mm[df$rr_mm > wet_day_mm],
  probs = p_extreme,
  type = 8,
  na.rm = TRUE
)

df <- df %>%
  mutate(
    exc  = as.integer(rr_mm > thr),
    time = as.numeric(DATE - min(DATE)),
    doy  = yday(DATE)
  )

cat("Threshold:", round(thr, 2), "mm\n")
cat("Exceedance rate:", round(mean(df$exc), 4), "\n")

df_lag <- df %>%
  arrange(DATE) %>%
  mutate(lag1 = dplyr::lag(exc)) %>%
  filter(!is.na(lag1))


# Read NAO file (daily)
nao_daily <- read_excel("NAO_cleaned.xlsx") %>%
  mutate(
    DATE = as.Date(DATE),
    nao  = nao_index_cdas
  ) %>%
  select(DATE, nao) %>%
  filter(!is.na(nao))

# Join NAO onto your lagged station dataframe (df_lag)
df_lag_nao <- df_lag %>%
  left_join(nao_daily, by = "DATE") %>%
  filter(!is.na(nao))

cat("Rows after NAO join:", nrow(df_lag_nao), "\n")

# lineaire nao
m_lag_nao <- gam(
  exc ~ s(time, k = 30) +
    s(doy, bs = "cc", k = 20) +
    lag1 +
    nao,
  data   = df_lag_nao,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

summary(m_lag_nao)
gam.check(m_lag_nao)

cat("Lag OR:", round(exp(coef(m_lag_nao)["lag1"]), 2), "\n")
cat("NAO OR:", round(exp(coef(m_lag_nao)["nao"]), 2), "\n")






# nao x seizoen
m_test <- gam(
  exc ~ s(time) + s(doy, bs="cc") + lag1 + nao +
    ti(nao, doy, bs=c("tp","cc")),
  family = binomial,
  data = df_lag_nao,
  method = "REML"
)
summary(m_test)


# niet lineaire nao
# deze niet doen, mogelijk uitleggen waarom
m_lag_nao_s <- gam(
  exc ~ s(time, k = 30) +
    s(doy, bs = "cc", k = 20) +
    lag1 +
    s(nao, k = 10),
  data   = df_lag_nao,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

summary(m_lag_nao_s)
gam.check(m_lag_nao_s)


# validation
library(mgcv)
library(dplyr)
library(pROC)
library(ggplot2)

eval_station_split_pred <- function(df, formula, train_end_date) {
  
  df <- df %>% arrange(DATE)
  
  train <- df %>% filter(DATE <= as.Date(train_end_date))
  test  <- df %>% filter(DATE >  as.Date(train_end_date))
  
  mod <- gam(
    formula,
    data = train,
    family = binomial(link="logit"),
    method = "REML",
    knots = list(doy = c(0.5, 366.5))
  )
  
  p <- predict(mod, newdata = test, type="response")
  y <- test$exc
  
  pi_hat <- mean(train$exc)
  brier_clim <- pi_hat * (1 - pi_hat)
  
  brier <- mean((y - p)^2)
  logloss <- -mean(y*log(p + 1e-15) + (1-y)*log(1-p + 1e-15))
  auc <- tryCatch(as.numeric(pROC::auc(y, p)), error = function(e) NA_real_)
  
  list(
    model = mod,
    scores = tibble(
      train_end = as.Date(train_end_date),
      n_train = nrow(train), n_test = nrow(test),
      event_rate_train = pi_hat,
      brier = brier,
      brier_climatology = brier_clim,
      logloss = logloss,
      auc = auc
    ),
    y = y,
    p = p
  )
}

out_base <- eval_station_split_pred(
  df,
  exc ~ s(time, k=30) + s(doy, bs="cc", k=20),
  "2010-12-31"
)
out_base$scores

out_lag <- eval_station_split_pred(
  df_lag,
  exc ~ s(time, k=30) + s(doy, bs="cc", k=20) + lag1,
  "2010-12-31"
)
out_lag$scores

calibration_plot_bins(out_lag$y, out_lag$p, n_bins = 10)
calibration_gam_smooth(out_lag$y, out_lag$p)



# recalibreren
# Zorg dat df gesorteerd is
df <- df %>% arrange(DATE)

train <- df %>% filter(DATE <= as.Date("2010-12-31"))
test  <- df %>% filter(DATE >  as.Date("2010-12-31"))

m_train <- gam(
  exc ~ s(time, k=30) + s(doy, bs="cc", k=20),
  data = train,
  family = binomial(link="logit"),
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

p_train <- predict(m_train, newdata = train, type = "response")
p_test  <- predict(m_train, newdata = test,  type = "response")

y_train <- train$exc
y_test  <- test$exc

p_test_cal <- recalibrate_platt(
  y_train = y_train,
  p_train = p_train,
  p_test  = p_test
)

brier_raw <- mean((y_test - p_test)^2)
logloss_raw <- -mean(y_test*log(p_test + 1e-15) +
                       (1-y_test)*log(1-p_test + 1e-15))
brier_cal <- mean((y_test - p_test_cal)^2)
logloss_cal <- -mean(y_test*log(p_test_cal + 1e-15) +
                       (1-y_test)*log(1-p_test_cal + 1e-15))
tibble(
  model = c("raw", "recalibrated"),
  brier = c(brier_raw, brier_cal),
  logloss = c(logloss_raw, logloss_cal)
)
calibration_plot_bins(y_test, p_test, n_bins = 10)
calibration_gam_smooth(y_test, p_test)
calibration_plot_bins(y_test, p_test_cal, n_bins = 10)
calibration_gam_smooth(y_test, p_test_cal)

