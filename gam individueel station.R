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

station_name <- "Alicante"

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
if (!is.na(script_dir) && nzchar(script_dir)) {
  setwd(script_dir)
}

# read individual station + filter on date
read_ecad_rr <- function(path, station_name = "alicante",
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

df <- read_ecad_rr("alicante.xlsx", station_name)

cat("Rows:", nrow(df), "\n")
cat("Date range:", format(min(df$DATE)), "to", format(max(df$DATE)), "\n")
print(summary(df$rr_mm))

# Define exceedance threshold and wet day
wet_day_mm <- 1
p_extreme  <- 0.95   

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

# fit GAM, k can be altered
m_indiv <- gam(
  exc ~ s(time, k = 30) + s(doy, bs = "cc", k = 20),
  data   = df,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

print(summary(m_indiv))

# gam check to png in folder
png("individual_gamcheck.png", width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_indiv)
dev.off()

# plot time trend using baseline GAM keeping seasonality fixed
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
    title = paste0(station_name, " — baseline GAM: time trend (doy fixed at 180)\n",
                   "threshold=", round(thr, 1), " mm, p=", p_extreme)
  ) +
  theme_minimal() -> p_time

print(p_time)

# plot seasonality holding time fixed (baseline GAM)
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
    title = paste0(station_name, " — baseline GAM: seasonal pattern (time fixed)\n",
                   "threshold=", round(thr, 1), " mm, p=", p_extreme)
  ) +
  theme_minimal() -> p_season

print(p_season)



# gam with one lag persistence
df_lag <- df %>%
  arrange(DATE) %>%
  mutate(lag1 = dplyr::lag(exc)) %>%
  filter(!is.na(lag1))

# fit gam with lag
m_lag <- gam(
  exc ~ s(time, k = 30) + s(doy, bs = "cc", k = 20) + lag1,
  data   = df_lag,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = list(doy = c(0.5, 366.5))
)

print(summary(m_lag))

# odds ratio for lag effect
b_lag <- coef(m_lag)["lag1"]
cat("\nLag1 odds ratio (exp(beta)) =", round(exp(b_lag), 3), "\n")

png(paste0("gamcheck_", tolower(station_name), "_lag1.png"),
    width = 1400, height = 1000, res = 150)
par(mar = c(4, 4, 2, 1))
gam.check(m_lag)
dev.off()

# plot time using GAM with lag
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
    title = paste0(station_name, " — GAM + lag1: time trend (doy=180, lag1=0)\n",
                   "thr=", round(thr, 1), " mm (p=", p_extreme, ")")
  ) +
  theme_minimal()

print(p_time)

# plot seasonality using GAM with lag
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
    title = paste0(station_name, " — GAM + lag1: seasonal pattern (time fixed, lag1=0)\n",
                   "thr=", round(thr, 1), " mm (p=", p_extreme, ")")
  ) +
  theme_minimal()

print(p_season)



# read NAO file
nao_daily <- read_excel("NAO_cleaned.xlsx") %>%
  mutate(
    DATE = as.Date(DATE),
    nao  = nao_index_cdas
  ) %>%
  select(DATE, nao) %>%
  filter(!is.na(nao))

# join nao to model with lag
df_lag_nao <- df_lag %>%
  left_join(nao_daily, by = "DATE") %>%
  filter(!is.na(nao))

cat("Rows after NAO join:", nrow(df_lag_nao), "\n")

# fit gam with lag and nao
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

# plot time using GAM with lag + NAO (NAO fixed to 0)
nd_time <- data.frame(
  time = seq(min(df_lag_nao$time), max(df_lag_nao$time), length.out = 800),
  doy  = 180,
  lag1 = 0,
  nao  = 0
)

pred_t <- predict(m_lag_nao, newdata = nd_time, type = "link", se.fit = TRUE)
nd_time$p  <- plogis(pred_t$fit)
nd_time$lo <- plogis(pred_t$fit - 1.96 * pred_t$se.fit)
nd_time$hi <- plogis(pred_t$fit + 1.96 * pred_t$se.fit)
nd_time$date <- min(df_lag_nao$DATE) + round(nd_time$time)

p_time <- ggplot(nd_time, aes(x = date, y = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "P(exceedance)",
    title = paste0(station_name,
                   " — GAM + lag1 + NAO: time trend (doy=180, lag1=0, nao=0)\n",
                   "thr=", round(thr, 1), " mm (p=", p_extreme, ")")
  ) +
  theme_minimal()

print(p_time)

# plot seasonality using GAM with lag + NAO (NAO fixed to 0)
nd_season <- data.frame(
  doy  = 1:366,
  time = median(df_lag_nao$time),
  lag1 = 0,
  nao  = 0
)

pred_s <- predict(m_lag_nao, newdata = nd_season, type = "link", se.fit = TRUE)
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
    title = paste0(
      station_name,
      " — GAM + lag1 + NAO: seasonal pattern (time fixed, lag1=0, nao=0)\n",
      "thr=", round(thr, 1), " mm (p=", p_extreme, ")"
    )
  ) +
  theme_minimal()

print(p_season)



# nao x seasonality
m_test <- gam(
  exc ~ s(time) + s(doy, bs="cc") + lag1 + nao +
    ti(nao, doy, bs=c("tp","cc")),
  family = binomial,
  data = df_lag_nao,
  method = "REML"
)
summary(m_test)
gam.check(m_test)

