############################################################
# CLEAN JOINT GAM SCRIPT (weekly joint exceedances)
# - Read & clean ECAD station files (Excel with raw CSV rows)
# - Create station-specific exceedance (p=0.95 wet-day > 1mm)
# - Aggregate to ISO-week "any exceedance"
# - Build weekly panel across stations
# - Construct joint exceedance outcomes (>=2, >=3)
# - Fit baseline GAM + k-sensitivity table
# - Serial correlation diagnostics (ACF + Ljung-Box)
############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(ggplot2)
library(purrr)

setwd("C:/Users/Sara/Downloads/Case study")

# ============================================================
# 1) Read + clean one station (daily)
# ============================================================
read_ecad_rr <- function(path, station_name) {
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
    filter(DATE >= as.Date("1970-01-01"), DATE <= as.Date("2020-12-31")) %>%
    arrange(DATE) %>%
    select(DATE, station, rr_mm)
}

# ============================================================
# 2) Daily exceedance indicator (station-specific threshold)
# ============================================================
make_daily_exceedance <- function(df_station, wet_day_mm = 1, p_extreme = 0.95) {
  thr <- quantile(df_station$rr_mm[df_station$rr_mm > wet_day_mm],
                  probs = p_extreme, type = 8, na.rm = TRUE)
  out <- df_station %>%
    mutate(exc = as.integer(rr_mm > thr)) %>%
    select(DATE, exc)
  
  list(df = out, thr = as.numeric(thr))
}

# ============================================================
# 3) Weekly aggregation: any exceedance in ISO week
# ============================================================
to_weekly_any <- function(df_daily_exc) {
  df_daily_exc %>%
    mutate(
      iso_year = isoyear(DATE),
      iso_week = isoweek(DATE)
    ) %>%
    group_by(iso_year, iso_week) %>%
    summarise(week_exc = as.integer(any(exc == 1)), .groups = "drop") %>%
    arrange(iso_year, iso_week)
}

# ============================================================
# 4) Build aligned weekly panel (inner join on iso_year/week)
# ============================================================
build_weekly_panel <- function(stations, wet_day_mm = 1, p_extreme = 0.95) {
  
  thr_tbl <- tibble(station = names(stations), thr_mm = NA_real_)
  weekly_list <- list()
  
  for (st in names(stations)) {
    df_raw <- read_ecad_rr(stations[[st]], st)
    tmp <- make_daily_exceedance(df_raw, wet_day_mm = wet_day_mm, p_extreme = p_extreme)
    thr_tbl$thr_mm[thr_tbl$station == st] <- tmp$thr
    
    w <- to_weekly_any(tmp$df) %>%
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
      t   = row_number(),        # time index
      woy = as.integer(iso_week) # week of year (1..53)
    )
  
  list(df_week = df_week, thresholds = thr_tbl)
}

# ============================================================
# 5) Joint outcomes
# ============================================================
add_joint_outcomes <- function(df_week) {
  exc_cols <- grep("^week_exc_", names(df_week), value = TRUE)
  
  df_week %>%
    mutate(
      n_exc = rowSums(across(all_of(exc_cols))),
      joint_2p_w = as.integer(n_exc >= 2),
      joint_3p_w = as.integer(n_exc >= 3)
    )
}

# ============================================================
# 6) Fit weekly GAM
# ============================================================
fit_weekly_gam <- function(data, yvar, k_time = 30, k_woy = 20) {
  fml <- as.formula(paste0(
    yvar, " ~ s(t, k = ", k_time, ") + s(woy, bs = 'cc', k = ", k_woy, ")"
  ))
  
  gam(
    fml,
    family = binomial(link = "logit"),
    method = "REML",
    data = data,
    knots = list(woy = c(0.5, 53.5))
  )
}

# ============================================================
# 7) Diagnostics helpers
# ============================================================
serial_diag <- function(mod, lags = c(10, 20)) {
  rP <- residuals(mod, type = "pearson")
  rD <- residuals(mod, type = "deviance")
  
  cat("\n--- Ljung-Box tests (Pearson residuals) ---\n")
  for (L in lags) print(Box.test(rP, lag = L, type = "Ljung-Box"))
  
  cat("\n--- Ljung-Box tests (Deviance residuals) ---\n")
  for (L in lags) print(Box.test(rD, lag = L, type = "Ljung-Box"))
  
  # ACF plots
  acf(rP, na.action = na.pass, main = "ACF Pearson residuals")
  acf(rD, na.action = na.pass, main = "ACF Deviance residuals")
}

# Extract k-check table from gam.check
get_kcheck_table <- function(mod) {
  gc <- gam.check(mod, rep = 0)  # rep=0 avoids simulation checks
  ktab <- gc$k.check
  
  tibble(
    term    = rownames(ktab),
    k_prime = ktab[, "k'"],
    edf     = ktab[, "edf"],
    k_index = ktab[, "k-index"],
    p_value = ktab[, "p-value"]
  )
}

# ============================================================
# FIXED: k-check extraction + k-sensitivity table
# ============================================================

get_kcheck_table <- function(mod) {
  # mgcv::k.check returns a matrix with rows = smooths
  ktab <- mgcv::k.check(mod)
  
  tibble::tibble(
    term    = rownames(ktab),
    k_prime = ktab[, "k'"],
    edf     = ktab[, "edf"],
    k_index = ktab[, "k-index"],
    p_value = ktab[, "p-value"]
  )
}

k_sensitivity <- function(data, yvar, grid) {
  out_rows <- list()
  mods <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    kt <- grid$k_time[i]
    kw <- grid$k_woy[i]
    
    m <- fit_weekly_gam(data, yvar, k_time = kt, k_woy = kw)
    mods[[i]] <- m
    
    kc_wide <- get_kcheck_table(m) %>%
      dplyr::select(term, k_index, p_value) %>%
      tidyr::pivot_wider(
        names_from = term,
        values_from = c(k_index, p_value),
        names_glue = "{.value}_{term}"
      )
    
    out_rows[[i]] <- tibble::tibble(
      yvar = yvar,
      spec = paste0("k_time=", kt, ", k_woy=", kw),
      AIC  = AIC(m),
      REML = m$gcv.ubre
    ) %>% dplyr::bind_cols(kc_wide)
  }
  
  list(models = mods, table = dplyr::bind_rows(out_rows))
}


# ============================================================
# RUN SECTION
# ============================================================

# ---- (A) Define station files (edit filenames if needed)
stations <- c(
  Madrid    = "madrid.xlsx",
  Barcelona = "barcelona.xlsx",
  Asturias  = "asturias.xlsx",
  Sevilla   = "sevilla.xlsx",
  Vigo      = "vigo.xlsx",
  Alicante  = "alicante.xlsx",
  Valencia  = "valencia.xlsx"
)

# ---- (B) Build weekly panel
panel <- build_weekly_panel(stations, wet_day_mm = 1, p_extreme = 0.95)
df_week <- add_joint_outcomes(panel$df_week)

cat("\n--- Station thresholds (mm) ---\n")
print(panel$thresholds)

cat("\n--- Event counts ---\n")
print(df_week %>% summarise(
  n_weeks = n(),
  events_2p = sum(joint_2p_w),
  events_3p = sum(joint_3p_w),
  rate_2p = mean(joint_2p_w),
  rate_3p = mean(joint_3p_w)
))

# ---- (C) Baseline models
m_2p <- fit_weekly_gam(df_week, "joint_2p_w", k_time = 30, k_woy = 20)
m_3p <- fit_weekly_gam(df_week, "joint_3p_w", k_time = 30, k_woy = 20)

cat("\n--- Baseline summary: joint_2p_w ---\n")
print(summary(m_2p))
cat("\n--- k-check: joint_2p_w ---\n")
print(get_kcheck_table(m_2p))

cat("\n--- Baseline summary: joint_3p_w ---\n")
print(summary(m_3p))
cat("\n--- k-check: joint_3p_w ---\n")
print(get_kcheck_table(m_3p))

# ---- (D) Serial correlation diagnostics (baseline)
cat("\n================ Serial correlation diagnostics: m_2p ================\n")
serial_diag(m_2p)

# ---- (E) Optional lag robustness (only if you want)
df_week_lag <- df_week %>%
  arrange(week_start) %>%
  mutate(lag1 = dplyr::lag(joint_2p_w)) %>%
  filter(!is.na(lag1))

m_2p_lag <- gam(
  joint_2p_w ~ s(t, k = 60) + s(woy, bs = "cc", k = 20) + lag1,
  family = binomial(link = "logit"),
  method = "REML",
  data = df_week_lag,
  knots = list(woy = c(0.5, 53.5))
)

cat("\n--- Lag model summary (joint_2p_w) ---\n")
print(summary(m_2p_lag))

# ---- (F) k-sensitivity grid (recommended)
k_grid <- tibble(
  k_time = c(30, 60, 80),
  k_woy  = c(20, 20, 30)
)

sens_2p <- k_sensitivity(df_week, "joint_2p_w", k_grid)
sens_3p <- k_sensitivity(df_week, "joint_3p_w", k_grid)

cat("\n================ k-sensitivity table: joint_2p_w ================\n")
print(sens_2p$table)

cat("\n================ k-sensitivity table: joint_3p_w ================\n")
print(sens_3p$table)

# ---- (G) Quick plots (optional)
plot_trend <- function(data, model, title, woy_fixed = 26) {
  nd <- data.frame(
    t = seq(min(data$t), max(data$t), length.out = 400),
    woy = woy_fixed
  )
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  invlogit <- function(x) 1 / (1 + exp(-x))
  nd$fit <- invlogit(pr$fit)
  nd$lo  <- invlogit(pr$fit - 1.96 * pr$se.fit)
  nd$hi  <- invlogit(pr$fit + 1.96 * pr$se.fit)
  nd$week_start <- data$week_start[1] + (nd$t - 1) * 7
  
  ggplot(nd, aes(x = week_start, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 1) +
    theme_minimal() +
    labs(x = NULL, y = "Probability", title = title)
}

print(plot_trend(df_week, m_2p, "Baseline joint exceedance (>=2)"))
print(plot_trend(df_week, m_3p, "Baseline joint exceedance (>=3)"))

# flexiility test
# 1) model is al gefit
m_2p_fx <- gam(
  joint_2p_w ~ s(t, k = 60, fx = TRUE) + s(woy, bs = "cc", k = 20),
  family = binomial,
  method = "REML",
  data = df_week,
  knots = list(woy = c(0.5, 53.5))
)

# 2) check
gam.check(m_2p_fx)
          
          
m_2p_reml <- gam(
  joint_2p_w ~ s(t, k = 30) + s(woy, bs = "cc", k = 20),
  family = binomial,
  method = "REML",
  data = df_week,
  knots = list(woy = c(0.5, 53.5))
)
          
m_2p_gcv  <- gam(joint_2p_w ~ s(t, k=30) + s(woy, bs="cc", k=20),
                           family=binomial, method="GCV.Cp",
                           data=df_week, knots=list(woy=c(0.5,53.5)))
          
summary(m_2p_reml)
summary(m_2p_gcv)

# lamdas bekijken
m_2p_reml$sp      # smoothness parameters (lambda's)
m_2p_reml$edf     # edf per smooth

set.seed(1)
n <- nrow(df_week)
cut <- floor(0.8*n)

train <- df_week[1:cut, ]
test  <- df_week[(cut+1):n, ]

m_reml <- gam(joint_2p_w ~ s(t, k=30) + s(woy, bs="cc", k=20),
              family=binomial, method="REML",
              data=train, knots=list(woy=c(0.5,53.5)))

m_gcv  <- gam(joint_2p_w ~ s(t, k=30) + s(woy, bs="cc", k=20),
              family=binomial, method="GCV.Cp",
              data=train, knots=list(woy=c(0.5,53.5)))

p_reml <- predict(m_reml, newdata=test, type="response")
p_gcv  <- predict(m_gcv,  newdata=test, type="response")

logloss <- function(y, p) {
  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1-eps)
  -mean(y*log(p) + (1-y)*log(1-p))
}

logloss(test$joint_2p_w, p_reml)
logloss(test$joint_2p_w, p_gcv)



# lamda backtesten ofzo
m_reml <- gam(
  joint_2p_w ~ s(t, k=30) + s(woy, bs="cc", k=20),
  family = binomial, method = "REML",
  data = df_week,
  knots = list(woy = c(0.5, 53.5))
)

m_reml$sp
summary(m_reml)$s.table

fit_fixed_sp <- function(sp_t, sp_woy, data) {
  gam(
    joint_2p_w ~ s(t, k=30) + s(woy, bs="cc", k=20),
    family = binomial, method = "REML",
    sp = c(sp_t, sp_woy),
    data = data,
    knots = list(woy = c(0.5, 53.5))
  )
}

library(dplyr)
library(tidyr)
library(purrr)

sp_grid <- expand.grid(
  sp_t   = c(1, 10, 50, 100, 500, 2000),
  sp_woy = c(1, 10, 50, 100, 500)
)

res_sp <- pmap_dfr(sp_grid, function(sp_t, sp_woy) {
  m <- fit_fixed_sp(sp_t, sp_woy, df_week)
  
  tibble(
    sp_t   = sp_t,
    sp_woy = sp_woy,
    edf_t  = summary(m)$s.table["s(t)", "edf"],
    edf_woy= summary(m)$s.table["s(woy)", "edf"],
    p_t    = summary(m)$s.table["s(t)", "p-value"],
    AIC    = AIC(m)
  )
})

res_sp %>%
  arrange(AIC) %>%
  print(n = 20)

set.seed(1)
n <- nrow(df_week)
cut <- floor(0.8*n)

train <- df_week[1:cut, ]
test  <- df_week[(cut+1):n, ]

logloss <- function(y, p) {
  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1-eps)
  -mean(y*log(p) + (1-y)*log(1-p))
}

res_sp_cv <- pmap_dfr(sp_grid, function(sp_t, sp_woy) {
  
  m <- fit_fixed_sp(sp_t, sp_woy, train)
  p <- predict(m, newdata = test, type = "response")
  
  tibble(
    sp_t   = sp_t,
    sp_woy = sp_woy,
    edf_t  = summary(m)$s.table["s(t)", "edf"],
    edf_woy= summary(m)$s.table["s(woy)", "edf"],
    logloss= logloss(test$joint_2p_w, p)
  )
})

res_sp_cv %>%
  arrange(logloss) %>%
  print(n = 20)




