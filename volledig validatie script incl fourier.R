# ============================================================
# FULL ONE-SHOT WORKFLOW (Alicante) — GAM extremes + diagnostics + validation
# Order (story + supervisor-proof):
#  1) Data + split + threshold on TRAIN (no leakage)
#  2) Baseline GAM (REML)
#  3) k-check loop (basis adequacy)
#  4) Residual checks + ACF
#  5) Add lag (misspec fix) + re-run k-check + ACF
#  6) Optional: BAM AR(1) for serial correlation robustness (on lag+nao spec)
#  7) Add NAO + NAO×season (hypothesis) + k-check + diagnostics
#  8) Validation (TEST): Brier/LogLoss/AUC + calibration plots
#  9) Robustness:
#      - REML vs GCV for baseline
#      - spline robustness: doy cc vs tp; time tp vs cr
#      - k-sensitivity (bigger k) for baseline
#      - GAM vs BAM (AR1) overlays
# 10) Final selection: choose best by primary score (LogLoss) then Brier, then AUC
# Outputs: CSV tables + PNG diagnostics in OUT_DIR
# ============================================================

# ----------------------------
# 0) Packages + settings
# ----------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(mgcv)
  library(ggplot2)
  library(pROC)
  library(tibble)
  library(stringr)
})

# >>>> CHANGE THIS <<<<
DATA_DIR <- "C:/Users/Sara/Downloads/Case study"
setwd(DATA_DIR)

STATION_NAME <- "Madrid"
STATION_FILE <- "madrid.xlsx"
NAO_FILE     <- "NAO_cleaned.xlsx"   # must contain DATE + nao column

DATE_MIN <- as.Date("1970-01-01")
DATE_MAX <- as.Date("2020-12-31")

# Exceedance definition
WET_DAY_MM <- 1
P_EXTREME  <- 0.95

# Validation split (time-block)
TRAIN_END <- as.Date("2010-12-31")

# Cyclic seasonality settings
KNOTS_DOY <- list(doy = c(0.5, 366.5))

# Initial k (upper bound) for auto-tune
K_TIME_START <- 30
K_DOY_START  <- 20
K_TIME_MAX   <- 120
K_DOY_MAX    <- 60
K_STEP_TIME  <- 10
K_STEP_DOY   <- 5
MAX_K_ITER   <- 6
KCHECK_P_CUT <- 0.05
KCHECK_KI_CUT <- 0.98   # set 0.99 if stricter

# Optional k-sensitivity for baseline
K_TIME_BIG <- 60
K_DOY_BIG  <- 30

# Output directory
OUT_DIR <- file.path(DATA_DIR, "outputs_madrid_ALL")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Generic helpers
# ----------------------------
invlogit <- function(x) 1/(1+exp(-x))
clip01 <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1 - eps)

logloss <- function(y, p, eps = 1e-15) {
  p <- clip01(p, eps)
  -mean(y*log(p) + (1-y)*log(1-p))
}
brier <- function(y, p) mean((y - p)^2)

safe_auc <- function(y, p) {
  tryCatch(as.numeric(pROC::auc(y, p)), error=function(e) NA_real_)
}

save_gamcheck_png <- function(mod, filename) {
  png(filename, width = 1600, height = 1100, res = 150)
  par(mar = c(4,4,2,1))
  gam.check(mod)
  dev.off()
}

save_acf_png <- function(resid_vec, filename, main = "ACF (Pearson residuals)") {
  png(filename, width = 1400, height = 900, res = 150)
  acf(resid_vec, na.action = na.pass, main = main)
  dev.off()
}

run_kcheck_df <- function(mod) {
  as.data.frame(mgcv::k.check(mod)) %>%
    tibble::rownames_to_column("term")
}

k_fail_flags <- function(kc, p_cut = 0.05, kindex_cut = 0.98) {
  cn <- names(kc)
  kindex_col <- cn[grepl("k-index", cn, ignore.case = TRUE)][1]
  pval_col   <- cn[grepl("p.value|p-value|p\\.value", cn, ignore.case = TRUE)][1]
  edf_col    <- cn[grepl("^edf$", cn, ignore.case = TRUE)][1]
  if (is.na(kindex_col) || is.na(pval_col) || is.na(edf_col)) return(rep(FALSE, nrow(kc)))
  (kc[[kindex_col]] < kindex_cut) & (kc[[pval_col]] < p_cut) & (kc[[edf_col]] > 1)
}

# Compact model summary row
model_report_row <- function(mod, model_name) {
  st <- summary(mod)
  get_edf <- function(term) {
    if (!is.null(st$s.table) && term %in% rownames(st$s.table)) st$s.table[term,"edf"] else NA_real_
  }
  data.frame(
    model = model_name,
    method = mod$method,
    n = nrow(mod$model),
    aic = tryCatch(AIC(mod), error=function(e) NA_real_),
    edf_time = get_edf("s(time)"),
    edf_doy  = get_edf("s(doy)"),
    sp = paste(round(mod$sp, 6), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

# Calibration plot helpers
calibration_bins <- function(y, p, n_bins = 10) {
  d <- data.frame(y=y, p=p)
  
  # probeer quantile bins
  qs <- quantile(d$p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  qs <- unique(as.numeric(qs))
  
  # als te weinig unieke breaks: gebruik gelijke afstand bins
  if (length(qs) < 3) {
    qs <- seq(min(d$p, na.rm=TRUE), max(d$p, na.rm=TRUE), length.out = min(n_bins + 1, length(unique(d$p)) + 1))
    qs <- unique(as.numeric(qs))
  }
  
  # als nog steeds te weinig (extreem degenerate): jitter klein beetje
  if (length(qs) < 3) {
    d$p <- d$p + rnorm(nrow(d), sd = 1e-8)
    qs <- quantile(d$p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
    qs <- unique(as.numeric(qs))
  }
  
  d$bin <- cut(d$p, breaks = qs, include.lowest = TRUE)
  
  out <- d %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      p_mean = mean(p),
      obs = mean(y),
      .groups="drop"
    ) %>%
    mutate(
      se = sqrt(obs*(1-obs)/pmax(n,1)),
      lo = pmax(obs - 1.96*se, 0),
      hi = pmin(obs + 1.96*se, 1)
    )
  out
}
plot_calibration <- function(y, p, title, n_bins = 10, out_png = NULL) {
  bins <- calibration_bins(y, p, n_bins)
  g <- ggplot(bins, aes(x = p_mean, y = obs)) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    geom_point() +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=0) +
    labs(x="Mean predicted probability", y="Observed frequency", title=title) +
    coord_equal(
      xlim=c(0, max(bins$p_mean, na.rm=TRUE)),
      ylim=c(0, max(bins$hi, na.rm=TRUE))
    ) +
    theme_minimal()
  if (!is.null(out_png)) ggsave(out_png, g, width=7, height=5, dpi=150)
  g
}

# Curves for overlays (time trend + seasonality)
make_curve_ci <- function(mod, df_ref, doy_fixed=180, lag1_fixed=0, nao_fixed=0,
                          n=600, label="model") {
  nd <- data.frame(
    time = seq(min(df_ref$time), max(df_ref$time), length.out=n),
    doy  = doy_fixed
  )
  if ("lag1" %in% all.vars(formula(mod))) nd$lag1 <- lag1_fixed
  if ("nao"  %in% all.vars(formula(mod))) nd$nao  <- nao_fixed
  
  pr <- predict(mod, newdata=nd, type="link", se.fit=TRUE)
  out <- nd
  out$date <- min(df_ref$DATE) + round(out$time)
  out$fit  <- invlogit(pr$fit)
  out$lo   <- invlogit(pr$fit - 1.96*pr$se.fit)
  out$hi   <- invlogit(pr$fit + 1.96*pr$se.fit)
  out$which <- label
  out
}
make_season_ci <- function(mod, df_ref, time_fixed=NULL, lag1_fixed=0, nao_fixed=0,
                           n=366, label="model") {
  if (is.null(time_fixed)) time_fixed <- median(df_ref$time)
  nd <- data.frame(doy=1:366, time=time_fixed)
  if ("lag1" %in% all.vars(formula(mod))) nd$lag1 <- lag1_fixed
  if ("nao"  %in% all.vars(formula(mod))) nd$nao  <- nao_fixed
  
  pr <- predict(mod, newdata=nd, type="link", se.fit=TRUE)
  out <- nd
  out$fit  <- invlogit(pr$fit)
  out$lo   <- invlogit(pr$fit - 1.96*pr$se.fit)
  out$hi   <- invlogit(pr$fit + 1.96*pr$se.fit)
  out$which <- label
  out
}

# ----------------------------
# 2) Data read + clean + split + threshold on TRAIN
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
read_nao_daily <- function(path) {
  # change nao_index_cdas if your NAO file uses a different name
  read_excel(path) %>%
    mutate(DATE = as.Date(DATE)) %>%
    transmute(DATE, nao = nao_index_cdas) %>%
    filter(!is.na(DATE), !is.na(nao))
}
add_time_vars <- function(df) {
  df %>% arrange(DATE) %>%
    mutate(
      time = as.numeric(DATE - min(DATE)),
      doy  = yday(DATE)
    )
}
make_threshold <- function(rr_mm, wet_day_mm = 1, p_extreme = 0.95) {
  quantile(rr_mm[rr_mm > wet_day_mm], probs = p_extreme, type = 8, na.rm = TRUE)
}
make_exc <- function(df, thr) df %>% mutate(exc = as.integer(rr_mm > thr))
add_lag1 <- function(df) {
  df %>% arrange(DATE) %>%
    mutate(lag1 = dplyr::lag(exc, 1)) %>%
    filter(!is.na(lag1))
}

df0 <- read_ecad_rr(STATION_FILE, STATION_NAME, DATE_MIN, DATE_MAX) %>% add_time_vars()

train0 <- df0 %>% filter(DATE <= TRAIN_END)
test0  <- df0 %>% filter(DATE >  TRAIN_END)

thr_train <- make_threshold(train0$rr_mm, WET_DAY_MM, P_EXTREME)
train0 <- make_exc(train0, thr_train)
test0  <- make_exc(test0,  thr_train)

df_full <- bind_rows(train0, test0) %>% arrange(DATE) %>% add_time_vars()
df_ref  <- df_full

train_lag <- add_lag1(train0)
# For lag1 in test, use last train day to avoid dropping first test day
test_lag <- bind_rows(train0, test0) %>%
  add_lag1() %>%
  filter(DATE > TRAIN_END)

nao_daily <- read_nao_daily(NAO_FILE)
train_lag_nao <- train_lag %>% left_join(nao_daily, by="DATE") %>% filter(!is.na(nao))
test_lag_nao  <- test_lag  %>% left_join(nao_daily, by="DATE") %>% filter(!is.na(nao))
standardize_nao_main <- function(train_df, test_df) {
  mu <- mean(train_df$nao, na.rm=TRUE)
  sdv <- sd(train_df$nao, na.rm=TRUE)
  if (is.na(sdv) || sdv == 0) sdv <- 1
  train_df$nao <- (train_df$nao - mu) / sdv
  test_df$nao  <- (test_df$nao  - mu) / sdv
  list(train=train_df, test=test_df)
}
nao_std <- standardize_nao_main(train_lag_nao, test_lag_nao)
train_lag_nao <- nao_std$train
test_lag_nao  <- nao_std$test

cat("Station:", STATION_NAME, "\n")
cat("Train rows:", nrow(train0), " Test rows:", nrow(test0), "\n")
cat("Train exceedance rate:", round(mean(train0$exc), 4), "\n")
cat("Test exceedance rate :", round(mean(test0$exc),  4), "\n")
cat("Threshold (train p):", round(thr_train,2), "mm\n")

# ----------------------------
# 3) Model builders (so we can refit in loops)
# ----------------------------
fit_baseline <- function(train_df, k_time, k_doy, bs_time="tp", bs_doy="cc", method="REML") {
  gam(exc ~ s(time, bs=bs_time, k=k_time) + s(doy, bs=bs_doy, k=k_doy),
      data=train_df, family=binomial("logit"),
      method=method,
      knots=if (bs_doy=="cc") KNOTS_DOY else NULL)
}

add_fourier_terms <- function(dat, K, period = 365.25, prefix = "harm") {
  if (K < 1) stop("K must be >= 1 for Fourier terms")
  out <- dat
  for (j in seq_len(K)) {
    ang <- 2 * pi * j * out$doy / period
    out[[paste0(prefix, "_sin", j)]] <- sin(ang)
    out[[paste0(prefix, "_cos", j)]] <- cos(ang)
  }
  out
}

fit_fourier <- function(train_df, K, k_time, period = 365.25, method = "REML") {
  train_df <- add_fourier_terms(train_df, K = K, period = period, prefix = "harm")
  seasonal_terms <- as.vector(rbind(
    paste0("harm_sin", seq_len(K)),
    paste0("harm_cos", seq_len(K))
  ))
  rhs <- c(sprintf("s(time, k = %d)", k_time), seasonal_terms)
  form <- as.formula(paste("exc ~", paste(rhs, collapse = " + ")))
  mod <- gam(form,
             data=train_df, family=binomial("logit"),
             method=method)
  list(model = mod, data = train_df)
}

fit_lag <- function(train_df, k_time, k_doy, bs_time="tp", bs_doy="cc", method="REML") {
  gam(exc ~ s(time, bs=bs_time, k=k_time) + s(doy, bs=bs_doy, k=k_doy) + lag1,
      data=train_df, family=binomial("logit"),
      method=method,
      knots=if (bs_doy=="cc") KNOTS_DOY else NULL)
}

fit_lag_nao <- function(train_df, k_time, k_doy, bs_time="tp", bs_doy="cc", method="REML") {
  gam(exc ~ s(time, bs=bs_time, k=k_time) + s(doy, bs=bs_doy, k=k_doy) + lag1 + nao,
      data=train_df, family=binomial("logit"),
      method=method,
      knots=if (bs_doy=="cc") KNOTS_DOY else NULL)
}

fit_int <- function(train_df, k_time, k_doy, k_ti_nao=10, bs_time="tp", method="REML") {
  # doy side of ti uses cyclic cc; k uses k_doy
  gam(exc ~ s(time, bs=bs_time, k=k_time) +
        s(doy, bs="cc", k=k_doy) +
        lag1 + nao +
        ti(nao, doy, bs=c("tp","cc"), k=c(k_ti_nao, k_doy)),
      data=train_df, family=binomial("logit"),
      method=method, knots=KNOTS_DOY)
}

auto_tune_k_time_doy <- function(model_name, fit_fun, fit_args,
                                 k_time_start, k_doy_start,
                                 k_step_time=10, k_step_doy=5,
                                 k_time_max=120, k_doy_max=60,
                                 max_iter=6, p_cut=0.05, kindex_cut=0.98) {
  
  k_time <- k_time_start
  k_doy  <- k_doy_start
  logs <- list()
  
  for (it in 1:max_iter) {
    mod <- do.call(fit_fun, c(fit_args, list(k_time=k_time, k_doy=k_doy)))
    kc  <- run_kcheck_df(mod)
    kc$model <- model_name
    kc$iter  <- it
    kc$k_time <- k_time
    kc$k_doy  <- k_doy
    logs[[it]] <- kc
    
    fail <- k_fail_flags(kc, p_cut=p_cut, kindex_cut=kindex_cut)
    need_time <- any(fail & grepl("time", kc$term))
    need_doy  <- any(fail & grepl("doy",  kc$term))
    
    if (!need_time && !need_doy) {
      return(list(final_model=mod,
                  log=bind_rows(logs),
                  final_k=data.frame(model=model_name, k_time=k_time, k_doy=k_doy, iter=it)))
    }
    
    if (need_time) k_time <- min(k_time + k_step_time, k_time_max)
    if (need_doy)  k_doy  <- min(k_doy  + k_step_doy,  k_doy_max)
    
    if ((need_time && k_time >= k_time_max) || (need_doy && k_doy >= k_doy_max)) break
  }
  
  # return last fit if max_iter/cap reached
  mod <- do.call(fit_fun, c(fit_args, list(k_time=k_time, k_doy=k_doy)))
  return(list(final_model=mod,
              log=bind_rows(logs),
              final_k=data.frame(model=model_name, k_time=k_time, k_doy=k_doy, iter=max_iter,
                                 max_iter_or_cap=TRUE)))
}

# ----------------------------
# 4) Step 2-3: Baseline REML + k-check loop
# ----------------------------
t_base <- auto_tune_k_time_doy(
  model_name="baseline_REML",
  fit_fun=fit_baseline,
  fit_args=list(train_df=train0, bs_time="tp", bs_doy="cc", method="REML"),
  k_time_start=K_TIME_START, k_doy_start=K_DOY_START,
  k_step_time=K_STEP_TIME, k_step_doy=K_STEP_DOY,
  k_time_max=K_TIME_MAX, k_doy_max=K_DOY_MAX,
  max_iter=MAX_K_ITER, p_cut=KCHECK_P_CUT, kindex_cut=KCHECK_KI_CUT
)
m_base_reml <- t_base$final_model
FOURIER_K <- 1:4
fourier_fits <- lapply(FOURIER_K, function(K) {
  fit_fourier(train0, K = K, k_time = t_base$final_k$k_time, period = 365.25, method = "REML")
})

# Diagnostics: baseline
save_gamcheck_png(m_base_reml, file.path(OUT_DIR, "gamcheck_baseline_REML.png"))
save_acf_png(residuals(m_base_reml, type="pearson"),
             file.path(OUT_DIR, "acf_baseline_REML.png"),
             main="ACF Pearson residuals — baseline REML")

# ----------------------------
# 5) Step 5-6: Add lag + re-run k-check
# ----------------------------
t_lag <- auto_tune_k_time_doy(
  model_name="lag_REML",
  fit_fun=fit_lag,
  fit_args=list(train_df=train_lag, bs_time="tp", bs_doy="cc", method="REML"),
  k_time_start=t_base$final_k$k_time, k_doy_start=t_base$final_k$k_doy,
  k_step_time=K_STEP_TIME, k_step_doy=K_STEP_DOY,
  k_time_max=K_TIME_MAX, k_doy_max=K_DOY_MAX,
  max_iter=MAX_K_ITER, p_cut=KCHECK_P_CUT, kindex_cut=KCHECK_KI_CUT
)
m_lag_reml <- t_lag$final_model

save_gamcheck_png(m_lag_reml, file.path(OUT_DIR, "gamcheck_lag_REML.png"))
save_acf_png(residuals(m_lag_reml, type="pearson"),
             file.path(OUT_DIR, "acf_lag_REML.png"),
             main="ACF Pearson residuals — lag REML")

# ----------------------------
# 6) Step 7: Add NAO + k-check
# ----------------------------
t_lagnao <- auto_tune_k_time_doy(
  model_name="lag_nao_REML",
  fit_fun=fit_lag_nao,
  fit_args=list(train_df=train_lag_nao, bs_time="tp", bs_doy="cc", method="REML"),
  k_time_start=t_lag$final_k$k_time, k_doy_start=t_lag$final_k$k_doy,
  k_step_time=K_STEP_TIME, k_step_doy=K_STEP_DOY,
  k_time_max=K_TIME_MAX, k_doy_max=K_DOY_MAX,
  max_iter=MAX_K_ITER, p_cut=KCHECK_P_CUT, kindex_cut=KCHECK_KI_CUT
)
m_lag_nao_reml <- t_lagnao$final_model

save_gamcheck_png(m_lag_nao_reml, file.path(OUT_DIR, "gamcheck_lag_nao_REML.png"))
save_acf_png(residuals(m_lag_nao_reml, type="pearson"),
             file.path(OUT_DIR, "acf_lag_nao_REML.png"),
             main="ACF Pearson residuals — lag+nao REML")

# ----------------------------
# 7) Step 7 (cont.): Add NAO×season interaction + k-check
# ----------------------------
t_int <- auto_tune_k_time_doy(
  model_name="int_REML",
  fit_fun=fit_int,
  fit_args=list(train_df=train_lag_nao, k_ti_nao=10, bs_time="tp", method="REML"),
  k_time_start=t_lagnao$final_k$k_time, k_doy_start=t_lagnao$final_k$k_doy,
  k_step_time=K_STEP_TIME, k_step_doy=K_STEP_DOY,
  k_time_max=K_TIME_MAX, k_doy_max=K_DOY_MAX,
  max_iter=MAX_K_ITER, p_cut=KCHECK_P_CUT, kindex_cut=KCHECK_KI_CUT
)
m_int_reml <- t_int$final_model

save_gamcheck_png(m_int_reml, file.path(OUT_DIR, "gamcheck_int_REML.png"))
save_acf_png(residuals(m_int_reml, type="pearson"),
             file.path(OUT_DIR, "acf_int_REML.png"),
             main="ACF Pearson residuals — lag+nao+int REML")

# ----------------------------
# 8) Save k-tuning logs (appendix proof)
# ----------------------------
k_log   <- bind_rows(t_base$log, t_lag$log, t_lagnao$log, t_int$log)
k_final <- bind_rows(t_base$final_k, t_lag$final_k, t_lagnao$final_k, t_int$final_k)

write.csv(k_log,   file.path(OUT_DIR, "k_tuning_log.csv"), row.names=FALSE)
write.csv(k_final, file.path(OUT_DIR, "k_final_by_step.csv"), row.names=FALSE)

# ----------------------------
# 9) Step 9: Robustness checks
# ----------------------------

# 9a) REML vs GCV for baseline (using final baseline k's)
m_base_gcv <- fit_baseline(train0,
                           k_time=t_base$final_k$k_time,
                           k_doy=t_base$final_k$k_doy,
                           bs_time="tp", bs_doy="cc", method="GCV.Cp")
save_gamcheck_png(m_base_gcv, file.path(OUT_DIR, "gamcheck_baseline_GCV.png"))
save_acf_png(residuals(m_base_gcv, type="pearson"),
             file.path(OUT_DIR, "acf_baseline_GCV.png"),
             main="ACF Pearson residuals — baseline GCV")

# 9b) k-sensitivity for baseline (bigger k, still REML)
m_base_kbig <- fit_baseline(train0, k_time=K_TIME_BIG, k_doy=K_DOY_BIG, bs_time="tp", bs_doy="cc", method="REML")
save_gamcheck_png(m_base_kbig, file.path(OUT_DIR, "gamcheck_baseline_REML_kbig.png"))

# 9c) spline robustness (baseline): doy cc vs tp; time tp vs cr
m_base_doy_tp <- fit_baseline(train0, k_time=t_base$final_k$k_time, k_doy=t_base$final_k$k_doy,
                              bs_time="tp", bs_doy="tp", method="REML")
save_gamcheck_png(m_base_doy_tp, file.path(OUT_DIR, "gamcheck_baseline_doyTP.png"))

m_base_time_cr <- fit_baseline(train0, k_time=t_base$final_k$k_time, k_doy=t_base$final_k$k_doy,
                               bs_time="cr", bs_doy="cc", method="REML")
save_gamcheck_png(m_base_time_cr, file.path(OUT_DIR, "gamcheck_baseline_timeCR.png"))

# 9d) BAM AR(1) robustness for lag+nao spec (serial correlation)
r_lagnao <- residuals(m_lag_nao_reml, type="pearson")
rho_hat  <- acf(r_lagnao, plot=FALSE)$acf[2]
df_bam <- train_lag_nao %>%
  arrange(DATE) %>%
  mutate(AR.start = c(TRUE, diff(DATE) != 1))

m_lag_nao_bam <- bam(
  exc ~ s(time, k=t_lagnao$final_k$k_time) +
    s(doy, bs="cc", k=t_lagnao$final_k$k_doy) +
    lag1 + nao,
  data=df_bam,
  family=binomial("logit"),
  method="REML",
  rho=rho_hat,
  AR.start=df_bam$AR.start,
  knots=KNOTS_DOY
)
save_gamcheck_png(m_lag_nao_bam, file.path(OUT_DIR, "gamcheck_lag_nao_BAM_AR1.png"))
save_acf_png(residuals(m_lag_nao_bam, type="pearson"),
             file.path(OUT_DIR, "acf_lag_nao_BAM_AR1.png"),
             main=paste0("ACF Pearson residuals — BAM AR1 rho=", round(rho_hat,3)))

# ----------------------------
# 10) Step 8: Validation on TEST (predictive)
# ----------------------------
score_model <- function(model_name, mod, test_df) {
  p <- predict(mod, newdata=test_df, type="response")
  data.frame(
    model=model_name,
    n_test=nrow(test_df),
    brier=brier(test_df$exc, p),
    logloss=logloss(test_df$exc, p),
    auc=safe_auc(test_df$exc, p),
    stringsAsFactors = FALSE
  )
}

score_fourier <- function(model_name, fit_obj, test_df, K, period = 365.25) {
  test_df2 <- add_fourier_terms(test_df, K = K, period = period, prefix = "harm")
  p <- predict(fit_obj$model, newdata = test_df2, type = "response")
  data.frame(
    model  = model_name,
    n_test = nrow(test_df2),
    brier  = brier(test_df2$exc, p),
    logloss= logloss(test_df2$exc, p),
    auc    = safe_auc(test_df2$exc, p),
    stringsAsFactors = FALSE
  )
}

val_rows <- list()

# Base comparison set: full test (baseline-only models)
p_clim_full <- mean(train0$exc)
val_rows$climatology_full <- data.frame(model="climatology_full", n_test=nrow(test0),
                                        brier=brier(test0$exc, rep(p_clim_full, nrow(test0))),
                                        logloss=logloss(test0$exc, rep(p_clim_full, nrow(test0))),
                                        auc=safe_auc(test0$exc, rep(p_clim_full, nrow(test0))))
val_rows$baseline_REML <- score_model("baseline_REML", m_base_reml, test0)
val_rows$baseline_GCV  <- score_model("baseline_GCV",  m_base_gcv,  test0)
val_rows$baseline_REML_kbig <- score_model("baseline_REML_kbig", m_base_kbig, test0)
val_rows$baseline_doyTP <- score_model("baseline_doyTP", m_base_doy_tp, test0)
val_rows$baseline_timeCR <- score_model("baseline_timeCR", m_base_time_cr, test0)
for (i in seq_along(FOURIER_K)) {
  K <- FOURIER_K[[i]]
  model_name <- paste0("baseline_fourier_K", K)
  val_rows[[model_name]] <- score_fourier(
    model_name = model_name,
    fit_obj = fourier_fits[[i]],
    test_df = test0,
    K = K,
    period = 365.25
  )
}

# Lag comparison set: lag-capable models (needs lag1)
val_rows$lag_REML <- score_model("lag_REML", m_lag_reml, test_lag)

# NAO comparison set: all models evaluated on NAO subset for fair ranking
test_nao <- test0 %>% left_join(nao_daily, by="DATE") %>% filter(!is.na(nao))
if (nrow(test_nao) > 0) {
  test_nao <- standardize_nao_main(train_lag_nao, test_nao)$test
  test_nao_lag <- bind_rows(train_lag_nao, test_nao) %>%
    add_lag1() %>%
    filter(DATE > TRAIN_END)
  
  for (i in seq_along(FOURIER_K)) {
    K <- FOURIER_K[[i]]
    model_name <- paste0("baseline_fourier_K", K, "_NAOsubset")
    val_rows[[model_name]] <- score_fourier(
      model_name = model_name,
      fit_obj = fourier_fits[[i]],
      test_df = test_nao,
      K = K,
      period = 365.25
    )
  }
  
  p_clim_nao <- rep(p_clim_full, nrow(test_nao))
  val_rows$climatology_NAOsubset <- data.frame(model="climatology_NAOsubset",
                                               n_test=nrow(test_nao),
                                               brier=brier(test_nao$exc, p_clim_nao),
                                               logloss=logloss(test_nao$exc, p_clim_nao),
                                               auc=safe_auc(test_nao$exc, p_clim_nao))
  val_rows$baseline_REML_NAOsubset <- score_model("baseline_REML_NAOsubset", m_base_reml, test_nao)
  val_rows$lag_REML_NAOsubset <- score_model("lag_REML_NAOsubset", m_lag_reml, test_nao_lag)
  val_rows$lag_nao_REML <- score_model("lag_nao_REML", m_lag_nao_reml, test_lag_nao)
  val_rows$int_REML <- score_model("int_REML", m_int_reml, test_lag_nao)
  val_rows$lag_nao_BAM_AR1 <- score_model("lag_nao_BAM_AR1", m_lag_nao_bam, test_lag_nao)
}

val_table <- bind_rows(val_rows)
write.csv(val_table, file.path(OUT_DIR, "validation_scores_ALL.csv"), row.names=FALSE)
cat("\nSaved:", file.path(OUT_DIR, "validation_scores_ALL.csv"), "\n")
print(val_table)

# Fair ranking (NAO subset only)
val_rank <- val_table %>%
  filter(grepl("NAOsubset|lag_nao|int_REML", model)) %>%
  arrange(logloss, brier, desc(auc))
write.csv(val_rank, file.path(OUT_DIR, "validation_scores_NAOsubset_ranked.csv"), row.names=FALSE)
cat("\nSaved:", file.path(OUT_DIR, "validation_scores_NAOsubset_ranked.csv"), "\n")

# Calibration plots for top-3 by LogLoss on NAO subset + baseline climatology
top3_nao <- head(val_rank$model, 3)
if (exists("test_nao") && nrow(test_nao) > 0) {
  for (mname in unique(c("climatology_full", top3_nao))) {
    if (mname == "climatology_full") {
      p <- rep(mean(train0$exc), nrow(test0))
      plot_calibration(test0$exc, p,
                       title=paste0(STATION_NAME, " — climatology (TEST full)"),
                       out_png=file.path(OUT_DIR, "calibration_climatology_full.png"))
    } else if (mname %in% c("baseline_REML","baseline_GCV","baseline_REML_kbig","baseline_doyTP","baseline_timeCR",
                            "baseline_REML_NAOsubset")) {
      mod <- if (mname == "baseline_REML_NAOsubset") m_base_reml else get(mname, inherits=TRUE)
      p <- predict(mod, newdata=test_nao, type="response")
      plot_calibration(test_nao$exc, p,
                       title=paste0(STATION_NAME, " — ", mname, " calibration (TEST NAO subset)"),
                       out_png=file.path(OUT_DIR, paste0("calibration_", mname, ".png")))
    } else if (mname %in% c("lag_REML_NAOsubset")) {
      p <- predict(m_lag_reml, newdata=test_nao_lag, type="response")
      plot_calibration(test_nao_lag$exc, p,
                       title=paste0(STATION_NAME, " — ", mname, " calibration (TEST NAO subset)"),
                       out_png=file.path(OUT_DIR, paste0("calibration_", mname, ".png")))
    } else {
      # lag+nao / int / bam: use test_lag_nao
      mod <- switch(mname,
                    "lag_nao_REML" = m_lag_nao_reml,
                    "int_REML" = m_int_reml,
                    "lag_nao_BAM_AR1" = m_lag_nao_bam,
                    NULL)
      if (!is.null(mod)) {
        p <- predict(mod, newdata=test_lag_nao, type="response")
        plot_calibration(test_lag_nao$exc, p,
                         title=paste0(STATION_NAME, " — ", mname, " calibration (TEST NAO subset)"),
                         out_png=file.path(OUT_DIR, paste0("calibration_", mname, ".png")))
      }
    }
  }
}

# ----------------------------
# 11) Overlays: base REML vs GCV, and lag+nao GAM vs BAM
# ----------------------------
curve_base_reml <- make_curve_ci(m_base_reml, df_ref, label="Baseline GAM (REML)")
curve_base_gcv  <- make_curve_ci(m_base_gcv,  df_ref, label="Baseline GAM (GCV)")
curves_base <- bind_rows(curve_base_reml, curve_base_gcv)

p_base_overlay <- ggplot(curves_base, aes(x=date, y=fit, linetype=which)) +
  geom_line(linewidth=1) +
  labs(x=NULL, y="P(exceedance)",
       title=paste0(STATION_NAME, " — Baseline: REML vs GCV (doy=180)")) +
  theme_minimal()
ggsave(file.path(OUT_DIR, "overlay_base_REML_vs_GCV.png"),
       p_base_overlay, width=10, height=5, dpi=150)

curve_gam <- make_curve_ci(m_lag_nao_reml, df_ref, label="Lag+NAO GAM (REML)")
curve_bam <- make_curve_ci(m_lag_nao_bam,  df_ref, label=paste0("Lag+NAO BAM AR1 rho=", round(rho_hat,3)))
curves_gb <- bind_rows(curve_gam, curve_bam)

p_gam_bam <- ggplot(curves_gb, aes(x=date, y=fit, linetype=which)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=which), alpha=0.12, color=NA) +
  geom_line(linewidth=1) +
  labs(x=NULL, y="P(exceedance)",
       title=paste0(STATION_NAME, " — Lag+NAO: GAM vs BAM AR(1) (doy=180)")) +
  theme_minimal()
ggsave(file.path(OUT_DIR, "overlay_lag_nao_GAM_vs_BAM.png"),
       p_gam_bam, width=11, height=5.5, dpi=150)

# ----------------------------
# 12) Concurvity reports (mainly for interaction model)
# ----------------------------
capture.output(concurvity(m_base_reml, full=TRUE),
               file=file.path(OUT_DIR, "concurvity_baseline_REML.txt"))
capture.output(concurvity(m_int_reml, full=TRUE),
               file=file.path(OUT_DIR, "concurvity_int_REML.txt"))

# ----------------------------
# 13) Final "best model" selection + write report (NAO subset ranking)
# ----------------------------
if (nrow(val_rank) == 0) {
  best <- val_table %>% slice(1)
  best_name <- best$model[1]
  selection_note <- "Selection rule (full test): minimize LogLoss; tie-breaker Brier; then maximize AUC."
} else {
  best <- val_rank %>% slice(1)
  best_name <- best$model[1]
  selection_note <- "Selection rule (NAO subset): minimize LogLoss; tie-breaker Brier; then maximize AUC."
}

cat("\nBEST on NAO subset (by LogLoss then Brier then AUC):\n")
print(best)

writeLines(
  c(
    paste0("Station: ", STATION_NAME),
    paste0("Train end: ", TRAIN_END),
    paste0("Threshold (train): ", round(thr_train,2), " mm (wet-day p=", P_EXTREME, ")"),
    "",
    selection_note,
    paste0("Best model: ", best_name),
    paste0("Scores (TEST): logloss=", round(best$logloss,6),
           "  brier=", round(best$brier,6),
           "  auc=", round(best$auc,6)),
    "",
    "k_final_by_step.csv and k_tuning_log.csv provide basis-dimension adequacy evidence.",
    "Diagnostics: gamcheck_*.png and acf_*.png.",
    "Calibration: calibration_*.png.",
    "Robustness overlays: overlay_*.png."
  ),
  con=file.path(OUT_DIR, "FINAL_REPORT.txt")
)

cat("\nDONE. All outputs in:\n", OUT_DIR, "\n")


# ============================================================
# 14) 3-FOLD TIME-BLOCKED CV (same specs) — scores + optional calibration
# ============================================================

FOLDS <- list(
  list(name="fold1", train_end=as.Date("1990-12-31"),
       test_start=as.Date("1991-01-01"), test_end=as.Date("2000-12-31")),
  list(name="fold2", train_end=as.Date("2000-12-31"),
       test_start=as.Date("2001-01-01"), test_end=as.Date("2010-12-31")),
  list(name="fold3", train_end=as.Date("2010-12-31"),
       test_start=as.Date("2011-01-01"), test_end=as.Date("2020-12-31"))
)

RUN_3FOLD_CV <- TRUE            # zet TRUE als je de langzame 3-fold CV wilt draaien
SAVE_FOLD_CALIBRATION <- TRUE     # zet FALSE als je geen extra plots wilt
CALIBRATION_BINS <- 10
CAL_X_MAX <- 0.015                # zoom x-as (anders wordt alles smal)
CAL_Y_MAX <- 0.03                 # zoom y-as

# --- fold-safe time/doy: per fold opnieuw (belangrijk!) ---
add_time_doy_fold <- function(df) {
  df %>% arrange(DATE) %>%
    mutate(
      time = as.numeric(DATE - min(DATE)),
      doy  = lubridate::yday(DATE)
    )
}

# --- threshold on TRAIN only ---
threshold_train <- function(train_df, wet_day_mm=WET_DAY_MM, p_extreme=P_EXTREME) {
  as.numeric(quantile(train_df$rr_mm[train_df$rr_mm > wet_day_mm],
                      probs=p_extreme, type=8, na.rm=TRUE))
}

# --- add exceedance ---
add_exc_thr <- function(df, thr) df %>% mutate(exc = as.integer(rr_mm > thr))

# --- lag1 (drop first day) ---
add_lag1_fold <- function(df){
  df %>% arrange(DATE) %>%
    mutate(lag1 = dplyr::lag(exc, 1)) %>%
    filter(!is.na(lag1))
}

# --- standardize NAO using TRAIN stats only (no leakage) ---
standardize_nao <- function(train_df, test_df) {
  mu <- mean(train_df$nao, na.rm=TRUE)
  sdv <- sd(train_df$nao, na.rm=TRUE)
  if (is.na(sdv) || sdv == 0) sdv <- 1
  train_df$nao <- (train_df$nao - mu) / sdv
  test_df$nao  <- (test_df$nao  - mu) / sdv
  list(train=train_df, test=test_df)
}

# --- calibration plot with zoom (fix "smal") ---
plot_calibration_zoom <- function(y, p, title, out_png=NULL,
                                  n_bins=CALIBRATION_BINS,
                                  x_max=CAL_X_MAX, y_max=CAL_Y_MAX) {
  bins <- calibration_bins(y, p, n_bins=n_bins)
  
  g <- ggplot(bins, aes(x=p_mean, y=obs)) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    geom_point() +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=0) +
    labs(x="Mean predicted probability", y="Observed frequency", title=title) +
    coord_cartesian(xlim=c(0, x_max), ylim=c(0, y_max)) +
    theme_minimal()
  
  if (!is.null(out_png)) ggsave(out_png, g, width=7, height=5, dpi=200)
  g
}

# --- scoring helper ---
score_df <- function(model_name, fold_name, thr, y, p) {
  data.frame(
    fold = fold_name,
    model = model_name,
    thr_mm = thr,
    n_test = length(y),
    event_rate = mean(y),
    brier = brier(y,p),
    logloss = logloss(y,p),
    auc = safe_auc(y,p),
    stringsAsFactors = FALSE
  )
}

cv_scores <- list()
cv_preds <- list()

if (RUN_3FOLD_CV) {
  for (sp in FOLDS) {
    
    cat("\n--- Running", sp$name, "train_end =", as.character(sp$train_end),
        "test =", as.character(sp$test_start), "to", as.character(sp$test_end), "\n")
    
    # split
    train <- df0 %>% filter(DATE <= sp$train_end)
    test  <- df0 %>% filter(DATE >= sp$test_start, DATE <= sp$test_end)
    
    # threshold on TRAIN only
    thr <- threshold_train(train, WET_DAY_MM, P_EXTREME)
    
    train <- train %>% add_exc_thr(thr) %>% add_time_doy_fold()
    test  <- test  %>% add_exc_thr(thr) %>% add_time_doy_fold()
    
    # climatology (from train)
    p_clim <- mean(train$exc)
    
    # --- baseline ---
    m_base <- fit_baseline(train, k_time=t_base$final_k$k_time, k_doy=t_base$final_k$k_doy,
                           bs_time="tp", bs_doy="cc", method="REML")
    p_base <- predict(m_base, newdata=test, type="response")
    
    cv_scores[[length(cv_scores)+1]] <- score_df("climatology", sp$name, thr, test$exc, rep(p_clim, nrow(test)))
    cv_scores[[length(cv_scores)+1]] <- score_df("baseline_REML", sp$name, thr, test$exc, p_base)
    
    # --- lag ---
    train_l <- add_lag1_fold(train)
    test_l  <- add_lag1_fold(test)
    
    m_lag <- fit_lag(train_l, k_time=t_lag$final_k$k_time, k_doy=t_lag$final_k$k_doy,
                     bs_time="tp", bs_doy="cc", method="REML")
    p_lag <- predict(m_lag, newdata=test_l, type="response")
    cv_scores[[length(cv_scores)+1]] <- score_df("lag_REML", sp$name, thr, test_l$exc, p_lag)
    
    # --- NAO merge (evaluate on NAO-available subset for fair comparison) ---
    train_nao <- train %>% left_join(nao_daily, by="DATE") %>% filter(!is.na(nao))
    test_nao  <- test  %>% left_join(nao_daily, by="DATE") %>% filter(!is.na(nao))
    
    # if NAO missing in this fold, skip NAO models
    if (nrow(train_nao) > 0 && nrow(test_nao) > 0 && sum(train_nao$exc) >= 5) {
      
      # standardize NAO on TRAIN stats only
      std <- standardize_nao(train_nao, test_nao)
      train_nao <- std$train
      test_nao  <- std$test
      
      # build lag versions on NAO subset
      train_nao_l <- add_lag1_fold(train_nao)
      test_nao_l  <- add_lag1_fold(test_nao)
      
      # IMPORTANT: for comparability, also score baseline/lag on the NAO-subset
      # (so every model sees the same evaluation days)
      p_clim_nao <- rep(p_clim, nrow(test_nao))
      p_base_nao <- predict(m_base, newdata=test_nao, type="response")
      p_lag_nao  <- predict(m_lag,  newdata=test_nao_l, type="response")
      
      cv_scores[[length(cv_scores)+1]] <- score_df("climatology_NAOsubset", sp$name, thr, test_nao$exc, p_clim_nao)
      cv_scores[[length(cv_scores)+1]] <- score_df("baseline_REML_NAOsubset", sp$name, thr, test_nao$exc, p_base_nao)
      cv_scores[[length(cv_scores)+1]] <- score_df("lag_REML_NAOsubset", sp$name, thr, test_nao_l$exc, p_lag_nao)
      
      # lag+NAO
      m_lagnao <- fit_lag_nao(train_nao_l,
                              k_time=t_lagnao$final_k$k_time, k_doy=t_lagnao$final_k$k_doy,
                              bs_time="tp", bs_doy="cc", method="REML")
      p_lagnao <- predict(m_lagnao, newdata=test_nao_l, type="response")
      cv_scores[[length(cv_scores)+1]] <- score_df("lag_nao_REML", sp$name, thr, test_nao_l$exc, p_lagnao)
      
      # interaction NAO × season
      m_int <- fit_int(train_nao_l,
                       k_time=t_int$final_k$k_time, k_doy=t_int$final_k$k_doy,
                       k_ti_nao=10, bs_time="tp", method="REML")
      p_int <- predict(m_int, newdata=test_nao_l, type="response")
      cv_scores[[length(cv_scores)+1]] <- score_df("int_REML", sp$name, thr, test_nao_l$exc, p_int)
      
      # --- calibration plots (zoomed) ---
      if (SAVE_FOLD_CALIBRATION) {
        cal_dir <- file.path(OUT_DIR, paste0("calibration_folds_", sp$name))
        dir.create(cal_dir, showWarnings=FALSE)
        
        plot_calibration_zoom(test_nao$exc, p_clim_nao,
                              title=paste0(STATION_NAME, " climatology (NAO subset) — ", sp$name),
                              out_png=file.path(cal_dir, "cal_climatology_NAOsubset.png"))
        
        plot_calibration_zoom(test_nao_l$exc, p_lagnao,
                              title=paste0(STATION_NAME, " lag+NAO — ", sp$name),
                              out_png=file.path(cal_dir, "cal_lag_nao_REML.png"))
        
        plot_calibration_zoom(test_nao_l$exc, p_int,
                              title=paste0(STATION_NAME, " int (NAO×season) — ", sp$name),
                              out_png=file.path(cal_dir, "cal_int_REML.png"))
      }
      
      cv_preds[[length(cv_preds)+1]] <- data.frame(
        fold = sp$name,
        model = "lag_nao_REML",
        y = test_nao_l$exc,
        p = p_lagnao
      )
      
    } else {
      cat("  (Skipping NAO models in", sp$name, "- insufficient NAO overlap or too few events.)\n")
    }
  }
  
  cv_table <- bind_rows(cv_scores)
  
  write.csv(cv_table, file.path(OUT_DIR, "cv3fold_scores_Madrid.csv"), row.names=FALSE)
  cat("\nSaved 3-fold CV scores to:", file.path(OUT_DIR, "cv3fold_scores_Madrid.csv"), "\n")
  
  # ---- mean over folds per model (overall) ----
  cv_mean <- cv_table %>%
    group_by(model) %>%
    summarise(
      mean_logloss = mean(logloss, na.rm=TRUE),
      mean_brier   = mean(brier, na.rm=TRUE),
      mean_auc     = mean(auc, na.rm=TRUE),
      n_folds      = n_distinct(fold),
      .groups="drop"
    ) %>%
    arrange(mean_logloss, mean_brier, desc(mean_auc))
  
  write.csv(cv_mean, file.path(OUT_DIR, "cv3fold_mean_scores_Madrid.csv"), row.names=FALSE)
  cat("Saved 3-fold mean scores to:", file.path(OUT_DIR, "cv3fold_mean_scores_Madrid.csv"), "\n")
  
  print(cv_mean)
  
  preds_df <- dplyr::bind_rows(cv_preds)
  
  if (nrow(preds_df) > 0) {
    df_lagnao <- preds_df %>%
      filter(model == "lag_nao_REML")
    
    plot_calibration_zoom(
      y = df_lagnao$y,
      p = df_lagnao$p,
      title = paste0(STATION_NAME, " — lag+NAO calibration (pooled CV)"),
      out_png = file.path(OUT_DIR, "calibration_lag_nao_pooled.png")
    )
  }
} else {
  cat("\nSkipping 3-fold CV. Set RUN_3FOLD_CV <- TRUE to enable.\n")
}

# ============================================================
# 15) JOINT WEEKLY VALIDATION (zoals validatie joint 2)
# ============================================================

RUN_JOINT_WEEKLY_VALIDATION <- TRUE

if (RUN_JOINT_WEEKLY_VALIDATION) {
  JOINT_STATIONS <- c(
    madrid    = "madrid.xlsx",
    barcelona = "barcelona.xlsx",
    asturias  = "asturias.xlsx",
    sevilla   = "sevilla.xlsx",
    valencia  = "valencia.xlsx",
    alicante  = "alicante.xlsx",
    vigo      = "vigo.xlsx"
  )

  JOINT_P_EXTREME <- 0.85
  JOINT_K <- 2

  JOINT_SPLITS <- list(
    list(train_end="1990-12-31", test_start="1991-01-01", test_end="2000-12-31"),
    list(train_end="2000-12-31", test_start="2001-01-01", test_end="2010-12-31"),
    list(train_end="2010-12-31", test_start="2011-01-01", test_end="2020-12-31")
  )

  JOINT_SAVE_PLOTS <- TRUE
  JOINT_PLOT_FOLD_TO_SAVE <- 3

  JOINT_OUT_DIR <- file.path(DATA_DIR, paste0("outputs_validation_joint_weekly_p", JOINT_P_EXTREME, "_K", JOINT_K))
  dir.create(JOINT_OUT_DIR, showWarnings = FALSE)

  joint_brier <- function(y, p) mean((y - p)^2)

  joint_logloss <- function(y, p, eps=1e-15){
    p <- pmin(pmax(p, eps), 1-eps)
    -mean(y*log(p) + (1-y)*log(1-p))
  }

  joint_score_all <- function(y, p){
    out <- c(
      brier = joint_brier(y, p),
      logloss = joint_logloss(y, p),
      auc = NA_real_
    )
    if (length(unique(y)) >= 2) {
      out["auc"] <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
    }
    out
  }

  joint_calibration_bins <- function(y, p, n_bins=10){
    df <- data.frame(y=y, p=p) %>%
      filter(!is.na(y), !is.na(p))

    if (nrow(df) == 0) {
      return(data.frame(n=integer(), p_mean=numeric(), obs=numeric(), lo=numeric(), hi=numeric()))
    }

    qs <- suppressWarnings(quantile(df$p, probs=seq(0,1,length.out=n_bins+1), na.rm=TRUE))
    br <- unique(as.numeric(qs))

    if (length(br) < 3) {
      br <- unique(seq(min(df$p), max(df$p), length.out = min(n_bins+1, 6)))
    }

    if (length(br) < 2) {
      dfc <- df %>%
        summarise(
          n = n(),
          p_mean = mean(p),
          obs = mean(y),
          .groups="drop"
        ) %>%
        mutate(
          se = sqrt(obs*(1-obs)/pmax(n,1)),
          lo = pmax(obs - 1.96*se, 0),
          hi = pmin(obs + 1.96*se, 1)
        )
      return(as.data.frame(dfc))
    }

    dfc <- df %>%
      mutate(bin = cut(p, breaks = br, include.lowest = TRUE)) %>%
      group_by(bin) %>%
      summarise(
        n = n(),
        p_mean = mean(p),
        obs = mean(y),
        .groups="drop"
      ) %>%
      mutate(
        se = sqrt(obs*(1-obs)/pmax(n,1)),
        lo = pmax(obs - 1.96*se, 0),
        hi = pmin(obs + 1.96*se, 1)
      )

    as.data.frame(dfc)
  }

  joint_plot_calibration <- function(y, p, title, n_bins=10){
    df <- data.frame(y=y, p=p) %>% filter(!is.na(y), !is.na(p))
    if (nrow(df) == 0) {
      return(ggplot() + ggtitle(paste0(title, " (no data)")) + theme_minimal())
    }
    n_bins2 <- min(n_bins, max(3, floor(nrow(df)/50)))
    bins <- joint_calibration_bins(df$y, df$p, n_bins=n_bins2)

    ggplot(bins, aes(x=p_mean, y=obs)) +
      geom_abline(slope=1, intercept=0, linetype=2) +
      geom_point() +
      geom_errorbar(aes(ymin=lo, ymax=hi), width=0) +
      labs(x="Mean predicted probability", y="Observed frequency", title=title) +
      theme_minimal()
  }

  joint_get_aic_bic <- function(model){
    c(
      AIC = AIC(model),
      BIC = BIC(model),
      logLik = as.numeric(logLik(model)),
      edf = sum(model$edf),
      n = nobs(model)
    )
  }

  all_daily_joint <- bind_rows(
    lapply(names(JOINT_STATIONS), function(st){
      read_ecad_rr(JOINT_STATIONS[[st]], st)
    })
  )

  nao_weekly <- read_excel(NAO_FILE) %>%
    mutate(
      DATE = as.Date(DATE),
      iso_year = isoyear(DATE),
      iso_week = isoweek(DATE),
      nao = suppressWarnings(as.numeric(nao_index_cdas))
    ) %>%
    filter(!is.na(nao)) %>%
    group_by(iso_year, iso_week) %>%
    summarise(nao = mean(nao), .groups = "drop")

  fit_joint_base <- function(train){
    gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15),
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  fit_joint_lag <- function(train){
    gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + lag1,
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  fit_joint_nao <- function(train){
    gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + nao,
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  fit_joint_lag_nao <- function(train){
    gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + lag1 + nao,
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  fit_joint_naoSeason <- function(train){
    gam(exc_joint ~ s(time, k=20) +
          s(woy, bs="cc", k=15) +
          s(woy, by=nao, bs="cc", k=15),
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  fit_joint_lag_naoSeason <- function(train){
    gam(exc_joint ~ s(time, k=20) +
          s(woy, bs="cc", k=15) +
          lag1 +
          s(woy, by=nao, bs="cc", k=15),
        data=train, family=binomial(), method="REML",
        knots=list(woy=c(0.5, 53.5)))
  }

  make_weekly_joint <- function(daily, thresholds, joint_k){
    daily_exc <- daily %>%
      left_join(thresholds, by="station") %>%
      mutate(exc = as.integer(rr_mm > thr))

    daily_exc %>%
      mutate(
        week_id = floor_date(DATE, unit = "week", week_start = 1),
        woy = isoweek(week_id),
        iso_year = isoyear(week_id),
        iso_week = isoweek(week_id)
      ) %>%
      group_by(week_id, station, iso_year, iso_week, woy) %>%
      summarise(exc_station = as.integer(any(exc == 1)), .groups="drop") %>%
      group_by(week_id, iso_year, iso_week, woy) %>%
      summarise(
        n_extreme = sum(exc_station),
        exc_joint = as.integer(n_extreme >= joint_k),
        .groups="drop"
      ) %>%
      arrange(week_id) %>%
      mutate(
        time = row_number(),
        lag1 = dplyr::lag(exc_joint)
      ) %>%
      filter(!is.na(lag1))
  }

  run_joint_split <- function(daily, split){
    train_daily <- daily %>% filter(DATE <= as.Date(split$train_end))
    test_daily  <- daily %>% filter(DATE >= as.Date(split$test_start),
                                   DATE <= as.Date(split$test_end))

    thresholds <- train_daily %>%
      filter(rr_mm > wet_day_mm) %>%
      group_by(station) %>%
      summarise(
        thr = quantile(rr_mm, probs=JOINT_P_EXTREME, type=8, na.rm=TRUE),
        .groups="drop"
      )

    weekly_all <- make_weekly_joint(daily, thresholds, JOINT_K) %>%
      left_join(nao_weekly, by = c("iso_year","iso_week")) %>%
      arrange(week_id)

    train <- weekly_all %>% filter(week_id <= as.Date(split$train_end))
    test  <- weekly_all %>% filter(week_id >= as.Date(split$test_start),
                                   week_id <= as.Date(split$test_end))

    train_nao <- train %>% filter(!is.na(nao))
    test_nao  <- test  %>% filter(!is.na(nao))

    if (nrow(test_nao) == 0) {
      p_clim <- mean(train$exc_joint)

      m_base <- fit_joint_base(train)
      p_base <- predict(m_base, test, type="response")

      m_lag  <- fit_joint_lag(train)
      p_lag  <- predict(m_lag, test, type="response")

      s_clim <- joint_score_all(test$exc_joint, rep(p_clim, nrow(test)))
      s_base <- joint_score_all(test$exc_joint, p_base)
      s_lag  <- joint_score_all(test$exc_joint, p_lag)

      scores <- data.frame(
        train_end = split$train_end,
        test_start = split$test_start,
        test_end = split$test_end,
        thr_mean = mean(thresholds$thr, na.rm = TRUE),
        model = c("climatology","base","base+lag","base+nao","base+lag+nao","base+naoSeason","base+lag+naoSeason"),
        brier   = c(s_clim["brier"], s_base["brier"], s_lag["brier"], rep(NA_real_, 4)),
        logloss = c(s_clim["logloss"], s_base["logloss"], s_lag["logloss"], rep(NA_real_, 4)),
        auc     = c(s_clim["auc"], s_base["auc"], s_lag["auc"], rep(NA_real_, 4))
      )

      diagnostics <- data.frame(
        train_end = split$train_end,
        test_start = split$test_start,
        test_end = split$test_end,
        model = c("base","base+lag"),
        rbind(joint_get_aic_bic(m_base), joint_get_aic_bic(m_lag))
      )

      pred_long <- tibble()
      return(list(scores=scores, diagnostics=diagnostics, pred_long=pred_long))
    }

    nao_mu <- mean(train_nao$nao)
    nao_sd <- sd(train_nao$nao)
    if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
    train_nao <- train_nao %>% mutate(nao = (nao - nao_mu)/nao_sd)
    test_nao  <- test_nao  %>% mutate(nao = (nao - nao_mu)/nao_sd)

    eval <- test_nao

    p_clim <- mean(train$exc_joint)
    p_clim_eval <- rep(p_clim, nrow(eval))

    m_base <- fit_joint_base(train)
    p_base_eval <- predict(m_base, eval, type="response")

    m_lag <- fit_joint_lag(train)
    p_lag_eval <- predict(m_lag, eval, type="response")

    if (nrow(train_nao) < 50 || sum(train_nao$exc_joint) < 5) {
      m_nao <- m_naoS <- NULL
      m_lagnao <- m_lagnaoS <- NULL
      p_nao <- p_naoS <- rep(NA_real_, nrow(eval))
      p_lagnao <- p_lagnaoS <- rep(NA_real_, nrow(eval))
    } else {
      m_nao <- fit_joint_nao(train_nao)
      p_nao <- predict(m_nao, eval, type="response")

      m_lagnao <- fit_joint_lag_nao(train_nao)
      p_lagnao <- predict(m_lagnao, eval, type="response")

      m_naoS <- fit_joint_naoSeason(train_nao)
      p_naoS <- predict(m_naoS, eval, type="response")

      m_lagnaoS <- fit_joint_lag_naoSeason(train_nao)
      p_lagnaoS <- predict(m_lagnaoS, eval, type="response")
    }

    s_clim <- joint_score_all(eval$exc_joint, p_clim_eval)
    s_base <- joint_score_all(eval$exc_joint, p_base_eval)
    s_lag  <- joint_score_all(eval$exc_joint, p_lag_eval)
    s_nao  <- joint_score_all(eval$exc_joint, p_nao)
    s_lagnao <- joint_score_all(eval$exc_joint, p_lagnao)
    s_naoS <- joint_score_all(eval$exc_joint, p_naoS)
    s_lagnaoS <- joint_score_all(eval$exc_joint, p_lagnaoS)

    scores <- data.frame(
      train_end = split$train_end,
      test_start = split$test_start,
      test_end = split$test_end,
      thr_mean = mean(thresholds$thr, na.rm = TRUE),
      model = c("climatology","base","base+lag","base+nao","base+lag+nao","base+naoSeason","base+lag+naoSeason"),
      brier   = c(s_clim["brier"], s_base["brier"], s_lag["brier"], s_nao["brier"], s_lagnao["brier"], s_naoS["brier"], s_lagnaoS["brier"]),
      logloss = c(s_clim["logloss"], s_base["logloss"], s_lag["logloss"], s_nao["logloss"], s_lagnao["logloss"], s_naoS["logloss"], s_lagnaoS["logloss"]),
      auc     = c(s_clim["auc"], s_base["auc"], s_lag["auc"], s_nao["auc"], s_lagnao["auc"], s_naoS["auc"], s_lagnaoS["auc"])
    )

    diag_list <- list(
      `base` = m_base,
      `base+lag` = m_lag
    )
    if (!is.null(m_nao)) {
      diag_list <- c(diag_list, list(
        `base+nao` = m_nao,
        `base+naoSeason` = m_naoS,
        `base+lag+nao` = m_lagnao,
        `base+lag+naoSeason` = m_lagnaoS
      ))
    }

    diagnostics <- bind_rows(lapply(names(diag_list), function(mname){
      mod <- diag_list[[mname]]
      out <- as.data.frame(as.list(joint_get_aic_bic(mod)))
      out$model <- mname
      out
    })) %>%
      mutate(
        train_end = split$train_end,
        test_start = split$test_start,
        test_end = split$test_end
      ) %>%
      select(train_end, test_start, test_end, model, AIC, BIC, logLik, edf, n)

    pred_long <- bind_rows(
      tibble(fold=NA_integer_, model="climatology", y=eval$exc_joint, p=p_clim_eval),
      tibble(fold=NA_integer_, model="base", y=eval$exc_joint, p=p_base_eval),
      tibble(fold=NA_integer_, model="base+nao", y=eval$exc_joint, p=p_nao),
      tibble(fold=NA_integer_, model="base+naoSeason", y=eval$exc_joint, p=p_naoS),
      tibble(fold=NA_integer_, model="base+lag", y=eval$exc_joint, p=p_lag_eval),
      tibble(fold=NA_integer_, model="base+lag+nao", y=eval$exc_joint, p=p_lagnao),
      tibble(fold=NA_integer_, model="base+lag+naoSeason", y=eval$exc_joint, p=p_lagnaoS)
    ) %>% filter(!is.na(p))

    list(scores=scores, diagnostics=diagnostics, pred_long=pred_long)
  }

  joint_scores <- list()
  joint_diagnostics <- list()
  joint_preds <- list()

  for (i in seq_along(JOINT_SPLITS)) {
    sp <- JOINT_SPLITS[[i]]

    res <- run_joint_split(all_daily_joint, sp)

    scores_i <- res$scores %>% mutate(fold = i)
    diagnostics_i <- res$diagnostics %>% mutate(fold = i)

    print(scores_i)

    joint_scores[[length(joint_scores) + 1]] <- scores_i
    joint_diagnostics[[length(joint_diagnostics) + 1]] <- diagnostics_i

    if (nrow(res$pred_long) > 0) {
      joint_preds[[length(joint_preds) + 1]] <- res$pred_long %>% mutate(fold = i)
    }

    if (JOINT_SAVE_PLOTS && i == JOINT_PLOT_FOLD_TO_SAVE && nrow(res$pred_long) > 0) {
      fold_dir <- file.path(JOINT_OUT_DIR, paste0("calibration_fold", i))
      dir.create(fold_dir, showWarnings = FALSE)

      fold_preds <- res$pred_long %>% mutate(fold=i)

      for (mname in sort(unique(fold_preds$model))) {
        dfm <- fold_preds %>% filter(model == mname)
        if (nrow(dfm) == 0) next

        g <- joint_plot_calibration(
          y = dfm$y,
          p = dfm$p,
          title = paste0("Joint weekly ", mname, " calibration (fold ", i, ")")
        )

        ggsave(
          filename = file.path(fold_dir, paste0("cal_joint_fold", i, "_", gsub("\\+","_", mname), ".png")),
          plot = g, width = 6, height = 4, dpi = 200
        )
      }
    }
  }

  joint_scores_df <- bind_rows(joint_scores)
  joint_diagnostics_df <- bind_rows(joint_diagnostics)
  joint_preds_df <- bind_rows(joint_preds)

  write.csv(joint_scores_df,       file.path(JOINT_OUT_DIR, "scores_joint_weekly.csv"), row.names = FALSE)
  write.csv(joint_diagnostics_df,  file.path(JOINT_OUT_DIR, "diagnostics_aic_bic_joint_weekly.csv"), row.names = FALSE)

  cat("\nDONE.\nSaved joint weekly outputs in:", JOINT_OUT_DIR, "\n")

  if (JOINT_SAVE_PLOTS && exists("joint_preds_df") && nrow(joint_preds_df) > 0) {
    pooled_dir <- file.path(JOINT_OUT_DIR, "calibration_pooled_allfolds")
    dir.create(pooled_dir, showWarnings = FALSE)

    for (mname in sort(unique(joint_preds_df$model))) {
      dfm <- joint_preds_df %>% filter(model == mname)
      if (nrow(dfm) == 0) next

      g <- joint_plot_calibration(
        y = dfm$y,
        p = dfm$p,
        title = paste0("Joint weekly ", mname, " calibration (POOLED folds)")
      )

      ggsave(
        filename = file.path(pooled_dir, paste0("cal_pooled_joint_", gsub("\\+","_", mname), ".png")),
        plot = g, width = 6, height = 4, dpi = 200
      )
    }

    cat("\nSaved pooled calibration plots in:\n", pooled_dir, "\n")
  }

  joint_summary_mean <- joint_scores_df %>%
    group_by(model) %>%
    summarise(
      mean_brier   = mean(brier, na.rm = TRUE),
      mean_logloss = mean(logloss, na.rm = TRUE),
      mean_auc     = mean(auc, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(joint_summary_mean,
            file = file.path(JOINT_OUT_DIR, "summary_mean_scores_joint_weekly.csv"),
            row.names = FALSE)

  joint_ranked_overall <- joint_summary_mean %>%
    mutate(
      rank_logloss = rank(mean_logloss, ties.method = "average"),
      rank_brier   = rank(mean_brier,   ties.method = "average"),
      rank_auc     = rank(-mean_auc,    ties.method = "average"),
      rank_total   = (rank_logloss + rank_brier + rank_auc) / 3
    ) %>%
    arrange(rank_total)

  write.csv(joint_ranked_overall,
            file.path(JOINT_OUT_DIR, "ranking_overall_models_joint_weekly.csv"),
            row.names = FALSE)

  if (nrow(joint_diagnostics_df) > 0) {
    joint_diag_summary <- joint_diagnostics_df %>%
      group_by(model) %>%
      summarise(
        mean_AIC = mean(AIC, na.rm=TRUE),
        mean_BIC = mean(BIC, na.rm=TRUE),
        mean_logLik = mean(logLik, na.rm=TRUE),
        mean_edf = mean(edf, na.rm=TRUE),
        .groups="drop"
      ) %>%
      arrange(mean_AIC)

    write.csv(joint_diag_summary,
              file.path(JOINT_OUT_DIR, "diagnostics_mean_aic_bic_joint_weekly.csv"),
              row.names = FALSE)

    cat("\nSaved diagnostics summary to:\n",
        file.path(JOINT_OUT_DIR, "diagnostics_mean_aic_bic_joint_weekly.csv"),
        "\n")
  }
}
