# ============================================================
# JOINT WEEKLY GAM — Validation & Calibration (incl. NAO models)
# Models:
#  - climatology
#  - base: s(time) + s(woy)
#  - base+lag: + lag1
#  - base+nao: + nao
#  - base+lag+nao: + lag1 + nao
#  - base+naoSeason: + s(woy, by=nao)   [NAO x season]
#  - base+lag+naoSeason: + lag1 + s(woy, by=nao)
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(pROC)
library(ggplot2)

DATA_DIR <- "C:/Users/Sara/Downloads/Case study"
setwd(DATA_DIR)

stations <- c(
  madrid    = "madrid.xlsx",
  barcelona = "barcelona.xlsx",
  asturias  = "asturias.xlsx",
  sevilla   = "sevilla.xlsx",
  valencia  = "valencia.xlsx",
  alicante  = "alicante.xlsx",
  vigo      = "vigo.xlsx"
)

NAO_FILE <- "NAO_cleaned.xlsx"  

wet_day_mm <- 1
p_extreme  <- 0.85
K_joint    <- 2   # stations extreme in a week

splits <- list(
  list(train_end="1990-12-31", test_start="1991-01-01", test_end="2000-12-31"),
  list(train_end="2000-12-31", test_start="2001-01-01", test_end="2010-12-31"),
  list(train_end="2010-12-31", test_start="2011-01-01", test_end="2020-12-31")
)

read_ecad_rr <- function(path, station_name,
                         date_min = as.Date("1970-01-01"),
                         date_max = as.Date("2020-12-31")){
  read_excel(path, col_names = FALSE) %>%
    select(1) %>%
    rename(raw = 1) %>%
    separate(raw, into=c("STAID","SOUID","DATE","RR","Q_RR"), sep=",", convert=TRUE) %>%
    mutate(
      station = station_name,
      DATE = as.Date(as.character(DATE), "%Y%m%d"),
      RR   = suppressWarnings(as.numeric(RR)),
      Q_RR = suppressWarnings(as.numeric(Q_RR)),
      rr_mm = RR / 10
    ) %>%
    filter(!is.na(DATE), !is.na(RR), !is.na(Q_RR), RR != -9999, Q_RR == 0) %>%
    filter(DATE >= date_min, DATE <= date_max) %>%
    select(DATE, station, rr_mm)
}

brier <- function(y, p) mean((y - p)^2)
logloss <- function(y, p, eps=1e-15){
  p <- pmin(pmax(p, eps), 1-eps)
  -mean(y*log(p) + (1-y)*log(1-p))
}

safe_auc <- function(y, p){
  # If test has only one class, AUC is undefined
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
}

# Calibration bins (for plotting)
calibration_bins <- function(y, p, n_bins=10){
  dfc <- data.frame(y=y, p=p) %>%
    mutate(bin = cut(
      p,
      breaks = quantile(p, probs=seq(0,1,length.out=n_bins+1), na.rm=TRUE),
      include.lowest=TRUE
    )) %>%
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
  dfc
}

plot_calibration <- function(y, p, title, n_bins=10){
  bins <- calibration_bins(y, p, n_bins=n_bins)
  ggplot(bins, aes(x=p_mean, y=obs)) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    geom_point() +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=0) +
    labs(x="Mean predicted probability", y="Observed frequency", title=title) +
    theme_minimal()
}

# read combine stations
all_daily <- bind_rows(
  lapply(names(stations), function(st){
    read_ecad_rr(stations[[st]], st)
  })
)

# thresholds
thresholds <- all_daily %>%
  filter(rr_mm > wet_day_mm) %>%
  group_by(station) %>%
  summarise(
    thr = quantile(rr_mm, probs=p_extreme, type=8, na.rm=TRUE),
    .groups="drop"
  )

daily_exc <- all_daily %>%
  left_join(thresholds, by="station") %>%
  mutate(exc = as.integer(rr_mm > thr))

# weekly joint panel
weekly <- daily_exc %>%
  mutate(
    week_id = floor_date(DATE, unit = "week", week_start = 1),
    woy = isoweek(week_id)
  ) %>%
  group_by(week_id, station) %>%
  summarise(exc_station = as.integer(any(exc == 1)), .groups="drop") %>%
  group_by(week_id) %>%
  summarise(
    n_extreme = sum(exc_station),
    exc_joint = as.integer(n_extreme >= K_joint),
    .groups="drop"
  ) %>%
  arrange(week_id) %>%
  mutate(
    time = row_number(),
    lag1 = dplyr::lag(exc_joint),
    woy  = isoweek(week_id),
    iso_year = isoyear(week_id),
    iso_week = isoweek(week_id)
  ) %>%
  filter(!is.na(lag1))

# nao
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

weekly <- weekly %>%
  select(-any_of(c("nao"))) %>%
  left_join(nao_weekly, by = c("iso_year","iso_week")) %>%
  arrange(week_id)

# modellen fitten
fit_base <- function(train){
  gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15),
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

fit_lag <- function(train){
  gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + lag1,
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

fit_nao <- function(train){
  gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + nao,
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

fit_lag_nao <- function(train){
  gam(exc_joint ~ s(time, k=20) + s(woy, bs="cc", k=15) + lag1 + nao,
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

# NAO x season
fit_naoSeason <- function(train){
  gam(exc_joint ~ s(time, k=20) +
        s(woy, bs="cc", k=15) +
        s(woy, by=nao, bs="cc", k=15),
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

fit_lag_naoSeason <- function(train){
  gam(exc_joint ~ s(time, k=20) +
        s(woy, bs="cc", k=15) +
        lag1 +
        s(woy, by=nao, bs="cc", k=15),
      data=train, family=binomial(), method="REML",
      knots=list(woy=c(0.5, 53.5)))
}

# validation loop
results <- list()

pred_store <- list()

for (i in seq_along(splits)) {
  sp <- splits[[i]]
  
  train <- weekly %>% filter(week_id <= as.Date(sp$train_end))
  test  <- weekly %>% filter(week_id >= as.Date(sp$test_start),
                             week_id <= as.Date(sp$test_end))
  
  train_nao <- train %>% filter(!is.na(nao))
  test_nao  <- test  %>% filter(!is.na(nao))

  nao_mu <- mean(train_nao$nao)
  nao_sd <- sd(train_nao$nao)
  if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
  
  train_nao <- train_nao %>% mutate(nao = (nao - nao_mu) / nao_sd)
  test_nao  <- test_nao  %>% mutate(nao = (nao - nao_mu) / nao_sd)
  
  # Baseline climatology (computed on train)
  p_clim <- mean(train$exc_joint)
  
  cat("\nFold", i,
      "train n =", nrow(train), "events =", sum(train$exc_joint), "rate =", mean(train$exc_joint),
      "| train_nao n =", nrow(train_nao), "\n")
  
  # Fit base models (no NAO needed)
  m_base <- fit_base(train)
  p_base <- predict(m_base, test, type="response")
  
  m_lag  <- fit_lag(train)
  p_lag  <- predict(m_lag, test, type="response")
  
  # Fit NAO models on NAO-complete subset
  # If too few rows, return NA scores for NAO models
  if (nrow(train_nao) < 50 || sum(train_nao$exc_joint) < 5) {
    p_nao <- p_lagnao <- p_naoS <- p_lagnaoS <- rep(NA_real_, nrow(test))
  } else {
    m_nao <- fit_nao(train_nao)
    p_nao <- predict(m_nao, test_nao, type="response")
    
    m_lagnao <- fit_lag_nao(train_nao)
    p_lagnao <- predict(m_lagnao, test_nao, type="response")
    
    m_naoS <- fit_naoSeason(train_nao)
    p_naoS <- predict(m_naoS, test_nao, type="response")
    
    m_lagnaoS <- fit_lag_naoSeason(train_nao)
    p_lagnaoS <- predict(m_lagnaoS, test_nao, type="response")
  }
  
  # Choose evaluation set:
  eval <- test_nao
  if (nrow(eval) == 0) stop("No NAO overlap in test fold; check NAO date coverage.")
  
  # Predictions
  p_clim_eval <- rep(p_clim, nrow(eval))
  p_base_eval <- predict(m_base, eval, type="response")
  p_lag_eval  <- predict(m_lag,  eval, type="response")
  
  pred_store[[i]] <- tibble(
    fold = i,
    y = eval$exc_joint,
    base = p_base_eval,
    `base+lag` = p_lag_eval,
    `base+nao` = p_nao,
    `base+lag+nao` = p_lagnao,
    `base+naoSeason` = p_naoS,
    `base+lag+naoSeason` = p_lagnaoS
  )
  
  scores <- data.frame(
    fold = i,
    model = c("climatology","base","base+lag","base+nao","base+lag+nao","base+naoSeason","base+lag+naoSeason"),
    brier = c(
      brier(eval$exc_joint, p_clim_eval),
      brier(eval$exc_joint, p_base_eval),
      brier(eval$exc_joint, p_lag_eval),
      brier(eval$exc_joint, p_nao),
      brier(eval$exc_joint, p_lagnao),
      brier(eval$exc_joint, p_naoS),
      brier(eval$exc_joint, p_lagnaoS)
    ),
    logloss = c(
      logloss(eval$exc_joint, p_clim_eval),
      logloss(eval$exc_joint, p_base_eval),
      logloss(eval$exc_joint, p_lag_eval),
      logloss(eval$exc_joint, p_nao),
      logloss(eval$exc_joint, p_lagnao),
      logloss(eval$exc_joint, p_naoS),
      logloss(eval$exc_joint, p_lagnaoS)
    ),
    auc = c(
      0.5,
      safe_auc(eval$exc_joint, p_base_eval),
      safe_auc(eval$exc_joint, p_lag_eval),
      safe_auc(eval$exc_joint, p_nao),
      safe_auc(eval$exc_joint, p_lagnao),
      safe_auc(eval$exc_joint, p_naoS),
      safe_auc(eval$exc_joint, p_lagnaoS)
    )
  )
  
  print(scores)
  results[[i]] <- scores
}

results_df <- bind_rows(results)
print(results_df)

# calibration plot joint
pred_all <- bind_rows(pred_store)

OUT_DIR_POOL <- file.path(DATA_DIR, paste0("calibration_pooled_K", K_joint))
dir.create(OUT_DIR_POOL, showWarnings = FALSE)

model_cols <- setdiff(names(pred_all), c("fold","y"))

for (mname in model_cols) {
  g <- plot_calibration(
    y = pred_all$y,
    p = pred_all[[mname]],
    title = paste0("Pooled calibration (K=", K_joint, ") — ", mname)
  )
  
  fname <- file.path(
    OUT_DIR_POOL,
    paste0("cal_pooled_K", K_joint, "_", gsub("\\+","_", mname), ".png")
  )
  
  ggsave(fname, plot=g, width=6, height=4, dpi=200)
}

cat("\nSaved pooled calibration plots to:\n", OUT_DIR_POOL, "\n")


# calibration plot fold 3
fold_to_plot <- 3
model_to_plot <- "base+lag+naoSeason"  # choose from scores$model

sp <- splits[[fold_to_plot]]
train <- weekly %>% filter(week_id <= as.Date(sp$train_end))
test  <- weekly %>% filter(week_id >= as.Date(sp$test_start),
                           week_id <= as.Date(sp$test_end))

train_nao <- train %>% filter(!is.na(nao))
test_nao  <- test  %>% filter(!is.na(nao))

nao_mu <- mean(train_nao$nao); nao_sd <- sd(train_nao$nao)
if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
train_nao <- train_nao %>% mutate(nao = (nao - nao_mu) / nao_sd)
test_nao  <- test_nao  %>% mutate(nao = (nao - nao_mu) / nao_sd)

# Fit models needed for plotting
m_base <- fit_base(train);           p_base <- predict(m_base, test_nao, type="response")
m_lag  <- fit_lag(train);            p_lag  <- predict(m_lag,  test_nao, type="response")

# NAO models
m_nao  <- fit_nao(train_nao);        p_nao  <- predict(m_nao, test_nao, type="response")
m_lagnao <- fit_lag_nao(train_nao);  p_lagnao <- predict(m_lagnao, test_nao, type="response")
m_naoS <- fit_naoSeason(train_nao);  p_naoS <- predict(m_naoS, test_nao, type="response")
m_lagnaoS <- fit_lag_naoSeason(train_nao); p_lagnaoS <- predict(m_lagnaoS, test_nao, type="response")

pred_map <- list(
  "base" = p_base,
  "base+lag" = p_lag,
  "base+nao" = p_nao,
  "base+lag+nao" = p_lagnao,
  "base+naoSeason" = p_naoS,
  "base+lag+naoSeason" = p_lagnaoS
)

p_plot <- pred_map[[model_to_plot]]
g <- plot_calibration(test_nao$exc_joint, p_plot,
                      title=paste0("Joint weekly calibration (K=", K_joint,
                                   "), fold ", fold_to_plot, " — ", model_to_plot))

# joint calibration plots
OUT_DIR <- file.path(DATA_DIR, paste0("calibration_plots_K", K_joint, "_fold", fold_to_plot))
dir.create(OUT_DIR, showWarnings = FALSE)

for (mname in names(pred_map)) {
  p_plot <- pred_map[[mname]]
  
  g <- plot_calibration(
    test_nao$exc_joint,
    p_plot,
    title = paste0("Joint weekly calibration (K=", K_joint,
                   "), fold ", fold_to_plot, " — ", mname)
  )
  
  fname <- file.path(
    OUT_DIR,
    paste0("cal_joint_K", K_joint,
           "_fold", fold_to_plot,
           "_", gsub("\\+","_", mname),
           ".png")
  )
  
  ggsave(fname, plot = g, width = 6, height = 4, dpi = 200)
}

cat("\nSaved all calibration plots to folder:\n", OUT_DIR, "\n")

# brier logloss etc ranken
ranked <- results_df %>%
  group_by(fold) %>%
  mutate(
    r_logloss = rank(logloss, ties.method = "average"),          # lager beter
    r_brier   = rank(brier,   ties.method = "average"),          # lager beter
    r_auc     = rank(-auc,    ties.method = "average")           # hoger beter -> min rank
  ) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(
    mean_logloss = mean(logloss, na.rm=TRUE),
    mean_brier   = mean(brier,   na.rm=TRUE),
    mean_auc     = mean(auc,     na.rm=TRUE),
    mean_rank_logloss = mean(r_logloss, na.rm=TRUE),
    mean_rank_brier   = mean(r_brier,   na.rm=TRUE),
    mean_rank_auc     = mean(r_auc,     na.rm=TRUE),
    # samengestelde rank (gewichten kun je aanpassen)
    rank_total = 0.5*mean_rank_logloss + 0.3*mean_rank_brier + 0.2*mean_rank_auc,
    .groups="drop"
  ) %>%
  arrange(rank_total)

print(ranked)

# ranking aic bic etc
aic_bic <- function(model){
  c(
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    edf = sum(model$edf)
  )
}

models_train <- list(
  base = m_base,
  base_lag = m_lag,
  base_nao = m_nao,
  base_lag_nao = m_lagnao,
  base_naoSeason = m_naoS,
  base_lag_naoSeason = m_lagnaoS
)

aic_table <- do.call(
  rbind,
  lapply(models_train, aic_bic)
)

aic_table <- as.data.frame(aic_table)
aic_table$model <- rownames(aic_table)
aic_table

