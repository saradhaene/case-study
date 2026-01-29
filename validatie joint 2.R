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

SAVE_PLOTS <- TRUE
PLOT_FOLD_TO_SAVE <- 3

OUT_DIR <- file.path(DATA_DIR, paste0("outputs_validation_joint_weekly_p", p_extreme, "_K", K_joint))
dir.create(OUT_DIR, showWarnings = FALSE)

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

score_all <- function(y, p){
  out <- c(
    brier = brier(y, p),
    logloss = logloss(y, p),
    auc = NA_real_
  )
  if (length(unique(y)) >= 2) {
    out["auc"] <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
  }
  out
}

# calibration bins (robust to few unique probs)
calibration_bins <- function(y, p, n_bins=10){
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

plot_calibration <- function(y, p, title, n_bins=10){
  df <- data.frame(y=y, p=p) %>% filter(!is.na(y), !is.na(p))
  if (nrow(df) == 0) {
    return(ggplot() + ggtitle(paste0(title, " (no data)")) + theme_minimal())
  }
  n_bins2 <- min(n_bins, max(3, floor(nrow(df)/50)))
  bins <- calibration_bins(df$y, df$p, n_bins=n_bins2)

  ggplot(bins, aes(x=p_mean, y=obs)) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    geom_point() +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=0) +
    labs(x="Mean predicted probability", y="Observed frequency", title=title) +
    theme_minimal()
}

# AIC/BIC helper
get_aic_bic <- function(model){
  c(
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    edf = sum(model$edf),
    n = nobs(model)
  )
}

# read combine stations
all_daily <- bind_rows(
  lapply(names(stations), function(st){
    read_ecad_rr(stations[[st]], st)
  })
)

# NAO weekly aggregation
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

make_weekly_joint <- function(daily, thresholds, K_joint){
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
      exc_joint = as.integer(n_extreme >= K_joint),
      .groups="drop"
    ) %>%
    arrange(week_id) %>%
    mutate(
      time = row_number(),
      lag1 = dplyr::lag(exc_joint)
    ) %>%
    filter(!is.na(lag1))
}

run_one_split <- function(daily, split){
  train_daily <- daily %>% filter(DATE <= as.Date(split$train_end))
  test_daily  <- daily %>% filter(DATE >= as.Date(split$test_start),
                                 DATE <= as.Date(split$test_end))

  thresholds <- train_daily %>%
    filter(rr_mm > wet_day_mm) %>%
    group_by(station) %>%
    summarise(
      thr = quantile(rr_mm, probs=p_extreme, type=8, na.rm=TRUE),
      .groups="drop"
    )

  weekly_all <- make_weekly_joint(daily, thresholds, K_joint) %>%
    left_join(nao_weekly, by = c("iso_year","iso_week")) %>%
    arrange(week_id)

  train <- weekly_all %>% filter(week_id <= as.Date(split$train_end))
  test  <- weekly_all %>% filter(week_id >= as.Date(split$test_start),
                                 week_id <= as.Date(split$test_end))

  # NAO subsets
  train_nao <- train %>% filter(!is.na(nao))
  test_nao  <- test  %>% filter(!is.na(nao))

  # If no NAO overlap in test, return base/lag scores and NA for NAO models
  if (nrow(test_nao) == 0) {
    p_clim <- mean(train$exc_joint)

    m_base <- fit_base(train)
    p_base <- predict(m_base, test, type="response")

    m_lag  <- fit_lag(train)
    p_lag  <- predict(m_lag, test, type="response")

    s_clim <- score_all(test$exc_joint, rep(p_clim, nrow(test)))
    s_base <- score_all(test$exc_joint, p_base)
    s_lag  <- score_all(test$exc_joint, p_lag)

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
      rbind(get_aic_bic(m_base), get_aic_bic(m_lag))
    )

    pred_long <- tibble()
    return(list(scores=scores, diagnostics=diagnostics, pred_long=pred_long))
  }

  # Standardize NAO on train stats
  nao_mu <- mean(train_nao$nao)
  nao_sd <- sd(train_nao$nao)
  if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
  train_nao <- train_nao %>% mutate(nao = (nao - nao_mu)/nao_sd)
  test_nao  <- test_nao  %>% mutate(nao = (nao - nao_mu)/nao_sd)

  eval <- test_nao

  # climatology baseline
  p_clim <- mean(train$exc_joint)
  p_clim_eval <- rep(p_clim, nrow(eval))

  # base + lag
  m_base <- fit_base(train)
  p_base_eval <- predict(m_base, eval, type="response")

  m_lag <- fit_lag(train)
  p_lag_eval <- predict(m_lag, eval, type="response")

  # NAO models on NAO-complete subset
  if (nrow(train_nao) < 50 || sum(train_nao$exc_joint) < 5) {
    m_nao <- m_naoS <- NULL
    m_lagnao <- m_lagnaoS <- NULL
    p_nao <- p_naoS <- rep(NA_real_, nrow(eval))
    p_lagnao <- p_lagnaoS <- rep(NA_real_, nrow(eval))
  } else {
    m_nao <- fit_nao(train_nao)
    p_nao <- predict(m_nao, eval, type="response")

    m_lagnao <- fit_lag_nao(train_nao)
    p_lagnao <- predict(m_lagnao, eval, type="response")

    m_naoS <- fit_naoSeason(train_nao)
    p_naoS <- predict(m_naoS, eval, type="response")

    m_lagnaoS <- fit_lag_naoSeason(train_nao)
    p_lagnaoS <- predict(m_lagnaoS, eval, type="response")
  }

  # Scores
  s_clim <- score_all(eval$exc_joint, p_clim_eval)
  s_base <- score_all(eval$exc_joint, p_base_eval)
  s_lag  <- score_all(eval$exc_joint, p_lag_eval)
  s_nao  <- score_all(eval$exc_joint, p_nao)
  s_lagnao <- score_all(eval$exc_joint, p_lagnao)
  s_naoS <- score_all(eval$exc_joint, p_naoS)
  s_lagnaoS <- score_all(eval$exc_joint, p_lagnaoS)

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

  # Diagnostics
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
    out <- as.data.frame(as.list(get_aic_bic(mod)))
    out$model <- mname
    out
  })) %>%
    mutate(
      train_end = split$train_end,
      test_start = split$test_start,
      test_end = split$test_end
    ) %>%
    select(train_end, test_start, test_end, model, AIC, BIC, logLik, edf, n)

  # Long predictions for pooled calibration
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

all_scores <- list()
all_diagnostics <- list()
all_preds <- list()

for (i in seq_along(splits)) {
  sp <- splits[[i]]

  res <- run_one_split(all_daily, sp)

  scores_i <- res$scores %>% mutate(fold = i)
  diagnostics_i <- res$diagnostics %>% mutate(fold = i)

  print(scores_i)

  all_scores[[length(all_scores) + 1]] <- scores_i
  all_diagnostics[[length(all_diagnostics) + 1]] <- diagnostics_i

  if (nrow(res$pred_long) > 0) {
    all_preds[[length(all_preds) + 1]] <- res$pred_long %>% mutate(fold = i)
  }

  # Save fold-specific calibration plots
  if (SAVE_PLOTS && i == PLOT_FOLD_TO_SAVE && nrow(res$pred_long) > 0) {
    fold_dir <- file.path(OUT_DIR, paste0("calibration_fold", i))
    dir.create(fold_dir, showWarnings = FALSE)

    fold_preds <- res$pred_long %>% mutate(fold=i)

    for (mname in sort(unique(fold_preds$model))) {
      dfm <- fold_preds %>% filter(model == mname)
      if (nrow(dfm) == 0) next

      g <- plot_calibration(
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

scores_df <- bind_rows(all_scores)
diagnostics_df <- bind_rows(all_diagnostics)
preds_df <- bind_rows(all_preds)

# Save tables
write.csv(scores_df,       file.path(OUT_DIR, "scores_joint_weekly.csv"), row.names = FALSE)
write.csv(diagnostics_df,  file.path(OUT_DIR, "diagnostics_aic_bic_joint_weekly.csv"), row.names = FALSE)

cat("\nDONE.\nSaved outputs in:", OUT_DIR, "\n")

# pooled calibration plots
if (SAVE_PLOTS && exists("preds_df") && nrow(preds_df) > 0) {
  pooled_dir <- file.path(OUT_DIR, "calibration_pooled_allfolds")
  dir.create(pooled_dir, showWarnings = FALSE)

  for (mname in sort(unique(preds_df$model))) {
    dfm <- preds_df %>% filter(model == mname)
    if (nrow(dfm) == 0) next

    g <- plot_calibration(
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

# mean scores per model
summary_mean <- scores_df %>%
  group_by(model) %>%
  summarise(
    mean_brier   = mean(brier, na.rm = TRUE),
    mean_logloss = mean(logloss, na.rm = TRUE),
    mean_auc     = mean(auc, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_mean,
          file = file.path(OUT_DIR, "summary_mean_scores_joint_weekly.csv"),
          row.names = FALSE)

# ranking
ranked_overall <- summary_mean %>%
  mutate(
    rank_logloss = rank(mean_logloss, ties.method = "average"),
    rank_brier   = rank(mean_brier,   ties.method = "average"),
    rank_auc     = rank(-mean_auc,    ties.method = "average"),
    rank_total   = (rank_logloss + rank_brier + rank_auc) / 3
  ) %>%
  arrange(rank_total)

write.csv(ranked_overall,
          file.path(OUT_DIR, "ranking_overall_models_joint_weekly.csv"),
          row.names = FALSE)

# diagnostics summary (AIC/BIC) averaged per model
if (nrow(diagnostics_df) > 0) {
  diag_summary <- diagnostics_df %>%
    group_by(model) %>%
    summarise(
      mean_AIC = mean(AIC, na.rm=TRUE),
      mean_BIC = mean(BIC, na.rm=TRUE),
      mean_logLik = mean(logLik, na.rm=TRUE),
      mean_edf = mean(edf, na.rm=TRUE),
      .groups="drop"
    ) %>%
    arrange(mean_AIC)

  write.csv(diag_summary,
            file.path(OUT_DIR, "diagnostics_mean_aic_bic_joint_weekly.csv"),
            row.names = FALSE)

  cat("\nSaved diagnostics summary to:\n",
      file.path(OUT_DIR, "diagnostics_mean_aic_bic_joint_weekly.csv"),
      "\n")
}
