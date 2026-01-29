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
  vigo      = "vigo.xlsx",
  valencia  = "valencia.xlsx",
  alicante  = "alicante.xlsx"
)

NAO_FILE <- "NAO_cleaned.xlsx"

wet_day_mm <- 1
p_extreme  <- 0.95

# fodls
splits <- list(
  list(train_end="1990-12-31", test_start="1991-01-01", test_end="2000-12-31"),
  list(train_end="2000-12-31", test_start="2001-01-01", test_end="2010-12-31"),
  list(train_end="2010-12-31", test_start="2011-01-01", test_end="2020-12-31")
)

# Output
SAVE_PLOTS <- TRUE
PLOT_FOLD_TO_SAVE <- 3

OUT_DIR <- file.path(DATA_DIR, paste0("outputs_validation_daily_p", p_extreme))
dir.create(OUT_DIR, showWarnings = FALSE)

read_ecad_rr <- function(path, station_name,
                         date_min = as.Date("1970-01-01"),
                         date_max = as.Date("2020-12-31")) {
  
  read_excel(path, col_names = FALSE) %>%
    select(1) %>%
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
    filter(!is.na(DATE), !is.na(RR), !is.na(Q_RR), RR != -9999, Q_RR == 0) %>%
    filter(DATE >= date_min, DATE <= date_max) %>%
    arrange(DATE) %>%
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

add_lag1 <- function(df){
  df %>%
    arrange(DATE) %>%
    mutate(lag1 = dplyr::lag(exc, 1)) %>%
    filter(!is.na(lag1))
}

# aic bic
aic_bic <- function(model){
  c(
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    edf = sum(model$edf),
    n = nobs(model)
  )
}

# Join NAO + add time/doy
add_features <- function(df, nao_daily){
  df %>%
    arrange(DATE) %>%
    left_join(nao_daily, by = "DATE") %>%
    mutate(
      time = as.numeric(DATE - min(DATE)),
      doy  = yday(DATE)
    )
}

make_train_test <- function(df, train_end, test_start, test_end=NULL,
                            wet_day_mm=1, p_extreme=0.90){
  
  train <- df %>% filter(DATE <= as.Date(train_end))
  test  <- df %>% filter(DATE >= as.Date(test_start))
  if (!is.null(test_end)) test <- test %>% filter(DATE <= as.Date(test_end))
  
  thr <- quantile(train$rr_mm[train$rr_mm > wet_day_mm],
                  probs = p_extreme, type = 8, na.rm = TRUE)
  
  train <- train %>% mutate(exc = as.integer(rr_mm > thr))
  test  <- test  %>% mutate(exc = as.integer(rr_mm > thr))
  
  list(train=train, test=test, thr=as.numeric(thr))
}

# modellen fitten
fit_base <- function(train){
  gam(exc ~ s(time, k=30) + s(doy, bs="cc", k=20),
      data=train, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

fit_lag <- function(train_l){
  gam(exc ~ s(time, k=30) + s(doy, bs="cc", k=20) + lag1,
      data=train_l, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

fit_nao <- function(train){
  gam(exc ~ s(time, k=30) + s(doy, bs="cc", k=20) + nao,
      data=train, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

fit_lag_nao <- function(train_l){
  gam(exc ~ s(time, k=30) + s(doy, bs="cc", k=20) + lag1 + nao,
      data=train_l, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

# NAO x season: varying coefficient over day-of-year
fit_naoSeason <- function(train){
  gam(exc ~ s(time, k=30) +
        s(doy, bs="cc", k=20) +
        s(doy, by=nao, bs="cc", k=20),
      data=train, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

fit_lag_naoSeason <- function(train_l){
  gam(exc ~ s(time, k=30) +
        s(doy, bs="cc", k=20) +
        lag1 +
        s(doy, by=nao, bs="cc", k=20),
      data=train_l, family=binomial(), method="REML",
      knots=list(doy=c(0.5, 366.5)))
}

# calibration plots
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

run_one_split <- function(df, station_name, train_end, test_start, test_end=NULL){
  
  spl <- make_train_test(df, train_end, test_start, test_end,
                         wet_day_mm=wet_day_mm, p_extreme=p_extreme)
  train <- spl$train
  test  <- spl$test
  
  # NAO subsets
  train_nao <- train %>% filter(!is.na(nao))
  test_nao  <- test  %>% filter(!is.na(nao))
  
  # If no NAO overlap in test, we can't compare NAO vs non-NAO fairly;
  # we still return base/lag scores on full test and NA for NAO models.
  if (nrow(test_nao) == 0) {
    
    # climatology on full train
    p_clim <- mean(train$exc)
    
    # fit base/lag
    m_base <- fit_base(train)
    p_base <- predict(m_base, newdata=test, type="response")
    
    train_l <- add_lag1(train)
    test_l  <- add_lag1(test)
    m_lag <- fit_lag(train_l)
    p_lag <- predict(m_lag, newdata=test_l, type="response")
    
    s_clim <- score_all(test$exc, rep(p_clim, nrow(test)))
    s_base <- score_all(test$exc, p_base)
    s_lag  <- score_all(test_l$exc, p_lag)
    
    scores <- data.frame(
      station = station_name,
      train_end = train_end,
      test_start = test_start,
      test_end = ifelse(is.null(test_end), NA, test_end),
      thr_mm = spl$thr,
      model = c("climatology","base","base+lag","base+nao","base+lag+nao","base+naoSeason","base+lag+naoSeason"),
      brier   = c(s_clim["brier"], s_base["brier"], s_lag["brier"], rep(NA_real_,4)),
      logloss = c(s_clim["logloss"], s_base["logloss"], s_lag["logloss"], rep(NA_real_,4)),
      auc     = c(s_clim["auc"], s_base["auc"], s_lag["auc"], rep(NA_real_,4))
    )
    
    diagnostics <- data.frame(
      station = station_name,
      train_end = train_end,
      test_start = test_start,
      test_end = ifelse(is.null(test_end), NA, test_end),
      model = c("base","base+lag"),
      rbind(aic_bic(m_base), aic_bic(m_lag))
    )
    
    pred_long <- tibble() # no pooled preds
    return(list(scores=scores, diagnostics=diagnostics, pred_long=pred_long))
  }
  
  # Standardize NAO using TRAIN stats (avoid leakage)
  nao_mu <- mean(train_nao$nao); nao_sd <- sd(train_nao$nao)
  if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
  train_nao <- train_nao %>% mutate(nao = (nao - nao_mu)/nao_sd)
  test_nao  <- test_nao  %>% mutate(nao = (nao - nao_mu)/nao_sd)
  
  # Evaluation set for fair comparison
  eval <- test_nao
  
  # climatology baseline (computed on FULL train)
  p_clim <- mean(train$exc)
  p_clim_eval <- rep(p_clim, nrow(eval))
  
  # Fit base + lag on FULL train (no NAO needed, but evaluated on eval)
  m_base <- fit_base(train)
  p_base_eval <- predict(m_base, newdata=eval, type="response")
  
  train_l <- add_lag1(train)
  m_lag <- fit_lag(train_l)
  
  eval_l <- add_lag1(eval)
  p_lag_eval <- predict(m_lag, newdata=eval_l, type="response")
  
  # Fit NAO models on train_nao
  # Guard against too few events
  if (sum(train_nao$exc) < 5) {
    m_nao <- m_naoS <- NULL
    m_lagnao <- m_lagnaoS <- NULL
    p_nao <- p_naoS <- rep(NA_real_, nrow(eval))
    p_lagnao <- p_lagnaoS <- rep(NA_real_, nrow(eval_l))
  } else {
    m_nao <- fit_nao(train_nao)
    p_nao <- predict(m_nao, newdata=eval, type="response")
    
    train_nao_l <- add_lag1(train_nao)
    m_lagnao <- fit_lag_nao(train_nao_l)
    
    eval_nao_l <- add_lag1(eval)
    p_lagnao <- predict(m_lagnao, newdata=eval_nao_l, type="response")
    
    m_naoS <- fit_naoSeason(train_nao)
    p_naoS <- predict(m_naoS, newdata=eval, type="response")
    
    m_lagnaoS <- fit_lag_naoSeason(train_nao_l)
    p_lagnaoS <- predict(m_lagnaoS, newdata=eval_nao_l, type="response")
  }
  
  # Scores
  s_clim <- score_all(eval$exc, p_clim_eval)
  s_base <- score_all(eval$exc, p_base_eval)
  s_lag  <- score_all(eval_l$exc, p_lag_eval)
  s_nao  <- score_all(eval$exc, p_nao)
  s_lagnao <- score_all(eval_l$exc, p_lagnao)
  s_naoS <- score_all(eval$exc, p_naoS)
  s_lagnaoS <- score_all(eval_l$exc, p_lagnaoS)
  
  scores <- data.frame(
    station = station_name,
    train_end = train_end,
    test_start = test_start,
    test_end = ifelse(is.null(test_end), NA, test_end),
    thr_mm = spl$thr,
    model = c("climatology","base","base+lag","base+nao","base+lag+nao","base+naoSeason","base+lag+naoSeason"),
    brier   = c(s_clim["brier"], s_base["brier"], s_lag["brier"], s_nao["brier"], s_lagnao["brier"], s_naoS["brier"], s_lagnaoS["brier"]),
    logloss = c(s_clim["logloss"], s_base["logloss"], s_lag["logloss"], s_nao["logloss"], s_lagnao["logloss"], s_naoS["logloss"], s_lagnaoS["logloss"]),
    auc     = c(s_clim["auc"], s_base["auc"], s_lag["auc"], s_nao["auc"], s_lagnao["auc"], s_naoS["auc"], s_lagnaoS["auc"])
  )
  
  # Diagnostics table
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
    out <- as.data.frame(as.list(aic_bic(mod)))
    out$model <- mname
    out
  })) %>%
    mutate(
      station = station_name,
      train_end = train_end,
      test_start = test_start,
      test_end = ifelse(is.null(test_end), NA, test_end)
    ) %>%
    select(station, train_end, test_start, test_end, model, AIC, BIC, logLik, edf, n)
  
  # Long predictions for pooled calibration plots
  pred_long <- bind_rows(
    tibble(station=station_name, fold=NA_integer_, model="climatology", y=eval$exc, p=p_clim_eval, group="eval"),
    tibble(station=station_name, fold=NA_integer_, model="base",        y=eval$exc, p=p_base_eval, group="eval"),
    tibble(station=station_name, fold=NA_integer_, model="base+nao",    y=eval$exc, p=p_nao, group="eval"),
    tibble(station=station_name, fold=NA_integer_, model="base+naoSeason", y=eval$exc, p=p_naoS, group="eval"),
    tibble(station=station_name, fold=NA_integer_, model="base+lag",    y=eval_l$exc, p=p_lag_eval, group="eval_l"),
    tibble(station=station_name, fold=NA_integer_, model="base+lag+nao", y=eval_l$exc, p=p_lagnao, group="eval_l"),
    tibble(station=station_name, fold=NA_integer_, model="base+lag+naoSeason", y=eval_l$exc, p=p_lagnaoS, group="eval_l")
  ) %>% filter(!is.na(p))
  
  list(scores=scores, diagnostics=diagnostics, pred_long=pred_long)
}

nao_daily <- read_excel(NAO_FILE) %>%
  mutate(
    DATE = as.Date(DATE),
    nao  = suppressWarnings(as.numeric(nao_index_cdas))
  ) %>%
  select(DATE, nao) %>%
  filter(!is.na(DATE), !is.na(nao)) %>%
  arrange(DATE)

# run all stations
all_scores <- list()
all_diagnostics <- list()
all_preds <- list()

for (st in names(stations)) {
  
  cat("\n==============================\n")
  cat("Station:", st, " file:", stations[[st]], "\n")
  cat("==============================\n")
  
  df <- read_ecad_rr(stations[[st]], st) %>%
    add_features(nao_daily)
  
  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    
    res <- run_one_split(df, st, sp$train_end, sp$test_start, sp$test_end)
    
    scores_i <- res$scores %>% mutate(fold = i)
    diagnostics_i <- res$diagnostics %>% mutate(fold = i)
    
    print(scores_i)
    
    all_scores[[length(all_scores) + 1]] <- scores_i
    all_diagnostics[[length(all_diagnostics) + 1]] <- diagnostics_i
    
    # attach fold label to preds
    if (nrow(res$pred_long) > 0) {
      all_preds[[length(all_preds) + 1]] <- res$pred_long %>% mutate(fold = i)
    }
    
    # Save fold-specific calibration plots (all models) for selected fold
    if (SAVE_PLOTS && i == PLOT_FOLD_TO_SAVE && nrow(res$pred_long) > 0) {
      st_dir <- file.path(OUT_DIR, paste0("calibration_", st, "_fold", i))
      dir.create(st_dir, showWarnings = FALSE)
      
      fold_preds <- res$pred_long %>% mutate(fold=i)
      
      for (mname in sort(unique(fold_preds$model))) {
        dfm <- fold_preds %>% filter(model == mname)
        if (nrow(dfm) == 0) next
        
        g <- plot_calibration(
          y = dfm$y,
          p = dfm$p,
          title = paste0(st, " ", mname, " calibration (fold ", i, ")")
        )
        
        ggsave(
          filename = file.path(st_dir, paste0("cal_", st, "_fold", i, "_", gsub("\\+","_", mname), ".png")),
          plot = g, width = 6, height = 4, dpi = 200
        )
      }
    }
  }
}

scores_df <- bind_rows(all_scores)
diagnostics_df <- bind_rows(all_diagnostics)
preds_df <- bind_rows(all_preds)

# Save tables
write.csv(scores_df,       file.path(OUT_DIR, "scores_all_stations_daily.csv"), row.names = FALSE)
write.csv(diagnostics_df,  file.path(OUT_DIR, "diagnostics_aic_bic_daily.csv"), row.names = FALSE)

cat("\nDONE.\nSaved outputs in:", OUT_DIR, "\n")

# pooled calibration plots
if (SAVE_PLOTS && exists("preds_df") && nrow(preds_df) > 0) {
  pooled_dir <- file.path(OUT_DIR, "calibration_pooled_allfolds")
  dir.create(pooled_dir, showWarnings = FALSE)
  
  for (st in unique(preds_df$station)) {
    stp <- preds_df %>% filter(station == st)
    
    for (mname in sort(unique(stp$model))) {
      dfm <- stp %>% filter(model == mname)
      if (nrow(dfm) == 0) next
      
      g <- plot_calibration(
        y = dfm$y,
        p = dfm$p,
        title = paste0(st, " ", mname, " calibration (POOLED folds)")
      )
      
      ggsave(
        filename = file.path(pooled_dir, paste0("cal_pooled_", st, "_", gsub("\\+","_", mname), ".png")),
        plot = g, width = 6, height = 4, dpi = 200
      )
    }
  }
  
  cat("\nSaved pooled calibration plots in:\n", pooled_dir, "\n")
}


# means scores per model
summary_mean <- scores_df %>%
  group_by(station, model) %>%
  summarise(
    mean_brier   = mean(brier, na.rm = TRUE),
    mean_logloss = mean(logloss, na.rm = TRUE),
    mean_auc     = mean(auc, na.rm = TRUE),
    .groups = "drop"
  )

# scores summary
summary_bss <- summary_mean %>%
  select(station, model, mean_brier) %>%
  pivot_wider(names_from = model, values_from = mean_brier) %>%
  mutate(
    BSS_base     = 1 - (base / climatology),
    BSS_base_lag = 1 - (`base+lag` / climatology)
  ) %>%
  select(station, BSS_base, BSS_base_lag)

# summary table
final_summary <- summary_mean %>%
  pivot_wider(
    names_from = model,
    values_from = c(mean_brier, mean_logloss, mean_auc),
    names_sep = "_"
  ) %>%
  left_join(summary_bss, by = "station") %>%
  arrange(station)

print(final_summary)
write.csv(final_summary,
          file = file.path(OUT_DIR, "summary_mean_scores_per_station_daily.csv"),
          row.names = FALSE)

cat("\nSummary table saved to:\n",
    file.path(OUT_DIR, "summary_mean_scores_per_station_daily.csv"),
    "\n")

# mean per model
overall_mean <- scores_df %>%
  group_by(model) %>%
  summarise(
    mean_brier   = mean(brier, na.rm = TRUE),
    mean_logloss = mean(logloss, na.rm = TRUE),
    mean_auc     = mean(auc, na.rm = TRUE),
    n_rows = n(),
    n_auc = sum(!is.na(auc)),
    .groups = "drop"
  )

# ranking
ranked_overall <- overall_mean %>%
  mutate(
    rank_logloss = rank(mean_logloss, ties.method = "average"),
    rank_brier   = rank(mean_brier,   ties.method = "average"),
    rank_auc     = rank(-mean_auc,    ties.method = "average"),
    rank_total   = (rank_logloss + rank_brier + rank_auc) / 3
  ) %>%
  arrange(rank_total)

print(ranked_overall)
write.csv(ranked_overall,
          file.path(OUT_DIR, "ranking_overall_models_daily.csv"),
          row.names = FALSE)

# ranking per station
ranked_per_station <- summary_mean %>%
  group_by(station) %>%
  mutate(
    rank_logloss = rank(mean_logloss, ties.method = "average"),
    rank_brier   = rank(mean_brier,   ties.method = "average"),
    rank_auc     = rank(-mean_auc,    ties.method = "average"),
    rank_total   = (rank_logloss + rank_brier + rank_auc) / 3
  ) %>%
  arrange(station, rank_total) %>%
  ungroup()

write.csv(ranked_per_station,
          file.path(OUT_DIR, "ranking_per_station_models_daily.csv"),
          row.names = FALSE)

cat("\nSaved rankings to:\n",
    file.path(OUT_DIR, "ranking_overall_models_daily.csv"), "\n",
    file.path(OUT_DIR, "ranking_per_station_models_daily.csv"), "\n")

# ============================================================
# OPTIONAL: diagnostics summary (AIC/BIC) averaged per station+model
# Note: AIC/BIC are in-sample TRAIN criteria; interpret separately from CV scores.
# ============================================================
diag_summary <- diagnostics_df %>%
  group_by(station, model) %>%
  summarise(
    mean_AIC = mean(AIC, na.rm=TRUE),
    mean_BIC = mean(BIC, na.rm=TRUE),
    mean_logLik = mean(logLik, na.rm=TRUE),
    mean_edf = mean(edf, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(station, mean_AIC)

write.csv(diag_summary,
          file.path(OUT_DIR, "diagnostics_mean_aic_bic_per_station.csv"),
          row.names = FALSE)

cat("\nSaved diagnostics summary to:\n",
    file.path(OUT_DIR, "diagnostics_mean_aic_bic_per_station.csv"),
    "\n")

# tabellen

cat("\n==============================\n")
cat("AIC / BIC summary over folds (TRAIN only)\n")
cat("==============================\n\n")

aic_all <- list()

for (st in names(stations)) {
  
  cat("\n--- Station:", st, "---\n")
  
  df <- read_ecad_rr(stations[[st]], st) %>%
    add_features(nao_daily)
  
  for (i in seq_along(splits)) {
    
    sp <- splits[[i]]
    
    spl <- make_train_test(
      df,
      train_end  = sp$train_end,
      test_start = sp$test_start,
      test_end   = sp$test_end,
      wet_day_mm = wet_day_mm,
      p_extreme  = p_extreme
    )
    
    train <- spl$train
    
    # NAO standardisatie op TRAIN
    train_nao <- train %>% filter(!is.na(nao))
    nao_mu <- mean(train_nao$nao)
    nao_sd <- sd(train_nao$nao)
    if (is.na(nao_sd) || nao_sd == 0) nao_sd <- 1
    train_nao <- train_nao %>% mutate(nao = (nao - nao_mu) / nao_sd)
    
    train_l     <- add_lag1(train)
    train_nao_l <- add_lag1(train_nao)
    
    # Fit models
    models <- list(
      base                = fit_base(train),
      base_lag            = fit_lag(train_l),
      base_nao            = fit_nao(train_nao),
      base_lag_nao        = fit_lag_nao(train_nao_l),
      base_naoSeason      = fit_naoSeason(train_nao),
      base_lag_naoSeason  = fit_lag_naoSeason(train_nao_l)
    )
    
    aic_fold <- do.call(
      rbind,
      lapply(names(models), function(mn){
        m <- models[[mn]]
        data.frame(
          station = st,
          fold    = i,
          model   = mn,
          AIC     = AIC(m),
          BIC     = BIC(m),
          logLik  = as.numeric(logLik(m)),
          edf     = sum(m$edf)
        )
      })
    )
    
    aic_all[[length(aic_all) + 1]] <- aic_fold
  }
}

aic_all_df <- bind_rows(aic_all)

# ---- Samenvatting over folds ----
aic_summary <- aic_all_df %>%
  group_by(station, model) %>%
  summarise(
    mean_AIC   = mean(AIC),
    sd_AIC     = sd(AIC),
    mean_BIC   = mean(BIC),
    mean_logLik= mean(logLik),
    mean_edf   = mean(edf),
    .groups = "drop"
  ) %>%
  arrange(station, mean_AIC)

print(aic_summary, n = 42)

write.csv(
  aic_summary,
  file.path(OUT_DIR, "diagnostics_aic_bic_mean_over_folds.csv"),
  row.names = FALSE
)

cat("\nSaved AIC/BIC fold-averaged diagnostics to:\n",
    file.path(OUT_DIR, "diagnostics_aic_bic_mean_over_folds.csv"),
    "\n")


