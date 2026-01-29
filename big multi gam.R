
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)
library(ggplot2)

setwd("C:/Users/Sara/Downloads/Case study")

# ----------------------------
# 1) Read + clean one station (daily)
# ----------------------------
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
    filter(DATE >= as.Date("1970-01-01"), DATE <= as.Date("2020-12-31")) %>%   # <-- add this
    arrange(DATE) %>%
    select(DATE, station, rr_mm)
}


# ----------------------------
# 2) DAILY exceedance indicator
# ----------------------------
make_daily_exceedance <- function(df_station, wet_day_mm = 1, p_extreme = 0.90) {
  thr <- quantile(df_station$rr_mm[df_station$rr_mm > wet_day_mm],
                  probs = p_extreme, type = 8, na.rm = TRUE)
  
  out <- df_station %>%
    mutate(exc = as.integer(rr_mm > thr)) %>%
    select(DATE, exc)
  
  list(df = out, thr = as.numeric(thr))
}

# ----------------------------
# 3) Weekly aggregation: any exceedance in ISO week
# ----------------------------
to_weekly_any <- function(df_daily_exc) {
  df_daily_exc %>%
    mutate(
      iso_year = isoyear(DATE),
      iso_week = isoweek(DATE)
    ) %>%
    group_by(iso_year, iso_week) %>%
    summarise(
      week_exc = as.integer(any(exc == 1)),
      .groups = "drop"
    ) %>%
    arrange(iso_year, iso_week)
}

# ----------------------------
# 4) Build aligned weekly panel (join on iso_year + iso_week)
# ----------------------------
build_weekly_panel <- function(stations, wet_day_mm = 1, p_extreme = 0.90) {
  
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
      t   = row_number(),           # time index (1..n_weeks)
      woy = as.integer(iso_week)    # week-of-year 1..53
    )
  
  list(df_week = df_week, thresholds = thr_tbl)
}

# ----------------------------
# 5) Joint outcomes for all stations
# ----------------------------
add_joint_outcomes <- function(df_week) {
  exc_cols <- grep("^week_exc_", names(df_week), value = TRUE)
  
  df_week %>%
    mutate(
      n_exc = rowSums(across(all_of(exc_cols))),
      joint_2p_w = as.integer(n_exc >= 2),
      joint_3p_w = as.integer(n_exc >= 3)
    )
}

# ----------------------------
# 6) Fit GAM with correct time variable (t)
# ----------------------------
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

# ----------------------------
# 7) Plot time trend (hold season fixed)
# ----------------------------
plot_weekly_time_trend <- function(data, model, title, woy_fixed = 26) {
  nd <- data.frame(
    t   = seq(min(data$t), max(data$t), length.out = 400),
    woy = woy_fixed
  )
  
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  
  invlogit <- function(x) 1 / (1 + exp(-x))
  nd$fit <- invlogit(pr$fit)
  nd$lo  <- invlogit(pr$fit - 1.96 * pr$se.fit)
  nd$hi  <- invlogit(pr$fit + 1.96 * pr$se.fit)
  
  # map t back to a date for plotting
  # approximate by linear mapping between first and last week_start
  nd$week_start <- data$week_start[1] + (nd$t - 1) * 7
  
  ggplot(nd, aes(x = week_start, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = "Probability", title = title) +
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
  Alicante  = "alicante.xlsx",
  Valencia  = "valencia.xlsx"
)

wet_day_mm <- 1
p_extreme  <- 0.90

panel <- build_weekly_panel(stations, wet_day_mm = wet_day_mm, p_extreme = p_extreme)
df_week <- panel$df_week
print(panel$thresholds)

df_week <- add_joint_outcomes(df_week)

summ_counts <- df_week %>%
  summarise(
    n_weeks   = n(),
    ones_2p   = sum(joint_2p_w),
    ones_3p   = sum(joint_3p_w),
    mean_2p   = mean(joint_2p_w),
    mean_3p   = mean(joint_3p_w),
    avg_n_exc = mean(n_exc)
  )
print(summ_counts)

m_2p <- fit_weekly_gam(df_week, "joint_2p_w")
m_3p <- fit_weekly_gam(df_week, "joint_3p_w")

print(summary(m_2p))
print(gam.check(m_2p))

print(summary(m_3p))
print(gam.check(m_3p))

p2 <- plot_weekly_time_trend(df_week, m_2p, "Weekly joint exceedance (>=2 of 5)")
p3 <- plot_weekly_time_trend(df_week, m_3p, "Weekly joint exceedance (>=3 of 5)")

print(p2)
print(p3)


# seizoensplot
plot_seasonality <- function(data, model, title) {
  
  nd <- data.frame(
    t   = median(data$t),     # fixeer tijd (langetermijntrend)
    woy = 1:53                # alle weken van het jaar
  )
  
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  
  invlogit <- function(x) 1 / (1 + exp(-x))
  nd$fit <- invlogit(pr$fit)
  nd$lo  <- invlogit(pr$fit - 1.96 * pr$se.fit)
  nd$hi  <- invlogit(pr$fit + 1.96 * pr$se.fit)
  
  ggplot(nd, aes(x = woy, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = seq(1, 53, by = 4)) +
    labs(
      x = "Week of year",
      y = "Probability",
      title = title
    ) +
    theme_minimal()
}

p_season_2p <- plot_seasonality(
  df_week,
  m_2p,
  "Seasonal pattern of weekly joint exceedance (≥2 of 7 stations)"
)

print(p_season_2p)



# tijd plot verschillende weken
plot_time_fixed_season_multi <- function(data, model, title, woy_values = c(10, 26, 45)) {
  
  nd <- expand.grid(
    t   = seq(min(data$t), max(data$t), length.out = 400),
    woy = as.integer(woy_values)
  )
  
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  
  invlogit <- function(x) 1 / (1 + exp(-x))
  nd$fit <- invlogit(pr$fit)
  nd$lo  <- invlogit(pr$fit - 1.96 * pr$se.fit)
  nd$hi  <- invlogit(pr$fit + 1.96 * pr$se.fit)
  
  nd$week_start <- data$week_start[1] + (nd$t - 1) * 7
  nd$woy <- factor(nd$woy, levels = woy_values, labels = paste("woy =", woy_values))
  
  ggplot(nd, aes(x = week_start, y = fit, group = woy)) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = "Probability", title = title) +
    theme_minimal()
}

p_time_multi <- plot_time_fixed_season_multi(
  df_week, m_2p,
  "Time trend (season fixed at selected weeks): weekly joint exceedance (≥2 of 7)",
  woy_values = c(10, 26, 45)
)
print(p_time_multi)

# tijdplot juiste labels/jaren
plot_time_fixed_season_index <- function(data, model, title, woy_fixed = 26) {
  
  nd <- data.frame(
    t   = seq(min(data$t), max(data$t), length.out = 400),
    woy = as.integer(woy_fixed)
  )
  
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  
  invlogit <- function(x) 1 / (1 + exp(-x))
  nd$fit <- invlogit(pr$fit)
  nd$lo  <- invlogit(pr$fit - 1.96 * pr$se.fit)
  nd$hi  <- invlogit(pr$fit + 1.96 * pr$se.fit)
  
  ggplot(nd, aes(x = t, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) +
    geom_line(linewidth = 1) +
    labs(
      x = "Time (weekly index)",
      y = "Probability",
      title = paste0(title, " (season fixed at week ", woy_fixed, ")")
    ) +
    theme_minimal()
}

print(plot_time_fixed_season_index(
  df_week, m_2p,
  "Time trend of weekly joint exceedance (≥2 of 7)"
))


# count exceedances
df_week %>%
  summarise(
    n_weeks = n(),
    joint_2p_events = sum(joint_2p_w),
    share_joint_2p = mean(joint_2p_w),
    joint_3p_events = sum(joint_3p_w),
    share_joint_3p = mean(joint_3p_w)
  )
df_week %>%
  mutate(season = case_when(
    woy %in% 1:8   ~ "Winter",
    woy %in% 9:21  ~ "Spring",
    woy %in% 22:34 ~ "Summer",
    TRUE           ~ "Autumn"
  )) %>%
  group_by(season) %>%
  summarise(
    n_weeks = n(),
    joint_2p_rate = mean(joint_2p_w),
    joint_3p_rate = mean(joint_3p_w)
  )
df_week %>%
  mutate(decade = floor(year(week_start) / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    joint_2p_rate = mean(joint_2p_w),
    joint_3p_rate = mean(joint_3p_w),
    n = n()
  )

exc_cols <- grep("^week_exc_", names(df_week), value = TRUE)
stations <- gsub("^week_exc_", "", exc_cols)

pairs <- combn(stations, 2, simplify = FALSE)

pairwise_counts <- lapply(pairs, function(p) {
  a <- paste0("week_exc_", p[1])
  b <- paste0("week_exc_", p[2])
  
  tibble(
    pair = paste(p[1], p[2], sep = "–"),
    joint_rate = mean(df_week[[a]] == 1 & df_week[[b]] == 1),
    joint_events = sum(df_week[[a]] == 1 & df_week[[b]] == 1)
  )
}) %>%
  bind_rows() %>%
  arrange(desc(joint_rate))

print(pairwise_counts)

df_week_season <- df_week %>%
  mutate(season = case_when(
    woy %in% 1:8   ~ "Winter",
    woy %in% 9:21  ~ "Spring",
    woy %in% 22:34 ~ "Summer",
    TRUE           ~ "Autumn"
  ))

pairwise_season <- lapply(pairs, function(p) {
  a <- paste0("week_exc_", p[1])
  b <- paste0("week_exc_", p[2])
  
  df_week_season %>%
    group_by(season) %>%
    summarise(
      pair = paste(p[1], p[2], sep = "–"),
      joint_rate = mean(.data[[a]] == 1 & .data[[b]] == 1),
      .groups = "drop"
    )
}) %>%
  bind_rows()

print(pairwise_season)

pairwise_season %>%
  arrange(season, desc(joint_rate)) %>%
  print(n = 60)





# lag
df_week2 <- df_week %>%
  arrange(week_start) %>%
  mutate(
    t = row_number(),              # CREATE time index
    lag1 = dplyr::lag(joint_2p_w)
  ) %>%
  filter(!is.na(lag1))

m_lag <- gam(
  joint_2p_w ~ s(t, k = 60) + s(woy, bs = "cc", k = 20) + lag1,
  family = binomial(link = "logit"),
  method = "REML",
  data = df_week2,
  knots = list(woy = c(0.5, 53.5))
)
summary(m_lag)




# nao
# ============================================================
# NAO -> weekly join + GAM for weekly joint exceedance (FIXED)
# ============================================================

library(readxl)
library(dplyr)
library(lubridate)
library(mgcv)

# 1) Read NAO
nao_df <- read_excel("NAO_cleaned.xlsx") %>%
  mutate(
    DATE = as.Date(DATE),
    iso_year = isoyear(DATE),
    iso_week = isoweek(DATE),
    nao = nao_index_cdas
  ) %>%
  select(iso_year, iso_week, nao) %>%
  filter(!is.na(nao))

# 2) Aggregate NAO to weekly
nao_week <- nao_df %>%
  group_by(iso_year, iso_week) %>%
  summarise(nao_week = mean(nao, na.rm = TRUE), .groups = "drop")

# 3) Join + rebuild time index safely (DON'T call it 't')
df_week_nao <- df_week %>%
  left_join(nao_week, by = c("iso_year", "iso_week")) %>%
  arrange(week_start) %>%
  mutate(
    time_idx = row_number()   # safe time index
  )

cat("NAO missing weeks:", sum(is.na(df_week_nao$nao_week)), "out of", nrow(df_week_nao), "\n")

# 4) Fit GAM with NAO
m_2p_nao <- gam(
  joint_2p_w ~ s(time_idx, k = 30) + s(woy, bs = "cc", k = 20) + nao_week,
  family = binomial(link = "logit"),
  method = "REML",
  data = df_week_nao,
  knots = list(woy = c(0.5, 53.5))
)

summary(m_2p_nao)
gam.check(m_2p_nao)




# specifieke storm / week zoeken overschrijding welke stations
df_week %>%
  filter(week_start == as.Date("2019-09-09")) %>%
  pivot_longer(
    cols = starts_with("week_exc_"),
    names_to = "station",
    values_to = "exceedance"
  ) %>%
  filter(exceedance == 1) %>%
  mutate(station = gsub("week_exc_", "", station)) %>%
  summarise(
    week_start = first(week_start),
    n_exc      = first(n_exc),
    joint_2p_w = first(joint_2p_w),
    joint_3p_w = first(joint_3p_w),
    stations   = paste(station, collapse = ", ")
  )


# topweken met stations
options(tibble.width = Inf)
options(pillar.width = Inf)

out_top <- df_week %>%
  arrange(desc(n_exc)) %>%
  slice_head(n = 10) %>%
  pivot_longer(
    cols = starts_with("week_exc_"),
    names_to = "station",
    values_to = "exceedance"
  ) %>%
  filter(exceedance == 1) %>%
  mutate(station = gsub("week_exc_", "", station)) %>%
  group_by(week_start, n_exc, joint_2p_w, joint_3p_w) %>%
  summarise(stations = paste(station, collapse = ", "), .groups = "drop") %>%
  arrange(desc(n_exc))

print(out_top, n = Inf)

# evt rolling nog bekijken