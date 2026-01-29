library(dplyr)
library(lubridate)

setwd("C:/Users/Sara/Downloads/Case study")
df <- read_excel("sevilla.xlsx")
df <- df %>%
  separate(
    col = 1,
    into = c("STAID", "SOUID", "DATE", "RR", "Q_RR"),
    sep = ",",
    convert = TRUE
  )
head(df)
str(df)

df$DATE <- as.Date(as.character(df$DATE), format = "%Y%m%d")
df$RR[df$RR == -9999] <- NA
df <- df[!is.na(df$RR), ]

head(df)
str(df)

df <- df %>%
  mutate(
    year  = year(DATE),
    month = month(DATE)
  )


monthly_mean <- df %>%
  group_by(month) %>%
  summarise(mean_rr = mean(RR, na.rm = TRUE))

plot(monthly_mean$month, monthly_mean$mean_rr,
     type = "l", lwd = 2,
     xlab = "Month",
     ylab = "Mean daily rainfall (mm)",
     main = "Mean monthly rainfall")

annual_mean <- df %>%
  group_by(year) %>%
  summarise(mean_rr = mean(RR, na.rm = TRUE))

plot(annual_mean$year, annual_mean$mean_rr,
     type = "l", lwd = 2,
     xlab = "Year",
     ylab = "Mean daily rainfall (mm)",
     main = "Annual mean rainfall over time")

lines(
  lowess(annual_mean$year, annual_mean$mean_rr, f = 0.2),
  col = "red", lwd = 2
)

# yearly mm
df <- df %>%
  mutate(year = year(DATE))

annual_total <- df %>%
  group_by(year) %>%
  summarise(
    total_mm = sum(RR, na.rm = TRUE),
    n_days   = n()
  )
summary(annual_total)









# beschrijvende statistiek stations samen
stations <- c(
  Madrid    = "madrid.xlsx",
  Barcelona = "barcelona.xlsx",
  Asturias  = "asturias.xlsx",
  Sevilla   = "sevilla.xlsx",
  Vigo      = "vigo.xlsx",
  Valencia  = "valencia.xlsx",
  Alicante  = "alicante.xlsx"
)

df_all <- bind_rows(lapply(names(stations), function(st) {
  read_ecad_rr(stations[[st]], st)
}))

desc_stats <- df_all %>%
  group_by(station) %>%
  summarise(
    mean_rr = mean(rr_mm),
    p95 = quantile(rr_mm[rr_mm > 1], 0.95, na.rm = TRUE),
    exc_rate = mean(rr_mm > p95, na.rm = TRUE),
    .groups = "drop"
  )

print(desc_stats)


