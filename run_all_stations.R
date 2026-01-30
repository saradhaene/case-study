# ============================================================
# Master runner: sequentially run validation workflow for all stations
# Uses environment overrides consumed by `validatie alicante.R`.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
})

DATA_DIR <- Sys.getenv("DATA_DIR", "C:/Users/Sara/Downloads/Case study")
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

NAO_FILE <- Sys.getenv("NAO_FILE", "NAO_cleaned.xlsx")

run_station <- function(st_name, st_file) {
  message("\n=== Running station: ", st_name, " (", st_file, ") ===")
  Sys.setenv(
    DATA_DIR = DATA_DIR,
    STATION_NAME = st_name,
    STATION_FILE = st_file,
    NAO_FILE = NAO_FILE
  )
  source(file.path(DATA_DIR, "validatie alicante.R"), local = TRUE)
}

for (st in names(stations)) {
  run_station(st, stations[[st]])
}
