# vaccine_trends_pull.R

# ── Packages ─────────────────────────────────────
pkgs <- c("gtrendsR", "tidyverse", "lubridate")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
lapply(pkgs, library, character.only = TRUE)

# ── Params ───────────────────────────────────────A
region_code <- "US-WA"
time_window <- "2022-01-01 2024-12-31"

vaccine_terms <- c(
  "whooping cough vaccine",
  "pertussis vaccine"
)

# ── Safe wrapper ────────────────────────────────
safe_gtrends <- function(...) tryCatch(gtrends(...), error = function(e) NULL)

# ── Pull trends in batches ──────────────────────
batches <- split(vaccine_terms, ceiling(seq_along(vaccine_terms) / 4))
trends_raw <- purrr::map_dfr(batches, \(b) {
  Sys.sleep(2)  # avoid rate-limiting
  safe_gtrends(b, geo = region_code, time = time_window)$interest_over_time
})

if (is.null(trends_raw) || nrow(trends_raw) == 0) stop("Google Trends failed.")

# ── Clean trends data ───────────────────────────
vaccine_trends <- trends_raw %>%
  transmute(date = floor_date(date, "week", week_start = 1),
            keyword,
            hits = as.numeric(replace(hits, hits == "<1", 0))) %>%
  pivot_wider(names_from = keyword, values_from = hits) %>%
  arrange(date)

# ── Save ────────────────────────────────────────
dir.create("data", showWarnings = FALSE)
readr::write_csv(vaccine_trends, "data/vaccine_trends.csv")
saveRDS(vaccine_trends, "data/vaccine_trends.rds")

message("✔ Vaccine Google Trends saved to data/vaccine_trends.csv and .rds")

