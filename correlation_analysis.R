###############################################
# WA Pertussis + Google Trends – full pipeline
#   • correlation scan
#   • ARIMA / ARIMAX fit
#   • publication tables
###############################################

# ── 0. PARAMETERS ─────────────────────────────────────────────────────────
region_code      <- "US-WA"
time_window      <- "2022-01-01 2024-12-31"
surge_split_date <- as.Date("2024-09-01")
keep_top_cor     <- 12                       # predictors to keep for model

terms <- c(                                   # full candidate list
  "whooping cough","pertussis","bordetella",
  "whooping cough symptoms","whooping cough adults","persistent cough",
  "violent cough","cough won't stop","long coughing fit","coughing fits",
  "coughing spasm","gasp for air when coughing","nonstop cough at night",
  "baby coughing fits","bad cough causes",
  "whooping cough treatment","pertussis treatment",
  "pertussis antibiotics","azithromycin whooping cough",
  "whooping cough vaccine","pertussis vaccine",
  "Tdap shot pregnant","Tdap pregnancy","Tdap booster",
  "pertusis","whooping cofh",
  "tos ferina","tos convulsiva","vacuna dtpa",
  "bordetella tos ferina","pertusis sintomas","bordetella sintomas"
)

# ── 1. PACKAGES ───────────────────────────────────────────────────────────
pkgs <- c("gtrendsR","tidyverse","lubridate","zoo",
          "forecast","glmnet","DT","knitr","kableExtra")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need))
  install.packages(need, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

# ── 2. GOOGLE-TRENDS (cached) ────────────────────────────────────────────
dir.create("data",      showWarnings = FALSE)
dir.create("artefacts", showWarnings = FALSE)

cache_file <- "data/trends_full_set.rds"
if (file.exists(cache_file)) {
  message("✔ Using cached Google-Trends pull: ", cache_file)
  trends_raw <- readRDS(cache_file)
} else {
  message("⟳ Downloading Google-Trends data …")
  trends_raw <- purrr::map_dfr(
    split(terms, ceiling(seq_along(terms) / 5)),
    \(batch) gtrends(batch, geo = region_code,
                     time = time_window, gprop = "web")$interest_over_time
  )
  if (!nrow(trends_raw)) stop("Google Trends download failed.")
  saveRDS(trends_raw, cache_file)
  message("✔ Saved to ", cache_file)
}

# ── 3. CLEAN GT (duplicate-safe) ─────────────────────────────────────────
trends_weekly <- trends_raw %>%
  transmute(date = floor_date(date, "week", week_start = 1),
            keyword,
            hits  = as.numeric(replace(hits, hits == "<1", 0))) %>%
  group_by(date, keyword) %>%                         # collapse dups
  summarise(hits = mean(hits), .groups = "drop") %>%
  pivot_wider(names_from = keyword, values_from = hits,
              names_prefix = "gt_", values_fill = 0) %>%
  janitor::clean_names() %>%
  arrange(date)

# ── 4. PERTUSSIS CASES (weekly) ───────────────────────────────────────────
pertussis_weekly <- readr::read_csv(
  "wa_formatted_doh.csv", show_col_types = FALSE,
  col_types = cols(date = col_date(), cases = col_double())
) %>%
  mutate(date = floor_date(date, "week", week_start = 1)) %>%
  group_by(date) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# ── 5. MERGE & LAG 0-2 ---------------------------------------------------
combined <- inner_join(trends_weekly, pertussis_weekly, by = "date") %>% arrange(date)
gt_cols  <- names(combined)[startsWith(names(combined), "gt_")]

lagged <- combined %>%
  mutate(across(all_of(gt_cols),
                list(l0 = identity,
                     l1 = \(x) lag(x, 1),
                     l2 = \(x) lag(x, 2)),
                .names = "{.col}_{.fn}")) %>%
  drop_na()

# ── 6. CORRELATION TABLE (all lags) ───────────────────────────────────────
lag_cols <- grep("^gt_.*_l[0-2]$", names(lagged), value = TRUE)
var_ok   <- lag_cols[apply(lagged[, lag_cols], 2, \(z) var(z, na.rm = TRUE) > 0)]
cor_vals <- sapply(var_ok, \(c) cor(lagged[[c]], lagged$cases, use = "complete.obs"))

cor_table_all <- tibble(
  variable    = names(cor_vals),
  correlation = cor_vals
) %>%
  arrange(desc(abs(correlation))) %>%
  mutate(
    `Search Term` = str_remove(variable, "^gt_") |> str_remove("_l[0-2]$") |> str_replace_all("_", " "),
    `Lag (weeks)` = str_extract(variable, "_l[0-2]$") |> str_remove("_l") |> as.integer(),
    Corr.         = correlation
  ) %>%
  select(`Search Term`, `Lag (weeks)`, Corr.)

write_csv(cor_table_all, "artefacts/all_gt_correlations_full.csv")
message("✔ Correlation table written to artefacts/all_gt_correlations_full.csv")

# ── 7. SELECT TOP-N PREDICTORS ────────────────────────────────────────────
keep_cols <- cor_table_all %>%
  slice_head(n = keep_top_cor) %>%
  transmute(colname = glue::glue("gt_{str_replace_all(`Search Term`, ' ', '_')}_l{`Lag (weeks)`}")) %>%
  pull()

# ensure column exists & has variance
keep_cols <- keep_cols[keep_cols %in% names(lagged)]
keep_cols <- keep_cols[apply(lagged[, keep_cols], 2, \(z) var(z, na.rm = TRUE) > 0)]

# ── 8. TRAIN / TEST SPLIT ────────────────────────────────────────────────
train <- lagged %>% filter(as.Date(date) <= surge_split_date)
test  <- lagged %>% filter(as.Date(date) >  surge_split_date)

x_train <- as.matrix(train[, keep_cols])
x_test  <- as.matrix(test [, keep_cols])

# ── 9. MODELS & ACCURACY ─────────────────────────────────────────────────
cv_fit <- cv.glmnet(x_train, train$cases, alpha = 0.5, family = "gaussian")
train$enet_score <- predict(cv_fit, x_train, s = "lambda.min")[,1]
test $enet_score <- predict(cv_fit, x_test,  s = "lambda.min")[,1]

y_train <- ts(train$cases, frequency = 52)

fit_arimax <- auto.arima(y_train, xreg = x_train, stepwise = FALSE, approximation = FALSE)
fc_arimax  <- forecast(fit_arimax, xreg = x_test)

fit_arima_base <- auto.arima(y_train)
fc_base        <- forecast(fit_arima_base, h = nrow(test))

rmse <- \(p, o) sqrt(mean((p - o)^2))
accuracy_table <- tibble(
  model = c("ARIMA_baseline", "ARIMAX"),
  RMSE  = c(rmse(fc_base$mean, test$cases),
            rmse(fc_arimax$mean, test$cases))
)

print(accuracy_table, digits = 3)

# ── 10. SAVE ARTEFACTS ────────────────────────────────────────────────────
saveRDS(list(
  cor_table_full  = cor_table_all,
  keep_cols       = keep_cols,
  arimax_fit      = fit_arimax,
  arima_baseline  = fit_arima_base,
  accuracy_table  = accuracy_table
), "artefacts/wa_pertussis_full_term_model.Rds")

message("✔ artefacts/wa_pertussis_full_term_model.Rds written")

# ── 11. PUBLICATION TABLES ────────────────────────────────────────────────
focus_terms <- c("whooping cough",
                 "whooping cough symptoms",
                 "pertussis")

model_terms <- cor_table_all %>%
  filter(`Search Term` %in% focus_terms, `Lag (weeks)` %in% 0:2) %>%
  arrange(`Search Term`, `Lag (weeks)`)

write_csv(model_terms, "artefacts/gt_terms_used_in_model.csv")

# Table A1 – interactive (HTML output)
datatable(
  cor_table_all,
  rownames = FALSE,
  options  = list(pageLength = 25, scrollX = TRUE),
  caption  = htmltools::HTML(
    "<b>Table A1.</b> All Google-Trends term × lag combinations evaluated (sorted by |ρ|)."
  )
)

# Table 1 – static (HTML/PDF via knitr)
knitr::kable(
  model_terms,
  booktabs = TRUE,
  digits   = 3,
  caption  = "Table 1. Google-Trends predictors retained for the ARIMAX model."
) %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Helvetica")