###############################################
# WA Pertussis + Google Trends — surge-aware v2
# weekly GT, robust fits, year-specific RMSE
###############################################

# ── 0. PARAMETERS ─────────────────────────────────────────────────────────
region_code      <- "US-WA"
time_window      <- "2014-01-01 2024-12-31"   # Google-Trends span
surge_split_date <- as.Date("2024-09-01")     # test-window start
keep_top_cor     <- 12                        # GT feature cap

terms <- c(
  "whooping cough","pertussis","bordetella",
  "whooping cough symptoms","whooping cough adults","persistent cough",
  "violent cough","cough won't stop","long coughing fit","coughing fits",
  "coughing spasm","gasp for air when coughing","nonstop cough at night",
  "baby coughing fits","bad cough causes",
  "whooping cough treatment","pertussis treatment",
  "pertussis antibiotics","azithromycin whooping cough",
  "whooping cough vaccine","pertussis vaccine",
  "Tdap shot pregnant","Tdap pregnancy","Tdap booster",
  "pertusis","bordatella","whooping cofh",
  "tos ferina","tos convulsiva","vacuna dtpa",
  "bordetella tos ferina","pertusis sintomas","bordatella sintomas"
)

# ── 1. PACKAGES ───────────────────────────────────────────────────────────
pkgs <- c(
  "gtrendsR","tidyverse","lubridate","zoo",
  "forecast","nlme","MuMIn","janitor","glue","glmnet"
)
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

# ── 2. GOOGLE-TRENDS (cached weekly) ──────────────────────────────────────
dir.create("data", showWarnings = FALSE)
cache_file <- "data/trends_raw_all_years.rds"

if (file.exists(cache_file)) {
  message("✔ Using cached Google-Trends pull: ", cache_file)
  trends_raw <- readRDS(cache_file)
} else {
  batches <- split(terms, ceiling(seq_along(terms) / 5))      # 5 keywords / request
  trends_raw <- purrr::map_dfr(
    batches,
    \(b) gtrends(b, geo = region_code, time = time_window, gprop = "web")$interest_over_time
  )
  if (!nrow(trends_raw)) stop("Google Trends download failed.")
  saveRDS(trends_raw, cache_file)
}

# ── 3. pull data from google trends ───────────────────────────────────────────────────────────
trends_weekly <- trends_raw %>%
  transmute(
    date   = floor_date(date, "week", week_start = 1),
    keyword,
    hits   = as.numeric(replace(hits, hits == "<1", 0))
  ) %>%
  pivot_wider(names_from = keyword, values_from = hits, names_prefix = "gt_") %>%
  janitor::clean_names() %>%
  arrange(date)

# ── 4. PERTUSSIS CASES (weekly) ───────────────────────────────────────────
pertussis_weekly <- readr::read_csv(
  "wa_all_years.csv",                        # ← combined 2014-19 + 2022-24 file
  show_col_types = FALSE,
  col_types = cols(
    date  = col_date(),
    cases = col_double()
  )
) %>%
  mutate(date = floor_date(date, "week", week_start = 1)) %>%
  filter(!year(date) %in% c(2020, 2021)) %>%   # drop 2020–21
  group_by(date) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# ── 5. MERGE & LAG (0-3 weeks) ────────────────────────────────────────────
combined <- inner_join(trends_weekly, pertussis_weekly, by = "date") %>% arrange(date)
gt_cols  <- names(combined)[startsWith(names(combined), "gt_")]

lagged <- combined %>%
  mutate(across(
    all_of(gt_cols),
    list(
      l0 = \(x) x,
      l1 = \(x) lag(x, 1),
      l2 = \(x) lag(x, 2),
      l3 = \(x) lag(x, 3)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  drop_na()

lagged$date <- as.Date(lagged$date)

# ── 6. TRAIN / TEST SPLIT ─────────────────────────────────────────────────
train <- lagged %>% filter(date <= surge_split_date)
test  <- lagged %>% filter(date >  surge_split_date)

# ── 7. FEATURE SCREENING ─────────────────────────────────────────────────
lag_cols <- grep("^gt_.*_l[0-3]$", names(lagged), value = TRUE)
var_ok   <- lag_cols[apply(train[, lag_cols], 2, \(z) var(z, na.rm = TRUE) > 0)]
cor_vals <- sapply(var_ok, \(c) cor(train[[c]], train$cases, use = "complete.obs"))
keep_cols <- names(sort(abs(cor_vals), decreasing = TRUE))[1:min(keep_top_cor, length(var_ok))]

x_pre    <- as.matrix(train[, keep_cols])
good_cols<- colnames(x_pre)[seq_len(qr(x_pre)$rank)]
x_train  <- as.matrix(train[, good_cols])
x_test   <- as.matrix(test [, good_cols])

# ── 8. ELASTIC-NET ────────────────────────────────────────────────────────
cv_fit <- cv.glmnet(x_train, train$cases, alpha = 0.5, family = "gaussian")
train$enet_score <- as.numeric(predict(cv_fit, x_train, s = "lambda.min"))
test$enet_score  <- as.numeric(predict(cv_fit, x_test,  s = "lambda.min"))

# ── 9. ARIMAX & BASELINE ARIMA ────────────────────────────────────────────
y_train <- ts(train$cases, frequency = 52)

fit_arimax <- auto.arima(y_train, xreg = x_train, stepwise = FALSE, approximation = FALSE)
fc_arimax  <- forecast(fit_arimax, xreg = x_test)

fit_arima_base <- auto.arima(y_train)
fc_base        <- forecast(fit_arima_base, h = nrow(test))

# ── 10. GLS (AR1)  fallback to LM ─────────────────────────────────────────
train$id <- test$id <- factor(1)
gls_formula <- as.formula(paste("cases ~", paste(good_cols, collapse = " + ")))
gls_fit <- tryCatch(
  gls(gls_formula, data = train, correlation = corAR1(~1 | id)),
  error = \(e) { warning("GLS failed – using LM: ", e$message); lm(gls_formula, data = train) }
)
gls_pred <- predict(gls_fit, newdata = test)

# ── 11. METRICS (test window) ─────────────────────────────────────────────
rmse <- \(p, o) sqrt(mean((p - o)^2))
accuracy_table <- tibble(
  model = c("ARIMA_baseline", "ARIMAX", if (inherits(gls_fit,"gls")) "GLS" else "LM"),
  RMSE  = c(
    rmse(fc_base$mean,   test$cases),
    rmse(fc_arimax$mean, test$cases),
    rmse(gls_pred,       test$cases)
  )
)

# ── 12. Correlation table (top GT predictors) ─────────────────────────────
cor_table <- lagged %>%
  select(cases, all_of(good_cols)) %>%
  summarise(across(-cases, ~ cor(.x, cases, use = "complete.obs"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

# helper: leave-one-year-out evaluation ------------------------------------
year_eval <- function(target_year) {
  train_idx <- lubridate::year(lagged$date) != target_year
  test_idx  <- !train_idx
  if (!any(test_idx)) return(tibble())

  # baseline ARIMA
  fit_base <- auto.arima(ts(lagged$cases[train_idx], frequency = 52))
  fc_base  <- forecast(fit_base, h = sum(test_idx))
  err_base <- fc_base$mean - lagged$cases[test_idx]

  # ARIMAX
  fit_x <- tryCatch(
    auto.arima(ts(lagged$cases[train_idx], frequency = 52),
               xreg = as.matrix(lagged[train_idx, good_cols])),
    error = \(e) NULL
  )
  if (is.null(fit_x)) return(tibble())

  fc_x   <- forecast(fit_x, xreg = as.matrix(lagged[test_idx, good_cols]))
  err_x  <- fc_x$mean - lagged$cases[test_idx]

  dm_p <- tryCatch({
    forecast::dm.test(err_base, err_x, alternative = "greater")$p.value
  }, error = function(e) NA_real_)

  tibble(
    year        = target_year,
    RMSE_ARIMA  = sqrt(mean(err_base^2)),
    RMSE_ARIMAX = sqrt(mean(err_x^2)),
    DM_pvalue   = dm_p,
    significant = if (!is.na(dm_p)) dm_p < 0.05 else NA
  )
}

# ── 13. LEAVE-ONE-YEAR-OUT RMSE (2022-24) ─────────────────────────────────
year_rmse_table <- purrr::map_dfr(2014:2024, year_eval)

# add mean weekly cases & normalised RMSE ----------------------------------
year_means <- readr::read_csv(
  "wa_all_years.csv", show_col_types = FALSE,
  col_types = cols(date = col_date(), cases = col_double())
) %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(mean_cases = mean(cases), .groups = "drop")

year_rmse_table <- year_rmse_table %>%
  left_join(year_means, by = "year") %>%
  mutate(
    nRMSE_ARIMA  = RMSE_ARIMA  / mean_cases,
    nRMSE_ARIMAX = RMSE_ARIMAX / mean_cases
  )

# Diebold-Mariano on surge-period test set --------------------------------
dm_p <- dm.test(
  e1 = fc_base$mean - test$cases,
  e2 = fc_arimax$mean - test$cases,
  alternative = "greater"
)$p.value

message("\nHeld-out RMSE table (plus normalised RMSE):")
print(year_rmse_table, digits = 3)

message(sprintf(
  "\nDiebold–Mariano test on 2024-surge test window:\n  p-value = %.4f  (smaller ⇒ ARIMAX significantly better)",
  dm_p
))

# ── 14. SAVE ARTEFACTS ────────────────────────────────────────────────────
dir.create("artefacts", showWarnings = FALSE)
saveRDS(list(
  df_full        = lagged,
  df_train       = train,
  df_test        = test,
  fc_arimax      = fc_arimax,
  fc_base        = fc_base,
  fit_arima_base = fit_arima_base,
  arimax_fit     = fit_arimax,
  gls_fit        = gls_fit,
  accuracy_table = accuracy_table,
  year_rmse_table= year_rmse_table,
  cor_table      = cor_table
), "artefacts/wa_pertussis_models.Rds")

message("✔ artefacts/wa_pertussis_models.Rds written")

library(dplyr)
library(stringr)

cor_table <- readRDS("artefacts/wa_pertussis_models.Rds")$cor_table

pretty_terms <- cor_table %>%
  mutate(
    term      = str_remove(variable, "^gt_") %>% str_remove("_l[0-3]$") %>% str_replace_all("_", " "),
    lag_weeks = str_extract(variable, "_l[0-3]$") %>% str_remove("_l") %>% as.integer()
  ) %>%
  select(`Search Term` = term, `Lag (Weeks)` = lag_weeks, `Correlation` = correlation) %>%
  arrange(desc(abs(Correlation)))

print(pretty_terms, n = 12)

# ── 15. ESTIMATE LEAD TIME TO DOH OUTBREAK ALERT ─────────────────────────

# Set official DOH alert date
doh_alert_date <- as.Date("2024-11-12")

# Add forecast dates to ARIMAX predictions
arimax_preds <- tibble(
  date      = test$date,
  forecast  = as.numeric(fc_arimax$mean)
)

# Estimate baseline using the final 4 weeks of training
baseline_mean <- tail(train$cases, 4) |> mean(na.rm = TRUE)
surge_threshold <- baseline_mean * 1.5  # Define surge as 50% above recent avg

# Identify first forecasted surge
predicted_surge <- arimax_preds %>%
  filter(forecast >= surge_threshold) %>%
  arrange(date) %>%
  slice(1)

if (nrow(predicted_surge) == 0) {
  message("⚠ No predicted surge detected above threshold.")
} else {
  lead_days <- as.integer(doh_alert_date - predicted_surge$date)
  message(glue::glue(
    "\n First predicted surge week: {predicted_surge$date}",
    "\n DOH official alert: {doh_alert_date}",
    "\n Model lead time: {lead_days} days before DOH alert."
  ))
}