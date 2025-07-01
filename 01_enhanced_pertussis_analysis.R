###############################################
# Cleaned WA Pertussis + Google Trends Script
# Works end-to-end (R ≥ 4.1)
###############################################

# 1. Load required packages -----------------------------------------------
pkgs <- c("gtrendsR", "tidyverse", "lubridate", "zoo",
          "forecast", "nlme", "MuMIn", "janitor", "glue")
install_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(install_missing))
  install.packages(install_missing, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# 2. Download Google Trends data ------------------------------------------
terms       <- c("pertussis", "whooping cough", "pertusis", "bordatella", "bordetella")
region_code <- "US-WA"
time_window <- "2022-01-01 2024-12-31"

trends_raw <- gtrends(
  keyword = terms,
  geo     = region_code,
  time    = time_window,
  gprop   = "web"
)$interest_over_time

if (nrow(trends_raw) == 0) stop("Google Trends download failed.")

# 3. Clean & reshape Trends data ------------------------------------------
trends_weekly <- trends_raw %>%
  transmute(
    date = floor_date(date, "week", week_start = 1),
    keyword,
    hits = as.numeric(replace(hits, hits == "<1", 0))
  ) %>%
  pivot_wider(
    names_from = keyword,
    values_from = hits,
    names_prefix = "gt_"
  ) %>%
  janitor::clean_names() %>%
  arrange(date)

# 4. Read pertussis case data --------------------------------------------
pertussis_weekly <- readr::read_csv(
  "wa_formatted.csv",
  show_col_types = FALSE,
  col_types = cols(date = col_date(), cases = col_double())
) %>%
  mutate(date = floor_date(date, "week", week_start = 1)) %>%
  group_by(date) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# 5. Merge & create GT lags ----------------------------------------------
combined <- trends_weekly %>%
  inner_join(pertussis_weekly, by = "date") %>%
  arrange(date)

lagged <- combined %>%
  mutate(across(starts_with("gt_pertussis"),
                list(l0 = \(x) x,
                     l1 = \(x) lag(x,1),
                     l2 = \(x) lag(x,2),
                     l3 = \(x) lag(x,3)),
                .names = "{.col}_{.fn}")) %>%
  drop_na()

# 6. Train/test split -----------------------------------------------------
split_idx <- floor(0.8 * nrow(lagged))
train <- lagged[1:split_idx, ]
test  <- lagged[(split_idx + 1):nrow(lagged), ]

x_cols <- grep("^gt_pertussis_l[0-3]$", names(lagged), value = TRUE)

# 7. ARIMAX model ---------------------------------------------------------
y_train <- ts(train$cases, frequency = 52)
x_train <- as.matrix(train[, x_cols])
x_test  <- as.matrix(test[,  x_cols])

fit_arimax <- auto.arima(y_train, xreg = x_train)
fc         <- forecast(fit_arimax, xreg = x_test)
print(accuracy(fc$mean, test$cases))

# 8. GLS model ------------------------------------------------------------
train$id <- factor(1)
test$id  <- factor(1)

gls_formula <- as.formula(paste("cases ~", paste(x_cols, collapse = " + ")))
gls_fit <- gls(
  gls_formula, data = train,
  correlation = corAR1(form = ~ 1 | id)
)

summary(gls_fit)

test$gls_pred <- predict(gls_fit, newdata = test)
rmse_gls <- sqrt(mean((test$gls_pred - test$cases)^2))
cat("GLS AR(1) RMSE =", round(rmse_gls, 2), "\n")

# 9. Model averaging ------------------------------------------------------
global_lm <- lm(gls_formula, data = train)
options(na.action = "na.fail")
dredge_set <- dredge(global_lm, rank = "AICc", trace = FALSE)
avg_mod <- model.avg(dredge_set, subset = delta < 4, fit = TRUE)

summary(avg_mod)
test$avg_pred <- predict(avg_mod, newdata = test)
rmse_avg <- sqrt(mean((test$avg_pred - test$cases)^2))
cat("Averaged-model RMSE:", round(rmse_avg, 2), "\n")

# 10. Forecast plot -------------------------------------------------------
autoplot(fc) +
  autolayer(
    ts(test$cases,
       start = end(time(fc$mean)) - length(test$cases) + 1,
       frequency = 52),
    series = "Actual cases", linetype = "dashed"
  ) +
  labs(
    title = "WA pertussis: ARIMAX forecast from Google-Trends lags",
    y = "Weekly confirmed / probable cases",
    x = "Year"
  )

  # ---------- 11. Persist objects for the plotting script ----------
dir.create("artefacts", showWarnings = FALSE)

saveRDS(
  list(
    df_full      = lagged,      # the weekly data with cases + GT lags
    df_train     = train,
    df_test      = test,
    fc_arimax    = fc,          # forecast::forecast() output
    arimax_fit   = fit_arimax,
    gls_fit      = gls_fit,
    avg_fit      = avg_mod,
    metrics      = list(
      arimax_RMSE = accuracy(fc$mean,  test$cases)["Test set","RMSE"],
      gls_RMSE    = rmse_gls,
      avg_RMSE    = rmse_avg
    )
  ),
  file = "artefacts/wa_pertussis_models.Rds"
)