#!/usr/bin/env Rscript
# ----------------------------------------------------
# 02_plot_metrics.R — Surge-era dashboard + GT overlay
# ----------------------------------------------------
# • Line plot: Actual, ARIMAX, ARIMA baseline (zoom 2024)
# • Bar chart: RMSE (ARIMA vs ARIMAX) for held-out 2023
# • NEW: two-line plot of cases vs Google-Trends (2014-24)
# ----------------------------------------------------

# 0. Packages ---------------------------------------------------------------
libs <- c("tidyverse", "lubridate", "glue", "patchwork",
          "scales", "forecast")
need <- setdiff(libs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(libs, library, character.only = TRUE)))

# 1. Load artefacts ---------------------------------------------------------
art <- readRDS("artefacts/wa_pertussis_models.Rds")
list2env(art, .GlobalEnv)            # df_*, arimax_fit, fit_arima_base, …

dir.create("artefacts", showWarnings = FALSE)

# 2. Long dataframe for model-comparison line plot -------------------------
arimax_train_df <- tibble(date = df_train$date,
                          series = "ARIMAX",
                          value  = as.numeric(fitted(arimax_fit)))
arimax_test_df  <- tibble(date = df_test$date,
                          series = "ARIMAX",
                          value  = as.numeric(fc_arimax$mean))

baseline_train_df <- tibble(date = df_train$date,
                            series = "ARIMA baseline",
                            value  = as.numeric(fitted(fit_arima_base)))
baseline_test_df  <- tibble(date = df_test$date,
                            series = "ARIMA baseline",
                            value  = as.numeric(fc_base$mean))

vis_long <- bind_rows(
  df_full |> mutate(series = "Actual", value = cases),
  arimax_train_df, arimax_test_df,
  baseline_train_df, baseline_test_df
)

vis_zoom <- vis_long %>%
  filter(date >= as.Date("2024-01-01") & date <= as.Date("2025-05-31")) %>%
  mutate(series = factor(series,
                         levels = c("Actual", "ARIMAX", "ARIMA baseline")))

# 3. Line plot: observed vs models (zoom) ----------------------------------
col_map <- c("Actual" = "black",
             "ARIMAX" = "#b2182b",
             "ARIMA baseline" = "gray40")
lt_map  <- c("Actual" = "solid",
             "ARIMAX" = "solid",
             "ARIMA baseline" = "dashed")

p_cases <- ggplot(vis_zoom,
                  aes(date, value, colour = series, linetype = series)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_linetype_manual(values = lt_map, guide = "none") +
  scale_y_continuous("Weekly pertussis cases", labels = comma) +
  labs(title    = "WA pertussis: observed vs. ARIMAX & baseline",
       subtitle = glue("Training window: {min(df_train$date)} – {max(df_train$date)}"),
       x        = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# 4. Bar chart — held-out 2023 RMSE ----------------------------------------
rmse_2023 <- year_rmse_table %>%
  filter(year == 2023) %>%
  select(RMSE_ARIMA, RMSE_ARIMAX) %>%
  pivot_longer(everything(), names_to = "Model", values_to = "RMSE") %>%
  mutate(Model = recode(Model,
                        RMSE_ARIMA  = "ARIMA baseline",
                        RMSE_ARIMAX = "ARIMAX (2023)"))

pad <- max(rmse_2023$RMSE) * 0.15

p_tbl <- ggplot(rmse_2023, aes(Model, RMSE, label = round(RMSE, 2))) +
  geom_col(fill = "#a6cee3", width = 0.6) +
  geom_text(vjust = -0.25, size = 4.5) +
  coord_cartesian(ylim = c(0, max(rmse_2023$RMSE) + pad)) +
  labs(title = "RMSE — leave-2023-out evaluation") +
  theme_minimal(base_size = 13) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1))

# 5. NEW two-line plot: cases vs Google-Trends ------------------------------
# Pick one unlagged GT column – strongest is whooping-cough search.
plot_df <- df_full %>%
  select(date,
         cases,
         search_hits = gt_whooping_cough_l0) %>%        # change term if desired
  pivot_longer(cols = -date, names_to = "series", values_to = "value")

# Rescale GT so both lines share a primary axis
scale_factor <- max(plot_df$value[plot_df$series == "cases"], na.rm = TRUE) /
                max(plot_df$value[plot_df$series == "search_hits"], na.rm = TRUE)

plot_df <- plot_df %>%
  mutate(value_scaled = ifelse(series == "search_hits",
                               value * scale_factor, value))

p_gt_cases <- ggplot(plot_df,
                     aes(date, value_scaled, colour = series)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(cases       = "black",
                                 search_hits = "#377eb8"),
                      labels = c(cases       = "Pertussis cases",
                                 search_hits = "Google Trends: whooping cough"),
                      name = NULL) +
  scale_y_continuous("Weekly pertussis cases",
                     labels = comma,
                     sec.axis = sec_axis(~ . / scale_factor,
                                         name = "GT hits (0–100)")) +
  labs(title = "Pertussis cases vs Google-Trends search interest (WA, 2014-2024)",
       x = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Save the two-line plot separately
ggsave("artefacts/wa_cases_vs_gt.pdf", p_gt_cases, width = 11, height = 4)
message("✔ Two-line plot written to artefacts/wa_cases_vs_gt.pdf")

# 6. Dashboard (zoom plot + bar chart) -------------------------------------
dashboard <- (p_cases / p_tbl) + plot_layout(heights = c(3, 1))
ggsave("artefacts/wa_pertussis_summary.pdf",
       dashboard, width = 11, height = 7)
message("✔ Dashboard written to artefacts/wa_pertussis_summary.pdf")