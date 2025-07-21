############################################################################
# Washington Pertussis Dashboard – full-year-2024 split view
#   • actual cases              = red
#   • ARIMAX (TS + GT)           = blue
#   • ARIMA baseline (TS only)   = orange dashed
############################################################################

library(shiny)
library(plotly)
library(tidyverse)
library(lubridate)
library(zoo)
library(DT)
library(scales)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────
cases <- read_csv("wa_doh_formatted_pertussis_case_data.csv", show_col_types = FALSE) %>%
  mutate(date = floor_date(date, "week", week_start = 1)) %>%
  group_by(date) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

vaccine_trends <- readRDS("data/vaccine_trends.rds") %>%
  pivot_longer(-date, names_to = "term", values_to = "hits") %>%
  mutate(date = as.Date(date), term = str_replace_all(term, "_", " "))

models <- readRDS("artefacts/wa_pertussis_models.Rds")

trends <- models$df_full %>%
  select(date, starts_with("gt_")) %>%
  pivot_longer(-date, names_to = "term", values_to = "hits") %>%
  mutate(term = str_remove(term, "^gt_") %>% str_replace_all("_", " "))

# ── 2. UI ─────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Washington Pertussis Dashboard"),
  tabsetPanel(
    tabPanel("Search Trends + Cases", plotlyOutput("casesTrendsPlot", height = 400)),
    tabPanel("Model Accuracy",       plotlyOutput("forecastPlot",   height = 450),
             DTOutput("accuracyTable")),
    tabPanel("Vaccine Perception",   plotlyOutput("vaccinePlot",    height = 400),
             DTOutput("vaccineCorTable"))
  )
)

# ── 3. SERVER ────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # A. Cases + Google-search plot -----------------------------------------
  output$casesTrendsPlot <- renderPlotly({
    trend_comb <- trends %>%
      filter(term %in% c("whooping cough","whooping cough symptoms","pertussis")) %>%
      group_by(date) %>%
      summarise(hits = mean(hits, na.rm = TRUE), .groups = "drop")
    
    gp <- ggplot() +
      geom_col(data = cases, aes(x = date, y = cases), fill = "#999999", alpha = 0.6) +
      geom_line(data = trend_comb,
                aes(x = date, y = hits * max(cases$cases, na.rm = TRUE) / 100),
                colour = "red", linewidth = 1) +
      labs(title = "Pertussis Cases + Google Search Interest",
           y = "Cases (gray) + Scaled Search Interest (red)", x = NULL) +
      theme_minimal()
    ggplotly(gp)
  })
  
  # B. Forecast plot -------------------------------------------------------
  # B. Forecast plot -------------------------------------------------------
  output$forecastPlot <- renderPlotly({
    full24   <- models$df_full %>% filter(year(date) == 2024)
    gt_cols  <- grep("^gt_", names(full24), value = TRUE)
    gt_mean  <- rowMeans(full24[gt_cols], na.rm = TRUE)
    gt_scaled <- rescale(gt_mean, to = c(0, max(full24$cases, na.rm = TRUE)))
    
    full24 <- full24 %>%
      mutate(gt_scaled = gt_scaled)
    
    pred_arimax <- c(fitted(models$arimax_fit), as.numeric(models$fc_arimax$mean))
    pred_arima  <- c(fitted(models$fit_arima_base), as.numeric(models$fc_base$mean))
    idx24       <- models$df_full$date >= as.Date("2024-01-01") & models$df_full$date <= as.Date("2024-12-31")
    
    p <- ggplot(full24) +
      geom_area(aes(x = date, y = gt_scaled), fill = "#6baed6", alpha = 0.25) +
      geom_line(aes(x = date, y = gt_scaled, colour = "Google Trends (scaled)"), linewidth = 0.5) +
      geom_line(aes(x = date, y = cases, colour = "Actual cases"), linewidth = 1.1) +
      geom_line(aes(x = models$df_full$date[idx24], y = pred_arimax[idx24], colour = "ARIMAX (TS + Google Trends)"), linewidth = 0.9) +
      geom_line(aes(x = models$df_full$date[idx24], y = pred_arima[idx24], colour = "ARIMA (TS-only baseline)"), linewidth = 0.9, linetype = "dashed") +
      scale_color_manual(values = c(
        "Actual cases" = "red",
        "ARIMAX (TS + Google Trends)" = "blue",
        "ARIMA (TS-only baseline)" = "orange",
        "Google Trends (scaled)" = "#6baed6"
      )) +
      scale_x_date(limits = as.Date(c("2024-01-01", "2024-12-31")),
                   date_breaks = "1 month", date_labels = "%b %Y", expand = c(0, 0)) +
      labs(
        title = "Pertussis Forecast vs Actual — 2024 (three 4-month blocks)",
        subtitle = "ARIMA uses only past case counts • ARIMAX blends time-series with Google-Trends signals",
        y = "Weekly cases", x = NULL
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
    
    ggplotly(p, tooltip = c("x", "y", "colour"))
  })
  
  # C. RMSE table ----------------------------------------------------------
  output$accuracyTable <- renderDT({
    full24 <- models$df_full %>% filter(year(date) == 2024)
    pred_arimax <- c(fitted(models$arimax_fit), as.numeric(models$fc_arimax$mean))
    pred_arima  <- c(fitted(models$fit_arima_base), as.numeric(models$fc_base$mean))
    idx24 <- models$df_full$date >= as.Date("2024-01-01") & models$df_full$date <= as.Date("2024-12-31")
    
    rmse <- function(from, to, label){
      m <- full24$date >= from & full24$date <= to
      tibble(Block = label,
             `RMSE ARIMA (TS)`  = sqrt(mean((pred_arima [idx24][m] - full24$cases[m])^2)),
             `RMSE ARIMAX (TS+GT)` = sqrt(mean((pred_arimax[idx24][m] - full24$cases[m])^2)))
    }
    
    tbl <- bind_rows(
      rmse(as.Date("2024-01-01"), as.Date("2024-04-30"), "Jan–Apr"),
      rmse(as.Date("2024-05-01"), as.Date("2024-08-31"), "May–Aug"),
      rmse(as.Date("2024-09-01"), as.Date("2024-12-31"), "Sep–Dec")
    ) %>% mutate(across(starts_with("RMSE"), round, 1))
    
    datatable(tbl, rownames = FALSE, options = list(dom = "t"),
              caption = htmltools::HTML(
                "<b>RMSE by 4-month block (lower = better accuracy)</b>"
              ))
  })
  
  # D. Vaccine-search plot -------------------------------------------------
  output$vaccinePlot <- renderPlotly({
    vax_plot_df <- vaccine_trends %>%
      mutate(term = paste("Search –", term)) %>%
      bind_rows(cases %>% rename(hits = cases) %>% mutate(term = "Actual cases"))
    
    gp <- ggplot(vax_plot_df,
                 aes(x = date, y = hits, colour = term)) +
      geom_line(linewidth = 1) +
      scale_x_date(breaks = seq(min(vax_plot_df$date), max(vax_plot_df$date), by = "3 months"),
                   date_labels = "%b '%y", expand = c(0.01,0)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Vaccine-related Searches & Pertussis Cases",
           y = "Google Trends index / Case count", x = NULL) +
      theme_minimal()
    
    ggplotly(gp, tooltip = c("y","colour")) %>%
      layout(xaxis = list(tickangle = -45, tickformat = "%b '%y"))
  })
  
  # E. Correlation table ---------------------------------------------------
  output$vaccineCorTable <- renderDT({
    wide <- vaccine_trends %>% pivot_wider(names_from = term, values_from = hits) %>% inner_join(cases, by = "date")
    cor_df <- wide %>% summarise(across(-c(date, cases), ~ cor(.x, wide$cases, use = "complete.obs"))) %>%
      pivot_longer(everything(), names_to = "Search Term", values_to = "Correlation") %>%
      arrange(desc(abs(Correlation))) %>% mutate(Correlation = round(Correlation, 3))
    
    datatable(cor_df, rownames = FALSE, options = list(dom = "tp"),
              caption = "Correlation of Vaccine Search Terms with Weekly Cases")
  })
}

shinyApp(ui, server)