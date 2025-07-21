![PRISM Banner](images/PRISM.png)

# PRISM: Pertussis Response via Internet Search Modeling.

*Forecasting pertussis outbreaks in Washington State using Google Trends and time-series modeling*

Research Report: https://docs.google.com/document/d/1bL3nd01IZx_GCCfaz8TaJ-O3YLieJIgtnP-QRV-Lp80/edit?usp=sharing 

This repository contains the code, data, and Shiny dashboard for the project “Pertussis Response via Internet Search Modeling”, which evaluates whether Google Trends data can improve the timeliness and accuracy of pertussis outbreak detection in Washington State. The project builds and evaluates ARIMA, ARIMAX, and elastic-net models by combining Google Trends search data with confirmed pertussis case reports. Tested so far on a statewide level, it includes a fully reproducible pipeline and an interactive R Shiny dashboard for public health use. To run the analysis or dashboard locally, clone the repository and install the required packages listed in each script (see below).

## Code - run in this order:

- **run_model_pipeline.R**: Main script that runs the full modeling pipeline. It pulls Google Trends data, processes pertussis case reports, computes lagged correlations, fits ARIMA/ARIMAX/elastic-net models, and outputs performance tables, lead-time calculations, and predictive artefacts used in the dashboard.

- **correlation_analysis.R**: Runs the full-model pipeline using an expanded set of Google Trends search terms. It performs correlation screening across all lag combinations, fits ARIMA and ARIMAX models using top-ranked predictors, and exports both evaluation artefacts and formatted publication tables. Outputs include a comprehensive correlation table (all_gt_correlations_full.csv), forecast accuracy results, and predictor summaries. Ideal for scanning additional GT features and reproducing full-term model results reported in the paper. 

- **vaccine_trends_pull.R**: For use in the dashboard. Downloads and processes weekly Google Trends data for vaccine-related search terms (e.g. "whooping cough vaccine," "pertussis vaccine") in Washington State. Cleans and aggregates the results by week; saves results in both .csv and .rds formats for use in the app.R dashboard visualization. Includes error-handling and rate-limit protection.

- **app.R**: Launches the interactive Pertussis Forecast Dashboard built in R Shiny. Displays weekly pertussis case counts alongside Google Trends interest in relevant symptoms, search-based forecasts (ARIMA vs. ARIMAX), vaccine-related search trends, and real-time model performance metrics (e.g., RMSE by block). Designed for public health officials to monitor trends, interpret forecasts, and track vaccine sentiment over time. To launch the dashboard locally, run `shiny::runApp("app.R")`

The dashboard is also viewable via this link: [https://annafelten.shinyapps.io/PRISM/](https://annafelten.shinyapps.io/PRISM/)

## Data Files

- **wa_doh_formatted_pertussis_data.csv** - Weekly confirmed and probable pertussis case counts from the Washington State Department of Health (2022–2024).
- **data/** - Directory for storing weekly updates, raw files, or preprocessed data that is used by the dashboard.

## Output & Misc

- **artefacts/** - Folder for model objects, saved plots, and serialized .Rds files used in the dashboard.

## Dependencies
This project requires R ≥ 4.2 and the following packages:

- `tidyverse`
- `forecast`
- `gtrendsR`
- `shiny`
- `lubridate`
- `zoo`
- `glmnet`
- `plotly`
- `DT`

## Quick Start

1. **Run the main modeling pipeline**  
   Execute `run_model_pipeline.R` to generate forecasts, train models, and output artefacts (model accuracy tables, saved predictors, forecast objects).

2. **Download vaccine-related Google Trends data**  
   Run `vaccine_trends_pull.R` to fetch and clean weekly search data for vaccine terms. Results are saved as both `.csv` and `.rds` files.

3. **Launch the interactive dashboard**  
   Start the Shiny app by running:
   ```r
   shiny::runApp("app.R")
