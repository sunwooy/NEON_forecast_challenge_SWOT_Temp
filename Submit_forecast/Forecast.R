install.packages(c("tidymodels", "doParallel", "ranger", "readr", "httr", "jsonlite"))

library(tidyverse)
library(lubridate)
library(slider)
library(dplyr)
library(purrr)
library(tidyr)
library(tidymodels)
library(doParallel)
library(ranger)
library(readr)
library(httr)
library(jsonlite)



# It's a random forest model that uses air temperature, a 3-day mean air temperature, and SWOT_WSE temperature as predictor for water temperature. 
my_model_id <- 'SWOT_AirTemp_RF'


targets <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")

# read in the sites data
aquatic_sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ci/refs/heads/main/neon4cast_field_site_metadata.csv") |>
  dplyr::filter(aquatics == 1)

focal_sites <- aquatic_sites |> 
  filter(field_site_subtype == 'Lake',
         field_site_id != 'TOOK') |> 
  pull(field_site_id)

# Filter the targets
targets <- targets %>%
  filter(site_id %in% focal_sites,
         variable == 'temperature')
targets %>%
  group_by(site_id) %>%
  summarise(max_val = max(datetime, na.rm = TRUE))


met_variables <- c("air_temperature")

# Past stacked weather -----
weather_past_s3 <- neon4cast::noaa_stage3()

weather_past <- weather_past_s3  |> 
  dplyr::filter(site_id %in% focal_sites,
                datetime >= ymd('2023-07-01'),
                variable %in% met_variables) |> 
  dplyr::collect()

# aggregate the past to mean values
weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |>
  group_by(site_id) |>
  mutate(
    airtemp_yday = slide_dbl(
      air_temperature,
      ~ mean(.x, na.rm = TRUE),
      .before = 2,
      .after = 0,
      .complete = TRUE
    )
  )


# # Future weather forecast --------
# # New forecast only available at 5am UTC the next day
forecast_date <- Sys.Date()
noaa_date <- forecast_date - days(1)

weather_future_s2 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))

weather_future <- weather_future_s2 |>
  dplyr::filter(datetime >= forecast_date,
                site_id %in% focal_sites,
                variable %in% met_variables) |>
  collect()

weather_future_daily <- weather_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))


forecast_df <- NULL
n_members = 310


# # Get horizon length by looking at available forecast dates that is not today. 
forecast_dates <- unique(weather_future_daily$datetime)
forecast_start_date <- forecast_date



# Function  to  get SWOT WSE
get_swot_wse <- function(site_ids=NULL,
                         start_date = "2023-07-11", # Science Orbit
                         end_date = Sys.Date()) {
  
  
  
  
  site_ids <- toupper(site_ids)
  
  # Mapping between NEON site IDs and PLD lake IDs
  PLD_add <- c(
    BARC = 7320214782,
    CRAM = 7250184002,
    LIRO = 7421110232,
    PRLA = 7420722382,
    PRPO = 7420733062,
    SUGG = 7320231232,
    TOOK = 8130554202
  )
  
  # If no site_ids provided, use all
  if (missing(site_ids) || is.null(site_ids) || length(site_ids) == 0) {
    site_ids <- names(PLD_add)
    message("No site_ids provided. Using all sites: ", paste(site_ids, collapse = ", "))
  } else {
    site_ids <- toupper(site_ids)
    
    # Validate inputs
    invalid_ids <- setdiff(site_ids, names(PLD_add))
    if (length(invalid_ids) > 0) {
      stop(
        "Invalid site ID(s): ",
        paste(invalid_ids, collapse = ", "),
        ". Valid options are: ",
        paste(names(PLD_add), collapse = ", ")
      )
    }
  }
  
  # Build site lookup table only for requested sites
  site_lookup <- tibble(
    field_site_id = site_ids,
    PLD_LakeID = unname(PLD_add[site_ids])
  )
  
  collection_name <- "SWOT_L2_HR_LakeSP_D"
  swot_vars <- "lake_id,time,wse,wse_u,quality_f,dark_frac,ice_clim_f,xtrk_dist,crid"
  
  base_time <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
  file_results <- list()
  
  for (i in seq_len(nrow(site_lookup))) {
    
    curr_site <- site_lookup$field_site_id[i]
    pld_id <- site_lookup$PLD_LakeID[i]
    
    api_url <- paste0(
      "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/",
      "timeseries?feature=PriorLake",
      "&feature_id=", pld_id,
      "&output=csv",
      "&start_time=", start_date, "T00:00:00Z",
      "&end_time=", as.character(end_date), "T23:59:59Z",
      "&collection_name=", collection_name,
      "&fields=", URLencode(swot_vars, reserved = TRUE)
    )
    
    result <- tryCatch({
      r <- GET(api_url, timeout(60))
      stop_for_status(r)
      
      response <- fromJSON(content(r, as = "text", encoding = "UTF-8"))
      
      if (!is.null(response$error)) {
        message("Error for site ", curr_site, ": ", response$error)
        return(NULL)
      }
      
      if (!is.null(response$status) &&
          response$status == "200 OK" &&
          !is.null(response$results$csv)) {
        
        csv_str <- response$results$csv
        
        if (!nzchar(trimws(csv_str))) {
          message("No data returned for site ", curr_site)
          return(NULL)
        }
        
        df <- read_csv(I(csv_str), show_col_types = FALSE)
        
        if ("time" %in% names(df)) {
          df <- df %>% filter(time >= 0)
        }
        
        if (nrow(df) == 0) {
          message("No valid rows for site ", curr_site)
          return(NULL)
        }
        
        unit_cols <- names(df)[grepl("units", names(df), ignore.case = TRUE)]
        if (length(unit_cols) > 0) {
          df <- df %>% select(-all_of(unit_cols))
        }
        
        if ("time" %in% names(df)) {
          df <- df %>%
            mutate(time = base_time + time)
        }
        
        df %>%
          mutate(site_id = curr_site) %>%
          rename(PLD_LakeID = lake_id, datetime = time) %>%
          distinct()
        
        
      } else {
        message("Error Message for site ", curr_site)
        NULL
      }
      
    }, error = function(e) {
      message("Failed to process site ", curr_site, ": ", e$message)
      NULL
    })
    
    if (!is.null(result)) {
      file_results[[length(file_results) + 1]] <- result
    }
  }
  
  if (length(file_results) == 0) {
    message("No valid data returned for requested site(s).")
    return(NULL)
  }
  
  bind_rows(file_results)
}




# For automated forecast, call WSE 

WSE <- get_swot_wse()

# Calculate baseline
baseline_wse <- WSE %>%
  mutate(datetime = as.Date(datetime)) %>%
  filter(
    site_id %in% focal_sites,
    quality_f < 2,
    datetime >= as.Date("2023-01-01"),
    datetime <= as.Date("2024-12-31")
  ) %>%
  group_by(site_id) %>%
  summarize(
    wse_ref_mean = mean(wse, na.rm = TRUE),
    .groups = "drop"
  )


WSE_filt <- WSE %>%
  mutate(datetime = as.Date(datetime)) %>%
  filter(
    site_id %in% focal_sites,
    quality_f < 2
  ) %>%
  left_join(baseline_wse, by = "site_id") %>%
  group_by(site_id) %>%
  arrange(datetime, .by_group = TRUE) %>%
  mutate(
    wse_anom = wse - wse_ref_mean
  ) %>%
  ungroup()





targets_model <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  arrange(site_id, datetime) |>
  group_by(site_id) |>
  left_join(weather_past_daily, by = c("datetime","site_id")) |>
  ungroup()|>
  drop_na() |>
  left_join(WSE_filt, by = c("datetime", "site_id"))



# To calculate maximum horizon, take longest forecast date possible from NOAA and subtract with the most recent SWOT observation. Take the max of that
last_forecast_date <- max(as.Date(forecast_dates))

site_recent_swot <- WSE_filt %>%
  filter(!is.na(wse),
         as.Date(datetime) < as.Date(forecast_date)) %>%
  group_by(site_id) %>%
  summarize(most_recent_swot = max(as.Date(datetime)), .groups = "drop") %>%
  mutate(max_h_needed = as.integer(last_forecast_date - most_recent_swot))

forecast_date_count <- max(site_recent_swot$max_h_needed, na.rm = TRUE)



targets_model_expanded <- targets_model %>%
  filter(!is.na(wse)) %>%
  mutate(row_id = row_number()) %>%
  rowwise() %>%
  mutate(
    block = list(
      tibble(
        site_id = site_id,
        SWOT_date = datetime,
        horizon = 0:forecast_date_count,
        datetime = SWOT_date + lubridate::days(0:forecast_date_count),
        wse_used = wse_anom,
        wse_u_used = wse_u
      ) %>%
        left_join(
          targets_model %>% select(-c(PLD_LakeID, wse, wse_u, quality_f, dark_frac, ice_clim_f, xtrk_dist, crid, wse_ref_mean, wse_anom)),
          by = c("site_id", "datetime")
        )
    )
  ) %>%
  ungroup() %>%
  select(row_id, block) %>%
  unnest(block) %>%
  drop_na()



df_model <- targets_model_expanded %>%
  select(horizon, wse_used, temperature, air_temperature, airtemp_yday, site_id)

df_model <- df_model %>%
  mutate(horizon = as.integer(horizon))

# Parallel setup
n_cores <- parallel::detectCores() - 1
cl <- makePSOCKcluster(max(1, n_cores))
registerDoParallel(cl)


# Function to tune + fit one horizon-specific model

fit_one_horizon_rf <- function(dat, horizon_value) {
  
  dat_h <- dat %>%
    filter(horizon == horizon_value) %>%
    select(-horizon)
  
  # Skip tiny groups
  if (nrow(dat_h) < 30) {
    return(list(
      horizon = horizon_value,
      status = "too_few_rows",
      tune_res = NULL,
      best_params = NULL,
      workflow_final = NULL
    ))
  }
  
  # Random split with 5 folds
  folds <- vfold_cv(dat_h, v = 5, strata=site_id)
  
  rec <- recipe(temperature ~ wse_used + air_temperature + airtemp_yday, data = dat_h)
  
  rf_spec <- rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = tune()
  ) %>%
    set_mode("regression") %>%
    set_engine("ranger")
  
  wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(rf_spec)
  
  # mtry depends on number of predictors
  param_set <- extract_parameter_set_dials(wf) %>%
    update(
      mtry = mtry(range = c(1L, 3L)),
      min_n = min_n(range = c(2L, 20L)),
      trees = trees(c(300L, 500L))
    )
  
  tune_res <- tune_grid(
    wf,
    resamples = folds,
    grid = 20,
    metrics = metric_set(rmse),
    param_info = param_set,
    control = control_grid(save_pred = TRUE)
  )
  
  best_params <- select_best(tune_res, metric = "rmse")
  
  wf_final <- finalize_workflow(wf, best_params) %>%
    fit(dat_h)
  
  list(
    horizon = horizon_value,
    status = "ok",
    tune_res = tune_res,
    best_params = best_params,
    workflow_final = wf_final
  )
}

# Train one model per horizon
horizons <- sort(unique(df_model$horizon))

rf_by_horizon <- map(horizons, ~ fit_one_horizon_rf(df_model, .x))
names(rf_by_horizon) <- paste0("h", horizons)

# Stop cluster
stopCluster(cl)
registerDoSEQ()



all_site_forecasts <- list()


for (i in 1:length(focal_sites)){
  curr_site <- focal_sites[i]
  
  site_target <- targets_model_expanded %>%
    filter(site_id == curr_site)
  
  
  
  # Predict on training data to calculate the sd of residual for process uncertanity
  train_preds <- list()
  
  for (h in horizons) {
    
    model_obj <- rf_by_horizon[[paste0("h", h)]]
    
    dat_h <- df_model %>%
      filter(horizon == h)
    
    if (is.null(model_obj) || model_obj$status != "ok") next
    
    mod <- predict(model_obj$workflow_final, new_data = dat_h) %>%
      pull(.pred)
    
    train_preds[[paste0("h", h)]] <- dat_h %>%
      mutate(prediction = mod,
             resid = temperature - prediction)
  }
  
  train_preds_df <- bind_rows(train_preds)
  
  horizon_resid_sd <- train_preds_df %>%
    group_by(horizon) %>%
    summarize(resid_sd = sd(resid, na.rm = TRUE), .groups = "drop")
  
  
  
  
  future_air <- weather_future_daily |>
    filter(site_id == curr_site,
           datetime >= forecast_date) |>
    arrange(parameter, datetime)
  
  
  # past + observed met used to create 3-day air temp lag 
  past_air <- weather_past_daily %>%
    mutate(datetime = as.Date(datetime)) %>%
    filter(
      site_id == curr_site,
      datetime >= forecast_date - days(2),
      datetime < forecast_date
    ) %>%
    select(datetime, site_id, air_temperature) %>%
    crossing(weather_future_daily %>% distinct(parameter)) %>%
    arrange(parameter, datetime)
  
  # combine past + future met and calc 3-day mean air temp
  air_for_lags <- bind_rows(
    past_air,
    future_air %>%
      select(datetime, site_id, air_temperature, parameter)
  ) %>%
    arrange(parameter, datetime) %>%
    group_by(parameter) %>%
    mutate(
      airtemp_yday = slide_dbl(
        air_temperature,
        ~ mean(.x, na.rm = TRUE),
        .before = 2,
        .after = 0,
        .complete = TRUE
      )
    ) %>%
    ungroup() %>%
    drop_na(airtemp_yday)
  
  
  
  
  # Now we're ready to forecast
  forecast_total_unc <- tibble(
    forecast_date = rep(forecast_dates, times = n_members),
    ensemble_member = rep(1:n_members, each = length(forecast_dates)),
    forecast_variable = "water_temperature",
    value = as.double(NA),
    uc_type = "total",
    site_id = curr_site)
  
  # We call WSE filt for last SWOT observation because training data is  only up  to date where we have w ater temp  measurement 
  last_wse_row <- WSE_filt %>%
    filter(site_id == curr_site,
           datetime < forecast_date, 
           !is.na(wse)) %>%
    arrange(desc(datetime)) %>%
    slice(1)
  
  wse_date <- last_wse_row$datetime
  wse_value <- last_wse_row$wse_anom
  
  
  
  for (j in seq_along(forecast_dates)){
    
    curr_date <- forecast_dates[j]
    
    h <- as.integer(as.Date(curr_date) - as.Date(wse_date)) 
    model_name <- paste0("h", h)
    
    
    # Skip if negative horizon OR if horizon expands beyond training done (too few row or too long since last SWOT)
    if (h < 0 || !(model_name %in% names(rf_by_horizon))) {
      message(curr_date, " skipped: no model name for h=", h)
      next
    }
    
    
    model_obj <- rf_by_horizon[[paste0("h", h)]]
    
    if (is.null(model_obj) || model_obj$status != "ok") {
      message(curr_date, " skipped: model status = ", model_obj$status, " for h=", h)
      next
    }
    
    air_lag <- air_for_lags %>%
      filter(datetime == curr_date)
    
    if (nrow(air_lag) == 0) {
      message(curr_date, " skipped: no air_lag rows")
      next
    }
    
    # build prediction frame expected by the RF
    pred_df <- air_lag %>%
      # rename(ensemble_member = parameter) %>%
      mutate(
        # We do this because weather future ensemble starts with 0 while ours start with 1
        ensemble_member = as.integer(parameter) + 1,
        wse_used = wse_value,
        horizon_since_last_SWOT = h
      ) %>%
      select(horizon_since_last_SWOT, wse_used, air_temperature, airtemp_yday, site_id, ensemble_member, datetime)
    
    
    pred_df <- pred_df %>%
      slice(rep(1:n(), length.out = n_members)) %>%
      mutate(
        ensemble_member = row_number()) 
    
    
    # RF prediction
    rf_pred <- predict(model_obj$workflow_final, new_data = pred_df) %>%
      pull(.pred)
    
    pred_df <- pred_df %>%
      mutate(value = rf_pred)
    
    
    # Add horizon specific residual noise
    resid_sd_h <- horizon_resid_sd %>%
      filter(horizon == h) %>%
      pull(resid_sd)
    
    if (length(resid_sd_h) == 1 && !is.na(resid_sd_h) && resid_sd_h > 0) {
      pred_df <- pred_df %>%
        mutate(value = value + rnorm(n(), mean = 0, sd = resid_sd_h))
    }
    
    
    
    forecast_total_unc <- forecast_total_unc %>%
      rows_update(
        pred_df %>%
          select(
            forecast_date = datetime,
            ensemble_member,
            site_id,
            value),
        by = c("forecast_date", "ensemble_member", "site_id"))
    
  }
  all_site_forecasts[[curr_site]] <- forecast_total_unc
}


forecast_df <- bind_rows(all_site_forecasts)


final_df <- forecast_df %>%
  transmute(
    datetime = forecast_date,
    parameter = ensemble_member,
    prediction= value,
    variable = "temperature",
    site_id = site_id
  )


forecast_df_EFI <- final_df %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = my_model_id,
         reference_datetime = forecast_date,
         family = 'ensemble',
         duration = 'P1D',
         parameter = as.character(parameter),
         project_id = 'neon4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction, model_id, project_id)



theme <- 'aquatics'
date <- forecast_df_EFI$reference_datetime[1]
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')

write_csv(forecast_df_EFI, forecast_file)

neon4cast::forecast_output_validator(forecast_file)


neon4cast::submit(forecast_file =  forecast_file, ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit



forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))

