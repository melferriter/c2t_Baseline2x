# ============================
# 0. Libraries
# ============================
library(lidR)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(caret)
library(tibble)
library(lmerTest)
library(purrr)
library(yardstick)
library(kableExtra)
library(ggplot2)

rm(list = ls(globalenv()))

# ============================
# 1. File definitions
# ============================
files <- list(
  Fusion = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv",
    idcol  = "SFMcombo"
  ),
  LiDAR = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_side3x.csv",
    idcol  = "Side3x"
  ),
  SfM = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv",
    idcol  = "SFMonly"
  )
)

dtm  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

# ============================
# 2. CV setup
# ============================
set.seed(444)
ctrl <- trainControl(method = "cv", number = 5, repeats = 10, savePredictions = "final")

# ============================
# 3. Summarizer
# ============================
summarize_rf <- function(model, label) {
  preds <- model$pred %>%
    mutate(obs_cm = exp(obs), pred_cm = exp(pred))

  folds <- preds %>%
    group_by(Resample) %>%
    summarise(
      R2   = ifelse(var(obs_cm, na.rm = TRUE) == 0 | var(pred_cm, na.rm = TRUE) == 0,
                    NA_real_, cor(obs_cm, pred_cm, use = "complete.obs")^2),
      RMSE = rmse_vec(obs_cm, pred_cm),
      MAE  = mae_vec(obs_cm, pred_cm),
      Bias = mean(pred_cm - obs_cm, na.rm = TRUE),
      .groups = "drop"
    )

  folds %>%
    summarise(
      R2_mean   = mean(R2, na.rm = TRUE),   R2_sd   = sd(R2, na.rm = TRUE),
      RMSE_mean = mean(RMSE, na.rm = TRUE), RMSE_sd = sd(RMSE, na.rm = TRUE),
      MAE_mean  = mean(MAE, na.rm = TRUE),  MAE_sd  = sd(MAE, na.rm = TRUE),
      Bias_mean = mean(Bias, na.rm = TRUE), Bias_sd = sd(Bias, na.rm = TRUE)
    ) %>%
    mutate(Source = label)
}



summarize_cv_lm <- function(data, formula, label, k = 5, repeats = 1) {
  all_metrics <- list()
  
  for (r in 1:repeats) {
    folds <- createFolds(data$dbh_cm, k = k)
    
    metrics <- map_dfr(folds, function(fold) {
      train_data <- data[-fold, ]
      test_data  <- data[fold, ]
      
      if (nrow(test_data) < 2) return(NULL)   # skip if too small
      
      model <- lm(formula, data = train_data)
      preds <- exp(predict(model, newdata = test_data))
      obs   <- test_data$dbh_cm
      
      if (length(unique(obs)) < 2) return(NULL) # skip degenerate folds
      
      tibble(
        R2   = cor(obs, preds, use = "complete.obs")^2,
        RMSE = sqrt(mean((obs - preds)^2, na.rm = TRUE)),
        MAE  = mean(abs(obs - preds), na.rm = TRUE),
        Bias = mean(preds - obs, na.rm = TRUE)
      )
    })
    all_metrics[[r]] <- metrics
  }
  
  bind_rows(all_metrics) %>%
    summarise(
      R2_mean   = mean(R2, na.rm = TRUE),   R2_sd   = sd(R2, na.rm = TRUE),
      RMSE_mean = mean(RMSE, na.rm = TRUE), RMSE_sd = sd(RMSE, na.rm = TRUE),
      MAE_mean  = mean(MAE, na.rm = TRUE),  MAE_sd  = sd(MAE, na.rm = TRUE),
      Bias_mean = mean(Bias, na.rm = TRUE), Bias_sd = sd(Bias, na.rm = TRUE)
    ) %>%
    mutate(Source = label)
}

# ---- Helpers to grab CV predictions ----
collect_rf_preds <- function(model, data, label, flight) {
  # caret saved predictions on log-scale; convert and join species by rowIndex
  stopifnot(!is.null(model$pred$rowIndex))
  p <- model$pred %>%
    mutate(
      obs_cm  = exp(obs),
      pred_cm = exp(pred),
      Species = data$Species[rowIndex],
      Model   = label,
      Flight  = flight
    ) %>%
    dplyr::select(Flight, Model, Species, obs_cm, pred_cm)
  p
}

evaluate_lmer_cv_with_preds <- function(data, k = 5, repeats = 1, label, flight) {
  all_metrics <- list()
  all_preds   <- list()

  for (r in 1:repeats) {
    folds <- createFolds(data$dbh_cm, k = k, list = TRUE, returnTrain = FALSE)

    fold_metrics <- purrr::map_dfr(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train <- data[-test_idx, ]
      test  <- data[test_idx, ]

      model <- lmer(
        log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) + (1 | Species),
        data = train
      )

      preds <- tibble(
        obs_cm  = test$dbh_cm,
        pred_cm = exp(predict(model, newdata = test, allow.new.levels = TRUE)),
        Species = test$Species,
        Model   = label,
        Flight  = flight
      )
      all_preds[[length(all_preds)+1]] <<- preds

      tibble(
        R2    = cor(preds$obs_cm, preds$pred_cm, use = "complete.obs")^2,
        RMSE  = sqrt(mean((preds$obs_cm - preds$pred_cm)^2, na.rm = TRUE)),
        MAE   = mean(abs(preds$obs_cm - preds$pred_cm), na.rm = TRUE),
        Bias  = mean(preds$pred_cm - preds$obs_cm, na.rm = TRUE)
      )
    })
    all_metrics[[r]] <- fold_metrics
  }

  metrics <- bind_rows(all_metrics) %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    ) %>% mutate(Source = label)

  list(metrics = metrics, preds = bind_rows(all_preds))
}


# ============================
# 4. Main loop
# ============================
all_results <- list()

for (flight in names(files)) {
  message("Processing: ", flight)

  f <- files[[flight]]

  # ---- Load data ----
  las     <- readLAS(f$las)
  crowns  <- st_read(f$crowns, quiet = TRUE) %>% mutate(treeID = as.character(treeID))
  field   <- readr::read_csv(f$field)

  id_map <- list(
  Fusion = "SFMcombo",
  LiDAR  = "Side3x",
  SfM    = "SFMonly"
)

field <- readr::read_csv(f$field) %>%
  # Rename the correct ID column to `treeID`
  rename(treeID = !!sym(id_map[[flight]])) %>%
  mutate(treeID = as.character(str_trim(treeID)))

# --- Standardize lidar height column name ---
if ("SFMcomb_ht" %in% names(field)) {
  field <- field %>% rename(lidar_height = SFMcomb_ht)
}
if ("Side3x_ht" %in% names(field)) {
  field <- field %>% rename(lidar_height = Side3x_ht)
}
if ("SFMonly_ht" %in% names(field)) {
  field <- field %>% rename(lidar_height = SFMonly_ht)
}

# --- Add species flags ---
field <- field %>%
  mutate(
    Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE), "Sequoia", "Other"),
    is_sequoia = ifelse(Species == "Sequoia", 1, 0)
  ) %>%
  dplyr::select(treeID, dbh_cm, treecbh, lidar_height, Species, is_sequoia)

message("Columns in field after cleanup: ", paste(names(field), collapse = ", "))


# ✅ Check that dbh_cm really exists here
message("Columns in field after cleanup: ", paste(names(field), collapse = ", "))

# ✅ Check if any dbh_cm values are missing
if (all(is.na(field$dbh_cm))) {
  stop(paste("All dbh_cm values are NA in flight", flight))
}

 

  # ---- Normalize LAS ----
  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(raster(dtm), opt = "slope", unit = "degrees")

  # ---- Crown metrics ----
  metrics_list <- map_dfr(1:nrow(crowns), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)
    if (is.null(las_clip) || npoints(las_clip) == 0) return(tibble(treeID = crown_poly$treeID))
    Z <- las_clip@data$Z
    tibble(
      treeID = crown_poly$treeID,
      point_density = npoints(las_clip) / as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope, crown_poly, "mean"),
      max_z  = max(Z, na.rm = TRUE),
      mean_z = mean(Z, na.rm = TRUE),
      p25    = quantile(Z, 0.25, na.rm = TRUE),
      p50    = quantile(Z, 0.50, na.rm = TRUE),
      p75    = quantile(Z, 0.75, na.rm = TRUE),
      p95    = quantile(Z, 0.95, na.rm = TRUE),
      sd_z   = sd(Z, na.rm = TRUE),
      cv_z   = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE),
      ndvi_mean    = exact_extract(ndvi, crown_poly, "mean"),
      crown_area_m2 = as.numeric(st_area(crown_poly))
    )
  })
  
# --- Step 1: Build crowns metrics (only LiDAR-derived stuff) ---
crowns_metrics <- crowns %>%
  dplyr::select(treeID, geom) %>%
  left_join(metrics_list, by = "treeID") %>%
  dplyr::select(
    treeID, crown_area_m2, point_density, slope_mean, ndvi_mean,
    max_z, mean_z, p25, p50, p75, p95, sd_z, cv_z, geom
  )
# Notice: NO dbh_cm here

# --- Step 2: Clean field dataset (only ground-truth stuff) ---
field_clean <- field %>%
  mutate(
    treeID = as.character(treeID),
    Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE), "Sequoia", "Other"),
    is_sequoia = ifelse(Species == "Sequoia", 1, 0)
  ) %>%
  dplyr::select(treeID, dbh_cm, treecbh, lidar_height , Species, is_sequoia)

# --- Step 3: Join clean crowns + clean field ---
crowns_final <- crowns_metrics %>%
  inner_join(field_clean, by = "treeID") %>%
  filter(Species != "intact snag", treecbh > 0, dbh_cm > 0)


# Ensure dbh_cm exists after join
stopifnot("dbh_cm" %in% names(crowns_final))

crown_data <- crowns_final %>%
  filter(!is.na(dbh_cm), dbh_cm > 0, dbh_cm < 350) %>%
  mutate(
    log_dbh        = log(dbh_cm),
    log_height     = log(lidar_height),
    log_crown_area = log1p(crown_area_m2),
    crown_depth    = pmax(lidar_height - treecbh, 0)   # <-- add this
  ) %>%
  drop_na()



  # ---- Run baseline RF ----
    # 1. Baseline RF
    rf_base <- train(
      log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
        p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
        point_density + slope_mean + treecbh + ndvi_mean,
      data = crown_data,
      method = "rf", trControl = ctrl, tuneLength = 5, ntree = 1000
    )
    res_rf_base <- summarize_rf(rf_base, paste(flight, "RF baseline"))

    # 2. Mixed-effects regression
    # ============================
# Function: Cross-validated Mixed-effects Regression
# ============================
evaluate_lmer_cv <- function(data, k = 5, repeats = 1) {
  all_metrics <- list()
  
  for (r in 1:repeats) {
    folds <- createFolds(data$dbh_cm, k = k, list = TRUE, returnTrain = FALSE)
    
    fold_metrics <- purrr::map_dfr(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train <- data[-test_idx, ]
      test  <- data[test_idx, ]
      
      model <- lmer(
        log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) + (1 | Species),
        data = train
      )
      
      preds <- tibble(
        obs_cm  = test$dbh_cm,
        pred_cm = exp(predict(model, newdata = test, allow.new.levels = TRUE))
      )
      
      tibble(
        R2    = cor(preds$obs_cm, preds$pred_cm, use = "complete.obs")^2,
        RMSE  = sqrt(mean((preds$obs_cm - preds$pred_cm)^2, na.rm = TRUE)),
        MAE   = mean(abs(preds$obs_cm - preds$pred_cm), na.rm = TRUE),
        Bias  = mean(preds$pred_cm - preds$obs_cm, na.rm = TRUE)
      )
    })
    
    all_metrics[[r]] <- fold_metrics
  }
  
  # Pool across repeats
  bind_rows(all_metrics) %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    )
}


    mix_results <- evaluate_lmer_cv(crown_data, k = 5, repeats = 10) %>%
  mutate(Source = paste(flight, "Mixed-effects regression"))


    # 3. Sequoia-only RF
    rf_seq <- train(
      log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
        p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
        point_density + slope_mean + treecbh + ndvi_mean,
      data = filter(crown_data, is_sequoia == 1),
      method = "rf", trControl = ctrl, tuneLength = 5, ntree = 1000
    )
    res_rf_seq <- summarize_rf(rf_seq, paste(flight, "Sequoia RF"))

    # 4. Non-sequoia-only RF
    rf_nonseq <- train(
      log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
        p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
        point_density + slope_mean + treecbh + ndvi_mean,
      data = filter(crown_data, is_sequoia == 0),
      method = "rf", trControl = ctrl, tuneLength = 5, ntree = 1000
    )
    res_rf_nonseq <- summarize_rf(rf_nonseq, paste(flight, "Non-sequoia RF"))

    # 5. Log–log allometry (Sequoia)
seq_data <- crown_data %>%
  filter(is_sequoia == 1, crown_depth > 0)

res_loglog_seq <- summarize_cv_lm(
  seq_data,
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2) + log(crown_depth),
  paste(flight, "Sequoia log–log"),
  k = 5, repeats = 10
)

res_loglog_nonseq <- summarize_cv_lm(
  filter(crown_data, is_sequoia == 0),
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2),
  paste(flight, "Non-sequoia log–log"),
  k = 5, repeats = 10
)

# After rf_base is trained:
preds_rf_base <- collect_rf_preds(rf_base, crown_data, "RF (pooled)", flight)

# Mixed-effects with predictions:
mix_cv <- evaluate_lmer_cv_with_preds(crown_data, k = 5, repeats = 10,
                                      label = "Mixed-effects (pooled)", flight = flight)
mix_results <- mix_cv$metrics %>% mutate(Source = paste(flight, "Mixed-effects regression"))
preds_mix   <- mix_cv$preds
preds_store <- bind_rows(preds_rf_base, preds_mix)
if (!exists("all_preds_all")) {
  all_preds_all <- list()
}
all_preds_all[[flight]] <- preds_store

    # ---- Collect ----
    all_results[[flight]] <- bind_rows(
      res_rf_base, mix_results,
      res_rf_seq, res_rf_nonseq,
      res_loglog_seq, res_loglog_nonseq
    )
# ---- Sequoia-only Random Forest ----
rf_seq <- train(
  log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + treecbh + ndvi_mean,
  data = filter(crown_data, is_sequoia == 1),
  method = "rf", trControl = ctrl, tuneLength = 5, ntree = 1000
)
res_rf_seq <- summarize_rf(rf_seq, paste(flight, "Sequoia RF"))

preds_rf_seq <- collect_rf_preds(rf_seq, filter(crown_data, is_sequoia == 1),
                                 "RF (Sequoia-only)", flight)


# ---- Sequoia log–log ----

seq_data <- crown_data %>% filter(is_sequoia == 1, crown_depth > 0)

# Cross-validated metrics (kept for reporting)
res_loglog_seq <- summarize_cv_lm(
  seq_data,
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2) + log(crown_depth),
  paste(flight, "Sequoia log–log"),
  k = 5, repeats = 10
)

# Predictions for plotting
model_loglog_seq <- lm(
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2) + log(crown_depth),
  data = seq_data
)

preds_loglog_seq <- seq_data %>%
  mutate(
    obs_cm = dbh_cm,
    pred_cm = exp(predict(model_loglog_seq, newdata = seq_data)),
    residual_cm = pred_cm - obs_cm,
    Model = "Log–log (Sequoia-only)",
    Flight = flight
  ) %>%
  dplyr::select(Flight, Model, Species, obs_cm, pred_cm, residual_cm)



# ---- Non-sequoia Random Forest ----
rf_nonseq <- train(
  log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + treecbh + ndvi_mean,
  data = filter(crown_data, is_sequoia == 0),
  method = "rf", trControl = ctrl, tuneLength = 5, ntree = 1000
)
res_rf_nonseq <- summarize_rf(rf_nonseq, paste(flight, "Non-sequoia RF"))

preds_rf_nonseq <- collect_rf_preds(rf_nonseq, filter(crown_data, is_sequoia == 0),
                                    "RF (Non-sequoia)", flight)


# ---- Non-sequoia log–log ----
nonseq_data <- crown_data %>% filter(is_sequoia == 0)

# Cross-validated metrics (kept for reporting)
res_loglog_nonseq <- summarize_cv_lm(
  nonseq_data,
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2),
  paste(flight, "Non-sequoia log–log"),
  k = 5, repeats = 10
)

# Predictions for plotting
model_loglog_nonseq <- lm(
  log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2),
  data = nonseq_data
)

preds_loglog_nonseq <- nonseq_data %>%
  mutate(
    obs_cm = dbh_cm,
    pred_cm = exp(predict(model_loglog_nonseq, newdata = nonseq_data)),
    residual_cm = pred_cm - obs_cm,
    Model = "Log–log (Non-sequoia)",
    Flight = flight
  ) %>%
  dplyr::select(Flight, Model, Species, obs_cm, pred_cm, residual_cm)


  preds_store <- bind_rows(
    preds_rf_base, preds_mix,
    preds_rf_seq, preds_loglog_seq,
    preds_rf_nonseq, preds_loglog_nonseq
  )
  
  # Store them in master list for this flight
  all_preds_all[[flight]] <- preds_store
  
  # Collect metrics
  all_results[[flight]] <- bind_rows(
    res_rf_base, mix_results,
    res_rf_seq, res_rf_nonseq,
    res_loglog_seq, res_loglog_nonseq
  )
}


# ============================
# 5. Final results
# ============================
final_results <- bind_rows(all_results)
print(final_results)


# === Assemble predictions once ===
preds_all <- bind_rows(all_preds_all) %>%
  filter(is.finite(obs_cm), is.finite(pred_cm)) %>%
  mutate(
    residual_cm = pred_cm - obs_cm,
    # Optional: control facet order
    Model = factor(
      Model,
      levels = c("Mixed-effects (pooled)",
                 "RF (pooled)",
                 "RF (Sequoia-only)",
                 "Log–log (Sequoia-only)",
                 "RF (Non-sequoia)",
                 "Log–log (Non-sequoia)")
    )
  )

# After you've combined all predictions into preds_all
# (make sure preds_all has treeID, obs_cm, pred_cm, Species, Model, Flight)
preds_all_collapsed <- preds_all %>%
  group_by(Flight, Model, Species, obs_cm) %>%
  summarise(
    pred_cm = mean(pred_cm, na.rm = TRUE),
    residual_cm = mean(residual_cm, na.rm = TRUE),
    .groups = "drop"
  )



---------
library(dplyr)
library(kableExtra)

# Function to format mean ± SD nicely
fmt_pm <- function(mean, sd, digits = 3) {
  ifelse(!is.na(mean),
         sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean, sd),
         NA)
}

# Build a compact summary table
make_summary_table <- function(results, flight_name) {
  results %>%
    filter(grepl(flight_name, Source)) %>%
    mutate(
      R2 = fmt_pm(R2_mean, R2_sd, 3),
      RMSE = fmt_pm(RMSE_mean, RMSE_sd, 1),
      MAE = fmt_pm(MAE_mean, MAE_sd, 1),
      Bias = fmt_pm(Bias_mean, Bias_sd, 2)
    ) %>%
    dplyr::select(Source, R2, RMSE, MAE, Bias) %>%
    rename(`Error Metrics` = Source,
           `R² (mean ± SD)` = R2,
           `RMSE (cm)` = RMSE,
           `MAE (cm)` = MAE,
           `Bias (cm)` = Bias) %>%
    kbl(caption = paste0(flight_name, " Models"),
        align = "lcccc") %>%
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
}

# Example: make tables for each flight
make_summary_table(final_results, "Fusion")
make_summary_table(final_results, "LiDAR")
make_summary_table(final_results, "SfM")


-----------
#figures 
  
library(ggplot2)

# Filter to LiDAR 3x flight (adjust if your label differs)
preds_lidar <- preds_all %>% dplyr::filter(Flight == "LiDAR")

# ---- Figure 1: Predicted vs Observed ----
fig1 <- ggplot(preds_lidar, aes(x = obs_cm, y = pred_cm, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 2) +
  coord_equal() +
  facet_wrap(~ Model, nrow = 1) +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Predicted vs. Observed DBH (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# ---- Figure 2: Residuals vs Predicted ----
preds_lidar <- preds_lidar %>%
  mutate(residual_cm = pred_cm - obs_cm)

fig2 <- ggplot(preds_lidar, aes(x = pred_cm, y = residual_cm, color = Species)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~ Model, nrow = 1) +
  labs(x = "Predicted DBH (cm)", y = "Residual (Pred − Obs, cm)",
       title = "Residuals vs. Predicted DBH (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# Residual calculation
preds_all <- preds_all %>%
  mutate(residual_cm = pred_cm - obs_cm)

# ---- Figure: Predicted vs Observed ----
fig_species_scatter <- ggplot(preds_all %>% filter(Flight == "LiDAR"),
                              aes(x = obs_cm, y = pred_cm, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 2) +
  coord_equal() +
  facet_wrap(~ Model, ncol = 2) +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Predicted vs. Observed DBH by Model (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# ---- Figure: Residuals vs Predicted ----
fig_species_resid <- ggplot(preds_all %>% filter(Flight == "LiDAR"),
                            aes(x = pred_cm, y = residual_cm, color = Species)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~ Model, ncol = 2) +
  labs(x = "Predicted DBH (cm)", y = "Residual (Pred − Obs, cm)",
       title = "Residuals vs. Predicted DBH by Model (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# Save
ggsave("fig_species_scatter.png", fig_species_scatter, width = 8, height = 6, dpi = 300)
ggsave("fig_species_resid.png",   fig_species_resid, width = 8, height = 6, dpi = 300)


# Residual calculation
preds_all <- preds_all %>%
  mutate(residual_cm = pred_cm - obs_cm)

# ---- Figure: Predicted vs Observed ----
fig_species_scatter <- ggplot(preds_all %>% filter(Flight == "LiDAR"),
                              aes(x = obs_cm, y = pred_cm, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 2) +
  coord_equal() +
  facet_wrap(~ Model, ncol = 2) +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Predicted vs. Observed DBH by Model (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# ---- Figure: Residuals vs Predicted ----
fig_species_resid <- ggplot(preds_all %>% filter(Flight == "LiDAR"),
                            aes(x = pred_cm, y = residual_cm, color = Species)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~ Model, ncol = 2) +
  labs(x = "Predicted DBH (cm)", y = "Residual (Pred − Obs, cm)",
       title = "Residuals vs. Predicted DBH by Model (LiDAR 3× Sidelap)",
       color = "Species") +
  theme_bw()

# Save
ggsave("fig_species_scatter.png", fig_species_scatter, width = 8, height = 6, dpi = 300)
ggsave("fig_species_resid.png",   fig_species_resid, width = 8, height = 6, dpi = 300)



# Print or save
print(fig1); print(fig2)
ggsave("fig1_pred_vs_obs_lidar.png", fig1, width = 8, height = 3.2, dpi = 300)
ggsave("fig2_residuals_lidar.png",  fig2, width = 8, height = 3.2, dpi = 300)

ggsave(
  filename = "fig2_residuals_lidar.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
) 
  
  
  
------------
# all at once
  
  
library(ggplot2)

make_scatter_fig <- function(flight_name) {
  ggplot(preds_all_collapsed %>% filter(Flight == flight_name),
         aes(x = obs_cm, y = pred_cm, color = Species)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6, size = 2) +
    coord_equal() +
    facet_wrap(~ Model, nrow = 2) +
    labs(
      x = "Observed DBH (cm)",
      y = "Predicted DBH (cm)",
      title = paste0("Predicted vs. Observed DBH — ", flight_name),
      color = "Species"
    ) +
    theme_bw()
}

make_residual_fig <- function(flight_name) {
  ggplot(preds_all_collapsed %>% filter(Flight == flight_name),
         aes(x = pred_cm, y = residual_cm, color = Species)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(alpha = 0.6, size = 2) +
    facet_wrap(~ Model, nrow = 2) +
    labs(
      x = "Predicted DBH (cm)",
      y = "Residual (Pred − Obs, cm)",
      title = paste0("Residuals vs. Predicted DBH — ", flight_name),
      color = "Species"
    ) +
    theme_bw()
}
  
make_scatter_fig("LiDAR", "Mixed-effects (pooled)")





make_scatter_fig <- function(flight1, model1, flight2, model2) {
  df <- preds_all_collapsed %>%
    filter((Flight == flight1 & Model == model1) |
             (Flight == flight2 & Model == model2)) %>%
    mutate(Panel = paste(Flight, Model, sep = " — "))
  
  ggplot(df, aes(x = obs_cm, y = pred_cm, color = Species)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6, size = 2) +
    coord_equal() +
    facet_wrap(~ Panel, nrow = 1) +
    labs(
      x = "Observed DBH (cm)", 
      y = "Predicted DBH (cm)",
      title = paste("Predicted vs. Observed DBH:", flight1, "vs", flight2),
      color = "Species"
    ) +
    theme_bw()
}


# Compare LiDAR and Fusion, both Mixed-effects
make_scatter_fig("LiDAR", "Log–log (Non-sequoia)", 
                 "LiDAR", "Log–log (Sequoia-only)")

# Compare LiDAR Mixed-effects vs LiDAR RF
make_scatter_fig("LiDAR", "Mixed-effects (pooled)", 
                 "LiDAR", "RF (pooled)")




make_residual_fig <- function(flight1, model1, flight2, model2) {
  df <- preds_all_collapsed %>%
    filter((Flight == flight1 & Model == model1) |
             (Flight == flight2 & Model == model2)) %>%
    mutate(Panel = paste(Flight, Model, sep = " — "))
  
  ggplot(df, aes(x = pred_cm, y = residual_cm, color = Species)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(alpha = 0.6, size = 2) +
    facet_wrap(~ Panel, nrow = 1) +
    labs(
      x = "Predicted DBH (cm)",
      y = "Residual (Pred − Obs, cm)",
      title = paste("Residuals vs. Predicted DBH:", flight1, "vs", flight2),
      color = "Species"
    ) +
    theme_bw()
}


# Compare LiDAR and Fusion residuals, Mixed-effects models
make_residual_fig("LiDAR", "Mixed-effects (pooled)", 
                  "LiDAR", "Mixed-effects (pooled)")

# Compare LiDAR RF vs LiDAR Mixed-effects residuals
make_residual_fig("LiDAR", "RF (pooled)", 
                  "LiDAR", "Mixed-effects (pooled)")

