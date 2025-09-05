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
set.seed(123)
ctrl <- trainControl(method = "cv", number = 5, savePredictions = "final")

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



summarize_cv_lm <- function(data, formula, label) {
  folds <- createFolds(data$dbh_cm, k = 5)   # 5-fold CV
  metrics <- map_dfr(folds, function(fold) {
    train_data <- data[-fold, ]
    test_data  <- data[fold, ]

    # Fit log–log LM
    model <- lm(formula, data = train_data)

    # Predictions
    preds <- exp(predict(model, newdata = test_data))
    obs   <- test_data$dbh_cm

    tibble(
      R2   = cor(obs, preds, use = "complete.obs")^2,
      RMSE = sqrt(mean((obs - preds)^2, na.rm = TRUE)),
      MAE  = mean(abs(obs - preds), na.rm = TRUE),
      Bias = mean(preds - obs, na.rm = TRUE)
    )
  })

  metrics %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    ) %>%
    mutate(Source = label)
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
evaluate_lmer_cv <- function(data, k = 5) {
  folds <- createFolds(data$dbh_cm, k = k, list = TRUE, returnTrain = FALSE)
  
  fold_metrics <- map_dfr(seq_along(folds), function(i) {
    test_idx <- folds[[i]]
    train <- data[-test_idx, ]
    test  <- data[test_idx, ]
    
    # Fit model
    model <- lmer(
      log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) + (1 | Species),
      data = train
    )
    
    # Predict
    preds <- tibble(
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(model, newdata = test, allow.new.levels = TRUE))
    )
    
    # Metrics
    tibble(
      R2    = cor(preds$obs_cm, preds$pred_cm)^2,
      RMSE  = sqrt(mean((preds$obs_cm - preds$pred_cm)^2)),
      MAE   = mean(abs(preds$obs_cm - preds$pred_cm)),
      Bias  = mean(preds$pred_cm - preds$obs_cm, na.rm = TRUE)
    )
  })
  
  # Average across folds
  fold_metrics %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    )
}

    mix_results <- evaluate_lmer_cv(crown_data, k = 5) %>%
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
  paste(flight, "Sequoia log–log")
)


    # 6. Log–log allometry (Non-sequoia)
    res_loglog_nonseq <- summarize_cv_lm(
      filter(crown_data, is_sequoia == 0),
      log(dbh_cm) ~ log(lidar_height) + log(crown_area_m2),
      paste(flight, "Non-sequoia log–log")
    )

    # ---- Collect ----
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