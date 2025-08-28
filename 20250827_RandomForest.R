# =========================
# Libraries
# =========================
library(lidR)
library(sf)
library(terra)
library(raster)
library(dplyr)
library(tidyverse)
library(purrr)
library(exactextractr)
library(tibble)
library(stringr)
library(randomForest)
library(caret)

# =========================
# 1. Load Data
# =========================
# Full LAS cloud (forest-level)
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las")

# Tree crown polygons (already delineated by cloud2trees)
crowns <- st_read("C:/Users/User/Desktop/RandomForest2/filtered_crowns_with_all_metrics.gpkg")

# Digital terrain model (for slope)
dtm <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")

# NDVI raster (separately generated from SfM or drone imagery)
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

# =========================
# 2. Normalize LAS
# =========================
# Classify ground
las_ground <- classify_ground(las, csf(
  sloop_smooth = FALSE, 
  class_threshold = 0.2, 
  cloth_resolution = 0.2, 
  rigidness = 2L, 
  iterations = 500L, 
  time_step = 0.65
))

# Normalize heights above ground
las_norm <- normalize_height(las_ground, knnidw())

# =========================
# 3. Pre-compute ancillary rasters
# =========================
# Convert DTM to slope raster
slope <- terrain(raster(dtm), opt = "slope", unit = "degrees")

# =========================
# 4. Loop over crowns and compute metrics
# =========================
metrics_list <- map_dfr(1:nrow(crowns), function(i) {
  crown_poly <- crowns[i, ]  # single crown polygon
  
  # Clip LAS points inside this crown
  las_clip <- clip_roi(las_norm, crown_poly)
  
  # Skip if empty
  if (is.null(las_clip) || npoints(las_clip) == 0) {
    return(tibble(treeID = crown_poly$treeID,
                  point_density = NA, slope_mean = NA,
                  max_z = NA, mean_z = NA, 
                  p25 = NA, p50 = NA, p75 = NA, p95 = NA, 
                  sd_z = NA, cv_z = NA, ndvi_mean = NA))
  }
  
  # Heights inside crown
  Z <- las_clip@data$Z
  area_m2 <- as.numeric(st_area(crown_poly))
  
  tibble(
    treeID = crown_poly$treeID,
    point_density = npoints(las_clip) / area_m2,
    slope_mean = exact_extract(slope, crown_poly, 'mean'),
    max_z = max(Z, na.rm = TRUE),
    mean_z = mean(Z, na.rm = TRUE),
    p25 = quantile(Z, 0.25, na.rm = TRUE),
    p50 = quantile(Z, 0.50, na.rm = TRUE),
    p75 = quantile(Z, 0.75, na.rm = TRUE),
    p95 = quantile(Z, 0.95, na.rm = TRUE),
    sd_z = sd(Z, na.rm = TRUE),
    cv_z = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE),
    ndvi_mean = exact_extract(ndvi, crown_poly, 'mean')
  )
})

# =========================
# 5. Join back to crowns and save
# =========================
summary(crowns)
crowns_metrics <- crowns %>%
  left_join(metrics_list, by = "treeID")

st_write(crowns_metrics, 
         "C:/Users/User/Desktop/RandomForest3/filtered_crowns_with_all_metrics.gpkg",
         delete_layer = TRUE)
--------------
crowns_metrics

# --- Read and clean fieldpoints ---
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/FieldPointsCombo.csv")


# Clean SFM column and rename
fieldpoints_clean <- fieldpoints %>%
  mutate(
    SFM = str_trim(SFM),                 # remove spaces
    SFM = str_replace_all(SFM, "[\r\n]", "")  # remove carriage returns / newlines
  ) %>%
  dplyr::rename(treeID = SFM) %>%
  as_tibble()   # force into tibble so select/join works

# Ensure crowns_metrics has treeID as character too
crowns_metrics$treeID <- as.character(crowns_metrics$treeID)
fieldpoints_clean$treeID <- as.character(fieldpoints_clean$treeID)

# Now join: keep only useful columns from fieldpoints
crowns_joined <- crowns_metrics %>%
  left_join(
    fieldpoints_clean %>%
      dplyr::select(treeID, dbh_cm, treecbh, treeheight),
    by = "treeID"
  )

summary(fieldpoints)

# Remove rows where dbh_cm is NA
crowns_final <- crowns_joined %>%
  filter(!is.na(dbh_cm.y))

# Add log(DBH) safely
crowns_final <- crowns_final %>%
  mutate(log_dbh = ifelse(dbh_cm.y > 0, log(dbh_cm.y), NA_real_)) %>%
  filter(!is.na(log_dbh))   # remove invalid rows


# Save back out
st_write(crowns_final, "C:/Users/User/Desktop/RandomForest3/crowns_with_fielddata.gpkg",
         delete_layer = TRUE)

----------------------------------------
  # Model DBH: Random Forest
----------------------------------------

crowns_final <- crowns_final %>%
  filter(dbh_cm.y > 0, dbh_cm.y < 350)
  
# Check for missing values
summary(crowns_final)

# Fit the model to predict DBH
set.seed(123)

rf_model <- randomForest(
  log_dbh ~ tree_height_m + is_sequoia + crown_area_m2 +
     p95.x + p75.x + p50.x + p25.x + mean_z.x + max_z.x +
    sd_z.x + cv_z.x + point_density + slope_mean.x + treecbh + ndvi_mean.x,
  data = crowns_final,
  importance = TRUE,
  ntree = 500
)

# View results
print(rf_model)
varImpPlot(rf_model)

# Predictions (already back-transformed from log to cm)
crowns_final$log_dbh_pred <- predict(rf_model, newdata = crowns_final)
crowns_final$dbh_pred_cm  <- exp(crowns_final$log_dbh_pred)

# Define observed and predicted
obs  <- crowns_final$dbh_cm.y        # observed DBH
pred <- crowns_final$dbh_pred_cm   # predicted DBH

# Compute metrics
rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
bias <- mean(pred - obs, na.rm = TRUE)
r2   <- cor(obs, pred, use = "complete.obs")^2

# Put results into a nice tibble
results <- tibble(
  R2   = round(r2, 3),
  RMSE = round(rmse, 2),
  Bias = round(bias, 2),
  N    = sum(!is.na(obs) & !is.na(pred))
)

print(results)

# Create a data frame with observed and predicted values
plot_df <- data.frame(
  obs  = obs,
  pred = pred
)

# Observed vs Predicted
ggplot(plot_df, aes(x = obs, y = pred)) +
  geom_point(color = "darkgreen", alpha = 0.7, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Observed DBH (cm)",
    y = "Predicted DBH (cm)",
    title = "Observed vs Predicted DBH"
  ) +
  theme_minimal(base_size = 14)

# Residuals vs Predicted
ggplot(results, aes(x = pred, y = residuals)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)",
       title = "Residuals vs Predicted DBH") +
  theme_minimal()


----
  
# ============================
# Libraries
# ============================

set.seed(123)  # reproducibility

# ============================
# 1. Prepare Data
# ============================
crown_data <- crowns_final %>%
  filter(!is.na(dbh_cm) & dbh_cm > 0) %>%
  mutate(log_dbh = log(dbh_cm))  # log-transform DBH for modeling

# ============================
# 2. CV Setup
# ============================
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"   # ✅ stores fold-wise predictions
)

# ============================
# 3. Fit Random Forest with 5-Fold CV
# ============================
rf_cv <- train(
  log_dbh ~ tree_height_m + is_sequoia + crown_area_m2 +
     p95.x + p75.x + p50.x + p25.x + mean_z.x + max_z.x +
    sd_z.x + cv_z.x + point_density + slope_mean.x + treecbh + ndvi_mean.x,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,     # try multiple mtry values
  importance = TRUE,
  ntree = 500
)

print(rf_cv)
varImp(rf_cv)   # variable importance

# ============================
# 4. Extract Out-of-Sample Predictions
# ============================
preds_log <- rf_cv$pred$pred   # CV predictions on log scale
obs_log   <- rf_cv$pred$obs

# Back-transform to cm
preds_dbh <- exp(preds_log)
obs_dbh   <- exp(obs_log)

# ============================
# 5. Compute Metrics
# ============================
r2   <- cor(preds_dbh, obs_dbh, use = "complete.obs")^2
rmse <- sqrt(mean((preds_dbh - obs_dbh)^2, na.rm = TRUE))
bias <- mean(preds_dbh - obs_dbh, na.rm = TRUE)

results <- tibble(R2 = r2, RMSE = rmse, Bias = bias, N = length(obs_dbh))
print(results)

# ============================
# 6. Plots
# ============================

# Observed vs Predicted
plot_df <- tibble(obs = obs_dbh, pred = preds_dbh)

ggplot(plot_df, aes(x = obs, y = pred)) +
  geom_point(color = "darkgreen", alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Observed vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

# Residuals vs Predicted
plot_df <- plot_df %>%
  mutate(resid = obs - pred)

ggplot(plot_df, aes(x = pred, y = resid)) +
  geom_point(color = "blue", alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)",
       title = "Residuals vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

-------------
  
  
  
crowns_final <- crowns_final %>%
  filter(dbh_cm.y > 0, dbh_cm.y < 350)

crown_data <- crowns_final

# ============================
# Libraries
# ============================
library(caret)
library(randomForest)
library(dplyr)
library(tibble)

set.seed(123)  # reproducibility

# ============================
# 1. Prepare Data
# ============================
crown_data <- crowns_final %>%
  filter(!is.na(dbh_cm.y) & dbh_cm.y > 0 & dbh_cm.y < 350) %>%   # remove bad DBH
  mutate(
    log_dbh       = log(dbh_cm.y),        # log-transform response
    log_height    = log(tree_height_m),   # log-transform height
    log_crown_area = log1p(crown_area_m2) # log1p handles zero crown area
  )

# ============================
# 2. Define CV setup
# ============================
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"   # keeps out-of-fold predictions
)

# ============================
# 3. Fit Random Forest Model
# ============================
rflog_model <- train(
  log_dbh ~ log_height + is_sequoia + log_crown_area +
    p95.x + p75.x + p50.x + p25.x + mean_z.x + max_z.x +
    sd_z.x + cv_z.x + point_density + slope_mean.x + treecbh + ndvi_mean.x,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,       # tests multiple mtry values
  importance = TRUE,
  ntree = 500
)

# View results
print(rflog_model)
varImp(rflog_model)


# Pull CV predictions (log scale)
preds_log <- rflog_model$pred$pred
obs_log   <- rflog_model$pred$obs

# Back-transform to cm
preds_dbh <- exp(preds_log)
obs_dbh   <- exp(obs_log)

# Metrics
r2   <- cor(preds_dbh, obs_dbh, use = "complete.obs")^2
rmse <- sqrt(mean((preds_dbh - obs_dbh)^2, na.rm = TRUE))
bias <- mean(preds_dbh - obs_dbh, na.rm = TRUE)

results <- tibble(R2 = r2, RMSE = rmse, Bias = bias, N = length(obs_dbh))
print(results)


# Observed vs Predicted
ggplot(data.frame(obs = obs_dbh, pred = preds_dbh), aes(x = obs, y = pred)) +
  geom_point(color = "forestgreen", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Observed vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

# Residuals vs Predicted
residuals <- preds_dbh - obs_dbh
ggplot(data.frame(pred = preds_dbh, resid = residuals), aes(x = pred, y = resid)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)",
       title = "Residuals vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

valid_idx <- which(!is.na(obs) & !is.na(preds_dbh))
obs <- obs[valid_idx]
preds_dbh <- preds_dbh[valid_idx]

residuals <- obs_dbh - preds_dbh
diagnostics_df <- crown_data %>%
  slice(valid_idx) %>%
  mutate(
    obs_dbh = obs_dbh,
    preds_dbh = preds_dbh,
    residuals = residuals
  )



# Example: residuals vs crown area
ggplot(diagnostics_df, aes(x = log_crown_area, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(x = "Crown Area (m²)", y = "Residuals (cm)",
       title = "Residuals vs Crown Area") +
  theme_minimal()


# Force numeric just for correlation
diagnostics_df_num <- diagnostics_df %>%
  select(residuals, tree_height_m, crown_area_m2, slope_mean.x, ndvi_mean.x) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix
cor_matrix <- cor(diagnostics_df_num, use = "complete.obs")
print(cor_matrix)

# Force numeric just for correlation
diagnostics_df_num <- diagnostics_df %>%
  select(residuals, tree_height_m, crown_area_m2, slope_mean.x, ndvi_mean.x) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix
cor_matrix <- cor(diagnostics_df_num, use = "complete.obs")
print(cor_matrix)

# Force numeric just for correlation
diagnostics_df_num <- diagnostics_df %>%
  select(residuals, tree_height_m, crown_area_m2, slope_mean.x, ndvi_mean.x) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix
cor_matrix <- cor(diagnostics_df_num, use = "complete.obs")
print(cor_matrix)







library(mgcv)
gam_resid <- gam(residuals ~ s(log_crown_area) + s(log_height), data = diagnostics_df)
summary(gam_resid)
plot(gam_resid)

