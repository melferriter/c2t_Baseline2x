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
rm(list = ls(globalenv()))
# =========================
# 1. Load Data
# =========================
# Full LAS cloud (forest-level)
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las")

# Tree crown polygons (already delineated by cloud2trees)
crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg")

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
crowns_metrics <- crowns %>%
  left_join(metrics_list, by = "treeID")

# --- Read and clean fieldpoints ---
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv")


# Clean SFM column and rename
fieldpoints_clean <- fieldpoints %>%
  mutate(
    SFMcombo = str_trim(SFMcombo),                 # remove spaces
    SFMcombo = str_replace_all(SFMcombo, "[\r\n]", "")  # remove carriage returns / newlines
  ) %>%
  dplyr::rename(treeID = SFMcombo) %>%
  as_tibble()   # force into tibble so select/join works

# Ensure crowns_metrics has treeID as character too
crowns_metrics$treeID <- as.character(crowns_metrics$treeID)
fieldpoints_clean$treeID <- as.character(fieldpoints_clean$treeID)


# Now join: keep only useful columns from fieldpoints
crowns_final <- crowns_metrics %>%
  left_join(
    fieldpoints_clean %>%
      dplyr::select(treeID, dbh_cm, treecbh, treeheight),
    by = "treeID"
  )


# ----------------------------------------
# Model DBH: Random Forest
# ----------------------------------------

# ============================
# 1. Prepare Data
# ============================
summary(crowns_final)

crowns_final %>%
  filter(dbh_cm.y == 0 | treecbh == 0) %>%
  select(treeID, dbh_cm.y, treecbh, treeheight)

crown_data <- crowns_final %>%
  filter(!is.na(dbh_cm.y), dbh_cm.y > 0, dbh_cm.y < 350) %>%   # clean DBH
  mutate(
    log_dbh        = log(dbh_cm.y),        # response
    log_height     = log(treeheight),    # field tree height
    log_crown_area = log1p(crown_area_m2) # predictor
  )

# ============================
# 2. Define CV setup
# ============================
set.seed(123)  # master seed
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"
)

# Ensure same randomization across folds
seeds <- vector("list", length = 5 + 1)  # 5 folds + final model
for(i in 1:5) seeds[[i]] <- sample.int(1000, 5)  # one per tuning parameter
seeds[[6]] <- sample.int(1000, 1)  # final model
ctrl$seeds <- seeds

rflog_model <- train(
  log_dbh ~ log_height + is_sequoia + log_crown_area +
    p95 + p75 + p50 + p25 + mean_z + max_z +
    sd_z + cv_z + point_density.x + slope_mean + treecbh + ndvi_mean,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  importance = TRUE,
  ntree = 1000
)


print(rflog_model)
varImp(rflog_model)

# ============================
# 4. Back-transform Predictions + Overall Metrics
# ============================

# Back-transform predictions and observations from log scale
preds_dbh <- exp(rflog_model$pred$pred)  # predicted DBH (cm)
obs_dbh   <- exp(rflog_model$pred$obs)   # observed DBH (cm)

# Filter valid cases
valid_idx <- which(!is.na(obs_dbh) & !is.na(preds_dbh))
obs_dbh   <- obs_dbh[valid_idx]
preds_dbh <- preds_dbh[valid_idx]

# ---- Overall model performance (all folds combined) ----
r2   <- cor(preds_dbh, obs_dbh, use = "complete.obs")^2   # coefficient of determination
rmse <- sqrt(mean((preds_dbh - obs_dbh)^2, na.rm = TRUE)) # root mean squared error
bias <- mean(preds_dbh - obs_dbh, na.rm = TRUE)           # average over/under prediction

results <- tibble(R2 = r2, RMSE = rmse, Bias = bias, N = length(obs_dbh))
print(results)

# ============================
# 5. caret’s per-fold metrics (RMSE, R², MAE)
# ============================

# caret automatically saves metrics per resample (fold)
cv_results <- rflog_model$resample
head(cv_results)  
# Contains: RMSE, Rsquared, MAE, Resample (e.g., Fold1...Fold5)

# Summarise mean ± SD across folds
cv_summary_caret <- cv_results %>%
  summarise(
    R2_mean   = mean(Rsquared, na.rm = TRUE),
    R2_sd     = sd(Rsquared, na.rm = TRUE),
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_sd   = sd(RMSE, na.rm = TRUE),
    MAE_mean  = mean(MAE, na.rm = TRUE),
    MAE_sd    = sd(MAE, na.rm = TRUE)
  )
print(cv_summary_caret)

# ============================
# 6. Custom per-fold metrics (with Bias added)
# ============================

# Extract CV predictions with fold IDs
preds <- rflog_model$pred %>%
  mutate(
    obs_dbh  = exp(obs),   # back-transform
    pred_dbh = exp(pred),
    resid    = pred_dbh - obs_dbh
  )

# Compute metrics per fold
cv_metrics <- preds %>%
  group_by(Resample) %>%
  summarise(
    R2   = cor(obs_dbh, pred_dbh, use = "complete.obs")^2,
    RMSE = sqrt(mean((pred_dbh - obs_dbh)^2, na.rm = TRUE)),
    MAE  = mean(abs(pred_dbh - obs_dbh), na.rm = TRUE),
    Bias = mean(resid, na.rm = TRUE),  # bias isn’t in caret’s output
    .groups = "drop"
  )
print(cv_metrics)

# Summarise mean ± SD across folds (including Bias)
cv_summary <- cv_metrics %>%
  summarise(
    R2_mean   = mean(R2),   R2_sd   = sd(R2),
    RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
    MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
    Bias_mean = mean(Bias), Bias_sd = sd(Bias)
  )
print(cv_summary)


# ============================
# 5. Diagnostics Plots
# ============================

# Scatter: observed vs predicted
ggplot(data.frame(obs = obs_dbh, pred = preds_dbh), aes(x = obs, y = pred)) +
  geom_point(color = "forestgreen", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "Observed vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

# Residuals vs predicted
residuals <- obs_dbh - preds_dbh
ggplot(data.frame(pred = preds_dbh, resid = residuals), aes(x = pred, y = resid)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)",
       title = "Residuals vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 14)

# ============================
# 6. Residual Diagnostics
# ============================
diagnostics_df <- crown_data %>%
  slice(valid_idx) %>%
  mutate(
    obs_dbh    = obs_dbh,
    preds_dbh  = preds_dbh,
    residuals  = obs_dbh - preds_dbh
  )

# Example: residuals vs crown area
ggplot(diagnostics_df, aes(x = log_crown_area, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(x = "Crown Area (m²)", y = "Residuals (cm)",
       title = "Residuals vs Crown Area") +
  theme_minimal()

# Correlation of residuals with predictors
diagnostics_df_num <- diagnostics_df %>%
  st_drop_geometry() %>%   # remove geometry if still present
  select(residuals, log_height, log_crown_area, slope_mean, ndvi_mean) %>%
  mutate(across(everything(), as.numeric))

cor_matrix <- cor(diagnostics_df_num, use = "complete.obs")
print(cor_matrix)


-------
N <- nrow(crown_data)
print(N)

------------
# ======================================
# Compare model accuracy side by side
# ======================================

# Example summaries (replace with your actual results)
triple_return <- tibble(
  Model = "Triple Return Sidelap",
  R2_mean = 0.89, R2_sd = 0.11,
  RMSE_mean = 29.8, RMSE_sd = 20.5,
  MAE_mean = 19.6, MAE_sd = 11.1,
  Bias_mean = -8.16, Bias_sd = 11.9,
  N = 51
)

sfm <- tibble(
  Model = "SfM",
  R2_mean = 0.89, R2_sd = 0.15,
  RMSE_mean = 30.6, RMSE_sd = 21.3,
  MAE_mean = 21.1, MAE_sd = 12.3,
  Bias_mean = -9.29, Bias_sd = 17.2,
  N = 52
)

fusion <- tibble(
  Model = "Triple Return + SfM Fusion",
  R2_mean = 0.84, R2_sd = 0.13,
  RMSE_mean = 28.0, RMSE_sd = 19.5,
  MAE_mean = 17.1, MAE_sd = 9.5,
  Bias_mean = -7.40, Bias_sd = 11.1,
  N = 50
)

# Combine all into one comparison table
model_comparison <- bind_rows(triple_return, sfm, fusion)

# View table
print(model_comparison)

# ======================================
# Identify "best" model by metric
# ======================================

# Highest R²
best_r2 <- model_comparison %>% arrange(desc(R2_mean)) %>% slice(1)

# Lowest RMSE
best_rmse <- model_comparison %>% arrange(RMSE_mean) %>% slice(1)

# Lowest MAE
best_mae <- model_comparison %>% arrange(MAE_mean) %>% slice(1)

# Lowest absolute Bias (closest to 0)
best_bias <- model_comparison %>% 
  arrange(abs(Bias_mean)) %>% slice(1)

cat("Best by R²:\n"); print(best_r2)
cat("Best by RMSE:\n"); print(best_rmse)
cat("Best by MAE:\n"); print(best_mae)
cat("Best by Bias:\n"); print(best_bias)

------------
# ============================
# Libraries
# ============================
library(tidyverse)

# ============================
# Example data (replace with your actual values)
# ============================
model_comparison <- tribble(
  ~Model, ~R2_mean, ~R2_sd, ~RMSE_mean, ~RMSE_sd, ~MAE_mean, ~MAE_sd, ~Bias_mean, ~Bias_sd, ~N,
  "Triple Return Sidelap",        0.89, 0.11, 29.8, 20.5, 19.6, 11.1, -8.16, 11.9, 51,
  "SfM",                          0.89, 0.15, 30.6, 21.3, 21.1, 12.3, -9.29, 17.2, 52,
  "Triple Return + SfM Fusion",   0.84, 0.13, 28.0, 19.5, 17.1,  9.5, -7.40, 11.1, 50
)

# ============================
# Reshape for plotting
# ============================
plot_data <- model_comparison %>%
  select(Model, R2_mean, R2_sd, RMSE_mean, RMSE_sd, MAE_mean, MAE_sd, Bias_mean, Bias_sd) %>%
  pivot_longer(
    cols = -Model,
    names_to = c("Metric", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  )

# ============================
# Plot (example: R², RMSE, MAE, Bias)
# ============================
ggplot(plot_data, aes(x = Model, y = mean, fill = Model)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, position = position_dodge(0.7)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Model Performance Comparison (5-fold CV)",
    y = "Metric Value (mean ± SD)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

