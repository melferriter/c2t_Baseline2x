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
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las")
# Tree crown polygons (already delineated by cloud2trees)
crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg")
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


# --- Read and clean fieldpoints ---
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv")


# Clean SFM column and rename
fieldpoints_clean <- fieldpoints %>%
  mutate(
    SFMonly = str_trim(SFMonly),                 # remove spaces
    Side3x = str_replace_all(SFMonly, "[\r\n]", "")  # remove carriage returns / newlines
  ) %>%
  dplyr::rename(treeID = SFMonly) %>%
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

summary(crown_data)

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

# Extract variable importance from your rf model
var_imp <- varImp(rflog_model, scale = TRUE)
var_imp_df <- as.data.frame(var_imp$importance) %>%
  tibble::rownames_to_column("Variable") %>%
  arrange(desc(Overall))

# ============================
# Variable Importance Plot (Publication Style)
# ============================
library(caret)
library(ggplot2)
library(dplyr)

# Extract and tidy variable importance
var_imp <- varImp(rflog_model, scale = TRUE)
var_imp_df <- as.data.frame(var_imp$importance) %>%
  tibble::rownames_to_column("Variable") %>%
  arrange(Overall)   # arrange ascending for nicer gradient

# Define nicer labels just for this plot
custom_labels <- c(
  p50            = "50th Percentile Height (m)",
  max_z          = "Max Height (m)",
  log_height     = "Log Tree Height",
  p75            = "75th Percentile Height (m)",
  p95            = "95th Percentile Height (m)",
  mean_z         = "Mean Height (m)",
  is_sequoia     = "Species",
  cv_z           = "Height CV",
  treecbh        = "Crown Base Height (m)",
  log_crown_area = "Log Crown Area (m²)",
  point_density.x = "Point Density",
  slope_mean     = "Mean Slope (°)",
  p25            = "25th Percentile Height (m)",
  sd_z           = "Height SD (m)",
  ndvi_mean      = "NDVI"
)

# Plot
ggplot(var_imp_df, aes(x = Overall, y = reorder(Variable, Overall), fill = Overall)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +  # gradient from light blue to dark
  geom_text(aes(label = round(Overall, 1)), 
            hjust = -0.2, size = 4, color = "black") +
  scale_y_discrete(labels = custom_labels) +# value labels
  labs(
    title = "Variable Importance",
    x = "Relative Importance",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major.y = element_blank(),   # cleaner y-axis
    panel.grid.minor = element_blank()
  ) +
  expand_limits(x = max(var_imp_df$Overall) * 1.1)

ggsave(
  filename = "VariableImportancsfmonly.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)

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

library(ggplot2)

# ============================
# Scatter: Observed vs Predicted DBH
# ============================
# Pull values from cv_summary
r2_text   <- paste0("R² = ", round(cv_summary$R2_mean, 2), " ± ", round(cv_summary$R2_sd, 2))
rmse_text <- paste0("RMSE = ", round(cv_summary$RMSE_mean, 1), " ± ", round(cv_summary$RMSE_sd, 1))
bias_text <- paste0("Bias = ", round(cv_summary$Bias_mean, 1), " ± ", round(cv_summary$Bias_sd, 1))

# Scatterplot with metrics box (mean ± SD)
ggplot(data.frame(obs = obs_dbh, pred = preds_dbh), aes(x = obs, y = pred)) +
  geom_point(color = "#FF10F0", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)",
       title = "SfM") +
  annotate("text", x = 50, y = 250, hjust = 0, vjust = 1,
           label = paste(r2_text, rmse_text, bias_text, sep = "\n"),
           size = 4.2, fontface = "bold") +
  theme_bw(base_size = 16)

ggsave(
  filename = "Scattersfmonly.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)

# ============================
# Residuals vs Predicted DBH
# ============================
residuals <- obs_dbh - preds_dbh
resid_df <- data.frame(pred = preds_dbh, resid = residuals)

ggplot(data.frame(pred = preds_dbh, resid = residuals), aes(x = pred, y = resid)) +
  geom_point(shape = 16, color = "black", fill = "black", alpha = 0.7, size = 2) +  # filled points
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "red", size = 1) +  # loess trend
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)",
       title = "Residuals vs Predicted DBH (5-Fold CV)") +
  theme_minimal(base_size = 16)


ggsave(
  filename = "Residualsfmonly.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)


------
library(patchwork)

p1 <- ggplot(data.frame(obs = obs_dbh, pred = preds_dbh), aes(x = obs, y = pred)) +
  geom_point(color = "forestgreen", alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 60, y = 240,
           label = paste0("R² = ", round(cv_summary$R2_mean, 2), " ± ", round(cv_summary$R2_sd, 2), "\n",
                          "RMSE = ", round(cv_summary$RMSE_mean, 1), " ± ", round(cv_summary$RMSE_sd, 1), "\n",
                          "Bias = ", round(cv_summary$Bias_mean, 1), " ± ", round(cv_summary$Bias_sd, 1)),
           hjust = 0, size = 5) +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)", title = "Observed vs Predicted DBH") +
  theme_minimal(base_size = 18)

p2 <- ggplot(data.frame(pred = preds_dbh, resid = residuals), aes(x = pred, y = resid)) +
  geom_point(shape = 16, color = "black", alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(x = "Predicted DBH (cm)", y = "Residuals (cm)", title = "Residuals vs Predicted DBH") +
  theme_minimal(base_size = 18)

# Combine side by side
combined <- p1 + p2  

ggsave(
  filename = "scattresid.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 12, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)




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


---------
  


library(ggplot2)
library(dplyr)

# Calculate mean tree height
mean_height <- mean(crowns_final$treeheight, na.rm = TRUE)

# Get histogram max for line end
hist_max <- ggplot_build(
  ggplot(crowns_final, aes(x = treeheight)) + 
    geom_histogram(bins = 20)
)$data[[1]]$count %>% max()

# Attractive histogram
ggplot(crowns_final, aes(x = treeheight)) +
  geom_histogram(bins = 25, 
                 fill = "steelblue", color = "white", alpha = 0.85) +
  
  # Mean line
  annotate("segment",
           x = mean_height, xend = mean_height,
           y = 0, yend = hist_max,
           color = "firebrick", linetype = "dashed", size = 1.2) +
  
  # Label for mean
  annotate("text", x = mean_height, y = hist_max * 0.95,
           label = paste0("Mean = ", round(mean_height, 1), " m"),
           color = "firebrick", fontface = "bold", angle = 0, vjust = -1.5,
           size = 5) +
  
  labs(x = "Tree Height (m)", y = "Number of Trees",
       title = "Distribution of Tree Heights") +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text  = element_text(size = 13, color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "distsfmonly.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)
