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

# =========================
# 1. Load Data
# =========================
# Full LAS cloud (forest-level)
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las")

# Tree crown polygons (already delineated by cloud2trees)
crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/Combo_crowns.gpkg")

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

st_write(crowns_metrics, 
         "C:/Users/User/Desktop/RandomForest2/filtered_crowns_with_all_metrics_3.gpkg",
         delete_layer = TRUE)
--------------
crowns_metrics

# read your CSV
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest2/fieldpointscombo2.csv")

# check structure
glimpse(fieldpoints)

# Clean SFM column and rename
fieldpoints_clean <- fieldpoints %>%
  mutate(
    SFM = str_trim(SFM),                 # remove spaces
    SFM = str_replace_all(SFM, "[\r\n]", "")  # remove carriage returns / newlines
  ) %>%
  dplyr::rename(treeID = SFM) %>%
  as_tibble()   # ✅ force into tibble so select/join works

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

# Remove rows where dbh_cm is NA
crowns_final <- crowns_joined %>%
  filter(!is.na(dbh_cm.y))

# Log-transform DBH
crowns_final$log_dbh <- log(crowns_final$dbh_cm.y)

----------------------------------------
  # Model DBH: Random Forest
----------------------------------------

  # Add log(DBH) safely
crowns_final <- crowns_final %>%
  mutate(log_dbh = ifelse(dbh_cm.y > 0, log(dbh_cm.y), NA_real_)) %>%
  filter(!is.na(log_dbh))

# Fit Random Forest on log(DBH)
set.seed(123)

rf_model <- randomForest(
  log_dbh ~ tree_height_m + is_sequoia + crown_area_m2 +
    point_density.y + p95 + p75 + p50 + mean_z + max_z +
    sd_z + cv_z + treecbh + slope_mean + ndvi_mean,
  data = crowns_final,
  importance = TRUE,
  ntree = 500
)




--------------------
# Seperate RF models for sequoia vs non sequoia
  
# Sequoia trees
sequoia_data <- RF_variables[RF_variables$is_sequoia == 1 & !is.na(RF_variables$dbh_cm), ]

# Non-Sequoia trees
non_sequoia_data <- RF_variables[RF_variables$is_sequoia == 0 & !is.na(RF_variables$dbh_cm), ]

rf_sequoia <- randomForest(
  dbh_cm ~ tree_height_m + crown_area_m2 + basal_area_m2 + radius_m +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z + point_density +
    slope_mean + aspect_circmean + ndvi_mean,
  data = sequoia_data,
  importance = TRUE,
  ntree = 500
)

print(rf_sequoia)

rf_nonsequoia <- randomForest(
  dbh_cm ~ tree_height_m_1 + crown_area_m2 + basal_area_m2 + radius_m +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z + point_density +
    slope_mean + aspect_circmean + ndvi_mean,
  data = non_sequoia_data,
  importance = TRUE,
  ntree = 500
)

print(rf_nonsequoia)

# Separate them first:
new_sequoia <- filtered_crowns[filtered_crowns$is_sequoia == 1, ]
new_nonsequoia <- filtered_crowns[filtered_crowns$is_sequoia == 0, ]

new_sequoia$predicted_dbh <- predict(rf_sequoia, new_sequoia)
new_nonsequoia$predicted_dbh <- predict(rf_nonsequoia, new_nonsequoia)

nrow(sequoia_data)
summary(sequoia_data$dbh_cm)

nrow(nonsequoia_data)
summary(nonsequoia_data$dbh_cm)
