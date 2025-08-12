
## pkgbuild helps us check for Rtools
install.packages("pkgbuild")
# check for Rtools which is required to build packages
pkgbuild::check_build_tools(debug = TRUE)
## remotes helps us get packages hosted on github
install.packages("remotes")
## install lasR from the r-univers
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
## install TreeLS from github
remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F, force=T)
## get cloud2trees
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)


# download the external data
cloud2trees::get_data()
# download the TreeMap data
cloud2trees::get_treemap()
# download the forest type data
cloud2trees::get_foresttype()
# download the landfire cbd data
cloud2trees::get_landfire()

# Libraries====
library(lidR)
library(sf)
library(terra)
library(stars)
library(raster)
library(alphashape3d)
library(plyr)
library(tidyverse)
library(devtools)
library(canopyLazR)
library(spdep)
library(sp)
library(geosphere)
library(rlas)
library(rgl)
library(pracma)
library(spatstat)
library(terra)
library(cloud2trees)
library(tidyverse)
library(sf)
library(purrr)
library(patchwork)
library(viridis)
library(dbplyr)
library(pacman)
library(minpack.lm)
library(FNN)
library(nngeo)
library(mgcv)
pacman::p_load(ggplot2, rgeos, propagate, dplyr, ggpubr, gridExtra)

rm(list = ls(globalenv()))
rm(list = ls())
gc()   
graphics.off() 
---------------------------------------------------
roi <- st_read("C:/Users/User/Desktop/LidarClip/LidarExtent.shp", crs = 32611)
roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las", filter = "-set_withheld_flag 0")
clipped <- clip_roi(las, roi_2d)
plot(clipped)
writeLAS(clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip_clip.las")

i = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip_clip.las"  

cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = i)

## Default ITD window size functions====
itd_tuning_ans <- itd_tuning(input_las_dir = i)
itd_tuning_ans$plot_samples
best_ws <- itd_tuning_ans$ws_fn_list$lin_fn


cloud2trees_ans_c <- cloud2trees::cloud2trees(
  output_dir = tempdir()
  , input_las_dir = i
  , dtm_res_m = 0.5
  , ws = best_ws
  , estimate_tree_dbh = TRUE
  , estimate_tree_type = F
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = TRUE
  , cbh_estimate_missing_cbh = F
)

paste(
  "Default trees extracted:"
  , cloud2trees_ans$crowns_sf %>% nrow()
  , "|| Custom trees extracted:"
  , cloud2trees_ans_c$crowns_sf %>% nrow()
)

# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$crowns_sf, mapping = ggplot2::aes(fill = dbh_m)) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")


# Convert LAS object to a dataframe
lidar_df <- as.data.frame(clipped@data)  # Extract X, Y, Z coordinates

# Convert to an sf object (ensure same CRS as tree crowns)
lidar_sf <- st_as_sf(lidar_df, coords = c("X", "Y"), crs = st_crs(cloud2trees_ans_c$crowns_sf))

# Find which LiDAR points are inside each crown area
intersections <- st_intersects(cloud2trees_ans_c$crowns_sf, lidar_sf, sparse = TRUE)

# Count LiDAR points per treeID
lidar_counts <- data.frame(
  treeID = cloud2trees_ans_c$crowns_sf$treeID,
  point_count = lengths(intersections)  # Number of LiDAR points per hull
)

# Merge with tree crowns and compute density
density <- cloud2trees_ans_c$crowns_sf %>%
  left_join(lidar_counts, by = "treeID") %>%
  mutate(point_density = point_count / as.numeric(st_area(geometry)))

density_filtered <- density %>%
  filter(treeID %in% filtered_crowns$treeID)

st_write(density_filtered,
         "C:/Users/User/Desktop/RandomForest/filtered_crowns_density.gpkg",
         layer = "density_filtered",
         driver = "GPKG",
         delete_dsn = TRUE)  # Overwrites existing file



ggplot() +
  geom_sf(data = density, aes(fill = point_density), color = "black", alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "LiDAR Density") +
  theme_minimal() +
  ggtitle("LiDAR Point Density Within Tree Crowns")




filtered <- filter(density, point_density > 300 & tree_height_m >25)  # Set empirically

ggplot() +
  geom_sf(data = filtered, aes(fill = point_density), color = "black", alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "LiDAR Density") +
  theme_minimal() +
  ggtitle("LiDAR Point Density Within Tree Crowns")

# Filter the crowns to only include those with treeIDs in your filtered list
filtered_crowns <- cloud2trees_ans_c$crowns_sf %>%
  filter(treeID %in% filtered$treeID)

ggplot2::ggplot() + 
  ggplot2::geom_sf(data = filtered_crowns, mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

st_write(filtered_crowns, "C:/Users/User/Desktop/RandomForest/filtered_crowns.gpkg")

filtered_crowns <- filtered_crowns %>%
  mutate(is_sequoia = ifelse(tree_height_m > 50, 1, 0))

--------------------
  # Load training data
  training_data <- read.csv("C:/Users/User/Desktop/20250614_TrainingData/202150614_Monarchs_Hwd.csv")

training_data <- training_data %>%
  filter(
    !(tree_height_m < 50 & is_sequoia == 1 & dbh_cm >200)  # Likely Sequoias misclassified as non-sequoia
  )

training_data <- training_data %>% drop_na()

ggplot(training_data, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia==1))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("orange", "forestgreen"), labels = c("Non-Sequoia", "Sequoia")) +
  labs(title = "Training Data: Sequoia vs Non-Sequoia Trees", color = "Type",
       x = "Tree Height (m)",
       y = "Tree Diameter at Breast Height (cm)") +
  theme_minimal() +
  theme(text=element_text(size=14))

# Split data by species
sequoia_data     <- training_data %>% filter(is_sequoia == 1)
nonsequoia_data  <- training_data %>% filter(is_sequoia == 0)

ggplot(training_data %>% filter(is_sequoia == 1), aes(x=tree_height_m)) +
  geom_histogram(binwidth = 5, fill = "darkgreen", color = "white") +
  labs(title = "Training Data Sequoia Tree Height Distribution",
       x = "Tree Height (m)",
       y = "Count") +
  theme_minimal()

ggplot(training_data %>% filter(is_sequoia == 0), aes(x = tree_height_m)) +
  geom_histogram(binwidth = 5, fill = "darkgreen", color = "white") +
  labs(title = "Training Data Sequoia Tree Height Distribution",
       x = "Tree Height (m)",
       y = "Count") +
  theme_minimal()

ggplot(sequoia_data, aes(x = dbh_cm)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "black") +
  theme_minimal()

ggplot(nonsequoia_data, aes(x = dbh_cm)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "black") +
  theme_minimal()

qqnorm(training_data$dbh_cm)
qqline(sequoia_data$tree_height_m, col = "red")
-------------
  
# Read and prepare control points
controlpoints <- st_read("C:/Users/User/Desktop/20250614_ControlPoints/ControlPoints_3.shp")
names(controlpoints)[names(controlpoints) == "DBH"] <- "dbh_cm"
names(controlpoints)[names(controlpoints) == "Hgt"] <- "tree_height_m"
st_crs(controlpoints) <- 3857
controlpoints <- st_transform(controlpoints, crs = st_crs(cloud2trees_ans_c$treetops_sf))

# Filter for Giant Sequoias
sequoias <- controlpoints %>% filter(Species == "Giant Sequoia")

# Get min and max values
summary_values <- sequoias %>%
  summarize(
    min_height = min(tree_height_m, na.rm = TRUE),
    max_height = max(tree_height_m, na.rm = TRUE),
    min_dbh = min(dbh_cm, na.rm = TRUE),
    max_dbh = max(dbh_cm, na.rm = TRUE)
  )

print(summary_values)

# Extract extent from LAS and turn into an sf polygon
las_extent <- st_as_sfc(st_bbox(clipped), crs = st_crs(controlpoints))
# Filter control points that fall within the LAS extent
controlpoints_in_las <- controlpoints[st_intersects(controlpoints, las_extent, sparse = FALSE), ]

matched_points <- st_join(controlpoints, filtered_crowns, join = st_within)

#
# Buffer crowns by 1 meter (adjust as needed)
buffered_crowns <- st_buffer(filtered_crowns, dist = 2)

matched_points <- st_join(controlpoints, filtered_crowns, join = st_within)
matched_points_buffered <- st_join(controlpoints, buffered_crowns, join = st_within)
matched_points_bufferedm <- st_join(controlpoints, buffered_crowns, join = st_intersects)

matched_only <- matched_points_bufferedm[!is.na(matched_points_bufferedm$treeID), ]
names(matched_only)[names(matched_only) == "dbh_cm.x"] <- "dbh_cm"
names(matched_only)[names(matched_only) == "tree_height_m.x"] <- "tree_height_m"


ggplot(matched_only, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("orange", "forestgreen"), labels = c("Non-Sequoia", "Sequoia")) +
  labs(title = "Control Data: Sequoia vs Non-Sequoia", color = "Type") +
  theme_minimal()



# Split control points
seq_control <- controlpoints%>% filter(is_sequoia == 1)
nonseq_control <- controlpoints%>% filter(is_sequoia == 0)
seq_controlid <- seq_control$c2t_id
nonseq_controlid <- nonseq_control$c2t_id
mean(seq_control$dbh_cm)
---------------
  
  training_subset <- training_data %>%
  filter(tree_height_m <= max(matched_only$tree_height_m, na.rm = TRUE))

training_sequoia <- training_subset %>% filter(is_sequoia == 1)
training_nonsequoia <- training_subset %>% filter(is_sequoia == 0)

model_gam_sequoia <- gam(dbh_cm ~ s(tree_height_m), data = training_sequoia, method = "REML")
model_gam_nonsequoia <- gam(dbh_cm ~ s(tree_height_m), data = training_nonsequoia, method = "REML")

filtered_crowns <- filtered_crowns %>%
  mutate(predicted_dbh_cm = ifelse(
    is_sequoia == 1,
    predict(model_gam_sequoia, newdata = .),
    predict(model_gam_nonsequoia, newdata = .)
  ))

comparison_df <- matched_only %>%
  mutate(treeID = as.character(treeID)) %>%
  inner_join(
    filtered_crowns %>%
      st_drop_geometry() %>%
      select(treeID, predicted_dbh_cm),
    by = "treeID"
  )

comparison_df <- comparison_df %>%
  mutate(residual = dbh_cm - predicted_dbh_cm)

summary(comparison_df$residual)

rmse <- sqrt(mean(comparison_df$residual^2))
mae <- mean(abs(comparison_df$residual))
cat("RMSE:", rmse, "\nMAE:", mae)
------------------------
  
  
  #GAM 
  library(mgcv)
# Fit GAMs separately for each group
model_gam_sequoia    <- gam(dbh_cm ~ s(tree_height_m), data = sequoia_data, method = "REML")
model_gam_nonsequoia <- gam(dbh_cm ~ s(tree_height_m), data = nonsequoia_data, method = "REML")


filtered_crowns$is_sequoia <- ifelse(
  filtered_crowns$tree_height_m > 50 & filtered_crowns$crown_area_m2 > 100, 1, 0
)


# Predict DBH for all treetops
filtered_crowns <- filtered_crowns %>%
  mutate(predicted_dbh_cm = ifelse(
    is_sequoia == 1,
    predict(model_gam_sequoia, newdata = .),
    predict(model_gam_nonsequoia, newdata = .)
  ))

# Residuals from training fits
plot(model_gam_sequoia, main = "GAM Fit - Sequoia")
plot(model_gam_nonsequoia, main = "GAM Fit - Non-Sequoia")

# Check residuals
hist(resid(model_gam_sequoia), main = "Residuals - Sequoia")
hist(resid(model_gam_nonsequoia), main = "Residuals - Non-Sequoia")

matched_only$treeID <- as.character(matched_only$treeID)
filtered_crowns$treeID <- as.character(filtered_crowns$treeID)

sum(matched_only$treeID %in% filtered_crowns$treeID)

filtered_crowns$predicted_dbh_cm <- predict(model_gam_sequoia, newdata = filtered_crowns)


comparison_df <- matched_only %>%
  inner_join(st_drop_geometry(filtered_crowns)[, c("treeID", "predicted_dbh_cm")],
             by = "treeID")

comparison_df$residual <- comparison_df$dbh_cm - comparison_df$predicted_dbh_cm
summary(comparison_df$residual)

ggplot(comparison_df, aes(x = dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    title = "Measured vs Predicted DBH",
    x = "Measured DBH (cm)",
    y = "Predicted DBH (cm)"
  ) +
  theme_minimal()

summary(comparison_df$residual)
summary(model_gam_sequoia)

hist(comparison_df$residual, main = "Residuals", xlab = "Measured - Predicted DBH (cm)", col = "gray")

r_squared <- cor(comparison_df$dbh_cm, comparison_df$predicted_dbh_cm)^2

rmse <- sqrt(mean(comparison_df$residual^2))

cat("R²:", round(r_squared, 3), "\nRMSE:", round(rmse, 2), "cm")



-----------------------
  # Fit log-log model for Sequoia
  model_log_sequoia <- lm(log(dbh_cm) ~ log(tree_height_m), data = sequoia_data)

# Fit log-log model for Non-Sequoia
model_log_nonsequoia <- lm(log(dbh_cm) ~ log(tree_height_m), data = nonsequoia_data)

filtered_crowns <- filtered_crowns %>%
  mutate(predicted_dbh_cm = ifelse(
    is_sequoia == 1,
    exp(predict(model_log_sequoia, newdata = filtered_crowns)),
    exp(predict(model_log_nonsequoia, newdata = filtered_crowns))
  ))

matched_only$treeID <- as.character(matched_only$treeID)
filtered_crowns$treeID <- as.character(filtered_crowns$treeID)

sum(matched_only$treeID %in% filtered_crowns$treeID)

comparison_df <- matched_only %>%
  inner_join(st_drop_geometry(filtered_crowns)[, c("treeID", "predicted_dbh_cm")],
             by = "treeID")

comparison_df$residual <- comparison_df$dbh_cm - comparison_df$predicted_dbh_cm
summary(comparison_df$residual)

ggplot(comparison_df, aes(x = dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(x = "Measured DBH (cm)", y = "Predicted DBH (cm)", title = "Measured vs Predicted DBH") +
  theme_minimal()

summary(model_log_sequoia)

filtered_crowns$predicted_dbh_cm <- exp(predict(model_log_sequoia, newdata = filtered_crowns))
hist(filtered_crowns$predicted_dbh_cm, main = "Predicted DBH", xlab = "DBH (cm)", col = "lightblue")
summary(filtered_crowns$predicted_dbh_cm)

# Compute residuals
comparison_df$residual <- comparison_df$dbh_cm - comparison_df$predicted_dbh_cm

# Summary of residuals
summary(comparison_df$residual)

# Plot predicted vs. observed
plot(comparison_df$dbh_cm, comparison_df$predicted_dbh_cm,
     xlab = "Measured DBH (cm)", ylab = "Predicted DBH (cm)",
     main = "Measured vs Predicted DBH", pch = 16)
abline(a = 0, b = 1, col = "blue", lty = 2)

# If filtered_crowns is an sf object
library(mapview)
mapview(filtered_crowns, zcol = "predicted_dbh_cm")

# Export to shapefile or GeoJSON
st_write(filtered_crowns, "predicted_crowns_dbh.geojson")

plot(fitted(model_log_sequoia), resid(model_log_sequoia),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted (Log-Linear)")
abline(h = 0, lty = 2, col = "red")

---------------------
  library(minpack.lm)  # for robust nonlinear fitting

# Sequoia model
model_nls_sequoia <- nlsLM(dbh_cm ~ a * tree_height_m^b,
                           data = sequoia_data,
                           start = list(a = 1, b = 1))

# Non-Sequoia model
model_nls_nonsequoia <- nlsLM(dbh_cm ~ a * tree_height_m^b,
                              data = nonsequoia_data,
                              start = list(a = 1, b = 1))

filtered_crowns <- filtered_crowns %>%
  mutate(predicted_dbh_cm = ifelse(
    is_sequoia == 1,
    predict(model_nls_sequoia, newdata = filtered_crowns),
    predict(model_nls_nonsequoia, newdata = filtered_crowns)
  ))

matched_only$treeID <- as.character(matched_only$treeID)
filtered_crowns$treeID <- as.character(filtered_crowns$treeID)

sum(matched_only$treeID %in% filtered_crowns$treeID)

comparison_df <- matched_only %>%
  inner_join(st_drop_geometry(filtered_crowns)[, c("treeID", "predicted_dbh_cm")],
             by = "treeID")

comparison_df$residual <- comparison_df$dbh_cm - comparison_df$predicted_dbh_cm
summary(comparison_df$residual)

ggplot(comparison_df, aes(x = dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(x = "Measured DBH (cm)", y = "Predicted DBH (cm)", title = "Measured vs Predicted DBH") +
  theme_minimal()

summary(model_nls_sequoia)

--------------
install.packages("randomForest")  # if not already installed
library(randomForest)

training_subset <- training_data %>%
  filter(tree_height_m <= max(matched_only$tree_height_m, na.rm = TRUE))

summary(training_data$tree_height_m)
summary(matched_only$tree_height_m)
summary(training_subset$tree_height_m)

training_sequoia <- training_subset %>% filter(is_sequoia == 1)
training_nonsequoia <- training_subset %>% filter(is_sequoia == 0)

# Sequoia model
rf_model_sequoia <- randomForest(
  dbh_cm ~ tree_height_m,
  data = training_sequoia,
  ntree = 500,
  importance = TRUE
)

# Non-Sequoia model
rf_model_nonsequoia <- randomForest(
  dbh_cm ~ tree_height_m,
  data = training_nonsequoia,
  ntree = 500,
  importance = TRUE
)

filtered_crowns <- filtered_crowns %>%
  mutate(predicted_dbh_cm_rf = ifelse(
    is_sequoia == 1,
    predict(rf_model_sequoia, newdata = .),
    predict(rf_model_nonsequoia, newdata = .)
  ))

comparison_rf <- matched_only %>%
  mutate(treeID = as.character(treeID)) %>%
  inner_join(
    filtered_crowns %>%
      st_drop_geometry() %>%
      select(treeID, predicted_dbh_cm_rf),
    by = "treeID"
  ) %>%
  mutate(residual_rf = dbh_cm - predicted_dbh_cm_rf)

# Calculate metrics
rmse_rf <- sqrt(mean(comparison_rf$residual_rf^2))
mae_rf <- mean(abs(comparison_rf$residual_rf))

cat("Random Forest RMSE:", rmse_rf, "\n")
cat("Random Forest MAE:", mae_rf, "\n")

# Measured vs Predicted
plot(comparison_rf$dbh_cm, comparison_rf$predicted_dbh_cm_rf,
     xlab = "Measured DBH (cm)", ylab = "Predicted DBH (cm)",
     main = "Measured vs Predicted DBH (Random Forest)", pch = 19)
abline(a = 0, b = 1, col = "blue", lty = 2)

# Residuals
hist(comparison_rf$residual_rf, breaks = 15, col = "gray",
     main = "Residuals (Random Forest)", xlab = "Residual (Measured - Predicted)")

------------------
  "C:/Users/User/Desktop/RandomForest/20250806_RandomForest_Variables.csv"
  
RF_variables <- read.csv("C:/Users/User/Desktop/RandomForest/20250806_RandomForest_Variables.csv")  # or use sf::st_read if spatial

RF_variables <- RF_variables %>%
  rename(treeID = TREEID)
# Filter RF_variables to include only rows with TREEID in matched_only
RF_filtered <- RF_variables[RF_variables$treeID %in% matched_only$treeID, ]

RF_filtered <- RF_filtered %>%
  left_join(matched_only %>% select(treeID, dbh_cm), by = "treeID")

# Drop unused or ID columns (e.g., TREEID, ZONE_CODE, COUNT)
RF_variables <- subset(RF_filtered, select = -c(treeID, ZONE_CODE, COUNT))

RF_variables <- RF_variables %>%
  rename(tree_height_m = tree_height_m_1)

matched_only
# Check for missing values
summary(RF_variables)

install.packages("randomForest")
library(randomForest)
# Fit the model to predict DBH
set.seed(123)

rf_model <- randomForest(
  dbh_cm ~ tree_height_m + is_sequoia + crown_area_m2 + basal_area_m2 +
    radius_m + p95 + p75 + p50 + mean_z + max_z +
    sd_z + cv_z + point_density + slope_mean + aspect_circmean,
  data = RF_variables,
  importance = TRUE,
  ntree = 500
)

# View results
print(rf_model)
varImpPlot(rf_model)

summary(RF_variables$dbh_cm)
hist(RF_variables$dbh_cm)
any(is.na(RF_variables$dbh_cm))

sum(is.na(RF_variables$dbh_cm))


rf_model <- randomForest(
  dbh_cm ~ tree_height_m_1 + crown_area_m2 + basal_area_m2 + radius_m + 
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z + point_density +
    slope_mean + aspect_circmean,
  data = RF_variables,
  importance = TRUE,
  ntree = 500
)

boxplot(RF_variables$dbh_cm)
boxplot(RF_variables$tree_height_m_1)


# Sequoia trees
sequoia_data <- RF_variables[RF_variables$is_sequoia == 1 & !is.na(RF_variables$dbh_cm), ]

# Non-Sequoia trees
non_sequoia_data <- RF_variables[RF_variables$is_sequoia == 0 & !is.na(RF_variables$dbh_cm), ]

rf_sequoia <- randomForest(
  dbh_cm ~ tree_height_m + crown_area_m2 + basal_area_m2 + radius_m +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z + point_density +
    slope_mean + aspect_circmean,
  data = sequoia_data,
  importance = TRUE,
  ntree = 500
)

print(rf_sequoia)
print(rf_nonsequoia)

rf_nonsequoia <- randomForest(
  dbh_cm ~ tree_height_m_1 + crown_area_m2 + basal_area_m2 + radius_m +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z + point_density +
    slope_mean + aspect_circmean,
  data = non_sequoia_data,
  importance = TRUE,
  ntree = 500
)



# Separate them first:
new_sequoia <- filtered_crowns[filtered_crowns$is_sequoia == 1, ]
new_nonsequoia <- filtered_crowns[filtered_crowns$is_sequoia == 0, ]

new_sequoia$predicted_dbh <- predict(rf_sequoia, new_sequoia)

nrow(sequoia_data)
summary(sequoia_data$dbh_cm)

