## pkgbuild helps us check for Rtools
install.packages("pkgbuild")
# check for Rtools which is required to build packages
pkgbuild::check_build_tools(debug = TRUE)
## remotes helps us get packages hosted on github
install.packages("remotes")
## install lasR from the r-univers
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
## install TreeLS from github
remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)
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
pacman::p_load(ggplot2, rgeos, propagate, dplyr, ggpubr, gridExtra)

rm(list = ls(globalenv()))
------------------
  
cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las")

## Default ITD window size functions====
itd_tuning_ans <- itd_tuning(input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las")
best_ws <- itd_tuning_ans$ws_fn_list$lin_fn
 

i = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las"  

cloud2trees_ans_c <- cloud2trees::cloud2trees(
  output_dir = tempdir()
  , input_las_dir = i
  , dtm_res_m = 0.5
  , ws = best_ws
  , estimate_tree_dbh = TRUE
  , estimate_tree_type = TRUE
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = TRUE
  , cbh_tree_sample_n = 555
  , cbh_estimate_missing_cbh = TRUE
)

cloud2trees_ans_c_nosnag <- cloud2trees_ans_c
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c$treetops_sf[cloud2trees_ans_c$treetops_sf$crown_area_m2 > 2, ]

-------------------------------
training_data_1 <- read.csv("C:/Users/User/Desktop/TrainingData_issequoia_2.csv")

cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  filter(tree_height_m > 10)
 
training_data_filtered <- training_data_1 %>% 
  filter(tree_height_m > 10)

model <- lm(dbh_cm ~ tree_height_m + I(pmax(tree_height_m - 50, 0)), data = training_data_filtered)
 
# Check for outliers and data quality
ggplot(training_data_1, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia))) +
  geom_point() +
  geom_smooth(method = "loess") +
  #facet_wrap(~is_sequoia) +
  theme_minimal()

# Simple but effective for your data pattern
piecewise_model <- lm(dbh_cm ~ tree_height_m + I(pmax(tree_height_m - 50, 0)), 
                      data = training_data_1)

# Plot the fitted line over your data to see if it looks reasonable
ggplot(training_data_filtered, aes(x = tree_height_m, y = dbh_cm)) +
  geom_point(aes(color = factor(is_sequoia))) +
  geom_smooth(method = "lm", formula = y ~ x + I(pmax(x - 50, 0)), se = TRUE) +
  labs(title = "Piecewise Linear Model Fit") +
  theme_minimal()

# Predict DBH for all trees using the piecewise model
cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm <- predict(
  piecewise_model, 
  newdata = cloud2trees_ans_c_nosnag$treetops_sf
)

# Check for any negative predictions (shouldn't happen but good to verify)
summary(cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm)

# Read control points and rename columns
controlpoints <- st_read("C:/Users/User/Desktop/ControlPoints_1.shp")
names(controlpoints)[names(controlpoints) == "DBH"] <- "dbh_cm"
names(controlpoints)[names(controlpoints) == "Hgt"] <- "tree_height_m"

st_crs(controlpoints) <- 3857

controlpoints_proj <- st_transform(controlpoints, 32611)

# Ensure both layers are in the same CRS
controlpoints <- st_transform(controlpoints, crs = st_crs(cloud2trees_ans_c_nosnag$treetops_sf))


treetops_xy <- st_coordinates(cloud2trees_ans_c_nosnag$treetops_sf)[, 1:2]
control_xy  <- st_coordinates(controlpoints)[, 1:2]

# Nearest neighbor search
knn_result <- get.knnx(data = treetops_xy, query = control_xy, k = 1)
nearest_ids <- knn_result$nn.index[, 1]
distances <- knn_result$nn.dist[, 1]

# Height difference filtering
height_diff <- abs(controlpoints$tree_height_m - cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids])

# Accept only matches within 5 meters distance and height diff
valid_match <- which(distances <= 10 & height_diff <= 30)

controlpoints$predicted_dbh_cm <- NA_real_
controlpoints$tree_height_m <- NA_real_

controlpoints$predicted_dbh_cm[valid_match] <- cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm[nearest_ids[valid_match]]
controlpoints$tree_height_m[valid_match]     <- cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids[valid_match]]

summary(distances)
summary(height_diff)

# Subset only trees that have both field DBH and predicted DBH
eval_subset <- subset(cloud2trees_ans_c_nosnag$treetops_sf, !is.na(dbh_cm) & !is.na(predicted_dbh_cm))

# Calculate residuals
eval_subset$residual_dbh <- eval_subset$dbh_cm - eval_subset$predicted_dbh_cm

# R-squared
r_squared <- 1 - sum(eval_subset$residual_dbh^2) / sum((eval_subset$dbh_cm - mean(eval_subset$dbh_cm))^2)
r_squared








# Try progressively looser criteria
distance_thresholds <- c(3, 5, 7, 10, 15)
height_thresholds <- c(5, 10, 15, 20, 30)

for(d_thresh in distance_thresholds) {
  for(h_thresh in height_thresholds) {
    valid_match_test <- which(distances <= d_thresh & height_diff <= h_thresh)
    cat("Distance ≤", d_thresh, "m, Height diff ≤", h_thresh, "m: ", 
        length(valid_match_test), "matches\n")
  }
}

# Look at the distribution of distances and height differences
summary(distances)
summary(height_diff)

# Plot to see the distribution
hist(distances, main = "Distance to nearest UAV tree", xlab = "Distance (m)")
hist(height_diff, main = "Height difference", xlab = "Height difference (m)")

# Are there systematic height differences?
plot(controlpoints_proj$tree_height_m, 
     cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids],
     xlab = "Field height (m)", ylab = "UAV height (m)")
abline(0, 1, col = "red")

# Are there positioning issues?
plot(distances, height_diff, 
     xlab = "Distance (m)", ylab = "Height diff (m)")


# Your existing matching code but with tighter criteria
knn_result <- get.knnx(data = treetops_xy, query = control_xy, k = 1)
nearest_ids <- knn_result$nn.index[, 1]
distances <- knn_result$nn.dist[, 1]

# Much stricter matching criteria
height_diff <- abs(controlpoints_proj$tree_height_m - 
                     cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids])
valid_match <- which(distances <= 15 & height_diff <= 30)  # Tighter thresholds

# Apply predictions only to valid matches
controlpoints_proj$predicted_dbh_cm <- NA_real_
controlpoints_proj$predicted_dbh_cm[valid_match] <- 
  cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm[nearest_ids[valid_match]]


# Calculate metrics for matched points only
valid_data <- controlpoints_proj[!is.na(controlpoints_proj$predicted_dbh_cm), ]
valid_data$residual_dbh <- valid_data$dbh_cm - valid_data$predicted_dbh_cm

# Compute metrics
rmse <- sqrt(mean(valid_data$residual_dbh^2))
bias <- mean(valid_data$residual_dbh)
r2 <- cor(valid_data$dbh_cm, valid_data$predicted_dbh_cm)^2

cat("RMSE:", rmse, "cm\n")
cat("Bias:", bias, "cm\n") 
cat("R²:", r2, "\n")
cat("n matched:", nrow(valid_data), "\n")

# Plot predicted vs observed
ggplot(valid_data, aes(x = dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Observed DBH (cm)", y = "Predicted DBH (cm)") +
  theme_minimal()
# 

# Simple bias correction
site_bias <- -89.14  # Your calculated bias
corrected_predictions <- cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm + site_bias

# Update your treetops data
cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm_corrected <- corrected_predictions

# Re-evaluate with corrected values
controlpoints_proj$predicted_dbh_cm_corrected <- 
  controlpoints_proj$predicted_dbh_cm + site_bias

# ADD THIS LINE - update valid_data with corrected predictions
valid_data$predicted_dbh_cm_corrected <- valid_data$predicted_dbh_cm + site_bias

# Check new metrics
new_residuals <- valid_data$dbh_cm - valid_data$predicted_dbh_cm_corrected
new_rmse <- sqrt(mean(new_residuals^2))
new_bias <- mean(new_residuals)

cat("Corrected RMSE:", new_rmse, "cm\n")
cat("Corrected Bias:", new_bias, "cm\n")

# Plot predicted vs observed WITH site correction
ggplot(valid_data, aes(x = dbh_cm, y = predicted_dbh_cm_corrected)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Observed DBH (cm)", y = "Corrected Predicted DBH (cm)", 
       title = "Site-Corrected DBH Predictions") +
  theme_minimal()