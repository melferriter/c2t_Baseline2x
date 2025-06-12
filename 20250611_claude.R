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

-------------------------------------------------------------------------
  
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

# Load and filter training data
training_data_raw <- read.csv("C:/Users/User/Desktop/TrainingData_issequoia_2.csv")
training_data <- training_data_raw %>%
  filter(
    tree_height_m >= 15,                 # Remove trees shorter than 15 m
    dbh_cm >= 10, dbh_cm <= 400,         # Remove DBH outliers
    !(tree_height_m > 50 & is_sequoia == 0)  # Likely Sequoias misclassified as non-sequoia
  )


# Split data by species
sequoia_data     <- training_data %>% filter(is_sequoia == 1)
nonsequoia_data  <- training_data %>% filter(is_sequoia == 0)

# Fit separate models for each species
model_seq <- lm(dbh_cm ~ tree_height_m + I(pmax(tree_height_m - 50, 0)), data = sequoia_data)
model_nonseq <- lm(dbh_cm ~ tree_height_m + I(pmax(tree_height_m - 50, 0)), data = nonsequoia_data)

# Visualize fits
ggplot(training_data, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia))) +
  geom_point() +
  geom_smooth(data = sequoia_data, method = "lm", formula = y ~ x + I(pmax(x - 50, 0)), se = TRUE, color = "darkgreen") +
  geom_smooth(data = nonsequoia_data, method = "lm", formula = y ~ x + I(pmax(x - 50, 0)), se = TRUE, color = "orange") +
  labs(title = "Separate Piecewise Models for Sequoia and Non-Sequoia") +
  theme_minimal()

# Step 1: Initial filtering to remove clearly invalid or tiny trees
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  filter(
    tree_height_m >= 20,           # Remove very small trees
    crown_area_m2 >= 20            # Remove tiny crowns (optional but helpful)
  )

# Step 2: Predict is_sequoia using reasonable biological thresholds
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  mutate(is_sequoia = ifelse(tree_height_m > 50 & crown_area_m2 > 100, 1, 0))

# Step 3: Filter out biologically implausible non-sequoias
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  filter(!(tree_height_m > 50 & is_sequoia == 0))


cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  mutate(predicted_dbh_cm = ifelse(is_sequoia == 1,
                                   predict(model_seq, newdata = .),
                                   predict(model_nonseq, newdata = .)))

# Read and prepare control points
controlpoints <- st_read("C:/Users/User/Desktop/ControlPoints_1.shp")
names(controlpoints)[names(controlpoints) == "DBH"] <- "dbh_cm"
names(controlpoints)[names(controlpoints) == "Hgt"] <- "tree_height_m"
st_crs(controlpoints) <- 3857
controlpoints <- st_transform(controlpoints, crs = st_crs(cloud2trees_ans_c_nosnag$treetops_sf))

# Match predictions to control points
treetops_xy <- st_coordinates(cloud2trees_ans_c_nosnag$treetops_sf)[, 1:2]
control_xy  <- st_coordinates(controlpoints)[, 1:2]
knn_result <- get.knnx(data = treetops_xy, query = control_xy, k = 1)
nearest_ids <- knn_result$nn.index[, 1]
distances <- knn_result$nn.dist[, 1]
height_diff <- abs(controlpoints$tree_height_m - cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids])

valid_match <- which(distances <= 10 & height_diff <= 30)

valid_data <- data.frame(
  observed_dbh_cm   = controlpoints$dbh_cm[valid_match],
  predicted_dbh_cm  = cloud2trees_ans_c_nosnag$treetops_sf$predicted_dbh_cm[nearest_ids[valid_match]],
  observed_height_m = controlpoints$tree_height_m[valid_match],
  predicted_height_m = cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m[nearest_ids[valid_match]],
  distance_m        = distances[valid_match],
  height_diff_m     = height_diff[valid_match]
)

# Evaluate model performance
valid_data$residual_dbh <- valid_data$observed_dbh_cm - valid_data$predicted_dbh_cm
rmse <- sqrt(mean(valid_data$residual_dbh^2))
bias <- mean(valid_data$residual_dbh)
r_squared <- 1 - sum(valid_data$residual_dbh^2) / sum((valid_data$observed_dbh_cm - mean(valid_data$observed_dbh_cm))^2)

cat("RMSE:", rmse, "cm\n")
cat("Bias:", bias, "cm\n")
cat("R²:", r_squared, "\n")

# Plot observed vs predicted
ggplot(valid_data, aes(x = observed_dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Observed vs. Predicted DBH", x = "Observed DBH (cm)", y = "Predicted DBH (cm)") +
  theme_minimal()

------------------------------------------------------------
  
summary(training_data)
summary(cloud2trees_ans_c_nosnag$treetops_sf$tree_height_m)


  
  