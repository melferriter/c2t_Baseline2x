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

i = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las"  

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
  , estimate_tree_cbh = FALSE
  , cbh_estimate_missing_cbh = F
)

# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$crowns_sf, mapping = ggplot2::aes(fill = dbh_m)) + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$treetops_sf, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

cloud2trees_ans_c_nosnag <- cloud2trees_ans_c
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c$treetops_sf[cloud2trees_ans_c$treetops_sf$crown_area_m2 > 2, ]

# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c_nosnag$crowns_sf, mapping = ggplot2::aes(fill = dbh_m)) + 
  ggplot2::geom_sf(data = cloud2trees_ans_c_nosnag$treetops_sf, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")


-------------------------------

# Load and filter training data
training_data_raw <- read.csv("C:/Users/User/Desktop/TrainingData_issequoia_2.csv")
training_data <- training_data_raw %>%
  filter(
    tree_height_m >= 15, tree_height_m <= 95,                 # Remove trees shorter than 15 m
    dbh_cm >= 10, dbh_cm <= 350,         # Remove DBH outliers
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
  geom_point(position = position_jitter(width = 0, height = 2)) +
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

# Split UAV treetops
treetops_seq     <- cloud2trees_ans_c_nosnag$treetops_sf %>% filter(is_sequoia == 1)
treetops_nonseq  <- cloud2trees_ans_c_nosnag$treetops_sf %>% filter(is_sequoia == 0)

# Split control points
control_seq      <- controlpoints %>% filter(is_sequoia == 1)
control_nonseq   <- controlpoints %>% filter(is_sequoia == 0)

# Coordinates
xy_treetops_seq    <- st_coordinates(treetops_seq)[, 1:2]
xy_treetops_nonseq <- st_coordinates(treetops_nonseq)[, 1:2]
xy_control_seq     <- st_coordinates(control_seq)[, 1:2]
xy_control_nonseq  <- st_coordinates(control_nonseq)[, 1:2]

# Nearest neighbors (within group)
knn_seq     <- get.knnx(data = xy_treetops_seq, query = xy_control_seq, k = 1)
knn_nonseq  <- get.knnx(data = xy_treetops_nonseq, query = xy_control_nonseq, k = 1)

# SEQUOIA
dist_seq <- knn_seq$nn.dist[, 1]
id_seq   <- knn_seq$nn.index[, 1]
height_diff_seq <- abs(control_seq$tree_height_m - treetops_seq$tree_height_m[id_seq])
valid_seq <- which(dist_seq <= 10 & height_diff_seq <= 30)

# NON-SEQUOIA
dist_nonseq <- knn_nonseq$nn.dist[, 1]
id_nonseq   <- knn_nonseq$nn.index[, 1]
height_diff_nonseq <- abs(control_nonseq$tree_height_m - treetops_nonseq$tree_height_m[id_nonseq])
valid_nonseq <- which(dist_nonseq <= 10 & height_diff_nonseq <= 30)

valid_data <- bind_rows(
  data.frame(
    is_sequoia        = 1,
    observed_dbh_cm   = control_seq$dbh_cm[valid_seq],
    predicted_dbh_cm  = treetops_seq$predicted_dbh_cm[id_seq[valid_seq]],
    observed_height_m = control_seq$tree_height_m[valid_seq],
    predicted_height_m = treetops_seq$tree_height_m[id_seq[valid_seq]],
    distance_m         = dist_seq[valid_seq],
    height_diff_m      = height_diff_seq[valid_seq]
  ),
  data.frame(
    is_sequoia        = 0,
    observed_dbh_cm   = control_nonseq$dbh_cm[valid_nonseq],
    predicted_dbh_cm  = treetops_nonseq$predicted_dbh_cm[id_nonseq[valid_nonseq]],
    observed_height_m = control_nonseq$tree_height_m[valid_nonseq],
    predicted_height_m = treetops_nonseq$tree_height_m[id_nonseq[valid_nonseq]],
    distance_m         = dist_nonseq[valid_nonseq],
    height_diff_m      = height_diff_nonseq[valid_nonseq]
  )
)

# Add residual column
valid_data$residual_dbh <- valid_data$observed_dbh_cm - valid_data$predicted_dbh_cm


model_metrics <- valid_data %>%
  group_by(is_sequoia) %>%
  summarise(
    RMSE  = sqrt(mean(residual_dbh^2)),
    Bias  = mean(residual_dbh),
    R2    = 1 - sum(residual_dbh^2) / sum((observed_dbh_cm - mean(observed_dbh_cm))^2),
    n     = n()
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

ggplot(training_data, aes(x = t, fill = factor(t$tree_height_m))) +
  geom_histogram(binwidth = 10) +
  theme_minimal()

t <- training_data%>% filter(is_sequoia == 1)
t$tree_height_m

--------------
valid_data$is_sequoia <- cloud2trees_ans_c_nosnag$treetops_sf$is_sequoia[nearest_ids[valid_match]]

valid_seq     <- valid_data %>% filter(is_sequoia == 1)
valid_nonseq  <- valid_data %>% filter(is_sequoia == 0)

evaluate_model <- function(data, label) {
  data$residual_dbh <- data$observed_dbh_cm - data$predicted_dbh_cm
  rmse <- sqrt(mean(data$residual_dbh^2))
  bias <- mean(data$residual_dbh)
  r2 <- 1 - sum(data$residual_dbh^2) / sum((data$observed_dbh_cm - mean(data$observed_dbh_cm))^2)
  
  cat("\n", label, "\n")
  cat("  RMSE:", round(rmse, 2), "cm\n")
  cat("  Bias:", round(bias, 2), "cm\n")
  cat("  R²:", round(r2, 4), "\n")
  
  # Optional plot
  ggplot(valid_seq, aes(x = observed_dbh_cm, y = predicted_dbh_cm)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(title = paste(label, "Observed vs. Predicted DBH"),
         x = "Observed DBH (cm)", y = "Predicted DBH (cm)") +
    theme_minimal()
}

# Plot observed vs predicted
ggplot(valid_nonseq, aes(x = observed_dbh_cm, y = predicted_dbh_cm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Observed vs. Predicted DBH", x = "Observed DBH (cm)", y = "Predicted DBH (cm)") +
  theme_minimal()





  
  