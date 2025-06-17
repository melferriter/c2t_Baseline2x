
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
library(mgcv)
pacman::p_load(ggplot2, rgeos, propagate, dplyr, ggpubr, gridExtra)

rm(list = ls(globalenv()))
---------------------------------------------------

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
  , estimate_tree_cbh = TRUE
  , cbh_estimate_missing_cbh = F
)

#cloud2trees_ans_c_nosnag <- cloud2trees_ans_c
#cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c$treetops_sf[cloud2trees_ans_c$treetops_sf$crown_area_m2 > 30, ]

# Step 2: Predict is_sequoia using reasonable biological thresholds
cloud2trees_ans_c$treetops_sf <- cloud2trees_ans_c$treetops_sf %>%
  mutate(is_sequoia = ifelse(crown_area_m2 > 115 & tree_height_m > 45, 1, 0))


#cloud2trees_ans_c$treetops_sf[2]

ggplot(cloud2trees_ans_c$treetops_sf %>% filter(is_sequoia == 1), aes(x = tree_height_m)) +
  geom_histogram(binwidth = 5, fill = "darkgreen", color = "white") +
  labs(title = "C2T Sequoia Height Distribution",
       x = "Height (m)",
       y = "Count") +
  theme_minimal()


ggplot(cloud2trees_ans_c$treetops_sf %>% filter(is_sequoia == 0), aes(x = tree_height_m)) +
  geom_histogram(binwidth = 5, fill = "darkgreen", color = "white") +
  labs(title = "C2T Non Sequoia Height Distribution w/ Snags",
       x = "Height (m)",
       y = "Count") +
  theme_minimal()

ggplot(cloud2trees_ans_c_nosnag$treetops_sf %>% filter(is_sequoia == 0), aes(x = tree_height_m)) +
  geom_histogram(binwidth = 5, fill = "darkgreen", color = "white") +
  labs(title = "C2T Non Sequoia Height Distribution w/out snags",
       x = "Height (m)",
       y = "Count") +
  theme_minimal()

cloud2trees_ans_c_nosnag <- cloud2trees_ans_c
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c$treetops_sf[cloud2trees_ans_c$treetops_sf$tree_height_m > 10, ]

# Extract coordinates and combine with attribute data
#treetops_df <- cloud2trees_ans_c$treetops_sf %>%
#  mutate(X = st_coordinates(.)[, 1],
#         Y = st_coordinates(.)[, 2]) %>%
#  st_drop_geometry()  # Remove geometry column for CSV export

# Export to CSV
#write.csv(treetops_df, "C:/Users/User/Desktop/treetops_export.csv", row.names = FALSE)

--------------------------
  
# Load and filter training data
#training_data <- read.csv("C:/Users/User/Desktop/TrainingData_issequoia_2.csv")
training_data <- read.csv("C:/Users/User/Desktop/202150614_Monarchs_Hwd.csv")

training_data <- training_data %>%
  filter(
    !(tree_height_m < 50 & is_sequoia == 1 & dbh_cm >200)  # Likely Sequoias misclassified as non-sequoia
  )

training_data <- training_data %>% drop_na()

ggplot(training_data, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia==1))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("forestgreen", "orange"), labels = c("Sequoia", "Not Sequoia")) +
  labs(title = "Training Data: Sequoia vs Non-Sequoia", color = "Type") +
  theme_minimal()

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

shapiro.test(sequoia_data$tree_height_m)
shapiro.test(sequoia_data$dbh_cm)
shapiro.test(nonsequoia_data$tree_height_m)
shapiro.test(nonsequoia_data$dbh_cm)
-------------------
# Read and prepare control points
controlpoints <- st_read("C:/Users/User/Desktop/ControlPoints_3.shp")
names(controlpoints)[names(controlpoints) == "DBH"] <- "dbh_cm"
names(controlpoints)[names(controlpoints) == "Hgt"] <- "tree_height_m"
st_crs(controlpoints) <- 3857
controlpoints <- st_transform(controlpoints, crs = st_crs(cloud2trees_ans_c$treetops_sf))

ggplot(controlpoints, aes(x = tree_height_m, y = dbh_cm, color = factor(is_sequoia))) +
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
#GAM 

# Fit GAMs separately for each group
model_gam_sequoia    <- gam(dbh_cm ~ s(tree_height_m), data = sequoia_data, method = "REML")
model_gam_nonsequoia <- gam(dbh_cm ~ s(tree_height_m), data = nonsequoia_data, method = "REML")

# Predict DBH for all treetops
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
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

# Extract numeric prefix from treeID
cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  mutate(tree_prefix = str_extract(treeID, "^[^_]+"))

# Filter by matching prefix with controlpoints$c2t_id
filtered_treetops <- cloud2trees_ans_c_nosnag$treetops_sf %>%
  filter(tree_prefix %in% controlpoints$c2t_id)


cloud2trees_ans_c_nosnag$treetops_sf %>% filter(controlpoints$c2t_id)
control_id <- controlpoints$c2t_id

# Join controlpoints with predicted DBH using matching ID
matched_data <- controlpoints %>%
  mutate(
    predicted_dbh_cm = cloud2trees_ans_c$treetops_sf$predicted_dbh_cm[control_id],
    residual_dbh = dbh_cm - predicted_dbh_cm,
    crown_area_m2 = cloud2trees_ans_c$treetops_sf$crown_area_m2[control_id],
    comp_dist_to_nearest_m = cloud2trees_ans_c$treetops_sf$comp_dist_to_nearest_m[control_id],
    slope_deg = cloud2trees_ans_c$treetops_sf$slope_deg[control_id],
    tree_cbh_m = cloud2trees_ans_c$treetops_sf$tree_cbh_m[control_id]
  )

--------------
matched_data$residual_dbh <- matched_data$dbh_cm - matched_data$predicted_dbh_cm
#resid_model <- lm(residual_dbh ~ crown_area_m2, data = matched_data)
resid_model <- lm(residual_dbh ~ crown_area_m2 +tree_cbh_m, data = matched_data)
#resid_model <- lm(residual_dbh ~ crown_area_m2 * nearest_tree_dist_m, data = matched_data)
#resid_gam <- gam(residual_dbh ~ s(crown_area_m2) + s(comp_dist_to_nearest_m), data = matched_data, method = "REML")




cloud2trees_ans_c$treetops_sf <- cloud2trees_ans_c$treetops_sf %>%
  mutate(predicted_dbh_cm_adjusted = predicted_dbh_cm + predict(resid_model, newdata = .))
#cloud2trees_ans_c$treetops_sf <- cloud2trees_ans_c$treetops_sf %>%
#  mutate(predicted_dbh_cm_adjusted = predicted_dbh_cm + predict(resid_gam, newdata = .))

matched_data$predicted_dbh_cm_adjusted <- cloud2trees_ans_c$treetops_sf$predicted_dbh_cm_adjusted[control_id]

matched_data$residual_adjusted <- matched_data$dbh_cm - matched_data$predicted_dbh_cm_adjusted

rmse <- sqrt(mean(matched_data$residual_adjusted^2))
bias <- mean(matched_data$predicted_dbh_cm_adjusted - matched_data$dbh_cm)
r2 <- 1 - sum((matched_data$dbh_cm - matched_data$predicted_dbh_cm_adjusted)^2) / 
  sum((matched_data$dbh_cm - mean(matched_data$dbh_cm))^2)

cat("Adjusted DBH Prediction:/n")
cat("  RMSE:", round(rmse, 2), "cm/n")
cat("  Bias:", round(bias, 2), "cm/n")
cat("  R²:", round(r2, 3), "/n")

ggplot(matched_data, aes(x = dbh_cm, y = predicted_dbh_cm_adjusted)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Observed vs Adjusted Predicted DBH",
    x = "Observed DBH (cm)",
    y = "Adjusted Predicted DBH (cm)"
  ) +
  theme_minimal()


--------------

controlpoints %>% cloud2trees_ans_c$treetops_sf[controlpoints$c2t_id]
cloud2trees_ans_c$treetops_sf[744,]

# Accuracy function
evaluate_accuracy <- function(matched_data, label = "") {
  rmse <- sqrt(mean(matched_data$residual_adjusted^2, na.rm = TRUE))
  bias <- mean(matched_data$residual_adjusted, na.rm = TRUE)
  r2 <- 1 - sum((matched_data$dbh_cm - matched_data$predicted_dbh_cm_adjusted)^2, na.rm = TRUE) /
    sum((matched_data$dbh_cm - mean(matched_data$dbh_cm))^2, na.rm = TRUE)

  
  cat("/n", label, "Accuracy Metrics:/n")
  cat("  RMSE:", round(rmse, 2), "cm/n")
  cat("  Bias:", round(bias, 2), "cm/n")
  cat("  R²:", round(r2, 3), "/n")
  
  return(invisible(data.frame(RMSE = rmse, Bias = bias, R2 = r2)))
}

# Evaluate overall accuracy
evaluate_accuracy(matched_data, "All Trees")


# Evaluate separately for Sequoia vs Non-Sequoia
evaluate_accuracy(filter(matched_data, is_sequoia == 1), "Sequoia")
evaluate_accuracy(filter(matched_data, is_sequoia == 0), "Non-Sequoia")

ggplot(matched_data, aes(x = predicted_dbh_cm_adjusted, y = residual_adjusted, color = factor(is_sequoia))) +
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted DBH (cm)", y = "Residual (cm)", color = "Sequoia?")


ggplot(matched_data, aes(x = dbh_cm, y = predicted_dbh_cm_adjusted, color = factor(is_sequoia))) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Observed DBH (cm)",
    y = "Predicted DBH (cm)",
    title = "Observed vs Predicted DBH",
    color = "Sequoia?"
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "darkgreen"), labels = c("Non-Sequoia", "Sequoia")) +
  theme_minimal()

plot(resid_gam, pages = 1, residuals = TRUE)

  
----------------------

# 1. Load your DEM (must be aligned and in the same CRS as your crown polygons)
dem <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xCross_20241022194753_DSM.tif")  # Replace with your DEM path

# 2. Calculate slope in degrees
slope_raster <- terrain(dem, v = "slope", unit = "degrees")

# 3. Extract mean slope for each crown polygon
# Assume crowns_sf is in cloud2trees_ans_c$crowns_sf and is an sf POLYGON
crowns_sf <- cloud2trees_ans_c$crowns_sf
crowns_spat <- vect(crowns_sf)  # convert sf to terra-compatible SpatVector

# Extract mean slope per polygon
mean_slope_df <- terra::extract(slope_raster, cloud2trees_ans_c$treetops_sf, fun = mean, na.rm = TRUE)
# This returns a dataframe with ID and mean slope

# 4. Add the mean slope to the original crowns data
cloud2trees_ans_c$treetops_sf$slope_deg <- mean_slope_df$slope


--------------

  
