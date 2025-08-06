
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

cloud2trees_ans <- cloud2trees::cloud 
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


# Extract extent from LAS and turn into an sf polygon
las_extent <- st_as_sfc(st_bbox(clipped), crs = st_crs(controlpoints))
# Filter control points that fall within the LAS extent
controlpoints_in_las <- controlpoints[st_intersects(controlpoints, las_extent, sparse = FALSE), ]

matched_points <- st_join(controlpoints, filtered_crowns, join = st_within)



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
filtered_crowns$treetops_sf <- filtered_crowns %>%
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

