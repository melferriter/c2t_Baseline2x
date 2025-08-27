###################Melissa Ferriter 8/9/2025#################


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
library(data.table)
library(randomForest)
pacman::p_load(ggplot2, rgeos, propagate, dplyr, ggpubr, gridExtra)

rm(list = ls(globalenv()))
rm(list = ls())
gc()   
graphics.off() 
------------------------------------
  
#read in the SFM and 3XCross Las files and clip to ROI
roi <- st_read("C:/Users/User/Desktop/LidarClip/LidarExtent.shp", crs = 32611)
roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")
Las3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las", filter = "-set_withheld_flag 0")
sfm <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Raw/20241022/Sam/FreemanCreekGrove_south_group1_densified_point_cloud.las")
Las3x_clipped <- clip_roi(las, roi_2d)
sfm_clipped <- clip_roi(sfm, roi_2d)
writeLAS(Las3x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross _roi_20241022194753_20250809.las")
writeLAS(sfm_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/SFM_roi_20241022194753_20250809.las")


sfm <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Raw/20241022/Sam/FreemanCreekGrove_south_group1_densified_point_cloud.las")
side3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las", filter = "-set_withheld_flag 0")


crs(Las3x_clipped)
crs(sfm_clipped)

sfm_clipped <- spTransform(sfm_clipped, crs(Las3x_clipped))

names(Las3x_clipped@data)
names(sfm_clipped@data)

Las3x_clipped@data$R <- NA
Las3x_clipped@data$G <- NA
Las3x_clipped@data$B <- NA

# Ensure column order matches
sfm_clipped@data <- as.data.table(sfm_clipped@data[, names(Las3x_clipped@data), with = FALSE])

sfm_clipped@data <- sfm_clipped@data[, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

header_template <- Las3x_clipped@header
sfm_clipped@header <- header_template

# Define all standard LAS integer columns
las_int_cols <- c("ReturnNumber", "NumberOfReturns", "Classification", "ScanDirectionFlag",
                  "EdgeOfFlightline", "Synthetic_flag", "Keypoint_flag", "Withheld_flag",
                  "ScanAngleRank", "UserData", "PointSourceID", "ScannerChannel")

# Ensure both objects are using correct types
for (col in las_int_cols) {
  if (col %in% names(Las3x_clipped@data)) {
    Las3x_clipped@data[[col]] <- as.integer(Las3x_clipped@data[[col]])
  }
  if (col %in% names(sfm_clipped@data)) {
    sfm_clipped@data[[col]] <- as.integer(sfm_clipped@data[[col]])
  }
}

# Convert flag fields to logical
flag_fields <- c("Synthetic_flag", "Keypoint_flag", "Withheld_flag", "Overlap_flag")

for (col in flag_fields) {
  if (col %in% names(Las3x_clipped@data)) {
    Las3x_clipped@data[[col]] <- as.logical(Las3x_clipped@data[[col]])
  }
  if (col %in% names(sfm_clipped@data)) {
    sfm_clipped@data[[col]] <- as.logical(sfm_clipped@data[[col]])
  }
}

Las3x_clipped@data$R <- NULL
Las3x_clipped@data$G <- NULL
Las3x_clipped@data$B <- NULL

sfm_clipped@data$R <- NULL
sfm_clipped@data$G <- NULL
sfm_clipped@data$B <- NULL


las_combined <- rbind(Las3x_clipped, sfm_clipped)

las_combined@data$ReturnNumber[las_combined@data$ReturnNumber == 0] <- 1
las_combined@data$NumberOfReturns[las_combined@data$NumberOfReturns == 0] <- 1

header <- las_combined@header
header@PHB$XScaleFactor <- 0.001
header@PHB$YScaleFactor <- 0.001
header@PHB$ZScaleFactor <- 0.001

las_combined@header <- header

plot(las_combined)
summary(las_combined)

writeLAS(las_combined, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/SFM_3xSide_20241022194753_20250809.las")
------------------------------------

i = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM_3xSide_20241022194753_20250825.las"

cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = i)

## Default ITD window size functions====
itd_tuning_ans <- itd_tuning(input_las_dir = i)
itd_tuning_ans$plot_samples
best_ws <- itd_tuning_ans$ws_fn_list$lin_fn


sequoia_ws_lo = function(x) {(x * 0.023) + 2}

cloud2trees_ans_c <- cloud2trees::cloud2trees(
  output_dir = tempdir()
  , input_las_dir = las_path
  , dtm_res_m = 0.25
  , ws = sequoia_ws_lo
  , estimate_tree_dbh = F
  , estimate_tree_type = F
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = T
  , cbh_estimate_missing_cbh = F
)

paste(
  "Default trees extracted:"
  , cloud2trees_ans$crowns_sf %>% nrow()
  , "|| Custom trees extracted:"
  , cloud2trees_ans_c$crowns_sf %>% nrow()
)


------------------------------------
  
# Convert LAS object to a dataframe
lidar_df <- as.data.frame(las_combined@data)  # Extract X, Y, Z coordinates

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

----------------------------------------
  # Model DBH: Random Forest
----------------------------------------
  
crowns_metrics

# read your CSV
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest2/FieldPointsCombo.csv")
crowncombo <- read_csv("C:/Users/User/Desktop/combocrowns.csv")

# check structure
glimpse(fieldpoints)


library(dplyr)
library(stringr)

# Clean SFM column
fieldpoints_clean <- fieldpoints %>%
  mutate(SFM = str_trim(SFM),       # remove spaces
         SFM = str_replace_all(SFM, "[\r\n]", ""))  # remove carriage returns / newlines


# Now do the join
matched <- crowncombo %>%
  inner_join(SFM, by = "treeID")


crowns$treeID <- as.numeric(crowns_metrics$treeID)
fieldpoints$SFM <- as.numeric(fieldpoints$SFM)
----------------------------------------
  # Model DBH: Random Forest
----------------------------------------


# Check for missing values
summary(crowns_metrics)

# Fit the model to predict DBH
set.seed(123)

rf_model <- randomForest(
  dbh_cm ~ tree_height_m + is_sequoia + crown_area_m2 + basal_area_m2 +
    radius_m + p95 + p75 + p50 + mean_z + max_z +
    sd_z + cv_z + point_density + slope_mean + aspect_circmean + ndvi_mean,
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

boxplot(RF_variables$dbh_cm)
boxplot(RF_variables$tree_height_m_1)


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
