
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
roi <- st_read("C:/Users/User/Desktop/LidarExtent.shp", crs = 32611)
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


ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c_nosnag$crowns_sf, mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf(data = cloud2trees_ans_c_nosnag$treetops_sf, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

cloud2trees_ans_c_nosnag <- cloud2trees_ans_c$crowns

cloud2trees_ans_c_nosnag$treetops_sf <- cloud2trees_ans_c$treetops_sf[cloud2trees_ans_c$treetops_sf$tree_height_m > 10, ]

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

# height plot
plt_ht <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# diameter plot
plt_dbh <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = dbh_cm)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Purples", name = "tree DBH (cm)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# CBH plot
plt_cbh <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_cbh_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Greens", name = "tree CBH (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# combine with patchwork
plt_ht + plt_dbh + plt_cbh + patchwork::plot_layout(ncol = 2) &
  ggplot2::theme(
    legend.title = ggplot2::element_text(size = 8)
    , legend.text = ggplot2::element_text(size = 7)
  )