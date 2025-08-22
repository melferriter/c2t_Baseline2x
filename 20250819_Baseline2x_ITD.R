
# ---- Packages ----
library(lidR)
library(cloud2trees)
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
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

# ---- USER PARAMS ----

las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las")
las_path <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las"


#read in the SFM and 3XCross Las files and clip to ROI
#roi <- st_read("C:/Users/User/Desktop/clip.shp", crs = 32611)
#roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")

#roi2 <- st_read("C:/Users/User/Desktop/clip2.shp", crs = 32611)
#roi_2d2 <- st_zm(roi, drop = TRUE, what = "ZM")
  
#las_clipped <- clip_roi(las, roi_2d)
#writeLAS(las_clipped, "C:/Users/User/Desktop/clip.las")
#las_path <- "C:/Users/User/Desktop/clip.las"


#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = las_path)

cloud2trees_ans$chm_rast %>% terra::plot()

## Default ITD window size functions====
itd_tuning_ans <- itd_tuning(input_las_dir = las_path)
itd_tuning_ans$plot_samples
best_ws <- itd_tuning_ans$ws_fn_list$lin_fn

# # a constant window size has to be defined as:
#  ## rep(constant, times = length(x))
#  ## x*0 + constant
# my_constant <- function(x){ 
#    return( rep(3, times = length(x)) ) ## will always return 3
#   } 
# # a custom linear function
#  my_linear <- function(x) {(x * 0.1) + 3}
#  my_curved <- function(x){
#   3 + (0.1 * x) + (0.002 * (x^2))  # quadratic term accelerates window growth at higher heights
#  }
#  
# # More aggressive scaling for tall trees
# my_sequoia <- function(x){
#   ifelse(x > 50, (x * 0.2) + 5,   # for giant sequoias: larger window growth
#                 (x * 0.1) + 3)   # for other conifers: standard growth
# }
# 
# # Custom nonlinear window size for mixed conifer + sequoia
# my_sequoia_ws <- function(x){
#   # x = tree height (m)
#   ws <- ifelse(x > 45, (x * 0.25) + 6,   # for giant sequoias
#          ifelse(x > 25, (x * 0.15) + 4,   # for mid-size conifers
#                         (x * 0.1)  + 3))  # for small trees
#   return(ws)
# }
# 
# my_sequoia_ws <- function(x){
#   (0.002 * (x^2)) + 3   # quadratic growth with height
# }
# 
# # tuned for sequoia
# my_sequoia_ws2 <- function(x) {
#   (x * 0.2) + 5   # larger window for tall trees
# }
# 
# my_sequoia_ws <- function(x){
#   (x * 0.15) + 5   # bigger windows for tall crowns
# }
# 
# y_sequoia_ws <- function(x){
#   2 + (x^0.8)/5    # curved relationship
# }
# 
# # make window size scale more strongly with height
# y_sequoia_ws2 <- function(x){
#   (x * 0.25) + 2   # steeper slope, bigger windows for tall crowns
# }
# 
# y_sequoia_ws3 <- function(x) {
#   ifelse(x > 50, (x * 0.01) + 3,   # smaller window for very tall trees
#          (x * 0.07) + 3)          # normal growth for smaller conifers
# }
# 
# 
# # let's put these in a list to test with the best default function we saved from above
#   my_fn_list <- list(
#     y_sequoia_ws3 = y_sequoia_ws3
#     , my_linear = my_linear
#     , y_sequoia_ws = y_sequoia_ws
#   )
#   
#my_fn_list <- list(
#  my_linear = my_linear,
#  sequoia_ws_lo = function(x) {(x * 0.024) + 2},
#  sequoia_ws_low2 = function(x) {(x * 0.04) + 3},
#  sequoia_ws_low3 = function(x) {(x * 0.02) + 2}
#)

my_fn_list <- list(
  ws_hi  = function(x) {(x * 0.08) + 4},   # steeper slope
  ws_med = function(x) {(x * 0.06) + 3},   # moderate
  ws_low = function(x) {(x * 0.024) + 2}    # conservative
)


# run it with custom functions
itd_tuning_ans2 <- itd_tuning(
 input_las_dir = las_path
 , ws_fn_list = my_fn_list
 , n_samples = 3
) 

# look at the tuning plot
itd_tuning_ans2$plot_samples

#sequoia_ws_lo = function(x) {(x * 0.04) + 2}
sequoia_ws_lo = function(x) {(x * 0.023) + 2}

cloud2trees_ans_c <- cloud2trees::cloud2trees(
  output_dir = tempdir()
  , input_las_dir = las_path
  , dtm_res_m = 0.25
  , ws = sequoia_ws_lo
  , estimate_tree_dbh = F
  , estimate_tree_type = F
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = F
  , cbh_estimate_missing_cbh = F
)

saveRDS(cloud2trees_ans_c, "Baseline2x_cloud2trees_ans_c.rds")

paste(
  "Default trees extracted:"
  , cloud2trees_ans$crowns_sf %>% nrow()
  , "|| Custom trees extracted:"
  , cloud2trees_ans_c$crowns_sf %>% nrow()
)

library(patchwork)
# height plot
# plt_ht <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "Tree Height (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# diameter plot
# plt_dbh <-
#   cloud2trees_ans_c$crowns_sf %>% 
#   ggplot2::ggplot(mapping = ggplot2::aes(fill = dbh_cm)) + 
#   ggplot2::geom_sf() + 
#   ggplot2::scale_fill_distiller(palette = "Purples", name = "tree DBH (cm)", direction = 1) +
#   ggplot2::theme_void() +
#   ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# # CBH plot
# plt_cbh <-
#   cloud2trees_ans_c$crowns_sf %>% 
#   ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_cbh_m)) + 
#   ggplot2::geom_sf() + 
#   ggplot2::scale_fill_distiller(palette = "Greens", name = "tree CBH (m)", direction = 1) +
#   ggplot2::theme_void() +
#   ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# # combine with patchwork
# plt_ht + plt_dbh + plt_cbh + patchwork::plot_layout(ncol = 2) &
#   ggplot2::theme(
#     legend.title = ggplot2::element_text(size = 8)
#     , legend.text = ggplot2::element_text(size = 7)
#  )

c2t <- cloud2trees_ans_c

# Layers
crowns   <- c2t$crowns_sf       # polygons; must include treeID, tree_height_m
treetops <- c2t$treetops_sf     # points; must include treeID; has crown_area_m2 per your note
stopifnot(inherits(crowns, "sf"), inherits(treetops, "sf"))

message("Crowns: ", nrow(crowns), " | Treetops: ", nrow(treetops))

# ---- 3) LiDAR points as sf (XY only) ----
lidar_xy <- as.data.frame(las@data)[, c("X","Y")]
lidar_sf <- st_as_sf(lidar_xy, coords = c("X","Y"), crs = st_crs(crowns))

thr_pd <- 300   # points / m^2
thr_ht <- 25    # meters (tree height)
fallback_radius_m <- 1.0  # used only if treetops$crown_area_m2 is NA/invalid

# =========================
# A) Crown-based density (points in polygon / polygon area)
# =========================
hits_crown <- st_intersects(crowns, lidar_sf, sparse = TRUE)
crown_counts <- tibble(treeID = crowns$treeID, point_count = lengths(hits_crown))

density_crowns <- crowns %>%
  left_join(crown_counts, by = "treeID") %>%
  mutate(
    point_count   = replace_na(point_count, 0L),
    area_m2       = as.numeric(st_area(geometry)),
    point_density = point_count / area_m2,                 # pts / m^2
    is_sequoia = case_when(
      tree_height_m >= 55 ~ 1,         # tall trees are almost always sequoia
      tree_height_m >= 40 & crown_area_m2 > 100 ~ 1,  # medium trees w/ large crowns are likely sequoia
      TRUE ~ 0
    ),
    sequoia_group = factor(if_else(is_sequoia == 1, "Giant sequoia", "Other species"),
                           levels = c("Other species", "Giant sequoia"))
  )


#thr_ht  <- 55   # height in meters
#thr_pd  <- quantile(density_crowns$point_density, 0.75)  # top quartile of density


density_crowns_filtered <- density_crowns %>%
  filter(point_density > thr_pd, tree_height_m > thr_ht)

message("Filtered crowns: ", nrow(density_crowns_filtered),
        "  (pd > ", thr_pd, " & height > ", thr_ht, ")")
plot(density_crowns_filtered)

density_crowns_filtered %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "Tree Height (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

saveRDS(density_crowns_filtered, "density_crowns_filtered.rds")

# ============================================================
# Parameters
# ============================================================
thr_pd <- 300   # points / m^2
thr_ht <- 25       # minimum tree height (m)
fallback_radius_m <- 1.0   # fallback buffer radius if crown area invalid


# ------------------------------
# Prepare treetops with buffer
# ------------------------------
treetops <- treetops %>%
  mutate(
    crown_area_m2 = as.numeric(crown_area_m2),
    radius_m = sqrt(crown_area_m2 / pi),
    radius_m = ifelse(is.na(radius_m) | !is.finite(radius_m) | radius_m <= 0,
                      fallback_radius_m, radius_m)
  )

treetops_buf <- st_buffer(treetops, dist = treetops$radius_m)

# ------------------------------
# Count lidar points per treetop buffer
# ------------------------------
hits_top <- st_intersects(treetops_buf, lidar_sf, sparse = TRUE)
top_counts <- tibble(treeID = treetops$treeID,
                     point_count = lengths(hits_top))

# ------------------------------
# Build density table
# ------------------------------
density_treetops <- treetops_buf %>%
  mutate(buffer_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  left_join(top_counts, by = "treeID") %>%
  mutate(
    point_count   = replace_na(point_count, 0L),
    point_density = point_count / buffer_area_m2
  ) %>%
  # bring in height (ONLY ONCE, from crowns)
  left_join(
    crowns %>% 
      st_drop_geometry() %>% 
      dplyr::select(treeID, tree_height_m),
    by = "treeID"
  ) %>%
  # sequoia classification (⚠️ NO height here!)
  left_join(
    density_crowns %>%
      st_drop_geometry() %>%
      dplyr::select(treeID, is_sequoia, sequoia_group),
    by = "treeID"
  ) %>%
  # restore geometry
  mutate(geom = st_geometry(treetops_buf)) %>%
  st_as_sf(sf_column_name = "geom", crs = st_crs(treetops_buf))

# ------------------------------
# Filter with SAME thresholds as crowns
# ------------------------------
density_treetops_filtered <- density_treetops %>%
  filter(point_density > thr_pd, tree_height_m.x > thr_ht)

message("Filtered treetops: ", nrow(density_treetops_filtered),
        " (pd > ", thr_pd, " & height > ", thr_ht, ")")
plot(density_treetops_filtered)

density_treetops_points_filtered %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m.x)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "Tree Height (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

# ------------------------------
# Create point version
# ------------------------------
density_treetops_points <- st_set_geometry(
  density_treetops,
  st_geometry(treetops)[match(density_treetops$treeID, treetops$treeID)]
)

density_treetops_points_filtered <- st_set_geometry(
  density_treetops_filtered,
  st_geometry(treetops)[match(density_treetops_filtered$treeID, treetops$treeID)]
)

saveRDS(density_treetops_points_filtered, "density_treetops_points_filtered.rds")


ggplot2::ggplot() + 
  ggplot2::geom_sf(data = density_crowns_filtered, mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf(data = density_treetops_points_filtered, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "Tree Height (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

message("Filtered treetops (points): ", nrow(density_treetops_points_filtered))

--------------------------

# Both in one .gpkg file, each layer overwritten if it exists
st_write(density_crowns_filtered,
         "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2x_Baseline/2x_Baseline_crowns.gpkg",
         delete_layer = TRUE)

#st_write(density_treetops_points_filtered,
#         "trees.gpkg",
#         layer = "treetops",
#         delete_layer = TRUE)

