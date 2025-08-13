

baseline2x_o <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las")

# 2xBaseline

i = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las" 

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
  , estimate_tree_dbh = F
  , estimate_tree_type = F
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = F
  , cbh_estimate_missing_cbh = F
)

paste(
  "Default trees extracted:"
  , cloud2trees_ans$crowns_sf %>% nrow()
  , "|| Custom trees extracted:"
  , cloud2trees_ans_c$crowns_sf %>% nrow()
)

-------------------
  
# Convert LAS object to a dataframe
lidar_df <- as.data.frame(baseline2x_o@data)  # Extract X, Y, Z coordinates

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

filtered <- filter(density, point_density > 300 & tree_height_m >25)  # Set empirically


# Filter the crowns to only include those with treeIDs in your filtered list
filtered_crowns <- cloud2trees_ans_c$crowns_sf %>%
  filter(treeID %in% filtered$treeID)

density_filtered <- density %>%
  filter(treeID %in% filtered_crowns$treeID)

st_write(density_filtered,
         "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/clipped_filtered/2xBaseline_filtered_crowns_density.gpkg",
         layer = "density_filtered",
         driver = "GPKG",
         delete_dsn = TRUE)  # Overwrites existing file

st_write(filtered_crowns, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/clipped_filtered/2xBaseline_filtered_crowns.gpkg")

------------------

ggplot() +
  geom_sf(data = density, aes(fill = point_density), color = "black", alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "LiDAR Density") +
  theme_minimal() +
  ggtitle("LiDAR Point Density Within Tree Crowns")


ggplot() +
  geom_sf(data = filtered, aes(fill = point_density), color = "black", alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "LiDAR Density") +
  theme_minimal() +
  ggtitle("LiDAR Point Density Within Tree Crowns")


ggplot2::ggplot() + 
  ggplot2::geom_sf(data = filtered_crowns, mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")





