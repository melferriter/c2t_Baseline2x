# Load LAS file
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap2x_clipped.las")
las_b <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xSpeed_20241021230048_clip.las")

flights <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2x_Side/2x_Side_traj_Clip.shp")


las <- classify_ground(las, csf())
las_b <- classify_ground(las_b, csf())

las_norm <- normalize_height(las, knnidw())
las_norm_b <- normalize_height(las_b, knnidw())


# Example: 0.5 m resolution canopy height raster
chm <- grid_canopy(las_norm, res = 0.5, p2r())
chm_b <- grid_canopy(las_norm_b, res = 0.5, p2r())

# Define a 3x3 moving window
w <- matrix(1, 3, 3)

# Fill small gaps (use mean of neighbors only where NA exists)
chm_filled <- terra::focal(
  chm, 
  w = w, 
  fun = mean, 
  na.policy = "only", 
  na.rm = TRUE
)

# Fill small gaps (use mean of neighbors only where NA exists)
chm_filled_b <- terra::focal(
  chm_b, 
  w = w, 
  fun = mean, 
  na.policy = "only", 
  na.rm = TRUE
)
plot(chm_filled_b)
plot(chm_filled, add = TRUE, col = viridis::viridis(100))


ggplot() +
  # Background CHM (grayscale)
  geom_raster(data = df_big, aes(x = x, y = y, fill = height_big)) +
  scale_fill_gradient(low = "black", high = "white", na.value = NA, guide = "none") +
  
  new_scale_fill() +

  # Overlay CHM (viridis)
  geom_raster(data = df_small, aes(x = x, y = y, fill = height_small)) +
  scale_fill_viridis(option = "viridis", na.value = NA, guide = "none") +

  # Overlay flight path (sf line)
  geom_sf(data = flights, color = "lightblue", linewidth = .8) +

  # Coordinate system (must use coord_sf for sf objects)
  coord_sf() +
  theme_void() +
  labs(title = "75% Sidelap Flight") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +

  # Scalebar and north arrow (NOW correctly added inside the ggplot chain)
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black","white")) 
  

#annotation_north_arrow(
  location = "tr", 
  which_north = "true",
  style = north_arrow_fancy_orienteering,
  height = unit(1, "cm"),   # smaller height
  width  = unit(1, "cm")    # smaller width
)

library(patchwork)

#p1 + p2 + p3 + p4 + 
  #plot_layout(ncol = 4) & 
  #theme(plot.title = element_text(size = 12))


#ggsave("Flightlines_Figure.png", width = 12, height = 6, dpi = 600)