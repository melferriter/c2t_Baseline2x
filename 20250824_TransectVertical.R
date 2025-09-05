library(lidR)
library(ggplot2)
library(dplyr)

# Load your LAS file
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap2x_clipped.las")


#read in the SFM and 3XCross Las files and clip to ROI
roi <- st_read("C:/Users/User/Desktop/LidarExtent.shp", crs = 32611)
roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")
plot(roi)
----------------
  
las <- clip_roi(las, roi_2d)

## classify ground and normalize without DTM (better results)
las <- classify_ground(las, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
las <- normalize_height(las, knnidw())

writeLAS(las, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap2x_transect.las")

# Define rectangle extent for transect (adjust xmin/xmax/ymin/ymax)
transect_extent <- extent(364334, 364450.7 , 4000829 , 4000946) # Example coords

# Clip the transect
transect <- clip_rectangle(las, 
                           xleft   = transect_extent@xmin,
                           ybottom = transect_extent@ymin,
                           xright  = transect_extent@xmax,
                           ytop    = transect_extent@ymax)

df <- transect@data %>%
  as.data.frame() %>%
  dplyr::select(X, Y, Z, Intensity, Classification)

transect_df <- df



library(ggplot2)

ggplot(df, aes(x = X, y = Z, color = Z)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Point Cloud Cross-section",
       x = "Northing", 
       y = "Height (m)")

library(ggExtra)

p <- ggplot(df, aes(x = X, y = Z, color = Z)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis_c() +
  theme_minimal()

ggMarginal(p, type = "histogram", margins = "y", bins = 50)



# View transect
plot(transect)

# Or make a vertical slice profile
transect_profile <- plot_profiles(transect, color = "Z")



library(patchwork)

p1 <- ggplot(transect_df, aes(x = X, y = Z, color = Z)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  theme_minimal() +
  labs(title = "Point Cloud Cross-section",
       x = "UTM Easting", 
       y = "Height (m)")

p2 <- ggplot(transect_df, aes(x = Z)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Point Cloud Cross-section",
       x = "UTM Easting", 
       y = "Height (m)")

# Combine: profile + density
p1 | p2

ggsave("flights_2x_over50.png", flights_2x, width = 10, height = 5, dpi = 300)

