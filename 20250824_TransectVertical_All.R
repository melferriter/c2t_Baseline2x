# Load your LAS file
side2xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap2x_clipped.las")
cross2xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/crosshatch2x_clipped.las")
speed2xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/speed2x_clipped.las")
base2xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las")

side3xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las")
cross3xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/crosshatch3x_clipped.las")
speed3xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/speed3x_clipped.las")
base3xlas <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline3x_clipped.las")

#read in the SFM and 3XCross Las files and clip to ROI
roi <- st_read("C:/Users/User/Desktop/LidarExtent.shp", crs = 32611)
roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")
plot(roi)
----------------
  
side2xlas <- clip_roi(side2xlas, roi_2d)
base2xlas <- clip_roi(base2xlas, roi_2d)
cross2xlas <- clip_roi(cross2xlas, roi_2d)
speed2xlas <- clip_roi(speed2xlas, roi_2d)

side3xlas <- clip_roi(side3xlas, roi_2d)
base3xlas <- clip_roi(base3xlas, roi_2d)
cross3xlas <- clip_roi(cross3xlas, roi_2d)
speed3xlas <- clip_roi(speed3xlas, roi_2d)

## classify ground and normalize without DTM (better results)
side2xlas_norm <- classify_ground(side2xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
side2xlas_norm  <- normalize_height(side2xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
speed2xlas_norm <- classify_ground(speed2xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
speed2xlas_norm  <- normalize_height(speed2xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
base2xlas_norm <- classify_ground(base2xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
base2xlas_norm  <- normalize_height(base2xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
cross2xlas_norm <- classify_ground(cross2xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
cross2xlas_norm  <- normalize_height(cross2xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
side3xlas_norm <- classify_ground(side3xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
side3xlas_norm  <- normalize_height(side3xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
speed3xlas_norm <- classify_ground(speed3xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
speed3xlas_norm  <- normalize_height(speed3xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
base3xlas_norm <- classify_ground(base3xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
base3xlas_norm  <- normalize_height(base3xlas_norm , knnidw())

## classify ground and normalize without DTM (better results)
cross3xlas_norm <- classify_ground(cross3xlas, csf(sloop_smooth = FALSE, 
                                class_threshold = 0.2, 
                                cloth_resolution = 0.2, 
                                rigidness = 2L, 
                                iterations = 500L, 
                                time_step = 0.65))
cross3xlas_norm  <- normalize_height(cross3xlas_norm , knnidw())

writeLAS(side2xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2x_side/sidelap2x_transect.las")
writeLAS(speed2xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2xspeed/speed2x_transect.las")
writeLAS(base2xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2x_baseline/baselap2x_transect.las")
writeLAS(cross2xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/2x_cross/crosslap2x_transect.las")

writeLAS(side3xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_side/sidelap3x_transect.las")
writeLAS(speed3xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_speed/speedlap3x_transect.las")
writeLAS(base3xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_baseline/baselap3x_transect.las")
writeLAS(cross3xlas_norm, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_cross/crosslap3x_transect.las")


# =========================
# Libraries
# =========================
library(lidR)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(sf)
library(dplyr)

# =========================
# 1. Setup
# =========================
# Folder with your LAS files
las_dir <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Transect"

# List of LAS files to process (adjust names if needed)
las_files <- c(
  "baselap2x_transect.las",
  "crosslap2x_transect.las",
  "sidelap2x_transect.las",
  "speed2x_transect.las",
  "baselap3x_transect.las",
  "crosslap3x_transect.las",
  "sidelap3x_transect.las",
  "speedlap3x_transect.las"
)


# Define rectangle extent for transect
transect_extent <- extent(364334, 364450.7 , 4000829 , 4000946)

# =========================
# 2. Function to process one LAS
# =========================
process_las <- function(file, roi, transect_extent, out_dir = las_dir) {
  # Load LAS
  las <- readLAS(file.path(out_dir, file))
  if (is.empty(las)) {
    warning(paste("LAS is empty:", file))
    return(NULL)
  }
  
  
  # Save processed LAS
  out_file <- gsub("_clipped.las", "_transect.las", file)
  writeLAS(las, file.path(out_dir, out_file))
  
  # Clip transect
  transect <- clip_rectangle(
    las,
    xleft   = transect_extent@xmin,
    ybottom = transect_extent@ymin,
    xright  = transect_extent@xmax,
    ytop    = transect_extent@ymax
  )
  
  # Convert to dataframe
  df <- as.data.frame(transect@data) %>%
    dplyr::select(X, Y, Z, Intensity, Classification)
  
  # Build plots
  p1 <- ggplot(df, aes(x = X, y = Z, color = Z)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_viridis_c(option = "viridis") +
    theme_minimal() +
    labs(title = paste("Point Cloud Cross-section —", file),
         x = "UTM Easting", y = "Height (m)")
  
  p2 <- ggplot(df, aes(x = Z)) +
    geom_density(fill = "skyblue", alpha = 0.6) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Vertical Density",
         x = "Density", y = "Height (m)")
  
  combined <- p1 | p2
  return(combined)
}

# =========================
# 3. Run on all LAS files
# =========================
plot <- process_las("crosslap3x_transect.las", transect_extent = transect_extent)

# =========================
# 4. Combine all plots
# =========================
library(patchwork)
all_plots <- wrap_plots(plots, ncol = 2)
all_plots

# Optionally save
ggsave("All_Flights_Transects.png", plot, width = 18, height = 6, dpi = 300)



