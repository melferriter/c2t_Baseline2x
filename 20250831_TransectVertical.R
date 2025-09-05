# =========================
# Libraries
# =========================
library(lidR)
library(ggplot2)
library(patchwork)
library(sf)
library(dplyr)

# =========================
# 1. Setup
# =========================
las_dir <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las"

# Input LAS files
las_files <- c(
  "sidelap2x_clipped.las",
  "baseline2x_clipped.las",
  "crosshatch2x_clipped.las",
  "speed2x_clipped.las",
  "sidelap3x_clipped.las",
  "baseline3x_clipped.las",
  "crosshatch3x_clipped.las",
  "speed3x_clipped.las",
  "SFM_3xSide_20241022194753_20250825.las"
)

# Output subfolders (match in same order)
out_dirs <- c(
  "2x_side", "2x_baseline", "2x_cross", "2x_speed",
  "3x_side", "3x_baseline", "3x_cross", "3x_speed",
  "Combo"
)

# Clip region
roi <- st_read("C:/Users/User/Desktop/LidarExtent.shp", crs = 32611) %>%
  st_zm(drop = TRUE, what = "ZM")

# =========================
# 2. Processing Function
# =========================
process_and_save <- function(file, out_subdir, roi, las_dir) {
  las <- readLAS(file.path(las_dir, file))
  if (is.empty(las)) {
    warning(paste("Skipping empty LAS:", file))
    return(NULL)
  }
  
  # Clip to ROI
  las <- clip_roi(las, roi)
  
  # Classify & normalize
  las <- classify_ground(las, csf(
    sloop_smooth    = FALSE, 
    class_threshold = 0.2, 
    cloth_resolution= 0.2, 
    rigidness       = 2L, 
    iterations      = 500L, 
    time_step       = 0.65
  ))
  las <- normalize_height(las, knnidw())
  
  # Save
  out_path <- file.path(las_dir, out_subdir,
                        gsub("_clipped", "_transect", file))
  writeLAS(las, out_path)
  
  return(las)
}

# =========================
# 3. Run Processing
# =========================
dir.create(file.path(las_dir, "2x_side"), showWarnings = FALSE)
dir.create(file.path(las_dir, "2x_baseline"), showWarnings = FALSE)
dir.create(file.path(las_dir, "2x_cross"), showWarnings = FALSE)
dir.create(file.path(las_dir, "2x_speed"), showWarnings = FALSE)
dir.create(file.path(las_dir, "3x_side"), showWarnings = FALSE)
dir.create(file.path(las_dir, "3x_baseline"), showWarnings = FALSE)
dir.create(file.path(las_dir, "3x_cross"), showWarnings = FALSE)
dir.create(file.path(las_dir, "3x_speed"), showWarnings = FALSE)
dir.create(file.path(las_dir, "Combo"), showWarnings = FALSE)

las_normals <- mapply(process_and_save,
                      file = las_files,
                      out_subdir = out_dirs,
                      MoreArgs = list(roi = roi, las_dir = las_dir),
                      SIMPLIFY = FALSE)

# =========================
# 4. Example Transect Plot
# =========================
transect_extent <- extent(364334, 364450.7, 4000829, 4000946)

process_las <- function(file, transect_extent, out_dir = las_dir) {
  las <- readLAS(file.path(out_dir, file))
  if (is.empty(las)) return(NULL)
  
  transect <- clip_rectangle(las,
                             xleft   = transect_extent@xmin,
                             ybottom = transect_extent@ymin,
                             xright  = transect_extent@xmax,
                             ytop    = transect_extent@ymax)
  
  df <- as.data.frame(transect@data) %>% dplyr::select(X, Y, Z)
  df$Z_bin <- cut(df$Z, breaks = seq(0, max(df$Z, na.rm = TRUE), by = 1))
  
  # Cross-section
  p1 <- ggplot(df, aes(X, Z, color = Z)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_viridis_c(name = "Height (m)") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Cross-section —", file),
         x = "UTM Easting (m)", y = "Height (m)")
  
  # Vertical profile
  height_dist <- df %>%
    group_by(Z_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(mid_height = as.numeric(sub("\\((.+),.*", "\\1", Z_bin)) + 0.5,
           rel_density = 100 * count / sum(count))
  
  p2 <- ggplot(height_dist, aes(x = rel_density, y = mid_height)) +
    geom_col(fill = "skyblue", alpha = 0.6) +
    geom_line(aes(y = mid_height, x = rel_density), color = "black") +
    labs(title = "Vertical Return Profile",
         x = "Relative Density (%)", y = "Height (m)") +
    theme_minimal(base_size = 14)
  
  p1 | p2
}

# Example run
plot <- process_las("3x_Side.las", transect_extent)
ggsave("sfm_transect.png", plot, width = 10, height = 5, dpi = 300)
