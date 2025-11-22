library(lidR)
library(dplyr)
library(ggplot2)
library(patchwork)

# Path to the file you want
#LiDAR
#las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_side/sidelap3x_transect.las")
#SfM
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/sfmonly_transect.las")
#Combo
#las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/combo_transect.las")

las <- filter_poi(las, Z >= 0)
plot(las)

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



# Read LAS
las <- readLAS(r)
if (is.empty(las)) stop("❌ LAS file is empty")

# Clip to ROI if needed
transect <- clip_rectangle(
  las,
  xleft   = 364334,
  ybottom = 4000829,
  xright  = 364450.7,
  ytop    = 4000946
)

# Convert to data.frame
df <- as.data.frame(transect@data) %>%
  dplyr::select(X, Y, Z)

# Plot cross-section
p1 <- ggplot(df, aes(x = X, y = Z, color = Z)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis_c(option = "viridis", name = "Height (m)") +
  theme_minimal(base_size = 14) +
  labs(title = "Point Cloud Cross-section",
       x = "UTM Easting (m)", y = "Height (m)")


# Remove negative or ground-level points
df_filtered <- df %>%
  filter(Z > 0)   # keep only points above 2 m


# # Correct vertical density profile
# p2 <- ggplot(df_filtered, aes(x = Z)) +
#   geom_density(
#     aes(y = after_stat(density)),
#     fill = "skyblue", alpha = 0.5, color = "black", size = 0.8
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(
#     title = "Vertical Density Profile",
#     subtitle = "Z > 2 m)",
#     x = "Density",
#     y = "Height (m)"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", size = 16),
#     plot.subtitle = element_text(size = 12, color = "gray30"),
#     axis.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank()
#   )



# Right plot: comparable smooth vertical density profile
p2 <- ggplot(df_filtered, aes(y = Z)) +
  geom_density(
    aes(x = after_stat(density / max(density) * 100)),  # normalize to % of peak
    fill = "skyblue", alpha = 0.6, color = "skyblue"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Vertical Density Profile",
    x = "Relative Density (%)",
    y = "Height (m)"
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )# all flights scaled 0–100%





# Combine side-by-side
plot <- p1 | p2


ggsave(
  filename = "sfmonly_transect.png",  # file name
  plot = plot,                          # saves the most recent ggplot
  path = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/VerticalTransect",      # <-- change to your folder
  width = 10, height = 5,                       # adjust for paper layout
  dpi = 300                                    # high resolution for publication
)
 











