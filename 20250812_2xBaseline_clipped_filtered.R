# =========================
# LiDAR density for crowns (polygons) AND treetops (points)
# =========================

# ---- Packages ----
library(lidR)
library(cloud2trees)
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(units)

# ---- USER PARAMS (edit here) ----
las_path   <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las"
out_dir    <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/clipped_filtered"
out_gpkg   <- file.path(out_dir, "2xBaseline_outputs.gpkg")

# thresholds
thr_pd <- 300  # points / m^2
thr_ht <- 25   # meters

# treetop buffer options
use_area_matched_buffers <- TRUE  # TRUE = buffer area matches crown area; FALSE = fixed radius below
fixed_buffer_radius_m    <- 1.0   # used only if use_area_matched_buffers == FALSE

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 1) Read LAS ----
las <- readLAS(las_path)
stopifnot(!is.null(las))

# ---- 2) Run cloud2trees once (use your tuned ws if you have it) ----
# If you already computed c2t earlier, skip this block and set c2t <- your_object
c2t <- cloud2trees::cloud2trees(
  output_dir                = tempdir(),
  input_las_dir             = las_path,
  dtm_res_m                 = 0.5,
  estimate_tree_dbh         = FALSE,
  estimate_tree_type        = FALSE,
  estimate_tree_competition = TRUE,
  estimate_tree_cbh         = FALSE,
  cbh_estimate_missing_cbh  = FALSE
)

# Extract layers we’ll use
crowns   <- c2t$crowns_sf       # polygons (must have treeID, tree_height_m)
treetops <- c2t$treetops_sf     # points   (should have treeID)

stopifnot(inherits(crowns, "sf"), inherits(treetops, "sf"))

message("Crowns: ", nrow(crowns), " | Treetops: ", nrow(treetops))

# ---- 3) Build LiDAR point sf (XY only) ----
lidar_xy <- as.data.frame(las@data)[, c("X", "Y")]
lidar_sf <- st_as_sf(lidar_xy, coords = c("X","Y"), crs = st_crs(crowns))
# If unsure, ensure all are the same CRS:
# treetops <- st_transform(treetops, st_crs(crowns))
# lidar_sf <- st_transform(lidar_sf, st_crs(crowns))

# =========================
# A) Crown-based density (points inside polygons / polygon area)
# =========================

# point-in-polygon
hits_crown <- st_intersects(crowns, lidar_sf, sparse = TRUE)

# counts and density
crown_counts <- tibble(
  treeID      = crowns$treeID,
  point_count = lengths(hits_crown)
)

density_crowns <- crowns %>%
  left_join(crown_counts, by = "treeID") %>%
  mutate(
    point_count   = replace_na(point_count, 0L),
    area_m2       = as.numeric(st_area(geometry)),
    point_density = point_count / area_m2,                    # pts / m^2
    is_sequoia    = as.integer(tree_height_m > 50),           # placeholder rule; swap if you have a real flag
    sequoia_group = factor(if_else(is_sequoia == 1, "Giant sequoia", "Other species"),
                           levels = c("Other species", "Giant sequoia"))
  )

density_crowns_filtered <- density_crowns %>%
  filter(point_density > thr_pd, tree_height_m > thr_ht)

message("Filtered crowns: ", nrow(density_crowns_filtered),
        "  (pd > ", thr_pd, " & height > ", thr_ht, ")")

# =========================
# B) Treetop-based density (points within buffers / buffer area)
# =========================

# Decide buffer for each treetop:
# - If area-matched: radius = sqrt(crown_area/pi)
# - Else: fixed radius for all points
if (use_area_matched_buffers) {
  # join crown area to treetops by treeID
  crown_area_tbl <- density_crowns %>%
    st_drop_geometry() %>%
    transmute(treeID, crown_area_m2 = area_m2)

  treetops_buf <- treetops %>%
    left_join(crown_area_tbl, by = "treeID") %>%
    mutate(
      crown_area_m2 = ifelse(is.na(crown_area_m2), pi * fixed_buffer_radius_m^2, crown_area_m2),
      radius_m      = sqrt(crown_area_m2 / pi)
    ) %>%
    st_buffer(dist = radius_m)
} else {
  treetops_buf <- st_buffer(treetops, dist = fixed_buffer_radius_m)
}

# count LiDAR points within each treetop buffer
hits_top <- st_intersects(treetops_buf, lidar_sf, sparse = TRUE)

top_counts <- tibble(
  treeID      = treetops$treeID,
  point_count = lengths(hits_top)
)

# compute buffer area and density
density_treetops <- treetops_buf %>%
  st_sf() %>%                       # ensure sf after mutate/buffer
  mutate(buffer_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  transmute(treeID, buffer_area_m2) %>%
  left_join(top_counts, by = "treeID") %>%
  mutate(
    point_count   = replace_na(point_count, 0L),
    point_density = point_count / buffer_area_m2
  ) %>%
  left_join(st_drop_geometry(treetops), by = "treeID") %>%
  st_as_sf()

# bring in crown attributes for convenience (height, sequoia label)
density_treetops <- density_treetops %>%
  left_join(
    density_crowns %>%
      st_drop_geometry() %>%
      select(treeID, tree_height_m, is_sequoia, sequoia_group),
    by = "treeID"
  )

density_treetops_filtered <- density_treetops %>%
  filter(point_density > thr_pd, tree_height_m > thr_ht)

message("Filtered treetops: ", nrow(density_treetops_filtered),
        "  (pd > ", thr_pd, " & height > ", thr_ht, ")")

# ---- 4) Save outputs to a single GeoPackage ----
st_write(density_crowns,            out_gpkg, layer = "crowns_density",            driver = "GPKG", delete_layer = TRUE)
st_write(density_crowns_filtered,   out_gpkg, layer = "crowns_density_filtered",   driver = "GPKG", delete_layer = TRUE)
st_write(density_treetops,          out_gpkg, layer = "treetops_density",          driver = "GPKG", delete_layer = TRUE)
st_write(density_treetops_filtered, out_gpkg, layer = "treetops_density_filtered", driver = "GPKG", delete_layer = TRUE)

# =========================
# Quick figures
# =========================
theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold"))

# Crowns: all vs filtered
p_c_all <- ggplot(density_crowns) +
  geom_sf(aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt", breaks = pretty_breaks(5)) +
  labs(title = "Crown-based LiDAR point density — All crowns") +
  theme_pub
p_c_filt <- ggplot(density_crowns_filtered) +
  geom_sf(aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt") +
  labs(title = paste0("Crown density — Filtered (>", thr_pd, " pts/m² & >", thr_ht, " m)")) +
  theme_pub

# Treetops: buffers (plot buffers, not just points, so area context is clear)
p_t_all <- ggplot() +
  geom_sf(data = density_crowns, fill = NA, color = "grey90", size = 0.1) +
  geom_sf(data = st_set_geometry(density_treetops, st_geometry(treetops_buf)),
          aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt") +
  labs(title = "Treetop-based LiDAR point density (buffer method) — All trees") +
  theme_pub

p_t_filt <- ggplot() +
  geom_sf(data = density_crowns, fill = NA, color = "grey90", size = 0.1) +
  geom_sf(data = st_set_geometry(density_treetops_filtered, st_geometry(treetops_buf)[treetops$treeID %in% density_treetops_filtered$treeID]),
          aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt") +
  labs(title = paste0("Treetop density — Filtered (>", thr_pd, " pts/m² & >", thr_ht, " m)")) +
  theme_pub

print(p_c_all);  print(p_c_filt)
print(p_t_all);  print(p_t_filt)
# ggsave(file.path(out_dir, "crowns_density_all.png"), p_c_all, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "crowns_density_filtered.png"), p_c_filt, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "treetops_density_all.png"), p_t_all, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "treetops_density_filtered.png"), p_t_filt, width=6.5, height=5, dpi=300)
