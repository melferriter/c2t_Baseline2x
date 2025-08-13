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

# ---- USER PARAMS ----
las_path <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las"
out_dir  <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/clipped_filtered"
out_gpkg <- file.path(out_dir, "2xBaseline_outputs.gpkg")

thr_pd <- 300   # points / m^2
thr_ht <- 25    # meters (tree height)
fallback_radius_m <- 1.0  # used only if treetops$crown_area_m2 is NA/invalid

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 1) Read LAS ----
las <- readLAS(las_path)
stopifnot(!is.null(las))

# ---- 2) Run cloud2trees once (or set c2t <- your precomputed object) ----
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

# Layers
crowns   <- c2t$crowns_sf       # polygons; must include treeID, tree_height_m
treetops <- c2t$treetops_sf     # points; must include treeID; has crown_area_m2 per your note
stopifnot(inherits(crowns, "sf"), inherits(treetops, "sf"))

message("Crowns: ", nrow(crowns), " | Treetops: ", nrow(treetops))

# ---- 3) LiDAR points as sf (XY only) ----
lidar_xy <- as.data.frame(las@data)[, c("X","Y")]
lidar_sf <- st_as_sf(lidar_xy, coords = c("X","Y"), crs = st_crs(crowns))
# If needed: treetops <- st_transform(treetops, st_crs(crowns)); lidar_sf <- st_transform(lidar_sf, st_crs(crowns))

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
    is_sequoia    = as.integer(tree_height_m > 50),        # placeholder; replace if you have a real flag
    sequoia_group = factor(if_else(is_sequoia == 1, "Giant sequoia", "Other species"),
                           levels = c("Other species", "Giant sequoia"))
  )

density_crowns_filtered <- density_crowns %>%
  filter(point_density > thr_pd, tree_height_m > thr_ht)

message("Filtered crowns: ", nrow(density_crowns_filtered),
        "  (pd > ", thr_pd, " & height > ", thr_ht, ")")

# =========================
# B) Treetop-based density (points in area-matched buffer / buffer area)
# =========================
stopifnot("crown_area_m2" %in% names(treetops))  # you said this exists

treetops <- treetops %>%
  mutate(
    crown_area_m2 = as.numeric(crown_area_m2),
    radius_m = sqrt(crown_area_m2 / pi),
    radius_m = ifelse(is.na(radius_m) | !is.finite(radius_m) | radius_m <= 0,
                      fallback_radius_m, radius_m)
  )

treetops_buf <- st_buffer(treetops, dist = treetops$radius_m)

hits_top <- st_intersects(treetops_buf, lidar_sf, sparse = TRUE)
top_counts <- tibble(treeID = treetops$treeID, point_count = lengths(hits_top))

density_treetops <- treetops_buf %>%
  mutate(buffer_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  left_join(top_counts, by = "treeID") %>%
  mutate(
    point_count   = replace_na(point_count, 0L),
    point_density = point_count / buffer_area_m2
  ) %>%
  left_join(
    crowns %>% st_drop_geometry() %>% select(treeID, tree_height_m),
    by = "treeID"
  ) %>%
  left_join(
    density_crowns %>% st_drop_geometry() %>% select(treeID, is_sequoia, sequoia_group),
    by = "treeID"
  ) %>%
  left_join(
    st_drop_geometry(treetops) %>% select(treeID),  # keep treetop attributes if you have more
    by = "treeID"
  ) %>%
  left_join(                                                  # restore geometry = buffers
    tibble(treeID = treetops$treeID,
           geom = st_geometry(treetops_buf)),
    by = "treeID"
  ) %>%
  st_as_sf(sf_column_name = "geom", crs = st_crs(treetops_buf))

density_treetops_filtered <- density_treetops %>%
  filter(point_density > thr_pd, tree_height_m.x > thr_ht)

message("Filtered treetops: ", nrow(density_treetops_filtered),
        "  (pd > ", thr_pd, " & height > ", thr_ht, ")")

# ---- 4) Save outputs ----
st_write(density_crowns,            out_gpkg, layer = "crowns_density",            driver = "GPKG", delete_layer = TRUE)
st_write(density_crowns_filtered,   out_gpkg, layer = "crowns_density_filtered",   driver = "GPKG", delete_layer = TRUE)
st_write(density_treetops,          out_gpkg, layer = "treetops_density",          driver = "GPKG", delete_layer = TRUE)
st_write(density_treetops_filtered, out_gpkg, layer = "treetops_density_filtered", driver = "GPKG", delete_layer = TRUE)

# =========================
# Quick figures (comment out if not needed)
# =========================
theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold"))

# Crowns — all vs filtered
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

# Treetops (plot buffers so the density area is visible)
p_t_all <- ggplot() +
  geom_sf(data = density_crowns, fill = NA, color = "grey92", size = 0.1) +
  geom_sf(data = density_treetops, aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt") +
  labs(title = "Treetop-based LiDAR point density (area-matched buffers)") +
  theme_pub

p_t_filt <- ggplot() +
  geom_sf(data = density_crowns, fill = NA, color = "grey92", size = 0.1) +
  geom_sf(data = density_treetops_filtered, aes(fill = point_density), color = NA) +
  scale_fill_viridis_c("Points / m²", trans = "sqrt") +
  labs(title = paste0("Treetop density — Filtered (>", thr_pd, " pts/m² & >", thr_ht, " m)")) +
  theme_pub

print(p_c_all); print(p_c_filt); print(p_t_all); print(p_t_filt)
# ggsave(file.path(out_dir, "crowns_density_all.png"), p_c_all, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "crowns_density_filtered.png"), p_c_filt, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "treetops_density_all.png"), p_t_all, width=6.5, height=5, dpi=300)
# ggsave(file.path(out_dir, "treetops_density_filtered.png"), p_t_filt, width=6.5, height=5, dpi=300)



# reattach POINT geometry to the filtered treetops (they currently have buffer geometry)
tops_filtered_pts <- st_set_geometry(
  density_treetops_filtered,
  st_geometry(treetops)[match(density_treetops_filtered$treeID, treetops$treeID)]
)

ggplot() +
  geom_sf(data = density_crowns_filtered, fill = NA, color = "grey85", size = 0.2) +
  geom_sf(data = tops_filtered_pts, aes(color = tree_height_m.x),
          size = 1.8, shape = 16, alpha = 0.95) +
  scale_color_distiller(palette = "Oranges", direction = 1, name = "tree ht. (m)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())


#not filtered
# crowns as polygons
ggplot() +
  geom_sf(data = density_crowns, fill = NA, color = "grey70", size = 0.2) +
  # treetops as points
  geom_sf(data = treetops, aes(color = tree_height_m),    # or aes(color = point_density) if you have it
          size = 1.6, shape = 16, alpha = 0.9) +          # small solid points
  scale_color_distiller(palette = "Oranges", direction = 1, name = "tree ht. (m)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())


------------------
# === A) All treetops as points ===
# (this is just the c2t$treetops_sf, but we can add the density info)
treetops_all_pts <- treetops %>%
  left_join(
    density_treetops %>% 
      st_drop_geometry() %>%
      select(treeID, point_density, tree_height_m.x, is_sequoia, sequoia_group),
    by = "treeID"
  )

# === B) Filtered treetops as points ===
# Replace buffer geometry in density_treetops_filtered with original point geometry
treetops_filtered_pts <- st_set_geometry(
  density_treetops_filtered %>%
    st_drop_geometry(),
  st_geometry(treetops)[match(density_treetops_filtered$treeID, treetops$treeID)]
)

# === C) Write to GeoPackage for ArcGIS Pro ===
st_write(treetops_all_pts, out_gpkg, layer = "treetops_all_points", driver = "GPKG", delete_layer = TRUE)
st_write(treetops_filtered_pts, out_gpkg, layer = "treetops_filtered_points", driver = "GPKG", delete_layer = TRUE)

