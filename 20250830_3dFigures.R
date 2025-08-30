# =========================
# 1. Packages
# =========================
library(lidR)
library(sf)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(terra)  # if you want to use a DTM

# =========================
# 2. Inputs
# =========================
las_path   <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las"
crowns_gpkg <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg"
dtm_path <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif"


# =========================
# 3. Load data
# =========================
las     <- readLAS(las_path)
crowns  <- st_read(crowns_gpkg)

stopifnot(!is.empty(las), inherits(crowns, "sf"))

# Make sure crowns has sequoia flag
if (!"is_sequoia" %in% names(crowns)) {
  if ("sequoia_group" %in% names(crowns)) {
    crowns <- crowns %>%
      mutate(is_sequoia = ifelse(as.character(sequoia_group) == "Giant sequoia", 1L, 0L))
  } else {
    crowns <- crowns %>% mutate(is_sequoia = as.integer(tree_height_m > 50))
  }
}

# =========================
# 4. Pick example crowns
# =========================
ex_sequoia_id <- crowns %>%
  filter(is_sequoia == 1) %>%
  arrange(desc(tree_height_m)) %>%
  slice(1) %>% pull(treeID)

ex_other_id <- crowns %>%
  filter(is_sequoia == 0) %>%
  arrange(desc(tree_height_m)) %>%
  slice(1) %>% pull(treeID)

# =========================
# 5. Helper functions
# =========================
clip_crown_las <- function(las, crowns_sf, tree_id, dtm_path = NULL) {
  ply <- crowns_sf %>% filter(treeID == tree_id)
  if (nrow(ply) != 1) {
    stop("TreeID not found or not unique: ", tree_id)
  }
  
  # ✅ Use clip_roi instead of lasclip
  clp <- clip_roi(las, ply)   
  if (is.empty(clp)) return(NULL)
  
  # Optional: normalize to DTM
  if (!is.null(dtm_path)) {
    dtm <- rast(dtm_path)
    clp <- normalize_height(clp, dtm)
    clp@data$Height <- clp@data$Z
  } else {
    clp@data$Height <- clp@data$Z
  }
  clp
}

las_to_plotly <- function(las_obj, title = "3D crown", zmax = 70) {
  pts <- as.data.frame(las_obj@data)[, c("X","Y","Z","Height")]
  
  # Shift X and Y to local crown coords
  pts$X <- pts$X - min(pts$X, na.rm = TRUE)
  pts$Y <- pts$Y - min(pts$Y, na.rm = TRUE)
  
  plot_ly(
    pts, x = ~X, y = ~Y, z = ~Z,
    type = "scatter3d", mode = "markers",
    marker = list(
      size = 2, opacity = 0.75,
      color = ~Height, colorscale = "Viridis",
      cmin = 0, cmax = zmax,   # same scale across trees
      colorbar = list(title = "Height (m)")
    )
  ) %>%
    layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "", showticklabels = FALSE, visible = FALSE),
        yaxis = list(title = "", showticklabels = FALSE, visible = FALSE),
        zaxis = list(title = "", showticklabels = FALSE, visible = FALSE)
      )
    )
}


# =========================
# 6. Clip + Plot
# =========================
# Get all Sequoia IDs
sequoia_ids <- crowns %>%
  filter(is_sequoia == 1) %>%
  pull(treeID)
view(sequoia_ids)
# View them
head(sequoia_ids, 20)   # first 20 IDs
length(sequoia_ids)  

nonsequoia_ids <- crowns %>%
  filter(is_sequoia == 0) %>%
  pull(treeID)
view(sequoia_ids)

#"140_364363.1_4000854.9", "137_364442.4_4000856.1", "32_364383.9_4000909.6", "137_364442.4_4000856.1", 154_364379.9_4000844.4,
#140_364363.1_4000854.9, 94_364443.4_4000875.4, 82_364342.4_4000880.4, 39_364383.1_4000900.9
ex_sequoia_id <- "140_364363.1_4000854.9"
las_seq  <- clip_crown_las(las, crowns, ex_sequoia_id, dtm_path)
p_seq   <- las_to_plotly(las_seq,  paste0("Giant sequoia (", ex_sequoia_id, ")"))


ex_nonsequoia_id <- "5_364352.9_4000937.1"
las_other <- clip_crown_las(las, crowns, ex_nonsequoia_id, dtm_path)
las_to_plotly(las_other,paste0("Non-sequoia (", ex_other_id, ")"))



p_other <- las_to_plotly(las_other,paste0("Non-sequoia (", ex_other_id, ")"))


---
#"94_364443.4_4000875.4"
ex_sequoia_id <- "94_364443.4_4000875.4"
las_seq  <- clip_crown_las(las, crowns, ex_sequoia_id, dtm_path)
p_seq <- las_to_plotly(las_seq,  paste0("Giant sequoia (", ex_sequoia_id, ")"))
library(plotly)
library(webshot2)
library(reticulate)
reticulate::py_install("kaleido")

plotly::save_image(p_seq,
           file = "sequoia.png",
           width  = 4200,
           height = 3000) 
htmlwidgets::saveWidget(p_seq, "sequoia_fixed_angle.html", selfcontained = TRUE)
getwd()

# =========================
# 7. Save outputs
# =========================
saveWidget(p_seq,   "sequoia_crown_3D.html", selfcontained = TRUE)
saveWidget(p_other, "nonsequoia_crown_3D.html", selfcontained = TRUE)

writeLAS(las_seq,   "sequoia_crown.las")
writeLAS(las_other, "nonsequoia_crown.las")


----------------
  library(lidR)
library(sf)
library(dplyr)

# crown polygons must have treeID
crowns <- crowns %>% mutate(treeID = as.character(treeID))

# Compute point counts per crown
crown_stats <- crowns %>%
  st_drop_geometry() %>%
  select(treeID, is_sequoia, tree_height_m) %>%
  mutate(point_count = NA_integer_, area_m2 = NA_real_, density = NA_real_)

for (i in seq_len(nrow(crowns))) {
  this_id <- crowns$treeID[i]
  clp <- clip_roi(las, crowns[i, ])   # clip LAS to this crown polygon
  if (!is.empty(clp)) {
    pts <- npoints(clp)                        # number of points in crown
    area <- st_area(crowns[i, ]) %>% as.numeric()
    crown_stats$point_count[i] <- pts
    crown_stats$area_m2[i] <- area
    crown_stats$density[i] <- pts / area       # pts per m²
  }
}

sequoia_density <- crown_stats %>%
  filter(is_sequoia == 1) %>%
  arrange(desc(density))

head(sequoia_density, 30)   # top 10 densest Sequoias


# all Sequoias sorted by trunk density
all_trunks <- top_trunks %>%
  select(treeID, tree_height_m, trunk_density) %>%
  arrange(desc(trunk_density))
arrange(desc(sequoia_density))
print(sequoia_density)


----------------
# helper to compute trunk density for one crown
compute_trunk_density <- function(las, crown_poly, zmax_trunk = 10) {
  clp <- clip_roi(las, crown_poly)
  if (is.empty(clp)) return(NA_real_)
  
  # filter to trunk zone (below zmax_trunk)
  trunk_pts <- filter_poi(clp, Z <= zmax_trunk)
  
  if (is.empty(trunk_pts)) return(NA_real_)
  
  # compute density (points per m² of crown footprint)
  pts <- npoints(trunk_pts)
  area <- st_area(crown_poly) %>% as.numeric()
  density <- pts / area
  return(density)
}

# ensure IDs are characters
crowns <- crowns %>% mutate(treeID = as.character(treeID))

# compute trunk densities
crown_stats <- crowns %>%
  filter(is_sequoia == 1) %>%
  rowwise() %>%
  mutate(trunk_density = compute_trunk_density(las, cur_data())) %>%
  ungroup()

# sort by densest trunks
top_trunks <- crown_stats %>%
  arrange(desc(trunk_density))

head(top_trunks, 10)   # densest 10 sequoias

top_seq_id <- top_trunks$treeID[1]
las_top_seq <- clip_roi(las, crowns %>% filter(treeID == top_seq_id))

p_trunk <- las_to_plotly(las_top_seq,
                         paste0("High trunk density Sequoia (", top_seq_id, ")"),
                         zmax = 70)
p_trunk



