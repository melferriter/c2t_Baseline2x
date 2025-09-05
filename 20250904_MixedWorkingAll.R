# ============================
# 0. Libraries
# ============================
library(lidR)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(caret)
library(tibble)
library(lmerTest)
library(purrr)
library(yardstick)
library(kableExtra)
library(ggplot2)

rm(list = ls(globalenv()))

# ============================
# 1. File definitions
# ============================
files <- list(
  Fusion = list(
    crowns = st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg"),
    field  = read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv"),
    las = readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las")
  ),
  LiDAR = list(
    crowns = st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg"),
    field  = read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_side3x.csv"),
    las = readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las")
  ),
  SfM = list(
    crowns = st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg"),
    field  = read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv"),
    las = readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las")
  )
)

dtm <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

-----------------



# ============================
# 2. Cross-validation setup
# ============================
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"
)

# ============================
# 3. Helper functions
# ============================
all_results <- list()

for (flight in names(files)) {
  message("Processing: ", flight)

  # ---- Load inputs ----
  crowns <- files[[flight]]$crowns
  fieldpoints <- files[[flight]]$field
  las <- files[[flight]]$las

  # ---- LAS ground + normalize ----
  las_ground <- classify_ground(las, csf())
  las_norm   <- normalize_height(las_ground, knnidw())

  # ---- Ancillary rasters ----
  slope <- terrain(raster(dtm), opt = "slope", unit = "degrees")

  # ---- Crown-level metrics ----
  metrics_list <- map_dfr(1:nrow(crowns), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)

    if (is.null(las_clip) || npoints(las_clip) == 0) {
      return(tibble(treeID = crown_poly$treeID))
    }

    Z <- las_clip@data$Z
    area_m2 <- as.numeric(st_area(crown_poly))

    tibble(
      treeID = crown_poly$treeID,
      point_density = npoints(las_clip) / area_m2,
      slope_mean    = exact_extract(slope, crown_poly, "mean"),
      max_z  = max(Z, na.rm = TRUE),
      mean_z = mean(Z, na.rm = TRUE),
      p25    = quantile(Z, 0.25, na.rm = TRUE),
      p50    = quantile(Z, 0.50, na.rm = TRUE),
      p75    = quantile(Z, 0.75, na.rm = TRUE),
      p95    = quantile(Z, 0.95, na.rm = TRUE),
      sd_z   = sd(Z, na.rm = TRUE),
      cv_z   = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE),
      ndvi_mean = exact_extract(ndvi, crown_poly, "mean"),
      crown_area_m2 = area_m2
    )
  })

  crowns_metrics <- crowns %>%
    select(treeID, geom) %>%
    left_join(metrics_list, by = "treeID")

  # ---- Field data cleanup ----
  fieldpoints_clean <- fieldpoints %>%
    mutate(
      treeID = !!sym(names(fieldpoints)[1]),   # auto-detect ID col (SfMcombo/Side3x/SFMonly)
      Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE),
                       "Sequoia", "Other"),
      is_sequoia = ifelse(Species == "Sequoia", 1, 0)
    ) %>%
    select(treeID, dbh_cm, treecbh, treeheight, Species, is_sequoia)

  # ---- Join crowns + field ----
  crowns_final <- crowns_metrics %>%
    inner_join(fieldpoints_clean, by = "treeID") %>%
    select(
      treeID, dbh_cm, treecbh, treeheight, Species, is_sequoia,
      crown_area_m2, point_density, slope_mean, ndvi_mean,
      max_z, mean_z, p25, p50, p75, p95, sd_z, cv_z, geom
    ) %>%
    filter(Species != "intact snag", treecbh > 0)

  # ---- Prepare modeling dataset ----
  crown_data <- crowns_final %>%
    filter(!is.na(dbh_cm), dbh_cm > 0, dbh_cm < 350) %>%
    mutate(
      log_dbh        = log(dbh_cm),
      log_height     = log(treeheight),
      log_crown_area = log1p(crown_area_m2),
      crown_depth    = pmax(treeheight - treecbh, 0)
    )
}

  # ---- Run all models (from cleaned script) ----
  res_rf_base       <- summarize_rf(train(...), paste(flight, "RF baseline"))
  res_mix           <- summarize_cv_lmer(crown_data, paste(flight, "Mixed-effects regression"))
  res_rf_seq        <- summarize_rf(train(...), paste(flight, "Sequoia RF"))
  res_rf_nonseq     <- summarize_rf(train(...), paste(flight, "Non-sequoia RF"))
  res_loglog_seq    <- summarize_cv_lm(filter(crown_data, is_sequoia == 1),
                                       log(dbh_cm) ~ log(treeheight) + log(crown_area_m2) + log(crown_depth),
                                       paste(flight, "Sequoia log–log"))
  res_loglog_nonseq <- summarize_cv_lm(filter(crown_data, is_sequoia == 0),
                                       log(dbh_cm) ~ log(treeheight) + log(crown_area_m2),
                                       paste(flight, "Non-sequoia log–log"))

  # ---- Collect ----
  all_results[[flight]] <- bind_rows(
    res_rf_base, res_mix, res_rf_seq,
    res_rf_nonseq, res_loglog_seq, res_loglog_nonseq
  )
}

# ============================
# 5. Final results table
# ============================
final_results <- bind_rows(all_results) %>%
  select(Source, everything())

final_results %>%
  kbl(caption = "📊 Model Performance (5-fold CV, all flights)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))
