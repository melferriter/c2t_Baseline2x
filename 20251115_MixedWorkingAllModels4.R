# ============================================================
# 0) Setup
# ============================================================
rm(list = ls())

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
library(readr)
library(exactextractr)
library(glue)
library(boot)
library(tidyr)
library(ggpp)
library(MetBrewer)
library(scico)

# ------------------------------------------------------------
# Global options
# ------------------------------------------------------------
set.seed(12345)        # master seed — fully deterministic runs
CV_K      <- 3         # k-fold CV for all models
N_TREES   <- 1000      # RF tree count
THEME_BW  <- TRUE
options(mc.cores = 1)  # ensure reproducibility

if (THEME_BW) theme_set(theme_bw())

# ============================================================
# 1) File definitions
# ============================================================
files <- list(
  Fusion = list(
    las         = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las",
    crowns      = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg",
    field       = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv",
    idcol       = "SFMcombo",
    height_col  = "SFMcomb_ht"
  ),
  
  LiDAR = list(
    las         = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las",
    crowns      = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg",
    field       = "C:/Users/User/Desktop/RandomForest3/fieldpoints_side3x.csv",
    idcol       = "Side3x",
    height_col  = "Side3x_ht"
  ),
  
  SfM = list(
    las         = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las",
    crowns      = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg",
    field       = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv",
    idcol       = "SFMonly",
    height_col  = "SFMonly_ht"
  )
)


id_map <- list(Fusion = "SFMcombo", LiDAR = "Side3x", SfM = "SFMonly")

dtm  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")


###############################################################################
# 1.2) CROSS-VALIDATION CONTROL
###############################################################################

ctrl_default <- trainControl(
  method = "cv",
  number = CV_K,
  savePredictions = "final",
  returnResamp    = "final"
)

###########################################################################
### 2) HELPER FUNCTIONS
###########################################################################

# --- 2.1: Compute metrics
compute_metrics <- function(df) {
  obs  <- df$obs_cm
  pred <- df$pred_cm

  rmse_val <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae_val  <- mean(abs(obs - pred), na.rm = TRUE)
  bias_val <- mean(pred - obs, na.rm = TRUE)

  tibble(
    R2_cor   = cor(obs, pred, use="complete.obs")^2,
    RMSE     = rmse_val,
    RMSE_pct = 100 * rmse_val / mean(obs, na.rm = TRUE),
    MAE      = mae_val,
    Bias     = bias_val
  )
}

# --- 2.2: RF out-of-sample (OOS)
collect_rf_oos <- function(model, data, label, flight) {
  best <- model$bestTune
  model$pred %>%
    semi_join(best, by = names(best)) %>%
    transmute(
      id = rowIndex,
      obs_cm  = exp(obs),
      pred_cm = exp(pred),
      Species = data$Species[rowIndex],
      Model   = label,
      Flight  = flight
    ) %>%
    group_by(id, Species, Model, Flight) %>%
    summarise(
      obs_cm  = first(obs_cm),
      pred_cm = mean(pred_cm),
      .groups = "drop"
    )
}

# --- 2.3: Log-log CV (deterministic)
collect_loglog_oos <- function(data, formula, label, flight, k = CV_K) {
  set.seed(12345)
  if (nrow(data) < 20) k <- 2

  folds <- createFolds(data$dbh_cm, k = k, list = TRUE)

  out <- lapply(folds, function(test_idx) {
    train <- data[-test_idx, ]
    test  <- data[test_idx, ]

    fit <- lm(formula, data=train)

    tibble(
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(fit, newdata=test)),
      Species = test$Species,
      Model   = label,
      Flight  = flight
    )
  })

  bind_rows(out)
}

# --- 2.4: Mixed-effects pooled CV (ME1)
evaluate_lmer_cv_pooled <- function(data, k = CV_K, label, flight) {
  set.seed(12345)
  if (nrow(data) < 20) k <- 2

  folds <- createFolds(data$dbh_cm, k = k, list=TRUE)

  out <- lapply(folds, function(test_idx) {
    train <- data[-test_idx, ]
    test  <- data[test_idx, ]

    fit <- lmer(
      log(dbh_cm) ~ log_height + log1p(crown_area_m2) + (1|Species),
      data=train
    )

    tibble(
      id      = test_idx,
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(fit, newdata=test, allow.new.levels=TRUE)),
      Species = test$Species,
      Model   = label,
      Flight  = flight
    )
  })

  oos <- bind_rows(out) %>%
    group_by(id, Species, Model, Flight) %>%
    summarise(
      obs_cm  = first(obs_cm),
      pred_cm = mean(pred_cm),
      .groups = "drop"
    )

  list(preds=oos, metrics=compute_metrics(oos))
}


###########################################################################
### 3) MAIN TRAINING LOOP (ALL 6 MODELS)
###########################################################################

all_preds_all <- list()
all_training_data <- list()

for (flight in names(files)) {

  message("Processing: ", flight)
  f <- files[[flight]]

  # --- Load LAS + crowns
  las_raw <- readLAS(f$las)
  crowns  <- st_read(f$crowns, quiet=TRUE) %>%
    mutate(treeID = as.character(treeID))

  # --- Load field
  field <- read_csv(f$field, show_col_types=FALSE) %>%
    rename(treeID = !!sym(f$idcol)) %>%
    mutate(
      treeID = as.character(str_trim(treeID)),
      Species = ifelse(grepl("Sequoia", Species, ignore.case=TRUE),
                       "Sequoia", "Other"),
      is_sequoia = as.integer(Species == "Sequoia")
    ) %>%
    rename(lidar_height = !!sym(f$height_col)) %>%
    select(treeID, dbh_cm, lidar_height, Species, is_sequoia)

  # --- Normalize LAS
  las_norm <- normalize_height(classify_ground(las_raw, csf()), knnidw())
  slope_rast <- terrain(dtm, v="slope", unit="degrees")

  # --- Crown metrics
  metrics_list <- map_dfr(seq_len(nrow(crowns)), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)

    if (is.null(las_clip) || npoints(las_clip)==0)
      return(tibble(
        treeID = crown_poly$treeID,
        crown_area_m2 = as.numeric(st_area(crown_poly)),
        point_density = NA, slope_mean = NA, ndvi_mean = NA,
        max_z = NA, mean_z = NA, p25=NA, p50=NA, p75=NA, p95=NA,
        sd_z = NA, cv_z = NA
      ))

    Z <- las_clip@data$Z

    tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(st_area(crown_poly)),
      point_density = npoints(las_clip)/as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope_rast, crown_poly, "mean"),
      ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
      max_z = max(Z),
      mean_z = mean(Z),
      p25 = quantile(Z,0.25), p50 = quantile(Z,0.50),
      p75 = quantile(Z,0.75), p95 = quantile(Z,0.95),
      sd_z = sd(Z),
      cv_z = sd(Z)/mean(Z)
    )
  })

  crowns_final <- crowns %>%
    select(treeID, geom) %>%
    left_join(metrics_list, by="treeID") %>%
    inner_join(field, by="treeID") %>%
    filter(dbh_cm > 0, dbh_cm < 350)

  crown_data <- crowns_final %>%
    mutate(
      log_height     = log(pmax(lidar_height,0.5)),
      log_crown_area = log1p(crown_area_m2),
      log_dbh        = log(dbh_cm),
      Flight = flight
    ) %>%
    drop_na()

  all_training_data[[flight]] <- crown_data


  #######################################################################
  ### 3.1  Train ALL 6 models (ME1 option)
  #######################################################################

  # 1) RF pooled
  rf_base <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z+
      point_density+slope_mean+ndvi_mean,
    data = crown_data,
    method="rf",
    trControl=ctrl_default,
    tuneLength=5,
    ntree=N_TREES
  )
  oos_rf_base <- collect_rf_oos(rf_base, crown_data, "RF (pooled)", flight)

  # 2) Mixed-effects pooled
  mix_cv <- evaluate_lmer_cv_pooled(
    crown_data, k = CV_K,
    label="Mixed-effects (pooled)", flight=flight
  )
  oos_mix <- mix_cv$preds

  # 3) RF Sequoia-only
  seq_data <- filter(crown_data, is_sequoia==1)
  seq_ctrl <- if (nrow(seq_data) < 20)
                trainControl(method="cv", number=2, savePredictions="final")
              else ctrl_default

  rf_seq <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z+
      point_density+slope_mean+ndvi_mean,
    data = seq_data,
    method="rf",
    trControl = seq_ctrl,
    tuneLength=5,
    ntree=N_TREES
  )
  oos_rf_seq <- collect_rf_oos(rf_seq, seq_data, "RF (Sequoia-only)", flight)

  # 4) RF Non-sequoia
  nonseq_data <- filter(crown_data, is_sequoia==0)

  rf_nonseq <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z+
      point_density+slope_mean+ndvi_mean,
    data = nonseq_data,
    method="rf",
    trControl=ctrl_default,
    tuneLength=5,
    ntree=N_TREES
  )
  oos_rf_nonseq <- collect_rf_oos(rf_nonseq, nonseq_data, "RF (Non-sequoia)", flight)

  # 5) Log–log Sequoia-only
  oos_loglog_seq <- collect_loglog_oos(
    seq_data,
    log(dbh_cm) ~ log_height + log(pmax(crown_area_m2,0.01)),
    "Log–log (Sequoia-only)", flight
  )

  # 6) Log–log Non-sequoia
  oos_loglog_nonseq <- collect_loglog_oos(
    nonseq_data,
    log(dbh_cm) ~ log_height + log(pmax(crown_area_m2,0.01)),
    "Log–log (Non-sequoia)", flight
  )

  # Store
  all_preds_all[[flight]] <- bind_rows(
    oos_rf_base,
    oos_mix,
    oos_rf_seq,
    oos_rf_nonseq,
    oos_loglog_seq,
    oos_loglog_nonseq
  )
}


###########################################################################
### 4) FINAL LOG–LOG MODELS USING ALL TRAINING DATA
###########################################################################

training_combined <- bind_rows(all_training_data)

final_loglog_seq <- lm(
  log(dbh_cm) ~ log_height + log(pmax(crown_area_m2,0.01)),
  data = filter(training_combined, is_sequoia==1)
)

final_loglog_nonseq <- lm(
  log(dbh_cm) ~ log_height + log(pmax(crown_area_m2,0.01)),
  data = filter(training_combined, is_sequoia==0)
)


###########################################################################
### 5) COLLAPSE OOS PREDICTIONS + METRICS
###########################################################################

preds_all <- bind_rows(all_preds_all) %>%
  filter(is.finite(obs_cm), is.finite(pred_cm)) %>%
  mutate(
    residual_cm = pred_cm - obs_cm,
    Model = factor(Model, levels = c(
      "Mixed-effects (pooled)",
      "RF (pooled)",
      "RF (Sequoia-only)",
      "Log–log (Sequoia-only)",
      "RF (Non-sequoia)",
      "Log–log (Non-sequoia)"
    ))
  )

preds_all_collapsed <- preds_all %>%
  group_by(Flight, Model, Species, obs_cm) %>%
  summarise(
    pred_cm     = mean(pred_cm),
    residual_cm = mean(residual_cm),
    .groups="drop"
  )

metrics_labels <- preds_all_collapsed %>%
  group_by(Flight, Model) %>%
  group_modify(~ compute_metrics(.x)) %>%
  ungroup() %>%
  mutate(
    label = glue(
      "R^2 = {round(R2_cor,2)}",
      "\nRMSE = {round(RMSE,1)} cm ({round(RMSE_pct,1)}%)",
      "\nMAE = {round(MAE,1)}",
      "\nBias = {round(Bias,1)}"
    )
  )



# ============================================================
# 6) Figures

# ---- 6.1 Predicted vs Observed (Fusion) ----
metrics_labels <- preds_all_collapsed %>%
  filter(Flight == "SfM") %>%
  group_by(Model) %>%
  summarise(
    r2       = cor(obs_cm, pred_cm)^2,
    rmse     = sqrt(mean((pred_cm - obs_cm)^2)),
    mae      = mean(abs(pred_cm - obs_cm)),
    bias     = mean(pred_cm - obs_cm),
    mean_obs = mean(obs_cm),
    rmse_pct = 100 * rmse / mean_obs,
    xmin     = min(obs_cm),
    xmax     = max(obs_cm),
    ymin     = min(pred_cm),
    ymax     = max(pred_cm),
    .groups = "drop"
  ) %>%
  mutate(
    # stable label positions (top-left)
    x_lab = xmin + 0.05 * (xmax - xmin),
    y_lab = ymax - 0.10 * (ymax - ymin),

    label = glue(
      "R² = {round(r2, 2)}\n",
      "RMSE = {round(rmse,1)} cm ({round(rmse_pct,1)}%)\n",
      "MAE = {round(mae,1)}\n",
      "Bias = {round(bias,1)}"
    )
  ) %>%
  mutate(npcx = 0.02, npcy = 0.98)   # for npc positioning

fig_species_scatter <- ggplot(
  preds_all_collapsed %>% filter(Flight == "SfM"),
  aes(obs_cm, pred_cm, color = Species)
) +
  # 1: 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.7, color = "gray40") +

  # 2: regression + gray shading
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "gray80", linewidth = 0.8) +

  # 3: points
  geom_point(alpha = 0.45, size = 1.8) +
  coord_cartesian(clip = "off") +

  # 4: facets
  facet_wrap(~ Model, ncol = 2) +

  # 5: clean overlay labels (top-left, stable)
  geom_text_npc(
    data = metrics_labels,
    aes(npcx = npcx, npcy = npcy, label = label),
    hjust = 0, vjust = 1,
    size = 3.5, lineheight = 1.0,
    color = "black",
    inherit.aes = FALSE
  ) +

  # 6: colors
  scale_color_manual(
    values = c("Other"="#E69F00", "Sequoia"="#56B4E9")
  ) +

  # 7: labels + theme
  labs(
    x = "Observed DBH (cm)",
    y = "Predicted DBH (cm)",
    title = "Predicted vs. Observed DBH (SfM)",
    color = "Species"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "sans"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

fig_species_scatter


ggsave(
  "20251116_Predicted_vs_Observed_DBH_SfM.png",
  fig_species_scatter,
  width = 6.7,
  height = 6,
  dpi = 600
)



# ---- 6.2 Residuals vs Predicted (Fusion, key models) ----
resid_labels <- preds_all_collapsed %>%
  filter(Flight == "SfM") %>%
  group_by(Model) %>%                      # <<< ONLY Model, NOT Species
  summarise(
    n        = n(),
    bias     = mean(residual_cm),
    rmse     = sqrt(mean(residual_cm^2)),
    mae      = mean(abs(residual_cm)),
    mean_obs = mean(obs_cm),
    rmse_pct = 100 * rmse / mean_obs,
    .groups = "drop"
  ) %>%
  mutate(
    label = glue(
      "N = {n}\n",
      "Bias = {round(bias,1)} cm\n",
      "RMSE = {round(rmse,1)} cm ({round(rmse_pct,1)}%)"
    ),
    npcx = 0.02,
    npcy = 0.02
  )


library(ggpp)



fig_residuals <- ggplot(
  preds_all_collapsed %>% filter(Flight == "SfM"),
  aes(x = pred_cm, y = residual_cm, color = Species)
) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.7, color = "gray40") +
geom_smooth(method="lm", se=TRUE,
            color="black", fill="gray80",
            linewidth=0.8, linetype="dotted")+
  geom_point(alpha = 0.45, size = 1.8) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ Model, ncol = 2) +

  geom_text_npc(
    data = resid_labels,
    aes(npcx = npcx, npcy = npcy, label = label),
    hjust = 0, vjust = 0,
    size = 3.5,
    color = "black",
    lineheight = 1.0,
    inherit.aes = FALSE
  ) +

  scale_color_manual(values = c("Other" = "#E69F00", "Sequoia" = "#56B4E9")) +
  labs(
    x = "Predicted DBH (cm)",
    y = "Residual (Predicted – Observed, cm)",
    title = "Residuals vs. Predicted DBH (SfM)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "sans"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

fig_residuals



ggsave(
  "20251116_Residuals_DBH_SfM.png",
  fig_residuals,
  width = 6.7,
  height = 6,
  dpi = 600
)

# ---- 6.3 Predict for Full Data and Map (LiDAR)-----
las_full   <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las")
crowns_full <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Cross/3x_Cross_crowns_all.gpkg")%>%
  rename(lidar_height = tree_height_m)
dtm_full  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi_full <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

compute_crown_metrics <- function(las, crowns, dtm, ndvi) {

  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(dtm, v = "slope", unit = "degrees")

  metrics_list <- purrr::map_dfr(seq_len(nrow(crowns)), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- lidR::clip_roi(las_norm, crown_poly)

    # If no LiDAR returns inside crown polygon
    if (is.null(las_clip) || lidR::npoints(las_clip) == 0) {
      return(tibble::tibble(
        treeID        = crown_poly$treeID,
        crown_area_m2 = as.numeric(sf::st_area(crown_poly)),
        point_density = NA,
        slope_mean    = NA,
        ndvi_mean     = NA,
        max_z         = NA,
        mean_z        = NA,
        p25           = NA,
        p50           = NA,
        p75           = NA,
        p95           = NA,
        sd_z          = NA,
        cv_z          = NA
      ))
    }

    Z <- las_clip@data$Z

    tibble::tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(sf::st_area(crown_poly)),
      point_density = lidR::npoints(las_clip) / as.numeric(sf::st_area(crown_poly)),
      slope_mean    = exactextractr::exact_extract(slope, crown_poly, "mean"),
      ndvi_mean     = exactextractr::exact_extract(ndvi, crown_poly, "mean"),
      max_z         = max(Z, na.rm = TRUE),
      mean_z        = mean(Z, na.rm = TRUE),
      p25           = quantile(Z, 0.25, na.rm = TRUE),
      p50           = quantile(Z, 0.50, na.rm = TRUE),
      p75           = quantile(Z, 0.75, na.rm = TRUE),
      p95           = quantile(Z, 0.95, na.rm = TRUE),
      sd_z          = sd(Z, na.rm = TRUE),
      cv_z          = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE)
    )
  })

  dplyr::left_join(
    crowns %>% dplyr::select(treeID, geom),
    metrics_list,
    by = "treeID"
  )
}


crown_metrics_full <- compute_crown_metrics(
    las = las_full,
    crowns = crowns_full,
    dtm = dtm_full,
    ndvi = ndvi_full
  )

crown_metrics_full <- crown_metrics_full %>%
  left_join(
    crowns_full %>% 
      st_drop_geometry() %>% 
      select(treeID, lidar_height),
    by = "treeID"
  )

crown_data_full <- crown_metrics_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height     = log(lidar_height),
    log_crown_area = log1p(crown_area_m2)
  )

good <- complete.cases(crown_data_full)
rf_predictions_full <- exp(predict(rf_base, newdata = crown_data_full[good, ]))


pred_table <- tibble(
  treeID = crown_data_full$treeID[good],
  pred_DBH_cm = rf_predictions_full
)

crowns_full <- crowns_full %>%
  left_join(pred_table, by = "treeID")





ggplot() +
  geom_sf(
    data = crowns_full,
    aes(fill = pred_DBH_cm),
    color = "grey15",   # thin crown outlines
    size  = 0.10        # subtle boundary
  ) +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = c(50, 300),       # adjust if needed
    breaks = seq(50, 300, 50)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    text = element_text(family = "sans"),
    legend.key.height = unit(0.55, "cm"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Predicted Diameter at Breast Height (DBH)",
    subtitle = "Cloud2Trees crowns + Random Forest model",
    #caption = "Grey polygons = No crown metrics / no prediction",
    x = NULL,
    y = NULL
  )
  
------------------
  

# =============================
# 1. Split Sequoias / Non-Sequoias
# =============================
sequoia_full <- crowns_full %>% filter(is_sequoia == 1)
nonseq_full  <- crowns_full %>% filter(is_sequoia == 0)

# =============================
# 2. Build prediction dataframes
# =============================
pred_df_seq <- sequoia_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log(pmax(crown_area_m2, 0.01))
  )

pred_df_nonseq <- nonseq_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log(pmax(crown_area_m2, 0.01))
  )

# =============================
# 3. Predict Sequoias (log–log sequoia)
# =============================
preds_seq <- exp(predict(final_loglog_seq, newdata = pred_df_seq))

sequoia_full <- sequoia_full %>%
  select(-matches("^pred_DBH_cm")) %>%
  mutate(pred_DBH_cm = preds_seq)

# =============================
# 4. Predict Non-Sequoias (log–log non-sequoia)
# =============================
preds_nonseq <- exp(predict(final_loglog_nonseq, newdata = pred_df_nonseq))

nonseq_full <- nonseq_full %>%
  select(-matches("^pred_DBH_cm")) %>%
  mutate(pred_DBH_cm = preds_nonseq)

# =============================
# 5. Combine all crowns back together
# =============================
crowns_loglog_all <- bind_rows(sequoia_full, nonseq_full)

ggplot() +
  geom_sf(
    data = crowns_loglog_all,
    aes(fill = pred_DBH_cm),
    color = "grey20",
    size  = 0.12
  ) +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    na.value = "grey90",
    limits = c(50, 320)
  ) +
  theme_bw() +
  labs(
    title = "DBH Predictions – Log–Log Model (Sequoia + Non-Sequoia)",
    subtitle = "Sequoias predicted with Sequoia log–log; others with Non-Sequoia log–log",
    x = NULL, y = NULL
  )


-------------

  
# Reproject to lon/lat
crowns_ll <- st_transform(crowns_full, 4326)

# Compute centroids (1 per feature)
centroids <- st_centroid(crowns_ll)

# Extract lon/lat from centroids
crowns_ll$lon <- st_coordinates(centroids)[, "X"]
crowns_ll$lat <- st_coordinates(centroids)[, "Y"]

# Compute limits
lon_min <- min(crowns_ll$lon, na.rm = TRUE)
lon_max <- max(crowns_ll$lon, na.rm = TRUE)
lat_min <- min(crowns_ll$lat, na.rm = TRUE)
lat_max <- max(crowns_ll$lat, na.rm = TRUE)


pad <- 0.0003   # ~30 m padding
lon_mid <- mean(c(lon_min, lon_max))


# ------------------------------------------------------------
# 1. Build the Log–Log model map
# ------------------------------------------------------------
map_rf <- ggplot() +
  geom_sf(
    data = crowns_ll,
    aes(fill = pred_DBH_cm),
    color = "grey15",
    linewidth = 0.10
  ) +

  scale_fill_scico(
    palette   = "vik",
    name      = "Predicted DBH (cm)",
    direction = 1,
    na.value  = "grey80",
    limits    = range(crowns_ll$pred_DBH_cm, na.rm = TRUE),
    breaks    = seq(50, 300, 50)
  ) +

  # ---- AXIS CONTROL: 3 longitude ticks ----
  scale_x_continuous(
    breaks = c(lon_min, lon_mid, lon_max),
    labels = scales::label_number(accuracy = 0.0001)
  ) +

  coord_sf(
    xlim   = c(lon_min - pad, lon_max + pad),
    ylim   = c(lat_min - pad, lat_max + pad),
    expand = FALSE
  ) +

  labs(
    title = "Random Forest (Pooled)",
    x = NULL,
    y = NULL
  ) +

  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    text = element_text(family = "sans"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(size = 10),

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    #legend.key.width = unit(1.4, "cm"),
    #legend.key.height = unit(0.35, "cm"),
    legend.margin = margin(t = -5),
    legend.key.width  = unit(0.7, "cm"),   # shorter
    legend.key.height = unit(0.28, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

map_rf


ggsave(
  "20251116_map_rf.png",
  map_rf,
  width = 4,
  height = 4,
  dpi = 300
)





ggsave(
  "20251116_LiDAR_RF.png",
  map_rf,
  width = 4,
  height = 4,
  dpi = 300
)

  
  

ggsave(
  "20251116_LiDAR_RF.png",
  map_rf,
  width = 7,
  height = 7,
  dpi = 300
)

st_write(
  crowns_full,
  "crowns_full.gpkg",
  layer = "crowns_full",
  delete_dsn = TRUE
)

# ============================================================
# 2. Make LOG–LOG FULL-FOOTPRINT MAP
# ============================================================


# Reproject to lon/lat
crowns_loglog_all_ll <- st_transform(crowns_loglog_all, 4326)

# Compute centroids (1 per feature)
centroids <- st_centroid(crowns_loglog_all_ll)

# Extract lon/lat from centroids
crowns_loglog_all_ll$lon <- st_coordinates(centroids)[, "X"]
crowns_loglog_all_ll$lat <- st_coordinates(centroids)[, "Y"]

# Compute limits
lon_min <- min(crowns_loglog_all_ll$lon, na.rm = TRUE)
lon_max <- max(crowns_loglog_all_ll$lon, na.rm = TRUE)
lat_min <- min(crowns_loglog_all_ll$lat, na.rm = TRUE)
lat_max <- max(crowns_loglog_all_ll$lat, na.rm = TRUE)


pad <- 0.0003   # ~30 m padding
lon_mid <- mean(c(lon_min, lon_max))


# ------------------------------------------------------------
# 1. Build the Log–Log model map
# ------------------------------------------------------------
map_loglog <- ggplot() +
  geom_sf(
    data = crowns_loglog_all_ll,
    aes(fill = pred_DBH_cm),
    color = "grey15",
    linewidth = 0.10
  ) +

  scale_fill_scico(
    palette   = "vik",
    name      = "Predicted DBH (cm)",
    direction = 1,
    na.value  = "grey80",
    limits    = range(crowns_loglog_all_ll$pred_DBH_cm, na.rm = TRUE),
    breaks    = seq(50, 300, 50)
  ) +

  # ---- AXIS CONTROL: 3 longitude ticks ----
  scale_x_continuous(
    breaks = c(lon_min, lon_mid, lon_max),
    labels = scales::label_number(accuracy = 0.0001)
  ) +

  coord_sf(
    xlim   = c(lon_min - pad, lon_max + pad),
    ylim   = c(lat_min - pad, lat_max + pad),
    expand = FALSE
  ) +

  labs(
    title = "Log–Log Model (Sequoia + Non-Sequoia)",
    x = NULL,
    y = NULL
  ) +

  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    text = element_text(family = "sans"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(size = 10),

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    #legend.key.width = unit(1.4, "cm"),
    #legend.key.height = unit(0.35, "cm"),
    legend.margin = margin(t = -5),
    legend.key.width  = unit(0.7, "cm"),   # shorter
    legend.key.height = unit(0.28, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

map_loglog


ggsave(
  "20251116_map_loglog.png",
  map_loglog,
  width = 4,
  height = 4,
  dpi = 300
)



st_write(
  crowns_loglog_all,
  "crowns_loglog_all.gpkg",
  layer = "crowns_loglog_all",
  delete_dsn = TRUE
)



# ============================================================
# 3. SIDE-BY-SIDE COMPARISON MAP
# ============================================================

library(patchwork)

# Adjust RF plot
map_rf <- map_rf +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.key.width  = unit(0.75, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.margin = margin(t = -5)
  )

# Adjust Log-Log plot
map_loglog <- map_loglog +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.key.width  = unit(0.75, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.margin = margin(t = -5)
  )

map_rf <- map_rf +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = range(crowns_ll$pred_DBH_cm, na.rm = TRUE),
    breaks = scales::pretty_breaks(5)
  )

map_loglog <- map_loglog +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = range(crowns_loglog_all_ll$pred_DBH_cm, na.rm = TRUE),
    breaks = scales::pretty_breaks(5)
  )


# ------------------------------------------------------------
# FINAL 2-PANEL SIDE-BY-SIDE COMPARISON
# ------------------------------------------------------------

library(patchwork)

final_compare <- 
  map_rf + map_loglog +
  plot_layout(
    ncol = 2,
    widths = c(1, 1),          # equal panel widths
    guides = "keep"            # keep legends separate
  ) &
  theme(
    plot.margin = margin(5, 5, 5, 5),
    panel.spacing = unit(0, "pt")   # remove the big center gap
  )

final_compare <- 
  map_rf + plot_spacer() + map_loglog +
  plot_layout(
    ncol = 3,
    widths = c(1, -0.1, 1),    # negative width squeezes the middle space
    guides = "keep"
  ) &
  theme(
    panel.spacing = unit(0, "pt"),
    plot.margin = margin(5, 5, 5, 5)
  )


final_compare



ggsave(
  "DBH_compare_google.png",
  plot   = final_compare,
  width  = 8,        # wide
  height = 3.2,       # SHORT!!
  dpi    = 600,
  bg     = "white"
)


ggsave(
  "DBH_Comparison_highres.png",
  plot   = final_compare,
  width  = 10,       # adjust to fit your layout
  height = 6,
  units  = "in",
  dpi    = 600,      # extremely sharp in Google Docs
  bg     = "white"
)



ggsave(
  "DBH_Comparison.tiff",
  plot   = final_compare,
  width  = 9,
  height = 5,
  units  = "in",
  dpi    = 800,
  compression = "lzw"
)

ggsave(
  "20251116_RFvsLog_LiDAR.png",
  final_compare,
  width = 8,
  height = 8,
  dpi = 300
)


ggsave(
  "20251116_RFvsLog_LiDAR.png",
  plot = final_compare,
  width = 9,        # much smaller
  height = 4.2,
  units = "in",
  dpi = 300
)


ggsave(
  "20251116_RFvsLog_LiDAR.png",
  plot = final_compare,
 width = 8,
  height = 4.2,
  dpi = 150
)
------------
  
  
library(scico)
library(patchwork)


  
  
  
  
las_full   <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las")

crowns_full <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Cross/3x_Cross_crowns_all.gpkg") %>%
  rename(lidar_height = tree_height_m)

dtm_full  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi_full <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")


compute_crown_metrics <- function(las, crowns, dtm, ndvi) {

  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(dtm, v = "slope", unit = "degrees")

  metrics_list <- map_dfr(seq_len(nrow(crowns)), function(i) {

    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)

    if (is.null(las_clip) || npoints(las_clip) == 0) {
      return(tibble(
        treeID        = crown_poly$treeID,
        crown_area_m2 = as.numeric(st_area(crown_poly)),
        point_density = NA,
        slope_mean    = NA,
        ndvi_mean     = NA,
        max_z = NA, mean_z = NA, p25 = NA, p50 = NA, p75 = NA, p95 = NA,
        sd_z = NA, cv_z = NA
      ))
    }

    Z <- las_clip@data$Z

    tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(st_area(crown_poly)),
      point_density = npoints(las_clip) / as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope, crown_poly, "mean"),
      ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
      max_z = max(Z), mean_z = mean(Z),
      p25 = quantile(Z, 0.25), p50 = quantile(Z, 0.50),
      p75 = quantile(Z, 0.75), p95 = quantile(Z, 0.95),
      sd_z = sd(Z), cv_z = sd(Z) / mean(Z)
    )
  })

  left_join(
    crowns %>% select(treeID, geom),
    metrics_list,
    by = "treeID"
  )
}

crown_metrics_full <- compute_crown_metrics(
  las = las_full, crowns = crowns_full, dtm = dtm_full, ndvi = ndvi_full
)

# attach lidar height
crown_metrics_full <- crown_metrics_full %>%
  left_join(
    crowns_full %>% st_drop_geometry() %>% select(treeID, lidar_height, is_sequoia),
    by = "treeID"
  )

# prediction dataframe
crown_data_full <- crown_metrics_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log1p(crown_area_m2)
  )

good <- complete.cases(crown_data_full)
rf_pred <- exp(predict(rf_base, newdata = crown_data_full[good, ]))

rf_table <- tibble(
  treeID = crown_data_full$treeID[good],
  pred_RF = rf_pred
)

crowns_full <- crowns_full %>%
  left_join(rf_table, by = "treeID")   # RF predictions stored as pred_RF

# Split
seq_full    <- crowns_full %>% filter(is_sequoia == 1)
nonseq_full <- crowns_full %>% filter(is_sequoia == 0)

# Prediction dataframes
pred_seq <- seq_full %>% st_drop_geometry() %>%
  mutate(log_height = log(lidar_height), log_area = log1p(crown_area_m2))

pred_nonseq <- nonseq_full %>% st_drop_geometry() %>%
  mutate(log_height = log(lidar_height), log_area = log1p(crown_area_m2))

# Predictions
seq_pred    <- exp(predict(final_loglog_seq,    newdata = pred_seq))
nonseq_pred <- exp(predict(final_loglog_nonseq, newdata = pred_nonseq))

# Attach predictions
seq_full    <- seq_full    %>% mutate(pred_LogLog = seq_pred)
nonseq_full <- nonseq_full %>% mutate(pred_LogLog = nonseq_pred)

# Combine both species
crowns_loglog_all <- bind_rows(seq_full, nonseq_full)

map_rf <- ggplot() +
  geom_sf(data = crowns_full, aes(fill = pred_RF), color = "grey10", size = 0.1) +
  scale_fill_scico(palette = "vik", name = "Predicted DBH (cm)",
                   na.value = "grey80", limits = c(50,300)) +
  theme_bw() +
  labs(title = "Random Forest (Pooled)", subtitle = "Full-Footprint DBH")

map_loglog <- ggplot() +
  geom_sf(data = crowns_loglog_all, aes(fill = pred_LogLog), color = "grey10", size = 0.1) +
  scale_fill_scico(palette = "vik", name = "Predicted DBH (cm)",
                   na.value = "grey80", limits = c(50,300)) +
  theme_bw() +
  labs(title = "Log–Log (Species-Specific)", subtitle = "Full-Footprint DBH")

final_compare <- map_rf | map_loglog +
  plot_annotation(
    title = "DBH Predictions Across Full Footprint",
    subtitle = "Random Forest (pooled) vs Log–Log (species-specific)"
  )

final_compare







#######################################3


# ---- 6.4 Predict for Full Data and Map (Fusion)-----
sfm_clipped <- st_intersection(sfm_crowns, st_union(lidar_crowns))


las_full   <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las")
crowns_full <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Cross/3x_Cross_crowns_all.gpkg")%>%
  rename(lidar_height = tree_height_m)
dtm_full  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi_full <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

compute_crown_metrics <- function(las, crowns, dtm, ndvi) {

  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(dtm, v = "slope", unit = "degrees")

  metrics_list <- purrr::map_dfr(seq_len(nrow(crowns)), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- lidR::clip_roi(las_norm, crown_poly)

    # If no LiDAR returns inside crown polygon
    if (is.null(las_clip) || lidR::npoints(las_clip) == 0) {
      return(tibble::tibble(
        treeID        = crown_poly$treeID,
        crown_area_m2 = as.numeric(sf::st_area(crown_poly)),
        point_density = NA,
        slope_mean    = NA,
        ndvi_mean     = NA,
        max_z         = NA,
        mean_z        = NA,
        p25           = NA,
        p50           = NA,
        p75           = NA,
        p95           = NA,
        sd_z          = NA,
        cv_z          = NA
      ))
    }

    Z <- las_clip@data$Z

    tibble::tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(sf::st_area(crown_poly)),
      point_density = lidR::npoints(las_clip) / as.numeric(sf::st_area(crown_poly)),
      slope_mean    = exactextractr::exact_extract(slope, crown_poly, "mean"),
      ndvi_mean     = exactextractr::exact_extract(ndvi, crown_poly, "mean"),
      max_z         = max(Z, na.rm = TRUE),
      mean_z        = mean(Z, na.rm = TRUE),
      p25           = quantile(Z, 0.25, na.rm = TRUE),
      p50           = quantile(Z, 0.50, na.rm = TRUE),
      p75           = quantile(Z, 0.75, na.rm = TRUE),
      p95           = quantile(Z, 0.95, na.rm = TRUE),
      sd_z          = sd(Z, na.rm = TRUE),
      cv_z          = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE)
    )
  })

  dplyr::left_join(
    crowns %>% dplyr::select(treeID, geom),
    metrics_list,
    by = "treeID"
  )
}


crown_metrics_full <- compute_crown_metrics(
    las = las_full,
    crowns = crowns_full,
    dtm = dtm_full,
    ndvi = ndvi_full
  )

crown_metrics_full <- crown_metrics_full %>%
  left_join(
    crowns_full %>% 
      st_drop_geometry() %>% 
      select(treeID, lidar_height),
    by = "treeID"
  )

crown_data_full <- crown_metrics_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height     = log(lidar_height),
    log_crown_area = log1p(crown_area_m2)
  )

good <- complete.cases(crown_data_full)
rf_predictions_full <- exp(predict(rf_base, newdata = crown_data_full[good, ]))


pred_table <- tibble(
  treeID = crown_data_full$treeID[good],
  pred_DBH_cm = rf_predictions_full
)

crowns_full <- crowns_full %>%
  left_join(pred_table, by = "treeID")





ggplot() +
  geom_sf(
    data = crowns_full,
    aes(fill = pred_DBH_cm),
    color = "grey15",   # thin crown outlines
    size  = 0.10        # subtle boundary
  ) +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = c(50, 300),       # adjust if needed
    breaks = seq(50, 300, 50)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    text = element_text(family = "sans"),
    legend.key.height = unit(0.55, "cm"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Predicted Diameter at Breast Height (DBH)",
    subtitle = "Cloud2Trees crowns + Random Forest model",
    #caption = "Grey polygons = No crown metrics / no prediction",
    x = NULL,
    y = NULL
  )
  
------------------
  

# =============================
# 1. Split Sequoias / Non-Sequoias
# =============================
sequoia_full <- crowns_full %>% filter(is_sequoia == 1)
nonseq_full  <- crowns_full %>% filter(is_sequoia == 0)

# =============================
# 2. Build prediction dataframes
# =============================
pred_df_seq <- sequoia_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log(pmax(crown_area_m2, 0.01))
  )

pred_df_nonseq <- nonseq_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log(pmax(crown_area_m2, 0.01))
  )

# =============================
# 3. Predict Sequoias (log–log sequoia)
# =============================
preds_seq <- exp(predict(final_loglog_seq, newdata = pred_df_seq))

sequoia_full <- sequoia_full %>%
  select(-matches("^pred_DBH_cm")) %>%
  mutate(pred_DBH_cm = preds_seq)

# =============================
# 4. Predict Non-Sequoias (log–log non-sequoia)
# =============================
preds_nonseq <- exp(predict(final_loglog_nonseq, newdata = pred_df_nonseq))

nonseq_full <- nonseq_full %>%
  select(-matches("^pred_DBH_cm")) %>%
  mutate(pred_DBH_cm = preds_nonseq)

# =============================
# 5. Combine all crowns back together
# =============================
crowns_loglog_all <- bind_rows(sequoia_full, nonseq_full)

ggplot() +
  geom_sf(
    data = crowns_loglog_all,
    aes(fill = pred_DBH_cm),
    color = "grey20",
    size  = 0.12
  ) +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    na.value = "grey90",
    limits = c(50, 320)
  ) +
  theme_bw() +
  labs(
    title = "DBH Predictions – Log–Log Model (Sequoia + Non-Sequoia)",
    subtitle = "Sequoias predicted with Sequoia log–log; others with Non-Sequoia log–log",
    x = NULL, y = NULL
  )


-------------

  
# Reproject to lon/lat
crowns_ll <- st_transform(crowns_full, 4326)

# Compute centroids (1 per feature)
centroids <- st_centroid(crowns_ll)

# Extract lon/lat from centroids
crowns_ll$lon <- st_coordinates(centroids)[, "X"]
crowns_ll$lat <- st_coordinates(centroids)[, "Y"]

# Compute limits
lon_min <- min(crowns_ll$lon, na.rm = TRUE)
lon_max <- max(crowns_ll$lon, na.rm = TRUE)
lat_min <- min(crowns_ll$lat, na.rm = TRUE)
lat_max <- max(crowns_ll$lat, na.rm = TRUE)


pad <- 0.0003   # ~30 m padding
lon_mid <- mean(c(lon_min, lon_max))


# ------------------------------------------------------------
# 1. Build the Log–Log model map
# ------------------------------------------------------------
map_rf <- ggplot() +
  geom_sf(
    data = crowns_ll,
    aes(fill = pred_DBH_cm),
    color = "grey15",
    linewidth = 0.10
  ) +

  scale_fill_scico(
    palette   = "vik",
    name      = "Predicted DBH (cm)",
    direction = 1,
    na.value  = "grey80",
    limits    = range(crowns_ll$pred_DBH_cm, na.rm = TRUE),
    breaks    = seq(50, 300, 50)
  ) +

  # ---- AXIS CONTROL: 3 longitude ticks ----
  scale_x_continuous(
    breaks = c(lon_min, lon_mid, lon_max),
    labels = scales::label_number(accuracy = 0.0001)
  ) +

  coord_sf(
    xlim   = c(lon_min - pad, lon_max + pad),
    ylim   = c(lat_min - pad, lat_max + pad),
    expand = FALSE
  ) +

  labs(
    title = "Random Forest (Pooled)",
    x = NULL,
    y = NULL
  ) +

  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    text = element_text(family = "sans"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(size = 10),

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    #legend.key.width = unit(1.4, "cm"),
    #legend.key.height = unit(0.35, "cm"),
    legend.margin = margin(t = -5),
    legend.key.width  = unit(0.7, "cm"),   # shorter
    legend.key.height = unit(0.28, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

map_rf


ggsave(
  "20251116_map_rf.png",
  map_rf,
  width = 4,
  height = 4,
  dpi = 300
)





ggsave(
  "20251116_LiDAR_RF.png",
  map_rf,
  width = 4,
  height = 4,
  dpi = 300
)

  
  

ggsave(
  "20251116_LiDAR_RF.png",
  map_rf,
  width = 7,
  height = 7,
  dpi = 300
)

st_write(
  crowns_full,
  "crowns_full.gpkg",
  layer = "crowns_full",
  delete_dsn = TRUE
)

# ============================================================
# 2. Make LOG–LOG FULL-FOOTPRINT MAP
# ============================================================


# Reproject to lon/lat
crowns_loglog_all_ll <- st_transform(crowns_loglog_all, 4326)

# Compute centroids (1 per feature)
centroids <- st_centroid(crowns_loglog_all_ll)

# Extract lon/lat from centroids
crowns_loglog_all_ll$lon <- st_coordinates(centroids)[, "X"]
crowns_loglog_all_ll$lat <- st_coordinates(centroids)[, "Y"]

# Compute limits
lon_min <- min(crowns_loglog_all_ll$lon, na.rm = TRUE)
lon_max <- max(crowns_loglog_all_ll$lon, na.rm = TRUE)
lat_min <- min(crowns_loglog_all_ll$lat, na.rm = TRUE)
lat_max <- max(crowns_loglog_all_ll$lat, na.rm = TRUE)


pad <- 0.0003   # ~30 m padding
lon_mid <- mean(c(lon_min, lon_max))


# ------------------------------------------------------------
# 1. Build the Log–Log model map
# ------------------------------------------------------------
map_loglog <- ggplot() +
  geom_sf(
    data = crowns_loglog_all_ll,
    aes(fill = pred_DBH_cm),
    color = "grey15",
    linewidth = 0.10
  ) +

  scale_fill_scico(
    palette   = "vik",
    name      = "Predicted DBH (cm)",
    direction = 1,
    na.value  = "grey80",
    limits    = range(crowns_loglog_all_ll$pred_DBH_cm, na.rm = TRUE),
    breaks    = seq(50, 300, 50)
  ) +

  # ---- AXIS CONTROL: 3 longitude ticks ----
  scale_x_continuous(
    breaks = c(lon_min, lon_mid, lon_max),
    labels = scales::label_number(accuracy = 0.0001)
  ) +

  coord_sf(
    xlim   = c(lon_min - pad, lon_max + pad),
    ylim   = c(lat_min - pad, lat_max + pad),
    expand = FALSE
  ) +

  labs(
    title = "Log–Log Model (Sequoia + Non-Sequoia)",
    x = NULL,
    y = NULL
  ) +

  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    text = element_text(family = "sans"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(size = 10),

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    #legend.key.width = unit(1.4, "cm"),
    #legend.key.height = unit(0.35, "cm"),
    legend.margin = margin(t = -5),
    legend.key.width  = unit(0.7, "cm"),   # shorter
    legend.key.height = unit(0.28, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

map_loglog


ggsave(
  "20251116_map_loglog.png",
  map_loglog,
  width = 4,
  height = 4,
  dpi = 300
)



st_write(
  crowns_loglog_all,
  "crowns_loglog_all.gpkg",
  layer = "crowns_loglog_all",
  delete_dsn = TRUE
)



# ============================================================
# 3. SIDE-BY-SIDE COMPARISON MAP
# ============================================================

library(patchwork)

# Adjust RF plot
map_rf <- map_rf +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.key.width  = unit(0.75, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.margin = margin(t = -5)
  )

# Adjust Log-Log plot
map_loglog <- map_loglog +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.key.width  = unit(0.75, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.margin = margin(t = -5)
  )

map_rf <- map_rf +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = range(crowns_ll$pred_DBH_cm, na.rm = TRUE),
    breaks = scales::pretty_breaks(5)
  )

map_loglog <- map_loglog +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)",
    direction = 1,
    na.value = "grey80",
    limits = range(crowns_loglog_all_ll$pred_DBH_cm, na.rm = TRUE),
    breaks = scales::pretty_breaks(5)
  )


# ------------------------------------------------------------
# FINAL 2-PANEL SIDE-BY-SIDE COMPARISON
# ------------------------------------------------------------

library(patchwork)

final_compare <- 
  map_rf + map_loglog +
  plot_layout(
    ncol = 2,
    widths = c(1, 1),          # equal panel widths
    guides = "keep"            # keep legends separate
  ) &
  theme(
    plot.margin = margin(5, 5, 5, 5),
    panel.spacing = unit(0, "pt")   # remove the big center gap
  )

final_compare <- 
  map_rf + plot_spacer() + map_loglog +
  plot_layout(
    ncol = 3,
    widths = c(1, -0.1, 1),    # negative width squeezes the middle space
    guides = "keep"
  ) &
  theme(
    panel.spacing = unit(0, "pt"),
    plot.margin = margin(5, 5, 5, 5)
  )


final_compare



ggsave(
  "DBH_compare_google.png",
  plot   = final_compare,
  width  = 8,        # wide
  height = 3.2,       # SHORT!!
  dpi    = 600,
  bg     = "white"
)


ggsave(
  "DBH_Comparison_highres.png",
  plot   = final_compare,
  width  = 10,       # adjust to fit your layout
  height = 6,
  units  = "in",
  dpi    = 600,      # extremely sharp in Google Docs
  bg     = "white"
)



ggsave(
  "DBH_Comparison.tiff",
  plot   = final_compare,
  width  = 9,
  height = 5,
  units  = "in",
  dpi    = 800,
  compression = "lzw"
)

ggsave(
  "20251116_RFvsLog_LiDAR.png",
  final_compare,
  width = 8,
  height = 8,
  dpi = 300
)


ggsave(
  "20251116_RFvsLog_LiDAR.png",
  plot = final_compare,
  width = 9,        # much smaller
  height = 4.2,
  units = "in",
  dpi = 300
)


ggsave(
  "20251116_RFvsLog_LiDAR.png",
  plot = final_compare,
 width = 8,
  height = 4.2,
  dpi = 150
)
------------
  
  
library(scico)
library(patchwork)


  
  
  
  
las_full   <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las")

crowns_full <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Cross/3x_Cross_crowns_all.gpkg") %>%
  rename(lidar_height = tree_height_m)

dtm_full  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi_full <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")


compute_crown_metrics <- function(las, crowns, dtm, ndvi) {

  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(dtm, v = "slope", unit = "degrees")

  metrics_list <- map_dfr(seq_len(nrow(crowns)), function(i) {

    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)

    if (is.null(las_clip) || npoints(las_clip) == 0) {
      return(tibble(
        treeID        = crown_poly$treeID,
        crown_area_m2 = as.numeric(st_area(crown_poly)),
        point_density = NA,
        slope_mean    = NA,
        ndvi_mean     = NA,
        max_z = NA, mean_z = NA, p25 = NA, p50 = NA, p75 = NA, p95 = NA,
        sd_z = NA, cv_z = NA
      ))
    }

    Z <- las_clip@data$Z

    tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(st_area(crown_poly)),
      point_density = npoints(las_clip) / as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope, crown_poly, "mean"),
      ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
      max_z = max(Z), mean_z = mean(Z),
      p25 = quantile(Z, 0.25), p50 = quantile(Z, 0.50),
      p75 = quantile(Z, 0.75), p95 = quantile(Z, 0.95),
      sd_z = sd(Z), cv_z = sd(Z) / mean(Z)
    )
  })

  left_join(
    crowns %>% select(treeID, geom),
    metrics_list,
    by = "treeID"
  )
}

crown_metrics_full <- compute_crown_metrics(
  las = las_full, crowns = crowns_full, dtm = dtm_full, ndvi = ndvi_full
)

# attach lidar height
crown_metrics_full <- crown_metrics_full %>%
  left_join(
    crowns_full %>% st_drop_geometry() %>% select(treeID, lidar_height, is_sequoia),
    by = "treeID"
  )

# prediction dataframe
crown_data_full <- crown_metrics_full %>%
  st_drop_geometry() %>%
  mutate(
    log_height = log(lidar_height),
    log_area   = log1p(crown_area_m2)
  )

good <- complete.cases(crown_data_full)
rf_pred <- exp(predict(rf_base, newdata = crown_data_full[good, ]))

rf_table <- tibble(
  treeID = crown_data_full$treeID[good],
  pred_RF = rf_pred
)

crowns_full <- crowns_full %>%
  left_join(rf_table, by = "treeID")   # RF predictions stored as pred_RF

# Split
seq_full    <- crowns_full %>% filter(is_sequoia == 1)
nonseq_full <- crowns_full %>% filter(is_sequoia == 0)

# Prediction dataframes
pred_seq <- seq_full %>% st_drop_geometry() %>%
  mutate(log_height = log(lidar_height), log_area = log1p(crown_area_m2))

pred_nonseq <- nonseq_full %>% st_drop_geometry() %>%
  mutate(log_height = log(lidar_height), log_area = log1p(crown_area_m2))

# Predictions
seq_pred    <- exp(predict(final_loglog_seq,    newdata = pred_seq))
nonseq_pred <- exp(predict(final_loglog_nonseq, newdata = pred_nonseq))

# Attach predictions
seq_full    <- seq_full    %>% mutate(pred_LogLog = seq_pred)
nonseq_full <- nonseq_full %>% mutate(pred_LogLog = nonseq_pred)

# Combine both species
crowns_loglog_all <- bind_rows(seq_full, nonseq_full)

map_rf <- ggplot() +
  geom_sf(data = crowns_full, aes(fill = pred_RF), color = "grey10", size = 0.1) +
  scale_fill_scico(palette = "vik", name = "Predicted DBH (cm)",
                   na.value = "grey80", limits = c(50,300)) +
  theme_bw() +
  labs(title = "Random Forest (Pooled)", subtitle = "Full-Footprint DBH")

map_loglog <- ggplot() +
  geom_sf(data = crowns_loglog_all, aes(fill = pred_LogLog), color = "grey10", size = 0.1) +
  scale_fill_scico(palette = "vik", name = "Predicted DBH (cm)",
                   na.value = "grey80", limits = c(50,300)) +
  theme_bw() +
  labs(title = "Log–Log (Species-Specific)", subtitle = "Full-Footprint DBH")

final_compare <- map_rf | map_loglog +
  plot_annotation(
    title = "DBH Predictions Across Full Footprint",
    subtitle = "Random Forest (pooled) vs Log–Log (species-specific)"
  )

final_compare