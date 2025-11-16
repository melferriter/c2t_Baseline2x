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

set.seed(444)

# ------------------------------------------------------------
# Global options
# ------------------------------------------------------------
CV_K      <- 3      # k-fold CV across all models
N_TREES   <- 1000   # RF trees
THEME_BW  <- TRUE

if (THEME_BW) theme_set(theme_bw())

# ============================================================
# 1) File definitions
# ============================================================
files <- list(
  Fusion = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv",
    idcol  = "SFMcombo"
  ),
  LiDAR = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_side3x.csv",
    idcol  = "Side3x"
  ),
  SfM = list(
    las    = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las",
    crowns = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg",
    field  = "C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv",
    idcol  = "SFMonly"
  )
)

id_map <- list(Fusion = "SFMcombo", LiDAR = "Side3x", SfM = "SFMonly")

dtm  <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")

# ============================================================
# 2) Cross-validation setup
# ============================================================
ctrl <- trainControl(
  method          = "cv",
  number          = CV_K,
  savePredictions = "final",
  returnResamp    = "final"
)

# ============================================================
# 3) Helpers
# ============================================================

# ---- 3.1 Basic metrics ----
compute_metrics <- function(df) {
  obs  <- df$obs_cm
  pred <- df$pred_cm

  rmse_val <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae_val  <- mean(abs(obs - pred), na.rm = TRUE)
  bias_val <- mean(pred - obs, na.rm = TRUE)

  tibble(
    R2_cor    = cor(obs, pred, use = "complete.obs")^2,
    R2_pseudo = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    RMSE      = rmse_val,
    MAE       = mae_val,
    Bias      = bias_val,
    RMSE_pct  = 100 * rmse_val / mean(obs, na.rm = TRUE)
  )
}

# ---- 3.2 Pooled metrics (yardstick version if needed) ----
pooled_metrics <- function(df, truth = obs_cm, estimate = pred_cm) {
  metrics <- metric_set(rmse, rsq_trad, mae)
  out <- metrics(df, truth = {{ truth }}, estimate = {{ estimate }})

  truth_vec <- df[[deparse(substitute(truth))]]
  est_vec   <- df[[deparse(substitute(estimate))]]

  pct_rmse <- 100 * rmse_vec(truth_vec, est_vec) /
    mean(truth_vec, na.rm = TRUE)

  list(tbl = out, pct_rmse = pct_rmse)
}

# ---- 3.3 RF: collect out-of-sample predictions ----
collect_rf_oos <- function(model, data, label, flight) {
  best <- model$bestTune

  model$pred %>%
    semi_join(best, by = names(best)) %>%
    transmute(
      id      = rowIndex,
      obs_cm  = exp(obs),
      pred_cm = exp(pred),
      Species = data$Species[rowIndex],
      Model   = label,
      Flight  = flight
    ) %>%
    group_by(id, Species, Model, Flight) %>%
    summarise(
      obs_cm  = first(obs_cm),
      pred_cm = mean(pred_cm, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- 3.4 Log–log models with k-fold CV (no repeats) ----
collect_loglog_oos <- function(data, formula, label, flight, k = CV_K) {
  all_preds <- vector("list", k)
  fold_ids  <- createFolds(data$dbh_cm, k = k)

  for (i in seq_along(fold_ids)) {
    test_idx <- fold_ids[[i]]
    train    <- data[-test_idx, ]
    test     <- data[ test_idx, ]

    if (nrow(test) < 2) next

    fit <- lm(formula, data = train)

    all_preds[[i]] <- tibble(
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(fit, newdata = test)),
      Species = test$Species,
      Model   = label,
      Flight  = flight
    )
  }

  bind_rows(all_preds)
}

# ---- 3.5 Mixed-effects CV (k-fold, no repeats) ----
evaluate_lmer_cv_pooled <- function(data, k = CV_K, label, flight) {

  # For tiny-N subsets, lower k
  if (nrow(data) < 15 && k > 3) k <- 3

  folds     <- createFolds(data$dbh_cm, k = k, list = TRUE)
  all_preds <- vector("list", k)

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train    <- data[-test_idx, ]
    test     <- data[ test_idx, ]

    fit <- lmer(
      log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) + (1 | Species),
      data = train
    )

    all_preds[[i]] <- tibble(
      id      = test_idx,
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(fit, newdata = test, allow.new.levels = TRUE)),
      Species = test$Species,
      Model   = label,
      Flight  = flight
    )
  }

  oos <- bind_rows(all_preds) %>%
    group_by(id, Species, Model, Flight) %>%
    summarise(
      obs_cm  = first(obs_cm),
      pred_cm = mean(pred_cm, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    preds   = oos,
    metrics = compute_metrics(oos)
  )
}

# ---- 3.6 Formatter for mean ± SD ----
fmt_pm <- function(mean, sd, digits = 3) {
  ifelse(
    !is.na(mean),
    sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean, sd),
    NA
  )
}

# ---- 3.7 Bootstrap CIs for bias / RMSE (used in residual plots) ----
boot_ci <- function(residuals, stat = c("bias", "rmse"), R = 2000, seed = 42) {
  stat <- match.arg(stat)
  set.seed(seed)

  f <- switch(
    stat,
    bias = function(x, i) mean(x[i], na.rm = TRUE),
    rmse = function(x, i) sqrt(mean(x[i]^2, na.rm = TRUE))
  )

  b  <- boot(residuals, statistic = f, R = R)
  ci <- boot.ci(b, type = "perc")$percent[4:5]
  ci
}

# ============================================================
# 4) Main loop across flights
# ============================================================
all_preds_all <- list()

for (flight in names(files)) {
  message("Processing: ", flight)
  f <- files[[flight]]

  # ---- 4.1 Load LAS & crowns ----
  las    <- readLAS(f$las)
  crowns <- st_read(f$crowns, quiet = TRUE) %>%
    mutate(treeID = as.character(treeID))

  # ---- 4.2 Load & clean field data ----
  field <- readr::read_csv(f$field, show_col_types = FALSE) %>%
    rename(treeID = !!sym(id_map[[flight]])) %>%
    mutate(
      treeID = as.character(str_trim(treeID)),
      Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE),
                       "Sequoia", "Other"),
      is_sequoia = as.integer(Species == "Sequoia")
    ) %>%
    rename_with(~ "lidar_height",
                .cols = matches("SFMcomb_ht|Side3x_ht|SFMonly_ht")) %>%
    dplyr::select(treeID, dbh_cm, lidar_height, Species, is_sequoia)

  # ---- 4.3 Normalize LAS & compute terrain ----
  las_norm <- normalize_height(classify_ground(las, csf()), knnidw())
  slope    <- terrain(dtm, v = "slope", unit = "degrees")

  # ---- 4.4 Crown metrics (LiDAR + NDVI + area) ----
  metrics_list <- map_dfr(seq_len(nrow(crowns)), function(i) {
    crown_poly <- crowns[i, ]
    las_clip   <- clip_roi(las_norm, crown_poly)

    if (is.null(las_clip) || npoints(las_clip) == 0) {
      return(tibble(
        treeID        = crown_poly$treeID,
        crown_area_m2 = as.numeric(st_area(crown_poly))
      ))
    }

    Z <- las_clip@data$Z

    tibble(
      treeID        = crown_poly$treeID,
      point_density = npoints(las_clip) / as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope, crown_poly, "mean"),
      max_z         = max(Z, na.rm = TRUE),
      mean_z        = mean(Z, na.rm = TRUE),
      p25           = quantile(Z, 0.25, na.rm = TRUE),
      p50           = quantile(Z, 0.50, na.rm = TRUE),
      p75           = quantile(Z, 0.75, na.rm = TRUE),
      p95           = quantile(Z, 0.95, na.rm = TRUE),
      sd_z          = sd(Z, na.rm = TRUE),
      cv_z          = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE),
      ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
      crown_area_m2 = as.numeric(st_area(crown_poly))
    )
  })

  crowns_metrics <- crowns %>%
    dplyr::select(treeID, geom) %>%
    left_join(metrics_list, by = "treeID") %>%
    dplyr::select(
      treeID, crown_area_m2, point_density, slope_mean, ndvi_mean,
      max_z, mean_z, p25, p50, p75, p95, sd_z, cv_z, geom
    )

  # ---- 4.5 Join field + crown metrics ----
  crowns_final <- crowns_metrics %>%
    inner_join(field, by = "treeID") %>%
    filter(Species != "intact snag", dbh_cm > 0)

  stopifnot("dbh_cm" %in% names(crowns_final))

  crown_data <- crowns_final %>%
    filter(!is.na(dbh_cm), dbh_cm > 0, dbh_cm < 350) %>%
    mutate(
      log_dbh        = log(dbh_cm),
      log_height     = log(lidar_height),
      log_crown_area = log1p(crown_area_m2)
    ) %>%
    drop_na()

  # ========================================================
  # 4.6 Models (all with pooled OOS predictions)
  # ========================================================

  # ---- 4.6.1 RF (pooled) ----
  rf_base <- train(
    log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
      p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
      point_density + slope_mean + ndvi_mean,
    data        = crown_data,
    method      = "rf",
    trControl   = ctrl,
    tuneLength  = 5,
    ntree       = N_TREES
  )
  oos_rf_base <- collect_rf_oos(rf_base, crown_data, "RF (pooled)", flight)

  # ---- 4.6.2 Mixed-effects (pooled) ----
  mix_cv  <- evaluate_lmer_cv_pooled(
    crown_data, k = CV_K,
    label = "Mixed-effects (pooled)", flight = flight
  )
  oos_mix <- mix_cv$preds

  # ---- 4.6.3 Species-specific RF ----
  rf_seq <- train(
    log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
      p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
      point_density + slope_mean + ndvi_mean,
    data        = filter(crown_data, is_sequoia == 1),
    method      = "rf",
    trControl   = ctrl,
    tuneLength  = 5,
    ntree       = N_TREES
  )
  oos_rf_seq <- collect_rf_oos(
    rf_seq, filter(crown_data, is_sequoia == 1),
    "RF (Sequoia-only)", flight
  )

  rf_nonseq <- train(
    log(dbh_cm) ~ log(lidar_height) + log1p(crown_area_m2) +
      p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
      point_density + slope_mean + ndvi_mean,
    data        = filter(crown_data, is_sequoia == 0),
    method      = "rf",
    trControl   = ctrl,
    tuneLength  = 5,
    ntree       = N_TREES
  )
  oos_rf_nonseq <- collect_rf_oos(
    rf_nonseq, filter(crown_data, is_sequoia == 0),
    "RF (Non-sequoia)", flight
  )

  # ---- 4.6.4 Log–log allometry ----
  oos_loglog_seq <- collect_loglog_oos(
    filter(crown_data, is_sequoia == 1),
    log(dbh_cm) ~ log(lidar_height) +
      log(pmax(crown_area_m2, 0.01)),
    "Log–log (Sequoia-only)",
    flight
  )

  oos_loglog_nonseq <- collect_loglog_oos(
    filter(crown_data, is_sequoia == 0),
    log(dbh_cm) ~ log(lidar_height) +
      log(pmax(crown_area_m2, 0.01)),
    "Log–log (Non-sequoia)",
    flight
  )

  # ---- 4.6.5 Combine predictions for this flight ----
  preds_store <- bind_rows(
    oos_rf_base,
    oos_mix,
    oos_rf_seq,
    oos_rf_nonseq,
    oos_loglog_seq,
    oos_loglog_nonseq
  )

  all_preds_all[[flight]] <- preds_store
}

# ============================================================
# 5) Combine predictions + metrics
# ============================================================
preds_all <- bind_rows(all_preds_all) %>%
  filter(is.finite(obs_cm), is.finite(pred_cm)) %>%
  mutate(
    residual_cm = pred_cm - obs_cm,
    Model = factor(
      Model,
      levels = c(
        "Mixed-effects (pooled)",
        "RF (pooled)",
        "RF (Sequoia-only)",
        "Log–log (Sequoia-only)",
        "RF (Non-sequoia)",
        "Log–log (Non-sequoia)"
      )
    )
  )

# Collapse any remaining duplicates (e.g., from multiple CV splits)
preds_all_collapsed <- preds_all %>%
  group_by(Flight, Model, Species, obs_cm) %>%
  summarise(
    pred_cm     = mean(pred_cm,     na.rm = TRUE),
    residual_cm = mean(residual_cm, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 5.1 Metrics per Flight × Model ----
metrics_labels <- preds_all_collapsed %>%
  group_by(Flight, Model) %>%
  group_modify(~ compute_metrics(.x)) %>%
  ungroup() %>%
  mutate(
    label = glue(
      "R^2 = {round(R2_cor, 2)}",
      "\nRMSE = {round(RMSE, 1)} cm ({round(RMSE_pct, 1)}%)",
      "\nMAE = {round(MAE, 1)}",
      "\nBias = {round(Bias, 1)}"
    )
  )

# ============================================================
# 6) Figures

# ---- 6.1 Predicted vs Observed (Fusion) ----
metrics_labels <- preds_all_collapsed %>%
  filter(Flight == "Fusion") %>%
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
  preds_all_collapsed %>% filter(Flight == "Fusion"),
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
    title = "Predicted vs. Observed DBH (Fusion)",
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
  "Predicted_vs_Observed_DBH_Fusion.png",
  fig_species_scatter,
  width = 6.7,
  height = 6,
  dpi = 600
)



# ---- 6.2 Residuals vs Predicted (Fusion, key models) ----
resid_labels <- preds_all_collapsed %>%
  filter(Flight == "Fusion") %>%
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
  preds_all_collapsed %>% filter(Flight == "Fusion"),
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
    title = "Residuals vs. Predicted DBH (Fusion)"
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
  "Residuals_DBH_Fusion.png",
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



library(ggplot2)
library(scico)

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
    caption = "Grey polygons = No crown metrics / no prediction",
    x = NULL,
    y = NULL
  )



----------
  





# --- Calculate slope & aspect (in radians, required for shade()) ---
slope  <- terrain(dtm_full, v = "slope",  unit = "radians")
aspect <- terrain(dtm_full, v = "aspect", unit = "radians")

# --- Hillshade (illumination) ---
hill <- shade(slope, aspect)

# --- Convert to data frames for ggplot ---
dtm_df  <- as.data.frame(dtm_full, xy = TRUE)
hill_df <- as.data.frame(hill,      xy = TRUE)


# Terrain color palette (forest style)
terrain_cols <- scico::scico(50, palette = "broc")  # earthy greens/browns


ggplot() +
  # ---- Hillshade on alpha only ----
  geom_raster(
    data = hill_df,
    aes(x = x, y = y, alpha = hillshade),
    fill = "grey20"
  ) +
  scale_alpha(range = c(0.3, 0.9), guide = "none") +

  # ---- Crowns with DBH colors ----
  geom_sf(
    data = crowns_full,
    aes(fill = pred_DBH_cm),
    color = "white",
    size = 0.08
  ) +
  scale_fill_scico(
    palette = "vik",
    name = "Predicted DBH (cm)"
  ) +

  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "Predicted Diameter at Breast Height (DBH)",
    subtitle = "Cloud2Trees crowns + Random Forest model + Hillshade"
  )
