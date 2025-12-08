# ============================================================
# 0) SETUP
# ============================================================

rm(list = ls())

library(lidR)
library(terra)
library(sf)
library(stringr)
library(caret)
library(tibble)
library(lmerTest)
library(purrr)
library(ggplot2)
library(readr)
library(exactextractr)
library(glue)
library(tidyr)
library(dplyr)
library(scico)
library(patchwork)
library(ggpp)

set.seed(12345)
options(mc.cores = 1)

# ============================================================
# 1) GLOBAL OPTIONS
# ============================================================

N_TREES <- 1000           # RF tree count
K_POOLED <- 5             # pooled CV folds
K_NONSEQ <- 5             # non-sequoia folds
K_SEQ <- 3      
CV_K <- 3     # default for pooled models
N_TREES <- 1000


THEME_BW <- TRUE
if (THEME_BW) theme_set(theme_bw())

# ============================================================
# File paths for all flights
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

# ============================================================
# 2) METRIC FUNCTIONS
# ============================================================

compute_metrics <- function(df) {
  obs  <- df$obs_cm
  pred <- df$pred_cm
  
  tibble(
    R2    = cor(obs, pred)^2,
    RMSE  = sqrt(mean((pred - obs)^2)),
    MAE   = mean(abs(pred - obs)),
    Bias  = mean(pred - obs)
  )
}

# ============================================================
# 3) UNIVERSAL PARK-APPROVED CV FUNCTION
# ============================================================

run_kfold_oos <- function(data, k, fit_fun, pred_fun, label, flight) {
  
  folds <- createFolds(data$dbh_cm, k = k, list = TRUE)
  
  oos <- map_dfr(seq_along(folds), function(i) {
    
    test_idx  <- folds[[i]]
    train_dat <- data[-test_idx, ]
    test_dat  <- data[test_idx, ]
    
    model_fit <- fit_fun(train_dat)
    preds <- pred_fun(model_fit, test_dat)
    
    tibble(
      id      = test_idx,
      obs_cm  = test_dat$dbh_cm,
      pred_cm = preds,
      Species = test_dat$Species,
      Model   = label,
      Flight  = flight
    )
  })
  
  list(
    preds = oos,
    metrics = compute_metrics(oos) %>%
      mutate(Model = label, Flight = flight)
  )
}

# ============================================================
# 4) MODEL DEFINITIONS
# ============================================================

### RANDOM FOREST (pooled + species-specific)
fit_rf <- function(dat) {
  randomForest::randomForest(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z +
      point_density+slope_mean+ndvi_mean,
    data = dat,
    ntree = N_TREES
  )
}
pred_rf <- function(model, dat) exp(predict(model, dat))

fit_rf_seq     <- fit_rf
fit_rf_nonseq  <- fit_rf


### MIXED-EFFECTS
fit_lmer <- function(dat) {
  lmer(
    log(dbh_cm) ~ log_height + log1p(crown_area_m2) + (1|Species),
    data = dat
  )
}
pred_lmer <- function(model, dat) exp(predict(model, dat, allow.new.levels = TRUE))


### LOG–LOG MODELS
fit_ll_seq <- function(dat) lm(log(dbh_cm) ~ log_height + log(crown_area_m2), dat)
pred_ll_seq <- function(model, dat) exp(predict(model, dat))

fit_ll_nonseq <- function(dat) lm(log(dbh_cm) ~ log_height + log(crown_area_m2), dat)
pred_ll_nonseq <- function(model, dat) exp(predict(model, dat))

# ============================================================
# 5) LOAD DATA FOR EACH FLIGHT
# ============================================================

all_training_data <- list()

for (flight in names(files)) {

  message("\n-------------------------")
  message("Processing flight: ", flight)
  message("-------------------------")

  f <- files[[flight]]

  # --- Load LAS
  las_raw <- readLAS(f$las)
  stopifnot(!is.null(las_raw))

  # --- Load crowns
  crowns <- st_read(f$crowns, quiet = TRUE) %>%
    mutate(treeID = as.character(treeID))

  # --- Load field
  field <- read_csv(f$field, show_col_types = FALSE) %>%
    rename(treeID = !!sym(f$idcol)) %>%
    mutate(
      treeID = as.character(str_trim(treeID)),
      Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE),
                       "Sequoia", "Other"),
      is_sequoia = as.integer(Species == "Sequoia")
    ) %>%
    rename(lidar_height = !!sym(f$height_col)) %>%
    select(treeID, dbh_cm, lidar_height, Species, is_sequoia)

  
  # ============================================================
  # 6) NORMALIZE HEIGHTS
  # ============================================================

  las_norm <- las_raw %>%
    classify_ground(csf()) %>%
    normalize_height(knnidw())

    # ============================================================
  # 7) CROWN METRICS
  # ============================================================

  slope_rast <- terrain(dtm, v = "slope", unit = "degrees")

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
        max_z         = NA, mean_z = NA,
        p25 = NA, p50 = NA, p75 = NA, p95 = NA,
        sd_z = NA, cv_z = NA
      ))
    }

    Z <- las_clip@data$Z

    tibble(
      treeID        = crown_poly$treeID,
      crown_area_m2 = as.numeric(st_area(crown_poly)),
      point_density = npoints(las_clip) / as.numeric(st_area(crown_poly)),
      slope_mean    = exact_extract(slope_rast, crown_poly, "mean"),
      ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
      max_z = max(Z),
      mean_z = mean(Z),
      p25    = quantile(Z, 0.25),
      p50    = quantile(Z, 0.50),
      p75    = quantile(Z, 0.75),
      p95    = quantile(Z, 0.95),
      sd_z   = sd(Z),
      cv_z   = sd(Z) / mean(Z)
    )
  })

  # ============================================================
  # 8) BUILD FINAL MODELING DATASET
  # ============================================================

  crowns_final <- crowns %>%
    select(treeID, geom) %>%
    left_join(metrics_list, by = "treeID") %>%
    inner_join(field, by = "treeID") %>%
    filter(dbh_cm > 0, dbh_cm < 350)   # remove invalid or huge outliers

  crown_data <- crowns_final %>%
    mutate(
      log_height     = log(pmax(lidar_height, 0.5)),
      log_crown_area = log1p(crown_area_m2),
      log_dbh        = log(dbh_cm),
      Flight = flight
    ) %>%
    drop_na()

  all_training_data[[flight]] <- crown_data
}

# ============================================================
# 9) CROSS-VALIDATION CONTROL
# ============================================================

ctrl_default <- trainControl(
  method          = "cv",
  number          = CV_K,         # normally 3 folds
  savePredictions = "final",
  returnResamp    = "final"
)


# ============================================================
# 10) HELPER FUNCTIONS
# ============================================================

compute_metrics <- function(df) {
  obs  <- df$obs_cm
  pred <- df$pred_cm

  tibble(
    R2_cor   = cor(obs, pred, use="complete.obs")^2,
    RMSE     = sqrt(mean((obs - pred)^2, na.rm=TRUE)),
    RMSE_pct = 100 * RMSE / mean(obs, na.rm=TRUE),
    MAE      = mean(abs(pred - obs), na.rm=TRUE),
    Bias     = mean(pred - obs, na.rm=TRUE)
  )
}

collect_rf_oos <- function(model, data, label, flight) {

  best <- model$bestTune

  model$pred %>%
    semi_join(best, by = names(best)) %>%     # keep rows for best mtry
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
      pred_cm = mean(pred_cm),
      .groups = "drop"
    )
}

collect_loglog_oos <- function(data, formula, label, flight, k = CV_K) {

  if (nrow(data) < 20) k <- 2
  set.seed(12345)

  folds <- createFolds(data$dbh_cm, k = k, list = TRUE)

  oos <- map_dfr(folds, function(test_idx) {

    fit  <- lm(formula, data = data[-test_idx, ])
    pred <- exp(predict(fit, newdata = data[test_idx, ]))

    tibble(
      obs_cm  = data$dbh_cm[test_idx],
      pred_cm = pred,
      Species = data$Species[test_idx],
      Model   = label,
      Flight  = flight
    )
  })

  oos
}

evaluate_lmer_cv_pooled <- function(data, label, flight, k = CV_K) {

  if (nrow(data) < 20) k <- 2
  set.seed(12345)

  folds <- createFolds(data$dbh_cm, k = k, list=TRUE)

  oos <- map_dfr(folds, function(test_idx) {

    fit <- lmer(
      log(dbh_cm) ~ log_height + log1p(crown_area_m2) + (1 | Species),
      data = data[-test_idx, ]
    )

    pred <- exp(predict(fit, newdata = data[test_idx, ], allow.new.levels = TRUE))

    tibble(
      id      = test_idx,
      obs_cm  = data$dbh_cm[test_idx],
      pred_cm = pred,
      Species = data$Species[test_idx],
      Model   = label,
      Flight  = flight
    )
  })

  oos %>%
    group_by(id, Species, Model, Flight) %>%
    summarise(
      obs_cm  = first(obs_cm),
      pred_cm = mean(pred_cm),
      .groups = "drop"
    )
}

# ============================================================
# 11) MAIN TRAINING LOOP — ALL 6 MODELS
# ============================================================

all_preds_all <- list()

for (flight in names(all_training_data)) {

  crown_data <- all_training_data[[flight]]
  message("\nFitting models for flight: ", flight)

  # 1) RF pooled ============================
  rf_pooled <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
      point_density + slope_mean + ndvi_mean,
    data      = crown_data,
    method    = "rf",
    trControl = ctrl_default,
    tuneLength = 5,
    ntree      = N_TREES
  )

  oos_rf_pooled <- collect_rf_oos(rf_pooled, crown_data,
                                  "RF (pooled)", flight)


  # 2) Mixed-effects pooled =================
  oos_mix <- evaluate_lmer_cv_pooled(
    data = crown_data,
    label = "Mixed-effects (pooled)",
    flight = flight
  )


  # 3) RF Sequoia-only ======================
  seq_data  <- filter(crown_data, is_sequoia == 1)
  seq_ctrl  <- if (nrow(seq_data) < 20)
                 trainControl(method="cv", number=2, savePredictions="final")
               else ctrl_default

  rf_seq <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z+
      point_density+slope_mean+ndvi_mean,
    data      = seq_data,
    method    = "rf",
    trControl = seq_ctrl,
    tuneLength = 5,
    ntree      = N_TREES
  )

  oos_rf_seq <- collect_rf_oos(rf_seq, seq_data,
                               "RF (Sequoia-only)", flight)


  # 4) RF Non-sequoia =======================
  nonseq_data <- filter(crown_data, is_sequoia == 0)

  rf_nonseq <- train(
    log(dbh_cm) ~ log_height + log_crown_area +
      p95+p75+p50+mean_z+max_z+sd_z+cv_z+
      point_density+slope_mean+ndvi_mean,
    data      = nonseq_data,
    method    = "rf",
    trControl = ctrl_default,
    tuneLength = 5,
    ntree      = N_TREES
  )

  oos_rf_nonseq <- collect_rf_oos(rf_nonseq, nonseq_data,
                                  "RF (Non-sequoia)", flight)


  # 5) Log–log Sequoia-only ================
  oos_log_seq <- collect_loglog_oos(
    seq_data,
    log(dbh_cm) ~ log_height + log(pmax(crown_area_m2, 0.01)),
    label = "Log–log (Sequoia-only)",
    flight = flight
  )

  # 6) Log–log Non-sequoia =================
  oos_log_nonseq <- collect_loglog_oos(
    nonseq_data,
    log(dbh_cm) ~ log_height + log(pmax(crown_area_m2, 0.01)),
    label = "Log–log (Non-sequoia)",
    flight = flight
  )

  # ---- Store predictions for this flight
  all_preds_all[[flight]] <- bind_rows(
    oos_rf_pooled,
    oos_mix,
    oos_rf_seq,
    oos_rf_nonseq,
    oos_log_seq,
    oos_log_nonseq
  )
}

# ============================================================
# 12) COLLAPSE ALL OUT-OF-SAMPLE PREDICTIONS
# ============================================================

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
    .groups = "drop"
  )

metrics_table <- preds_all_collapsed %>%
  group_by(Flight, Model) %>%
  group_modify(~ compute_metrics(.x)) %>%
  ungroup() %>%
  mutate(
    label = glue(
      "R² = {round(R2_cor, 2)}\n",
      "RMSE = {round(RMSE, 1)} cm ({round(RMSE_pct, 1)}%)\n",
      "MAE = {round(MAE, 1)} cm\n",
      "Bias = {round(Bias, 1)} cm"
    )
  )

library(ggpp)

metrics_labels_sfm <- preds_all_collapsed %>%
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
    npcx = 0.02,
    npcy = 0.98,
    label = glue(
      "R² = {round(r2,2)}\n",
      "RMSE = {round(rmse,1)} cm ({round(rmse_pct,1)}%)\n",
      "MAE = {round(mae,1)} cm\n",
      "Bias = {round(bias,1)} cm"
    )
  )

fig_scatter_sfm <- ggplot(
  preds_all_collapsed %>% filter(Flight == "SfM"),
  aes(x = obs_cm, y = pred_cm, color = Species)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.7, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "gray80", linewidth = 0.8) +
  geom_point(alpha = 0.45, size = 1.8) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ Model, ncol = 2) +
  geom_text_npc(
    data = metrics_labels_sfm,
    aes(npcx = npcx, npcy = npcy, label = label),
    hjust = 0, vjust = 1, size = 3.5,
    color = "black", lineheight = 1.0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("Other"="#E69F00", "Sequoia"="#56B4E9")) +
  labs(
    x = "Observed DBH (cm)",
    y = "Predicted DBH (cm)",
    title = "Predicted vs. Observed DBH — SfM",
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

fig_scatter_sfm

resid_labels <- preds_all_collapsed %>%
  filter(Flight == "LiDAR") %>%
  group_by(Model) %>%
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
      "N = {n}\nBias = {round(bias,1)} cm\n",
      "RMSE = {round(rmse,1)} cm ({round(rmse_pct,1)}%)"
    ),
    npcx = 0.02, npcy = 0.02
  )

