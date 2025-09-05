# ============================
# 0. Libraries
# ============================
library(lidR)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(caret)
library(dplyr)
library(tibble)
library(lmerTest)
library(purrr)
library(yardstick)
library(kableExtra)
library(ggplot2)
rm(list = ls(globalenv()))
# ============================
# 1. Read Data
# ============================

# Fusion
# las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/Combo/20250825_SFM_3xSide.las")
# crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMcombo_crowns.gpkg")
# fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmcombo_1.csv")

# # LiDAR
# las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las")
# crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/3x_Side/3x_Side_crowns.gpkg")
# fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_side3x.csv")

# SfM
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las")
crowns <- st_read("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM/SFMonly_crowns.gpkg")
fieldpoints <- read_csv("C:/Users/User/Desktop/RandomForest3/fieldpoints_sfmonly_1.csv")

dtm <- rast("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/DEMs/3xSide_20241022190851_DTM.tif")
ndvi <- rast("C:/Users/User/Desktop/RandomForest2/NDVI.tif")


dplyr::glimpse(crowns)
# ============================
# 2. Normalize LAS
# ============================
las_ground <- classify_ground(las, csf())
las_norm <- normalize_height(las_ground, knnidw())

# ============================
# 3. Ancillary rasters
# ============================
slope <- terrain(raster(dtm), opt = "slope", unit = "degrees")

# ============================
# 4. Crown-level metrics
# ============================
metrics_list <- map_dfr(1:nrow(crowns), function(i) {
  crown_poly <- crowns[i, ]
  las_clip <- clip_roi(las_norm, crown_poly)

  if (is.null(las_clip) || npoints(las_clip) == 0) {
    return(tibble(treeID = crown_poly$treeID))
  }

  Z <- las_clip@data$Z
  area_m2 <- as.numeric(st_area(crown_poly))

  tibble(
    treeID = crown_poly$treeID,
    point_density = npoints(las_clip) / area_m2,
    slope_mean = exact_extract(slope, crown_poly, 'mean'),
    max_z = max(Z, na.rm = TRUE),
    mean_z = mean(Z, na.rm = TRUE),
    p25 = quantile(Z, 0.25, na.rm = TRUE),
    p50 = quantile(Z, 0.50, na.rm = TRUE),
    p75 = quantile(Z, 0.75, na.rm = TRUE),
    p95 = quantile(Z, 0.95, na.rm = TRUE),
    sd_z = sd(Z, na.rm = TRUE),
    cv_z = sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE),
    ndvi_mean = exact_extract(ndvi, crown_poly, 'mean'),
    crown_area_m2 = area_m2
  )
})

crowns_metrics <- crowns %>%
  dplyr::select(treeID, geom) %>%   # keep only treeID + geometry
  left_join(metrics_list, by = "treeID")

crowns_metrics <- crowns_metrics %>%
  dplyr::select(treeID, crown_area_m2, point_density, slope_mean, ndvi_mean,
         max_z, mean_z, p25, p50, p75, p95, sd_z, cv_z, geom)


# ============================
# 5. Field data cleanup + join
# ============================
# Clean SFM column and rename
# --- Load and clean field data ---
fieldpoints_clean <- fieldpoints %>%
  mutate(
    treeID = SFMonly %>%
      str_trim() %>%                 # remove leading/trailing spaces
      str_replace_all("[\r\n]", ""), # strip carriage returns/newlines
    Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE),
                     "Sequoia", "Other"),
    is_sequoia = ifelse(Species == "Sequoia", 1, 0)
  ) %>%
  # Keep only the *essential* field data columns
  dplyr::select(treeID, dbh_cm, treecbh, treeheight, Species, is_sequoia)

# --- Join crowns (metrics) with field data ---
crowns_final <- crowns_metrics %>%
  inner_join(fieldpoints_clean, by = "treeID")

crowns_final <- crowns_final %>%
  dplyr::select(
    treeID, dbh_cm, treecbh, treeheight, Species, is_sequoia,
    crown_area_m2, point_density, slope_mean, ndvi_mean,
    max_z, mean_z, p25, p50, p75, p95, sd_z, cv_z, geom
  )

crowns_final_clean <- crowns_final %>%
  filter(Species != "intact snag")

# ============================
# 6. Prepare model data
# ============================

crown_data <- crowns_final_clean %>%
  filter(!is.na(dbh_cm), dbh_cm > 0, dbh_cm < 350) %>%
  mutate(
    log_dbh        = log(dbh_cm),
    log_height     = log(treeheight),
    log_crown_area = log1p(crown_area_m2)
  )

# ============================
# Cross-validation setup
# ============================
set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"
)

# ============================
# Helper functions
# ============================

# Summarize RF results
summarize_rf <- function(model, label) {
  preds <- model$pred %>%
    mutate(obs_cm = exp(obs), pred_cm = exp(pred))
  
  pooled <- preds %>%
    summarise(
      R2   = rsq_vec(obs_cm, pred_cm),
      RMSE = rmse_vec(obs_cm, pred_cm),
      MAE  = mae_vec(obs_cm, pred_cm),
      Bias = mean(pred_cm - obs_cm, na.rm = TRUE)
    ) %>%
    mutate(Source = paste(label, "Overall (all folds pooled)"))
  
  folds <- preds %>%
    group_by(Resample) %>%
    summarise(
      R2   = rsq_vec(obs_cm, pred_cm),
      RMSE = rmse_vec(obs_cm, pred_cm),
      MAE  = mae_vec(obs_cm, pred_cm),
      Bias = mean(pred_cm - obs_cm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    ) %>%
    mutate(Source = paste(label, "Average across folds (±SD)"))
  
  list(pooled = pooled, folds = folds)
}

# Summarize predictions (for LM or lmer)
summarize_preds <- function(preds, label) {
  preds %>%
    summarise(
      R2   = rsq_vec(obs_cm, pred_cm),
      RMSE = rmse_vec(obs_cm, pred_cm),
      MAE  = mae_vec(obs_cm, pred_cm),
      Bias = mean(pred_cm - obs_cm, na.rm = TRUE)
    ) %>%
    mutate(Source = label)
}

# ============================
# 1. Random Forest Variants
# ============================
set.seed(123)
# Baseline RF
rf_base <- train(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + treecbh + ndvi_mean,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  ntree = 1000
)
rf_base_results <- summarize_rf(rf_base, "RF baseline")

# RF + Crown Volume Proxy
# crown_data <- crown_data %>%
#   mutate(crown_volume_proxy = crown_area_m2 * treeheight)

# rf_crownvol <- train(
#   log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) + log1p(crown_volume_proxy) +
#     p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
#     point_density + slope_mean + treecbh + ndvi_mean,
#   data = crown_data,
#   method = "rf",
#   trControl = ctrl,
#   tuneLength = 5,
#   ntree = 1000
# )
# rf_crownvol_results <- summarize_rf(rf_crownvol, "RF + Crown Volume")

# ============================
# 2. Mixed-effects regression
# ============================
# Install if you haven't already
install.packages("lme4")

# Load the package
library(lme4)
mix_lmer <- lmer(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) + (1 | Species),
  data = crown_data
)

mix_preds <- tibble(
  obs_cm  = crown_data$dbh_cm,
  pred_cm = exp(predict(mix_lmer, newdata = crown_data))
)
mix_results_pooled <- summarize_preds(mix_preds, "Mixed-effects regression (all data pooled)")

# CV for lmer
evaluate_lmer_cv <- function(data, k = 5) {
  folds <- createFolds(data$dbh_cm, k = k, list = TRUE, returnTrain = FALSE)
  fold_metrics <- map_dfr(seq_along(folds), function(i) {
    test_idx <- folds[[i]]
    train <- data[-test_idx, ]
    test  <- data[test_idx, ]
    model <- lmer(
      log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) + (1 | Species),
      data = train
    )
    preds <- tibble(
      obs_cm  = test$dbh_cm,
      pred_cm = exp(predict(model, newdata = test, allow.new.levels = TRUE))
    )
    tibble(
      R2    = rsq_vec(preds$obs_cm, preds$pred_cm),
      RMSE  = rmse_vec(preds$obs_cm, preds$pred_cm),
      MAE   = mae_vec(preds$obs_cm, preds$pred_cm),
      Bias  = mean(preds$pred_cm - preds$obs_cm, na.rm = TRUE)
    )
  })
  fold_metrics %>%
    summarise(
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias)
    ) %>%
    mutate(Source = "Mixed-effects regression (CV folds ±SD)")
}
mix_results_cv <- evaluate_lmer_cv(crown_data, k = 5)

# ============================
# 3. Species-specific RF
# ============================
species_results <- crown_data %>%
  group_split(Species) %>%
  map_dfr(function(df) {
    rf_model <- train(
      log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) +
        p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
        point_density + slope_mean + treecbh + ndvi_mean,
      data = df,
      method = "rf",
      trControl = ctrl,
      tuneLength = 5,
      ntree = 1000
    )
    
    preds <- rf_model$pred %>%
      mutate(obs_cm = exp(obs), pred_cm = exp(pred)) %>%
      group_by(Resample) %>%              # group by CV fold
      summarise(
        R2   = cor(obs_cm, pred_cm)^2,
        RMSE = sqrt(mean((obs_cm - pred_cm)^2)),
        MAE  = mean(abs(obs_cm - pred_cm)),
        Bias = mean(pred_cm - obs_cm),
        .groups = "drop"
      ) %>%
      summarise(across(c(R2, RMSE, MAE, Bias),
                       list(mean = mean, sd = sd))) %>%
      mutate(Source = paste(unique(df$Species), "RF"))
    
    preds
  })

# ============================
# 4. Parametric log–log allometry
# ============================
crown_data <- crown_data %>%
  mutate(crown_depth = treeheight - treecbh,
         crown_depth = ifelse(crown_depth <= 0, NA, crown_depth))

# Sequoia model (Ht + Area + Depth)
seq_data <- filter(crown_data, is_sequoia == 1, !is.na(crown_depth))
lm_seq <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2) + log(crown_depth),
             data = seq_data)
seq_preds <- seq_data %>%
  mutate(obs_cm = dbh_cm,
         pred_cm = exp(predict(lm_seq, newdata = .)))
seq_results <- summarize_preds(seq_preds, "Sequoia Log–log allometry (Ht + Area + Depth)")

# Non-sequoia model (Ht + Area only)
nonseq_data <- filter(crown_data, is_sequoia == 0)
lm_nonseq <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2),
                data = nonseq_data)
nonseq_preds <- nonseq_data %>%
  mutate(obs_cm = dbh_cm,
         pred_cm = exp(predict(lm_nonseq, newdata = .)))
nonseq_results <- summarize_preds(nonseq_preds, "Non-sequoia Log–log allometry (Ht + Area)")

# ============================
# 5. Combine all results
# ============================
library(kableExtra)
all_results <- bind_rows(
  rf_base_results$pooled, rf_base_results$folds,
  mix_results_pooled, mix_results_cv,
  species_results,
  seq_results, nonseq_results
)

# Clean table
all_results_clean <- all_results %>%
  dplyr::select(Source, R2, R2_mean, R2_sd,
         RMSE, RMSE_mean, RMSE_sd,
         MAE, MAE_mean, MAE_sd,
         Bias, Bias_mean, Bias_sd)

# Pretty summary
all_results_clean %>%
  kbl(caption = "📊 Model Performance Comparison (RF variants, Mixed-effects, Sequoia & Non-sequoia allometry)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))

# ============================
# 6. Residual Diagnostics
# ============================
ggplot(bind_rows(
  mutate(seq_preds, Species = "Sequoia"),
  mutate(nonseq_preds, Species = "Non-sequoia")
), aes(x = obs_cm, y = pred_cm - obs_cm, color = Species)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  labs(x = "Observed DBH (cm)", y = "Residuals (Pred – Obs)",
       title = "Residuals vs Tree Size: Log–log Allometry") +
  theme_minimal()


--------
library(caret)
library(dplyr)

set.seed(123)
folds <- createFolds(seq_data$dbh_cm, k = 5)

seq_cv_results <- lapply(folds, function(fold) {
  train_data <- seq_data[-fold, ]
  test_data  <- seq_data[fold, ]
  
  mod <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2) + log(crown_depth),
            data = train_data)
  
  preds <- exp(predict(mod, newdata = test_data))
  obs   <- test_data$dbh_cm
  
  tibble(
    R2   = cor(obs, preds)^2,
    RMSE = sqrt(mean((obs - preds)^2)),
    MAE  = mean(abs(obs - preds)),
    Bias = mean(preds - obs)
  )
})

# Summarize mean ± SD
seq_cv_summary <- bind_rows(seq_cv_results) %>%
  summarise(across(everything(),
                   list(mean = mean, sd = sd)))

seq_cv_summary


set.seed(123)
folds <- createFolds(nonseq_data$dbh_cm, k = 5)

nonseq_cv_results <- lapply(folds, function(fold) {
  train_data <- nonseq_data[-fold, ]
  test_data  <- nonseq_data[fold, ]
  
  # log–log allometry for non-sequoias (no crown_depth term)
  mod <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2),
            data = train_data)
  
  preds <- exp(predict(mod, newdata = test_data))
  obs   <- test_data$dbh_cm
  
  tibble(
    R2   = cor(obs, preds)^2,
    RMSE = sqrt(mean((obs - preds)^2)),
    MAE  = mean(abs(obs - preds)),
    Bias = mean(preds - obs)
  )
})

# Summarize mean ± SD across folds
nonseq_cv_summary <- bind_rows(nonseq_cv_results) %>%
  summarise(across(everything(),
                   list(mean = mean, sd = sd)))

nonseq_cv_summary



-------------
library(ggplot2)
library(patchwork)

library(dplyr)
library(ggplot2)
library(patchwork)

# Step 1: Create log-transformed predictors
nonseq_data <- nonseq_data %>%
  mutate(
    log_height = log(treeheight),
    log_area   = log(crown_area_m2)
  )

# Step 2: Filter out potential sequoia outliers (dbh > 120 cm) and NAs
nonseq_data_clean <- nonseq_data %>%
  filter(!is.na(dbh_cm), dbh_cm <= 120)

# Step 3: Fit log–log allometry on cleaned dataset
lm_nonseq <- lm(log(dbh_cm) ~ log_height + log_area, data = nonseq_data_clean)

# Step 4: Predictions + residuals for log–log allometry
loglog_preds <- exp(predict(lm_nonseq, newdata = nonseq_data_clean))
nonseq_data_clean$loglog_resid <- loglog_preds - nonseq_data_clean$dbh_cm

p1 <- ggplot(nonseq_data_clean, aes(x = dbh_cm, y = loglog_resid, color = Species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs DBH: Log–log Allometry",
       x = "Observed DBH (cm)", y = "Residual (Pred − Obs)") +
  theme_minimal()

# Step 5: Predictions + residuals for Random Forest (same cleaned data)
rf_preds <- exp(predict(rf_nonseq, newdata = nonseq_data_clean))
nonseq_data_clean$rf_resid <- rf_preds - nonseq_data_clean$dbh_cm

p2 <- ggplot(nonseq_data_clean, aes(x = dbh_cm, y = rf_resid, color = Species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs DBH: Random Forest",
       x = "Observed DBH (cm)", y = "Residual (Pred − Obs)") +
  theme_minimal()

# Step 6: Combine side by side
p1 + p2

ggsave(
  filename = "heightvsdbh.png",  # file name
  plot = last_plot(),                          # saves the most recent ggplot
  path = "C:/Users/User/Desktop",      # <-- change to your folder
  width = 7, height = 5,                       # adjust for paper layout
  dpi = 1200                                    # high resolution for publication
)



ggplot(nonseq_data_clean, aes(x = treeheight, y = dbh_cm)) +
  geom_point(size = 2, alpha = 0.7) +                       # slightly larger, semi-transparent points
  geom_smooth(method = "lm", se = TRUE, color = "blue") +   # add SE shading for clarity
  labs(
    title = "Fusion Data: Tree Height vs DBH",
    x = "Tree Height (m)",
    y = "DBH (cm)"
  ) +
  theme_minimal(base_size = 14) 
