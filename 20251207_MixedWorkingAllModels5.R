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
library(randomForest)

set.seed(1234)
options(mc.cores = 1)

###############################################################################
# 1) GLOBAL OPTIONS
###############################################################################

N_TREES <- 1000

# folds for each model type
K_POOLED  <- 3    # pooled models
K_NONSEQ  <- 3    # non-sequoia
K_SEQ     <- 2    # sequoia-only (small sample)

# plot theme
theme_set(theme_bw())
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

###############################################################################
# 3) METRIC FUNCTION (Unified)
###############################################################################

compute_metrics <- function(df) {
obs  <- df$obs_cm
pred <- df$pred_cm

rmse_val <- sqrt(mean((pred - obs)^2))

tibble(
  R2       = cor(obs, pred)^2,
  RMSE     = rmse_val,
  RMSE_pct = 100 * rmse_val / mean(obs),
  MAE      = mean(abs(pred - obs)),
  Bias     = mean(pred - obs)
)
}

###############################################################################
# 4) UNIVERSAL PARK-APPROVED K-FOLD CV FUNCTION
###############################################################################

run_kfold_oos <- function(data, k, fit_fun, pred_fun, label, flight) {

folds <- createFolds(data$dbh_cm, k = k, list = TRUE)

fold_metrics <- list()

oos <- map_dfr(seq_along(folds), function(i) {
  test_idx  <- folds[[i]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[test_idx, ]
  
  model_fit <- fit_fun(train_dat)
  preds     <- pred_fun(model_fit, test_dat)
  
  # store fold-specific metrics
  fold_metrics[[i]] <<- compute_metrics(
    tibble(obs_cm = test_dat$dbh_cm, pred_cm = preds)
  ) %>% mutate(Fold = i, Model = label, Flight = flight)
  
  tibble(
    Fold    = i,
    id      = test_idx,
    obs_cm  = test_dat$dbh_cm,
    pred_cm = preds,
    Species = test_dat$Species,
    Model   = label,
    Flight  = flight
  )
})

list(
  preds   = oos,
  metrics = bind_rows(fold_metrics)
)
}

###############################################################################
# 5) MODEL DEFINITIONS
###############################################################################

# Random Forest
fit_rf <- function(dat) {
randomForest(
  log(dbh_cm) ~ log_height + log_crown_area +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + ndvi_mean,
  data = dat,
  ntree = N_TREES
)
}
pred_rf <- function(model, dat) exp(predict(model, dat))

# Mixed-effects
fit_lmer <- function(dat) {
lmer(
  log(dbh_cm) ~ log_height + log1p(crown_area_m2) + (1|Species),
  data = dat
)
}
pred_lmer <- function(model, dat) exp(predict(model, dat, allow.new.levels = TRUE))

# Log–log models
fit_ll <- function(dat) lm(log(dbh_cm) ~ log_height + log(crown_area_m2), dat)
pred_ll <- function(model, dat) exp(predict(model, dat))

###############################################################################
# 6) BUILD TRAINING DATA FOR ALL FLIGHTS
###############################################################################

all_training_data <- list()

for (flight in names(files)) {
message("Processing flight: ", flight)

f <- files[[flight]]

las_raw <- readLAS(f$las)
crowns  <- st_read(f$crowns, quiet = TRUE) %>%
  mutate(treeID = as.character(treeID))

field <- read_csv(f$field, show_col_types = FALSE) %>%
  rename(treeID = !!sym(f$idcol)) %>%
  mutate(
    treeID = as.character(str_trim(treeID)),
    Species = ifelse(grepl("Sequoia", Species, ignore.case = TRUE),
                     "Sequoia", "Other"),
    is_sequoia = as.integer(Species == "Sequoia")
  ) %>%
  rename(lidar_height = !!sym(f$height_col))

las_norm <- normalize_height(classify_ground(las_raw, csf()), knnidw())
slope_rast <- terrain(dtm, v="slope", unit="degrees")

# crown metrics
metrics_list <- map_dfr(seq_len(nrow(crowns)), function(i) {
  crown_poly <- crowns[i, ]
  las_clip   <- clip_roi(las_norm, crown_poly)
  
  if (is.null(las_clip) || npoints(las_clip) == 0)
    return(tibble(treeID=crown_poly$treeID, crown_area_m2=as.numeric(st_area(crown_poly)),
                  point_density=NA, slope_mean=NA, ndvi_mean=NA,
                  max_z=NA, mean_z=NA, p25=NA, p50=NA, p75=NA, p95=NA,
                  sd_z=NA, cv_z=NA))
  
  Z <- las_clip@data$Z
  
  tibble(
    treeID        = crown_poly$treeID,
    crown_area_m2 = as.numeric(st_area(crown_poly)),
    point_density = npoints(las_clip)/as.numeric(st_area(crown_poly)),
    slope_mean    = exact_extract(slope_rast, crown_poly, "mean"),
    ndvi_mean     = exact_extract(ndvi, crown_poly, "mean"),
    max_z = max(Z),
    mean_z = mean(Z),
    p25    = quantile(Z, 0.25),
    p50    = quantile(Z, 0.50),
    p75    = quantile(Z, 0.75),
    p95    = quantile(Z, 0.95),
    sd_z   = sd(Z),
    cv_z   = sd(Z)/mean(Z)
  )
})

crown_data <- crowns %>%
  select(treeID, geom) %>%
  left_join(metrics_list, by="treeID") %>%
  inner_join(field, by="treeID") %>%
  filter(dbh_cm > 0, dbh_cm < 350) %>%
  mutate(
    log_height     = log(pmax(lidar_height, 0.5)),
    log_crown_area = log1p(crown_area_m2),
    log_dbh        = log(dbh_cm),
    Flight         = flight
  ) %>%
  drop_na()

all_training_data[[flight]] <- crown_data
}

###############################################################################
# 7) CROSS-VALIDATION FOR ALL MODELS
###############################################################################

all_preds_all   <- list()
all_fold_metrics <- list()

for (flight in names(all_training_data)) {

dat <- all_training_data[[flight]]

rf_pool  <- run_kfold_oos(dat, K_POOLED, fit_rf,  pred_rf, "RF (pooled)", flight)
mix      <- run_kfold_oos(dat, K_POOLED, fit_lmer, pred_lmer, "Mixed-effects (pooled)", flight)
rf_seq   <- run_kfold_oos(filter(dat, is_sequoia==1), K_SEQ, fit_rf,  pred_rf, "RF (Sequoia-only)", flight)
rf_non   <- run_kfold_oos(filter(dat, is_sequoia==0), K_NONSEQ, fit_rf,  pred_rf, "RF (Non-sequoia)", flight)
ll_seq   <- run_kfold_oos(filter(dat, is_sequoia==1), K_SEQ, fit_ll, pred_ll, "Log–log (Sequoia-only)", flight)
ll_non   <- run_kfold_oos(filter(dat, is_sequoia==0), K_NONSEQ, fit_ll, pred_ll, "Log–log (Non-sequoia)", flight)

all_preds_all[[flight]] <- bind_rows(
  rf_pool$preds,
  mix$preds,
  rf_seq$preds,
  rf_non$preds,
  ll_seq$preds,
  ll_non$preds
)

all_fold_metrics[[flight]] <- bind_rows(
  rf_pool$metrics,
  mix$metrics,
  rf_seq$metrics,
  rf_non$metrics,
  ll_seq$metrics,
  ll_non$metrics
)
}

fold_metrics_all <- bind_rows(all_fold_metrics)

###############################################################################
# 8B) COLLAPSE ALL OUT-OF-SAMPLE PREDICTIONS INTO ONE TABLE
###############################################################################

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



###############################################################################
# 8) SPREAD TABLE (Park-required: mean, min, max)
###############################################################################

spread_table <- fold_metrics_all %>%
  group_by(Flight, Model) %>%
  summarise(
    mean_R2  = mean(R2),
    min_R2   = min(R2),
    max_R2   = max(R2),

    mean_RMSE = mean(RMSE),
    min_RMSE  = min(RMSE),
    max_RMSE  = max(RMSE),

    mean_MAE  = mean(MAE),
    min_MAE   = min(MAE),
    max_MAE   = max(MAE),

    mean_Bias = mean(Bias),
    min_Bias  = min(Bias),
    max_Bias  = max(Bias),

    .groups = "drop"
  )

print(spread_table)


###############################################################################
# 9) FINAL PLOT EXAMPLE
###############################################################################

library(ggpp)

# Compute metrics per model
metrics_labels <- preds_all %>%
filter(Flight == "LiDAR") %>%
group_by(Model) %>%
summarise(
  R2       = cor(obs_cm, pred_cm)^2,
  RMSE     = sqrt(mean((pred_cm - obs_cm)^2)),
  MAE      = mean(abs(pred_cm - obs_cm)),
  Bias     = mean(pred_cm - obs_cm),
  RMSE_pct = 100 * RMSE / mean(obs_cm),
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
    "R² = {round(R2,2)}\n",
    "RMSE = {round(RMSE,1)} cm ({round(RMSE_pct,1)}%)\n",
    "MAE = {round(MAE,1)} cm\n",
    "Bias = {round(Bias,1)} cm"
  )
)





fig_species_scatter <- ggplot(
preds_all %>% filter(Flight == "LiDAR"),
aes(x = obs_cm, y = pred_cm, color = Species)
) +
geom_abline(
  slope = 1, intercept = 0,
  linetype = "dashed",
  linewidth = 0.7, color = "gray40"
) +
geom_smooth(
  method = "lm", se = TRUE,
  color = "black", fill = "gray80",
  linewidth = 0.8
) +
geom_point(alpha = 0.45, size = 1.8) +
facet_wrap(~ Model, ncol = 2) +
geom_text_npc(
  data = metrics_labels,
  aes(npcx = npcx, npcy = npcy, label = label),
  hjust = 0, vjust = 1, size = 3.5,
  color = "black", lineheight = 1.0,
  inherit.aes = FALSE
) +
scale_color_manual(values = c("Other"="#E69F00", "Sequoia"="#56B4E9")) +
labs(
  x = "Observed DBH (cm)",
  y = "Predicted DBH (cm)",
  title = "Predicted vs. Observed DBH — LiDAR",
  color = "Species"
) +
theme_bw(base_size = 12) +
theme(
  strip.text = element_text(face = "bold"),
  panel.grid.major = element_line(color = "gray90"),
  panel.grid.minor = element_blank(),
  legend.position = "right"
)

fig_species_scatter





ggsave(
  "20251208_Predicted_vs_Observed_DBH_LiDAR.png",
  fig_species_scatter,
  width = 6.7,
  height = 6,
  dpi = 600
)







resid_labels <- preds_all %>%
filter(Flight == "LiDAR") %>%
group_by(Model) %>%
summarise(
  N        = n(),
  Bias     = mean(residual_cm),
  RMSE     = sqrt(mean(residual_cm^2)),
  MAE      = mean(abs(residual_cm)),
  RMSE_pct = 100 * RMSE / mean(obs_cm),
  .groups = "drop"
) %>%
mutate(
  label = glue(
    "N = {N}\n",
    "Bias = {round(Bias,1)} cm\n",
    "RMSE = {round(RMSE,1)} cm ({round(RMSE_pct,1)}%)\n",
    "MAE = {round(MAE,1)} cm"
  ),
  npcx = 0.02,
  npcy = 0.02
)



fig_residuals <- ggplot(
  preds_all %>% filter(Flight == "LiDAR"),
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
    title = "Residuals vs. Predicted DBH (LiDAR)"
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
  "20251208_Residuals_DBH_LiDAR.png",
  fig_residuals,
  width = 6.7,
  height = 6,
  dpi = 600
)



















write_csv(preds_all_collapsed, "preds_all_collapsed1.csv")
