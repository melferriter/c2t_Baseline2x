# ============================
# Libraries
# ============================
library(caret)
library(dplyr)
library(tibble)
library(lmerTest)
library(purrr)
library(yardstick)
library(kableExtra)

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
# Helper: summarize metrics
# ============================
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
# 1. Baseline RF
# ============================
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

# ============================
# 2. RF + Crown Volume Proxy
# ============================
crown_data <- crown_data %>%
  mutate(crown_volume_proxy = crown_area_m2 * treeheight)

rf_crownvol <- train(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) + log1p(crown_volume_proxy) +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + treecbh + ndvi_mean,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  ntree = 1000
)

rf_crownvol_results <- summarize_rf(rf_crownvol, "RF + Crown Volume")

# ============================
# 3. Weighted RF
# ============================
rf_weighted <- train(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) +
    p95 + p75 + p50 + mean_z + max_z + sd_z + cv_z +
    point_density + slope_mean + treecbh + ndvi_mean,
  data = crown_data,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  ntree = 1000,
  weights = crown_data$dbh_cm / mean(crown_data$dbh_cm, na.rm = TRUE)
)

rf_weighted_results <- summarize_rf(rf_weighted, "RF Weighted")

# ============================
# 4. Mixed-effects regression (Species random effect)
# ============================
mix_lmer <- lmer(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2) + (1 | Species),
  data = crown_data
)

# pooled predictions
mix_preds <- tibble(
  obs_cm  = crown_data$dbh_cm,
  pred_cm = exp(predict(mix_lmer, newdata = crown_data))
)

mix_results_pooled <- summarize_preds(mix_preds, "Mixed-effects regression (all data pooled)")

# CV function for lmer
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

----------
  

# ============================
# 5. Combine Results
# ============================
all_results <- bind_rows(
  rf_base_results$pooled, rf_base_results$folds,
  rf_crownvol_results$pooled, rf_crownvol_results$folds,
  rf_weighted_results$pooled, rf_weighted_results$folds,
  mix_results_pooled, mix_results_cv
)

# ============================
# 6. Format Results
# ============================
pooled_results <- all_results %>%
  select(Source, R2, RMSE, MAE, Bias) %>%
  filter(!is.na(R2))

fold_results <- all_results %>%
  select(Source, R2_mean, R2_sd, RMSE_mean, RMSE_sd,
         MAE_mean, MAE_sd, Bias_mean, Bias_sd) %>%
  filter(!is.na(R2_mean))

# Pretty tables
pooled_table <- pooled_results %>%
  kbl(caption = "📊 Pooled Metrics (all folds combined)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

fold_table <- fold_results %>%
  kbl(caption = "📊 Cross-Validation Metrics (average ± SD across folds)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

# Print
pooled_table
fold_table

# ============================
# 7. Species-specific RF
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
      mutate(obs_cm = exp(obs), pred_cm = exp(pred))
    
    summarize_preds(preds, paste(unique(df$Species), "RF"))
  })

species_table <- species_results %>%
  kbl(caption = "📊 Species-specific RF Performance") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

species_table

# ============================
# 8. Sequoia Parametric Allometry
# ============================

# ============================
# Sequoia Parametric Allometry (expanded metrics)
# ============================

seq_data <- crown_data %>% filter(Species == "Sequoia")

seq_lm <- lm(
  log(dbh_cm) ~ log(treeheight) + log1p(crown_area_m2),
  data = seq_data
)

seq_preds <- tibble(
  obs_cm  = seq_data$dbh_cm,
  pred_cm = exp(predict(seq_lm, newdata = seq_data))
)

seq_results <- seq_preds %>%
  summarise(
    R2   = rsq_vec(obs_cm, pred_cm),
    RMSE = rmse_vec(obs_cm, pred_cm),
    MAE  = mae_vec(obs_cm, pred_cm),
    Bias = mean(pred_cm - obs_cm, na.rm = TRUE)
  ) %>%
  mutate(Source = "Sequoia parametric allometry (log–log)")

# Add to the results table
all_results <- bind_rows(all_results, seq_results)

# Clean results table
all_results_clean <- all_results %>%
  select(Source, R2, R2_mean, R2_sd, RMSE, RMSE_mean, RMSE_sd,
         MAE, MAE_mean, MAE_sd, Bias, Bias_mean, Bias_sd)

# Pretty table
all_results_clean %>%
  kbl(caption = "📊 Model Performance Comparison (RF variants, Mixed-effects, Sequoia allometry)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))




-----------
# ============================
# Species-specific log–log allometry
# ============================
library(dplyr)
library(yardstick)

# Define crown depth = height - CBH
crown_data <- crown_data %>%
  mutate(crown_depth = treeheight - treecbh,
         crown_depth = ifelse(crown_depth <= 0, NA, crown_depth))

# --- Fit models ---
lm_seq <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2) + log(crown_depth),
             data = filter(crown_data, is_sequoia == 1, !is.na(crown_depth)))

lm_nonseq <- lm(log(dbh_cm) ~ log(treeheight) + log(crown_area_m2),
                data = filter(crown_data, is_sequoia == 0))

# --- Predictions ---
seq_preds <- filter(crown_data, is_sequoia == 1, !is.na(crown_depth)) %>%
  mutate(obs_cm = dbh_cm,
         pred_cm = exp(predict(lm_seq, newdata = .)),
         resid   = pred_cm - obs_cm)

nonseq_preds <- filter(crown_data, is_sequoia == 0) %>%
  mutate(obs_cm = dbh_cm,
         pred_cm = exp(predict(lm_nonseq, newdata = .)),
         resid   = pred_cm - obs_cm)

# --- Metrics helper ---
get_metrics <- function(df, label) {
  tibble(
    Model = label,
    R2    = yardstick::rsq_vec(df$obs_cm, df$pred_cm),
    RMSE  = yardstick::rmse_vec(df$obs_cm, df$pred_cm),
    MAE   = yardstick::mae_vec(df$obs_cm, df$pred_cm),
    Bias  = mean(df$pred_cm - df$obs_cm, na.rm = TRUE)
  )
}

# --- Evaluate ---
seq_results    <- get_metrics(seq_preds, "Sequoia Log–log allometry (Ht + Area + Depth)")
nonseq_results <- get_metrics(nonseq_preds, "Non-sequoia Log–log allometry (Ht + Area)")

allom_results <- bind_rows(seq_results, nonseq_results)

print(allom_results)

# ============================
# Residual Plots
# ============================
library(ggplot2)

ggplot(bind_rows(
  mutate(seq_preds, Species = "Sequoia"),
  mutate(nonseq_preds, Species = "Non-sequoia")
), aes(x = obs_cm, y = resid, color = Species)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  labs(x = "Observed DBH (cm)", y = "Residuals (Pred – Obs)",
       title = "Residuals vs Tree Size: Log–log Allometry") +
  theme_minimal()

