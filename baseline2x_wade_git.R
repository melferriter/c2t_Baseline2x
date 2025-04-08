## pkgbuild helps us check for Rtools
install.packages("pkgbuild")
# check for Rtools which is required to build packages
pkgbuild::check_build_tools(debug = TRUE)
## remotes helps us get packages hosted on github
install.packages("remotes")
## install lasR from the r-univers
install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
## install TreeLS from github
remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)
## get cloud2trees
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)


# download the external data
cloud2trees::get_data()
# download the TreeMap data
cloud2trees::get_treemap()
# download the forest type data
cloud2trees::get_foresttype()
# download the landfire cbd data
cloud2trees::get_landfire()


# Libraries====
library(lidR)
library(sf)
library(terra)
library(stars)
library(raster)
library(alphashape3d)
library(plyr)
library(tidyverse)
library(devtools)
library(canopyLazR)
library(spdep)
library(sp)
library(geosphere)
library(rlas)
library(rgl)
library(pracma)
library(spatstat)
library(terra)
library(cloud2trees)
library(tidyverse)
library(sf)
library(purrr)
library(patchwork)
library(viridis)
library(dbplyr)


rm(list = ls(globalenv()))

# working directory====
setwd("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/")


----------------------------------------------
  # Extract Trees from Point Cloud: Default====
----------------------------------------------
  
las <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las", filter = "-set_withheld_flag 0")
i <- "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las"

# run cloud2trees without customization
cloud2trees_ans <- cloud2trees::cloud2trees(output_dir = tempdir(), input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las")

# return from cloud2trees package
cloud2trees_ans %>% names()

# DTM
cloud2trees_ans$dtm_rast %>% terra::plot()

# CHM
cloud2trees_ans$chm_rast %>% terra::plot()

# tree crowns (spatial data frame)
cloud2trees_ans$crowns_sf %>% dplyr::glimpse()

# plot tree crown polygons
cloud2trees_ans$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

# tree top points
cloud2trees_ans$treetops_sf %>% dplyr::glimpse()

# plot tree top points
cloud2trees_ans$treetops_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(color = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_color_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans$crowns_sf, mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf(data = cloud2trees_ans$treetops_sf, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")


----------------------------------------------
# Individual Tree Detection (ITD) Tuning====
----------------------------------------------
## Default ITD window size functions====
itd_tuning_ans <- itd_tuning(input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las")

itd_tuning_ans %>% names()

itd_tuning_ans$plot_samples

best_ws <- itd_tuning_ans$ws_fn_list$lin_fn

ggplot2::ggplot() +
  ggplot2::geom_function(fun = best_ws, color = "brown", lwd = 1) +
  ggplot2::xlim(-5,60) +
  ggplot2::labs(x = "heights", y = "ws", color = "") +
  ggplot2::theme_light()

# a constant window size has to be defined as:
## x*0 + constant
my_constant <- function(x){(x * 0) + 3} ## will always return 3
# a custom linear function
my_linear <- function(x) {(x * 0.1) + 3}
# another custom
my_custom2 <- function(x) {0.15 * x^0.6 + 2}
# let's put these in a list to test with the best default function we saved from above
my_fn_list <- list(
  my_constant = my_constant
  , my_linear = my_custom2
  , best_default_ws = best_ws
)

# run it with custom functions
itd_tuning_ans2 <- itd_tuning(
  input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las"
  , ws_fn_list = my_fn_list
  , n_samples = 2
)

# look at the tuning plot
itd_tuning_ans2$plot_samples

#plot the custom function
ggplot2::ggplot() +
  ggplot2::geom_function(fun = my_custom2, color = "brown", lwd = 1) +
  ggplot2::xlim(-5,60) +
  ggplot2::labs(x = "heights", y = "ws", color = "") +
  ggplot2::theme_light()

----------------------------------------------
#Extract Trees from Point Cloud: Custom====
----------------------------------------------
  
input_las_dir = "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las"  

cloud2trees_ans_c <- cloud2trees::cloud2trees(
  output_dir = tempdir()
  , input_las_dir = i
  , dtm_res_m = 0.5
  , ws = my_custom2
  , estimate_tree_dbh = TRUE
  , estimate_tree_type = FALSE
  , estimate_tree_competition = TRUE
  , estimate_tree_cbh = TRUE
  , cbh_tree_sample_n = 555
  , cbh_estimate_missing_cbh = TRUE
)

paste(
  "Default DTM resolution:"
  , cloud2trees_ans$dtm_rast %>% terra::res() %>% paste(collapse = ",")
  , "|| Custom DTM resolution:"
  , cloud2trees_ans_c$dtm_rast %>% terra::res() %>% paste(collapse = ",")
)

cloud2trees_ans_c$crowns_sf %>% dplyr::glimpse()

paste(
  "Default trees extracted:"
  , cloud2trees_ans$crowns_sf %>% nrow()
  , "|| Custom trees extracted:"
  , cloud2trees_ans_c$crowns_sf %>% nrow()
)

# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$crowns_sf, mapping = ggplot2::aes(fill = dbh_m)) + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$treetops_sf, shape = 20) + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

# DTM
cloud2trees_ans_c$dtm_rast %>% terra::plot()

# CHM
cloud2trees_ans_c$chm_rast %>% terra::plot()

cloud2trees_ans_c$crowns_sf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = tree_height_m, y = dbh_cm)) + 
  ggplot2::geom_point(color = "navy", alpha = 0.6) +
  ggplot2::labs(x = "tree ht. (m)", y = "tree DBH (cm)") +
  ggplot2::scale_x_continuous(limits = c(0,NA)) +
  ggplot2::scale_y_continuous(limits = c(0,NA)) +
  ggplot2::theme_light()

cloud2trees_ans_c$crowns_sf %>%
  dplyr::arrange(is_training_cbh) %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = tree_height_m, y = tree_cbh_m, color=is_training_cbh)) + 
  ggplot2::geom_point() +
  ggplot2::labs(x = "tree ht. (m)", y = "tree CBH (m)") +
  ggplot2::scale_y_continuous(breaks = scales::extended_breaks(n=12)) +
  ggplot2::scale_x_continuous(breaks = scales::extended_breaks(n=14)) +
  ggplot2::scale_color_viridis_d(alpha = 0.8, name = "is CBH/nfrom cloud?") +
  ggplot2::theme_light()


# height plot
plt_ht <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_height_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Oranges", name = "tree ht. (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# diameter plot
plt_dbh <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = dbh_cm)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Purples", name = "tree DBH (cm)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# CBH plot
plt_cbh <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_cbh_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Greens", name = "tree CBH (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")
# combine with patchwork
plt_ht + plt_dbh + plt_cbh + patchwork::plot_layout(ncol = 2) &
  ggplot2::theme(
    legend.title = ggplot2::element_text(size = 8)
    , legend.text = ggplot2::element_text(size = 7)
  )


cloud2trees_ans_c$treetops_sf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(color = comp_dist_to_nearest_m)) + 
  ggplot2::geom_sf() +
  ggplot2::scale_color_distiller(palette = "Greys", name = "distance to nearest tree", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")


# plot tree top points on top of tree crowns 
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$crowns_sf, mapping = ggplot2::aes(fill = cloud2trees_ans_c$crowns_sf$predicted_dbh_cm)) + 
  ggplot2::scale_fill_distiller(palette = "Purples", name = "tree DBH (cm)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = tree_cbh_m)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Greens", name = "tree CBH (m)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")


---------------------------------------------
  # Spatial Stats====
---------------------------------------------

# Convert LAS object to a dataframe
lidar_df <- as.data.frame(las@data)  # Extract X, Y, Z coordinates

# Convert to an sf object (ensure same CRS as tree crowns)
lidar_sf <- st_as_sf(lidar_df, coords = c("X", "Y"), crs = st_crs(cloud2trees_ans_c$crowns_sf))

# Find which LiDAR points are inside each crown area
intersections <- st_intersects(cloud2trees_ans_c$crowns_sf, lidar_sf, sparse = TRUE)

# Count LiDAR points per treeID
lidar_counts <- data.frame(
  treeID = cloud2trees_ans_c$crowns_sf$treeID,
  point_count = lengths(intersections)  # Number of LiDAR points per hull
)

# Merge with tree crowns and compute density
density <- cloud2trees_ans_c$crowns_sf %>%
  left_join(lidar_counts, by = "treeID") %>%
  mutate(point_density = point_count / as.numeric(st_area(geometry)))


ggplot() +
  geom_sf(data = density, aes(fill = point_density), color = "black", alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "LiDAR Density") +
  theme_minimal() +
  ggtitle("LiDAR Point Density Within Tree Crowns")


# histograms
# Check available columns
colnames(cloud2trees_ans_c$crowns_sf)

# Summary of tree heights
summary(cloud2trees_ans_c$crowns_sf$tree_height_m)

ggplot(cloud2trees_ans_c$crowns_sf, aes(x = tree_height_m)) +
  geom_histogram(binwidth = 2, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Canopy Heights",
       x = "Tree Height (m)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(0, max(cloud2trees_ans_c$crowns_sf$tree_height_m, na.rm = TRUE), by = 5))


# Summary of crown area
summary(cloud2trees_ans_c$crowns_sf$crown_area_m2)

ggplot(cloud2trees_ans_c$crowns_sf, aes(x = crown_area_m2)) +
  geom_histogram(binwidth = 2, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Canopy Area",
       x = "Canopy Area (m²)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(0, max(cloud2trees_ans_c$crowns_sf$crown_area_m2, na.rm = TRUE), by = 30))

# Summary of crown area
summary(cloud2trees_ans_c$crowns_sf$crown_area_m2)

hist(cloud2trees_ans_c$crowns_sf$crown_area_m2, breaks=30)
boxplot(cloud2trees_ans_c$crowns_sf$crown_area_m2)

# Find the top 5 largest crown areas
head(sort(cloud2trees_ans_c$crowns_sf$crown_area_m2, decreasing=TRUE), 5)

# Identify which crown polygons have these large areas
large_crowns <- which(cloud2trees_ans_c$crowns_sf$crown_area_m2 > 10)

# Plot all crowns
plot(cloud2trees_ans_c$crowns_sf["crown_area_m2"])

# Highlight the problematic crowns in a different color
plot(cloud2trees_ans_c$crowns_sf[large_crowns,], col="red", add=TRUE)


hist(log10(cloud2trees_ans_c$crowns_sf$crown_area_m2), 
     breaks=30, 
     main="Histogram of log10(Crown Area)", 
     xlab="log10(Crown Area m²)")


boxplot(cloud2trees_ans_c$crowns_sf$crown_area_m2,
        main="Crown Area Distribution",
        ylab="Crown Area (m²)")


plot(density(cloud2trees_ans_c$crowns_sf$crown_area_m2),
     main="Density Plot of Crown Areas",
     xlab="Crown Area (m²)")


par(mfrow=c(1,2))
hist(cloud2trees_ans_c$crowns_sf$crown_area_m2[cloud2trees_ans_c$crowns_sf$crown_area_m2 <= 10],
     main="Crown Areas ≤ 10 m²", 
     xlab="Crown Area (m²)")
hist(cloud2trees_ans_c$crowns_sf$crown_area_m2[cloud2trees_ans_c$crowns_sf$crown_area_m2 > 10],
     main="Crown Areas > 10 m²", 
     xlab="Crown Area (m²)")
par(mfrow=c(1,1))


---------------------------------------------
# DBH Model====
---------------------------------------------

library(pacman)

pacman::p_load(ggplot2, rgeos, propagate, dplyr, ggpubr, gridExtra)
  
setwd(regional_input)
regional.trees <- read.csv("C:/Users/User/Desktop/wade_trainingdata.csv")  

dbh.cm <- regional.trees$dbh.cm
ht.m <- regional.trees$ht.m
  
dbh.mod <- nls(dbh.cm ~ b * ht.m^z, 
                 start = list(b = 2.2, z = 1))

dbh.mod.gs_mixed <- nls(dbh.cm ~ b * ht.m^z, 
                        start = list(b = 1.5, z = 0.8),
                        # Optional: if you have weights to give more importance to certain trees
                        # weights = weights_vector,
                        control = nls.control(maxiter = 100))

# More complex model with species-specific parameters
dbh.mod.mixed <- nls(dbh.cm ~ ifelse(species == "sequoia", 
                                     b_seq * ht.m^z_seq, 
                                     b_hw * ht.m^z_hw), 
                     start = list(b_seq = 1.5, z_seq = 0.8, 
                                  b_hw = 2.0, z_hw = 0.95))

summary(dbh.mod.gs_mixed)
regional.dbh.parameters <- dbh.mod.gs_mixed$m$getPars()


#H = 2.5 × (DBH)^.7

# Extract heights from cloud2trees
tree_heights <- cloud2trees_ans_c$crowns_sf$tree_height_m

# Example field data 
field_data <- data.frame(
  tree_id = c(1, 2, 3), 
  dbh_cm = c(200, 300, 150),
  height_m = c(75, 90, 50)
)

field_data_test <- read.csv("C:/Users/User/Desktop/wade_trainingdata.csv")

# Fit the model using your field measurements
dbh.mod.gs_mixed <- nls(dbh_cm ~ b * height_m^z, 
                        data = field_data_test,
                        start = list(b = 1.5, z = 0.8),
                        control = nls.control(maxiter = 100))

# View model summary
summary(dbh.mod.gs_mixed)


# Extract model coefficients
b_coef <- coef(dbh.mod.gs_mixed)["b"]
z_coef <- coef(dbh.mod.gs_mixed)["z"]

# Predict DBH for all trees in cloud2trees output
cloud2trees_ans_c$crowns_sf$predicted_dbh_cm <- b_coef * cloud2trees_ans_c$crowns_sf$tree_height_m^z_coef

plt_dbh <-
  cloud2trees_ans_c$crowns_sf %>% 
  ggplot2::ggplot(mapping = ggplot2::aes(fill = predicted_dbh_cm)) + 
  ggplot2::geom_sf() + 
  ggplot2::scale_fill_distiller(palette = "Purples", name = "tree DBH (cm)", direction = 1) +
  ggplot2::theme_void() 

ggplot2::ggplot() + 
  ggplot2::geom_sf(data = cloud2trees_ans_c$crowns_sf, mapping = ggplot2::aes(fill = predicted_dbh_cm)) + 
  ggplot2::scale_fill_distiller(palette = "Purples", name = "tree dbh (cm)", direction = 1) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

summary(cloud2trees_ans_c$crowns_sf$predicted_dbh_cm)


# If you have a validation dataset
validation_data <- field_data_test

# Calculate prediction errors
validation_data$predicted_dbh <- b_coef * validation_data$ht.m^z_coef
validation_data$error <- validation_data$predicted_dbh - validation_data$dbh.cm

# Calculate RMSE
rmse <- sqrt(mean(validation_data$error^2))
print(paste("RMSE:", rmse))

# Plot predicted vs observed
plot(validation_data$dbh.cm, validation_data$predicted_dbh,
     xlab = "Measured DBH (cm)", ylab = "Predicted DBH (cm)",
     main = "Validation of DBH-Height Model")
abline(0, 1, col = "red")











