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


ggplot2::ggplot() +
  ggplot2::geom_function(fun = my_custom2, color = "brown", lwd = 1) +
  ggplot2::xlim(-5,60) +
  ggplot2::labs(x = "heights", y = "ws", color = "") +
  ggplot2::theme_light()
