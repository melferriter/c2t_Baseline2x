###############################################################################
###############################################################################
###                                                                         ###
### Turnkey Photogrammetry Point Cloud Processing Script                    ###
### MethodsX paper script 1:                                                ###
### Written by Neal Swayze 3/26/2021                                        ###
### POC: nealswayze1@gmail.com                                              ###
###                                                                         ###
###############################################################################
###############################################################################

###############################################################################
#####################      Install/ Load Packages       #######################
###############################################################################

install.packages("pacman")
library(pacman)
pacman::p_load(lidR, ForestTools, raster, sp, devtools)

#devtools::install_github("tiagodc/TreeLS", force=T)
library(TreeLS)
library(ForestTools)
#library(rgdal)
library(beepr)
beep(sound = 1)

rm(list = ls(globalenv()))

###############################################################################
####################        PREPERATION STEPS         #########################
###############################################################################

# 1. Create a folder on your desktop with the name "engine"

# 2. Put your Structure From Motion point cloud in the folder

# 3. Point your root directory at the engine folder

###############################################################################
########################     Directory Structure      #########################
###############################################################################

### Starting Time #############################################################
start_time <- Sys.time()

### Establish Root Directory ##################################################
rootDir <- ("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/R Scripts/Wade_v3")

setwd(rootDir)

### Creating a set of processing folders for outputs ##########################
dir.create(file.path(rootDir, "/outputs"))
dir.create(file.path(rootDir, "/outputs/01_las"))
dir.create(file.path(rootDir, "/outputs/02_CHM"))
dir.create(file.path(rootDir, "/outputs/03_single_tree_csv"))
dir.create(file.path(rootDir, "/outputs/03_single_tree_shp"))

### Setting your output directories for processing steps ######################
las_output <- file.path(rootDir, "/outputs/01_las")
chm_output <- file.path(rootDir, "/outputs/02_CHM")
csv_output <- file.path(rootDir, "/outputs/03_single_tree_csv")
shp_output <- file.path(rootDir, "/outputs/03_single_tree_shp")


###############################################################################
#####################      PROCESSING STARTS BELOW      #######################
###############################################################################

### Reading in a single normalized .las file ##################################
las = readLAS("2xBaseline_20241021232243_clip.las")
las_check(las)
plot(las)

###############################################################################
#############      Classify Ground and Normalize Point Cloud      #############
###############################################################################
las <- classify_ground(las, csf(sloop_smooth = TRUE,  # Changed to TRUE
                                class_threshold = 0.5,  # Increased 
                                cloth_resolution = 0.5, # Increased
                                rigidness = 1L,        # Reduced
                                iterations = 700L,     # Increased
                                time_step = 0.5))
las <- normalize_height(las, knnidw())
beep(sound = 1)
plot(las)

### Save out height normalized las file #######################################
setwd(las_output)
writeLAS(las, "normalized_cloud.las", index = FALSE)


ground <- filter_poi(las, Classification == 2L)
plot(ground)
writeLAS(ground, "ground.las", index = FALSE)
beep(sound = 1)
###############################################################################
##########      Rasterizing point cloud to canopy height model      ###########
###############################################################################
chm_0.1 <- grid_canopy(las, res = 0.1, algorithm = p2r())
plot(chm_0.1)

###############################################################################
######################      Detecting Single Trees      #######################
###############################################################################

### Set exponential function for single tree detection and minimum tree ####### 
### height for detection.The below Expo function is based on work by    #######     ### Do we want descriptive annotation of why we are doing things?
### Creasy et al. 2021.                                                 #######
Expo <- function(x) (x <- (1 + x * 0.1))
height_min <- 1.37

### Extract tree locations and plot them ######################################
tree_locations <- vwf(chm_0.1, my_custom2, height_min)
beep(sound = 1)
n_height <- nrow(tree_locations)
cat(n_height," trees identified")
sp_summarise(tree_locations, variables = "height")                                
plot(tree_locations, add = TRUE, col = "red")

#
#################MELISSA TEST###############
f <- function(x) {x * 0.1 + 3}
my_custom2 <- function(x) {0.15 * x^0.6 + 2}
ttops <- locate_trees(chm_0.1, lmf())

############################################

###############################################################################
#######################      Extracting DBH Values      #######################
###############################################################################

### Remove points above 4 m to reduce memory requirements #####################
las <- filter_poi(las, Z < 40)
las <- filter_poi(las1, Z > 0.1)
las <- las2
x = plot(las)
las_check(las2)
beep(3)
### Extract tree map from a thinned point cloud ###############################
map <- treeMap(las, map.hough(min_h = .5,             # Slightly wider range
                              min_density = 0.005,    # Lower for sparse points
                              max_d = 15)) 

map1 <- treeMap(las1, map.eigen.voxel(
  max_curvature = 0.15,
  max_verticality = 15,
  voxel_spacing = 0.1,
  max_d = 9
))
beep(3)

add_treeMap(x, map, color = "yellow", size = 2)
add_treetops3d(x, map)
plot(map)

### Classify tree regions #####################################################
las <- treePoints(las1, map1, trp.crop())
add_treePoints(x, las, size = 9)
add_treeIDs(x, las, cex = 2, col = "yellow")

### Classify stem points ######################################################
las <- stemPoints(las, stm.hough(max_d = 15))
add_stemPoints(x, las, color = "red", size = 8)

### Make the plot's inventory #################################################
dbh_locations <- tlsInventory(las, 
                              dh = 1.37, 
                              d_method = shapeFit(shape = "circle", 
                                                  algorithm = "ransac")) 

dbh_locations <- subset(dbh_locations, Radius < 5)
add_tlsInventory(x, dbh_locations)
n_dbh <- nrow(dbh_locations)
cat(n_dbh," tree DBHs identified")
dbh_locations$DBH <- dbh_locations$Radius * 2 * 100
summary(dbh_locations)

xy <- dbh_locations[,c(2,3)]
dbh_shape <- SpatialPointsDataFrame(coords = xy, data = dbh_locations)

plot(dbh_shape, add = TRUE, col = "blue")

###############################################################################
###                Write out derived data products                          ###
###############################################################################

### Save CHM ##################################################################
setwd(chm_output)
writeRaster(chm_0.1, file = "chm_0.1m.tif", overwrite = TRUE)

### Save tree locations and heights to csv and shapefile ######################
setwd(csv_output)
write.csv(ttops, "ttops")
setwd(shp_output)
rgdal::writeOGR(ttops, dsn = '.', layer = "tree_heights", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)

### Save DBH locations and values to csv and shapefile #######################
setwd(csv_output)
write.csv(dbh_locations, "tree_dbh.csv")
setwd(shp_output)
writeOGR(dbh_shape, dsn = '.', layer = "dbh_locations", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)

### End Script ################################################################
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

