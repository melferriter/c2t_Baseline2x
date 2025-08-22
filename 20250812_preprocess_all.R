# 

----------------------------------------------
  # Double Return====
----------------------------------------------
baseline2x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xBaseline_20241021232243_clip.las")
sidelap2x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xSidelap_20241021223703_clip.las")
crosshatch2x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xCrosshatch_20241022000644_clip.las")
speed2x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/2x/2xSpeed_20241021230048_clip.las")

----------------------------------------------
  # Triple Return====
----------------------------------------------
baseline3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xBaseline_20241022184540_clip.las")
sidelap3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xSide_20241022190851_clip.las")
crosshatch3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xCross_20241022194753_clip_clip.las")
speed3x <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/3x/3xSpeed_20241022192741._cliplaz.las")

----------------------------------------------
  # Structure from Motion====
----------------------------------------------
sfm <- readLAS("E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Raw/20241022/Sam/FreemanCreekGrove_south_group1_densified_point_cloud.las")

----------------
#read in the SFM and 3XCross Las files and clip to ROI
roi <- st_read("C:/Users/User/Desktop/LidarClip/LidarExtent.shp", crs = 32611)
roi_2d <- st_zm(roi, drop = TRUE, what = "ZM")

----------------
  
baseline2x_clipped <- clip_roi(baseline2x, roi_2d)
sidelap2x_clipped <- clip_roi(sidelap2x, roi_2d)
crosshatch2x_clipped <- clip_roi(crosshatch2x, roi_2d)
speed2x_clipped <- clip_roi(speed2x, roi_2d)
 
baseline3x_clipped <- clip_roi(baseline3x, roi_2d)
sidelap3x_clipped <- clip_roi(sidelap3x, roi_2d)
crosshatch3x_clipped <- clip_roi(crosshatch3x, roi_2d)
speed3x_clipped <- clip_roi(speed3x, roi_2d)
 
sfm_clipped <- clip_roi(sfm, roi_2d)  


---------------

writeLAS(baseline2x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline2x_clipped.las")
writeLAS(sidelap2x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap2x_clipped.las")
writeLAS(crosshatch2x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/crosshatch2x_clipped.las")
writeLAS(speed2x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/speed2x_clipped.las")

writeLAS(baseline3x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/baseline3x_clipped.las")
writeLAS(sidelap3x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sidelap3x_clipped.las")
writeLAS(crosshatch3x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/crosshatch3x_clipped.las")
writeLAS(speed3x_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/speed3x_clipped.las")

writeLAS(sfm_clipped, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/sfm_clipped.las")


------------------------

crs(crosshatch3x_clipped)
crs(sfm_clipped)


names(crosshatch3x_clipped@data)
names(sfm_clipped@data)

crosshatch3x_clipped@data$R <- NA
crosshatch3x_clipped@data$G <- NA
crosshatch3x_clipped@data$B <- NA

library(data.table) 

# Find missing columns
missing_cols <- setdiff(names(crosshatch3x_clipped@data),
                        names(sfm_clipped@data))

# Add missing columns with NA values
for (col in missing_cols) {
  sfm_clipped@data[[col]] <- NA
}

# Now match the column order
sfm_clipped@data <- as.data.table(
    sfm_clipped@data[ , names(crosshatch3x_clipped@data), with = FALSE]
)

sfm_clipped@data <- sfm_clipped@data[, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

header_template <- crosshatch3x_clipped@header
sfm_clipped@header <- header_template

# Define all standard LAS integer columns
las_int_cols <- c("ReturnNumber", "NumberOfReturns", "Classification", "ScanDirectionFlag",
                  "EdgeOfFlightline", "Synthetic_flag", "Keypoint_flag", "Withheld_flag",
                  "ScanAngleRank", "UserData", "PointSourceID", "ScannerChannel")

# Ensure both objects are using correct types
for (col in las_int_cols) {
  if (col %in% names(crosshatch3x_clipped@data)) {
    crosshatch3x_clipped@data[[col]] <- as.integer(crosshatch3x_clipped@data[[col]])
  }
  if (col %in% names(sfm_clipped@data)) {
    sfm_clipped@data[[col]] <- as.integer(sfm_clipped@data[[col]])
  }
}

# Convert flag fields to logical
flag_fields <- c("Synthetic_flag", "Keypoint_flag", "Withheld_flag", "Overlap_flag")

for (col in flag_fields) {
  if (col %in% names(crosshatch3x_clipped@data)) {
    crosshatch3x_clipped@data[[col]] <- as.logical(crosshatch3x_clipped@data[[col]])
  }
  if (col %in% names(sfm_clipped@data)) {
    sfm_clipped@data[[col]] <- as.logical(sfm_clipped@data[[col]])
  }
}

crosshatch3x_clipped@data$R <- NULL
crosshatch3x_clipped@data$G <- NULL
crosshatch3x_clipped@data$B <- NULL

sfm_clipped@data$R <- NULL
sfm_clipped@data$G <- NULL
sfm_clipped@data$B <- NULL


las_combined <- rbind(crosshatch3x_clipped, sfm_clipped)

las_combined@data$ReturnNumber[las_combined@data$ReturnNumber == 0] <- 1
las_combined@data$NumberOfReturns[las_combined@data$NumberOfReturns == 0] <- 1

header <- las_combined@header
header@PHB$XScaleFactor <- 0.001
header@PHB$YScaleFactor <- 0.001
header@PHB$ZScaleFactor <- 0.001

las_combined@header <- header

# Convert ReturnNumber and NumberOfReturns to integer
las_combined@data$ReturnNumber <- as.integer(las_combined@data$ReturnNumber)
las_combined@data$NumberOfReturns <- as.integer(las_combined@data$NumberOfReturns)


plot(las_combined)
summary(las_combined)


writeLAS(las_combined, "E:/Grad School/Data/UAS/Sequoia_National_Forest/2024/2024101222_processed/Agisoft/ROI_Las/SFM_3xCross_20241022194753_20250809.las")

-----------------------------