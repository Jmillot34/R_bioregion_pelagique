# Environmental variable selection

# Library

library(terra)
library(dplyr)
library(here)
library(sf)
library(stars)
library(viridis)
library(ade4)
library(factoextra)
library(patchwork)
library(tmap)

## --------------- Climato done on meSU_climato_3D ---------------------------

## --------------- Final process on climato 3D ------------------------------- Not done 

# Conversion stars in raster
# --------------------------

scen <- "Gdl_0-200_month"

# Climato files

month <- c("jan", "feb", "mar","apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

# Output directories

processed_dir <- here("input", "fPCA", scen, "stack")

# List climatologie

clim_files <- list.files(path = here("input/fPCA", scen, "raster"), 
                         pattern = "\\.tif$", 
                         full.names = TRUE)

clim_names <- tools::file_path_sans_ext(basename(clim_files))

for (i in seq_along(month)) {
  
  clim_files_month <- clim_files[grepl(month[i], clim_names)]
  
  clim_list_month <- lapply(clim_files_month, terra::rast)

# Standardized file formatting
# -----------------------------

# 1. Identify common extent and resolution

# Common extent

all_extents <- lapply(clim_list_month, ext)

common_extent <- ext(
  max(sapply(all_extents, function(e) e[1])),  # max xmin
  min(sapply(all_extents, function(e) e[2])),  # min xmax
  max(sapply(all_extents, function(e) e[3])),  # max ymin
  min(sapply(all_extents, function(e) e[4]))   # min ymax
)
common_extent

# Most common resolution

## Get all resolutions
resolutions <- lapply(clim_list_month, res)

## Convert to data frame
res_df <- do.call(rbind, resolutions)
colnames(res_df) <- c("x_res", "y_res")

## Find most common resolution
res_table <- table(paste(res_df[,1], res_df[,2], sep = "_"))
most_common <- names(which.max(res_table))

## Extract x and y resolution
most_common_res <- as.numeric(strsplit(most_common, "_")[[1]])
cat("Most common resolution: x =", most_common_res[1], ", y =", most_common_res[2], "\n")

# 2.Resampling on the reference raster

# Reference raster
ref_rast <- rast(
  extent = common_extent,
  resolution = most_common_res,
  crs = "EPSG:4326")

# Resampling
resample_rast <- lapply(clim_list_month, function(r) resample(r, ref_rast, method = "bilinear"))

# 3. Stack of rasters + mask of the marine domain

# Stack
clim_stack <- rast(resample_rast)

# Mask
med_mask <- st_read(here::here("Input/GSA/GSAs_simplified.shp")) %>% 
  filter(!SECT_COD %in% c("GSA28", "GSA29", "GSA30")) %>% 
  st_transform(4326)

clim_stack_mask <- clim_stack %>% mask(med_mask)

# Save
writeRaster(clim_stack_mask, 
            filename = file.path(processed_dir, paste0("fPCA_3D_stack_", month[i], ".tif")),
            overwrite = TRUE)

}



