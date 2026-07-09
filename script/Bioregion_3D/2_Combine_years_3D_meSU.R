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

# 1- Read annual 3D files and combine along time
# ----------------------------------------------

# Output file
output_combine_dir <- here("input/copernicus/3D_med_raw/combine_year")

# Parameters
parameters <- list.files(path = here("input/copernicus/3D_med_raw/single_year"),
                         pattern = "\\.(rds)$", 
                         full.names = TRUE)

# Extract layer base names (remove last 14 characters)

names_layers <- sub("_.*", "", basename(parameters))

# Group files by layer name

group_layers <- split(parameters, names_layers)

# # Combine each group along time 

for (i in 1:length(group_layers)) {
  
  group <- group_layers[[i]]
  layer_name <- names(group_layers)[i]
  
  message("Processing layer: ", layer_name)
  
  # Read RDS files
  stars_list <- lapply(group, readRDS)
  
  # Combine along time
  combined <- do.call(c, c(stars_list, along = "time"))
  
  # Save combined object
  saveRDS(combined, 
          file.path(output_combine_dir, paste0(layer_name, "_2000_2024.rds")))
  
  # Clean memory
  
  rm(combined)
  
  gc()
  
}