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

# 1- Climatologies computation
# ----------------------------

# Output file

output_clim_dir <- here("input/climatologies/3D_med_clim_month")

# Function to calculate climatologies
# -----------------------------------

process_climato_month <- function(stars_obj, dir) {
  
  # Select years by index : otherwise issues with spatraster from 2000 and 2022 
  
  var_names <- names(stars_obj)
  
  #time_vals <- st_get_dimension_values(stars_obj, "time")
  
  #idx <- which(format(time_vals, "%Y") %in% c("2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009",
  # "2010", "2011", "2012","2013", "2014", "2015", "2016",
  # "2017","2018", "2019", "2020", "2021",  "2022", "2023", "2024"))
  
  #stars_obj <- stars_obj[,,,idx]
  
  # Boucle sur chaque variable
  
  for (var in var_names) {
    
    cat("Traitement de la variable:", var, "\n")
    
    # Extraire la variable
    
    var_data <- stars_obj[var]
    
    # 1. Moyenne par mois sur 2000-2024 (time dim = 12)
    
    month_mean <- aggregate(var_data, by = function(x) format(x, "%m"), FUN = mean, na.rm = TRUE)
    
    # Sauvegarder un fichier global
    
    #saveRDS(month_mean, 
    #file.path(dir, paste0(var, "_month_mean.rds")))
    
    # 2. Standard deviation par mois sur 2000-2024 (time dim = 12)
    
    month_sd <- aggregate(var_data, by = function(x) format(x, "%m"), FUN = sd, na.rm = TRUE)
    
    # Sauvegarder un fichier global
    
    #saveRDS(month_sd, 
    #file.path(dir, paste0(var, "_month_sd.rds")))
    
    # 3. Percentile 10 (minimum moyen)
    
    month_p10 <- aggregate(
      var_data,
      by = function(x) format(x, "%m"),
      FUN = function(x) quantile(x, probs = 0.10, na.rm = TRUE)
    )
    
    # 4. Percentile 90 (maximum moyen)
    
    month_p90 <- aggregate(
      var_data,
      by = function(x) format(x, "%m"),
      FUN = function(x) quantile(x, probs = 0.90, na.rm = TRUE)
    )
    
    # Sauvergarder un fichier par mois
    
    t <- time(var_data)
    
    months <- format(t, "%m")
    
    unique_months <- unique(months)
    
    month_name <- c("jan", "feb", "mar","apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
    
    # Mean
    
    for (m in 1:length(unique_months)) {
      
      obj_m <- month_mean[,m,,]
      
      saveRDS(obj_m, file.path(dir, "per_month", paste0(var,"_", month_name[m], "_mean.rds")))
    }
    
    # Sd
    
    for (m in 1:length(unique_months)) {
      
      obj_m <- month_sd[,m,,]
      
      saveRDS(obj_m, file.path(dir, "per_month", paste0(var,"_", month_name[m], "_sd.rds")))
    }
    
    # Minimum moyen
    
    for (m in 1:length(unique_months)) {
      
      obj_m <- month_p10[,m,,]
      
      saveRDS(obj_m, file.path(dir, "per_month", paste0(var,"_", month_name[m], "_min_moy.rds")))
    }
    
    # Maximum moyen
    
    for (m in 1:length(unique_months)) {
      
      obj_m <- month_p90[,m,,]
      
      saveRDS(obj_m, file.path(dir, "per_month", paste0(var,"_", month_name[m], "_max_moy.rds")))
    }
  }
}

# Apply the function

parameters_3D <- list.files(path = here("input/copernicus/3D_med_raw/combine_year"),
                            pattern = "\\.(rds)$", 
                            full.names = TRUE)

for (i in 1:length(parameters_3D)) {
  
  file <- readRDS(parameters_3D[i])
  
  process_climato_month(file, output_clim_dir)
  
  # Clean memory
  
  rm(file)
  
  gc()
  
}++++++++++++++++;