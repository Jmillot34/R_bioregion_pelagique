library(terra)
library(dplyr)
library(here)
library(sf)
library(stars)
library(viridis)
library(ade4)
library(factoextra)
library(tmap)
library(fdapace)

# Link toward explanations of fPCA
# -------------------------------

#https://github.com/EPauthenet/fda.oce : use by Valentin (fit 20 Bspline on a matrix level x stations x variables / then 
#apply a fPCA) : but cannot handle NA, and in our case we have incomplete profils
#https://rdrr.io/cran/fdapace/f/vignettes/fdapaceVig.Rmd : 

# Scenario
# --------

scen <- "srm_0-200_month"
#scen <- "med_0-200_month"

# Read Input data
# ---------------

# Crop on srm = sous-region marine

ext_srm <- st_bbox(c(
  xmin = 2,
  xmax = 12,
  ymin = 38.7, 
  ymax = 44.7
), crs = st_crs(4326))


# Month mean

data_3D_mean <- list.files(here("input/climatologies/3D_med_clim_month/per_month"), pattern = "mean", full.names = TRUE)

data_3D_name <- tools::file_path_sans_ext(basename(data_3D_mean))

# Read RDS

data_3D_files_raw <- lapply(data_3D_mean, function(x){
  
  data <- readRDS(x) %>% 
    st_crop(ext_srm) # Remove if Med
  
  # Only keep 3D files
  if ("elevation" %in% names(dim(data))) {
    
    data[,,, drop = TRUE] # remove an artifact of dimension (geometry)
  
    } 
  
  else {
    NULL
  }
})

keep <- !sapply(data_3D_files_raw, is.null)

data_3D_files <- data_3D_files_raw[keep]

data_3D_name  <- data_3D_name[keep]

# Output directory
# ----------------

fpca_dir <- here("input", "fPCA", scen)
dir.create(fpca_dir , recursive = TRUE, showWarnings = FALSE)

raster_dir <- here(fpca_dir, "raster")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

expl_var_dir <- here(fpca_dir, "expl_var")
dir.create(expl_var_dir, recursive = TRUE, showWarnings = FALSE)

map_dir <- here(fpca_dir, "map")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

stack_dir <- here(fpca_dir, "stack")
dir.create(stack_dir, recursive = TRUE, showWarnings = FALSE)

# Compute fPCA : first method
# ----------------------------

# 1. Initialisation

depth_range <- c(1:35)

pc1_ref <- NULL

pc2_ref <- NULL

# Table storing explained variance for each predictor
variance_table <- data.frame(
  predictor = character(),
  PC1 = numeric(),
  PC2 = numeric(),
  PC3 = numeric(),
  PC4 = numeric(),
  stringsAsFactors = FALSE
)

# 2. Function

for (i in seq_along(data_3D_files)) {
  
  cat("Computing fPCA on ", data_3D_name[i], "\n")

  # Crop by depth
  data_depth <- data_3D_files[[i]][,,,depth_range]
  
  # Convert stars object to long-format data.frame
  df <- as.data.frame(data_depth, long = TRUE)
  
  # Rename value column with parameter name
  colnames(df)[4] <- data_3D_name[i]
  
  # Create pixel ID = vertical profile
  df <- df %>%
    # Assign ID to each pixel
    group_by(longitude, latitude) %>%
    mutate(id = cur_group_id()) %>%
    # Remove profiles containing only NA
    filter(!all(is.na(.data[[data_3D_name[i]]])))
  # Remove profiles with at least one NA (uncomment if needed for complete profiles only)
  # filter(!any(is.na(.data[[data_3D_name[i]]])))
  
  df_meta <- df %>%
    distinct(longitude, latitude, id) %>%
    arrange(id)
  
  # Build FPCA inputs
  
  # Observations
  Ly <- split(df[[data_3D_name[i]]], df$id)
  
  # Depth (treated as time)
  Lt <- split(df$elevation, df$id)
  
  # Run fPCA
  fpca_res <- FPCA(
    Ly = Ly,
    Lt = Lt,
    optns = list(dataType = "Dense")  # Dense but truncated profiles (not sparse due to regular sampling)
  )
  
  # Explained variance
  # -----------------
  
  var_explained <- fpca_res$lambda / sum(fpca_res$lambda) * 100
  
  # Four principale components
  cat(
    "Variance explained:\n",
    "PC1 =", round(var_explained[1], 2), "%\n",
    "PC2 =", round(var_explained[2], 2), "%\n",
    "PC3 =", round(var_explained[3], 2), "%\n",
    "PC4 =", round(var_explained[4], 2), "%\n"
  )
  
  # Ensure at least 4 values
  var_tmp <- rep(NA, 4)
  var_tmp[1:min(4, length(var_explained))] <- var_explained[1:min(4, length(var_explained))]
  
  # Add row to summary table
  variance_table <- rbind(
    variance_table,
    data.frame(
      predictor = data_3D_name[i],
      PC1 = round(var_tmp[1], 2),
      PC2 = round(var_tmp[2], 2),
      PC3 = round(var_tmp[3], 2),
      PC4 = round(var_tmp[4], 2)
    )
  )
  
  # Save table
  write.csv(
    variance_table,
    file.path(expl_var_dir, paste0(data_3D_name[i], "expl_var.csv")),
    row.names = FALSE
  )

  # PC1 scores
  # ----------
  
  pc1 <- fpca_res$xiEst[,1]
  
  # Check PC1 sign consistency
  if (!is.null(pc1_ref)) {
    cor_val <- cor(pc1_ref, pc1, use = "complete.obs")
    
    # Align sign with reference PC1
    if (cor_val < 0) {
      cat("PC1 inversion yes \n")
      pc1 <- -pc1
      fpca_res$xiEst[,1] <- pc1
      fpca_res$phi[,1] <- -fpca_res$phi[,1]  # Important: eigenfunction must also be inverted
    }
  }
  
  # Update reference PC1
  pc1_ref <- pc1
  
  # Add PC1 scores to metadata
  df_meta$pc1 <- pc1
  
  # Convert to raster
  pc1_rast <- rast(df_meta[, c("longitude", "latitude", "pc1")], type = "xyz")
  names(pc1_rast) <- paste0(data_3D_name[i], "_pc1")
  
  # Save Map
  png(file.path(map_dir, paste0(data_3D_name[i], "_3D_pc1.png")), width = 1000, height = 800)
  plot(pc1_rast)
  dev.off()
  
  # Save Raster
  writeRaster(pc1_rast, file.path(raster_dir, paste0(data_3D_name[i], "_3D_pc1.tif")), overwrite = TRUE)
  
  # Check if a second axis exists
  # -----------------------------
  
  if (ncol(fpca_res$xiEst) >= 2) {
    
    # PC2 scores
    # ----------
    
    pc2 <- fpca_res$xiEst[,2]
    
    # Check PC2 sign consistency
    if (!is.null(pc2_ref)) {
      cor_val <- cor(pc2_ref, pc2, use = "complete.obs")
      
      # Align sign with reference PC2
      if (cor_val < 0) {
        cat("PC2 inversion yes \n")
        pc2 <- -pc2
        fpca_res$xiEst[,2] <- pc2
        fpca_res$phi[,2] <- -fpca_res$phi[,2]  # Important: eigenfunction must also be inverted
      }
    }
    
    # Update reference PC2
    pc2_ref <- pc2
    
    # Add PC2 scores to metadata
    df_meta$pc2 <- pc2
    
    # Convert to raster
    pc2_rast <- rast(df_meta[, c("longitude", "latitude", "pc2")], type = "xyz")
    names(pc2_rast) <- paste0(data_3D_name[i], "_pc2")
    
    # Save Map
    png(file.path(map_dir, paste0(data_3D_name[i], "_3D_pc2.png")), width = 1000, height = 800)
    plot(pc2_rast)
    dev.off()
    
    # Save Raster
    writeRaster(pc2_rast, file.path(raster_dir, paste0(data_3D_name[i], "_3D_pc2.tif")), overwrite = TRUE)
    
  }
  
}

# ============
# Stack fPCA
# ============

# Output directories

processed_dir <- here("input", "fPCA", scen, "stack")

# List climatologie

clim_files <- list.files(path = here("input/fPCA", scen, "raster"), 
                         pattern = "\\.tif$", 
                         full.names = TRUE)

clim_names <- tools::file_path_sans_ext(basename(clim_files))


# ====================================
# Standardized files format + stacking
# ====================================

# Apply the function per month 

month <- c("jan", "feb", "mar","apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

for (i in seq_along(month)) {
  
  clim_files_month <- clim_files[grepl(month[i], clim_names)]
  
  clim_list_month <- lapply(clim_files_month, terra::rast)
  
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


# # Compute fPCA : second method
# # ----------------------------
# 
# install.packages("devtools")
# devtools::install_github("Epauthenet/fda.oce")
# 
# chl_mean <- readRDS(here("input/climatologies/3D_med_clim/chl_mean.rds"))
# names(chl_mean) <- "chl"
# dissic_mean <- readRDS(here("input/climatologies/3D_med_clim/dissic_mean.rds"))
# names(dissic_mean) <- "dissic"
# nh4_mean <- readRDS(here("input/climatologies/3D_med_clim/nh4_mean.rds"))
# names(nh4_mean) <- "nh4"
# 
# # Converted in array
# 
# nh4_arr    <- nh4_mean[[1]]
# chl_arr    <- chl_mean[[1]]
# dissic_arr <- dissic_mean[[1]]
# 
# # Reshape to station x depth
# nlon <- dim(nh4_arr)[1]
# nlat <- dim(nh4_arr)[2]
# ndepth <- dim(nh4_arr)[3]
# 
# nstation <- nlon * nlat
# 
# nh4_mat <- matrix(nh4_arr, nrow = nstation, ncol = ndepth)
# chl_mat <- matrix(chl_arr, nrow = nstation, ncol = ndepth)
# dissic_mat <- matrix(dissic_arr, nrow = nstation, ncol = ndepth)
# 
# # Build final array
# arr <- array(NA, dim = c(ndepth, nstation, 3))
# 
# # Fill with the transposed matrix : because they were in station x depth and we want depth x station
# 
# arr[,,1] <- t(nh4_mat)
# arr[,,2] <- t(chl_mat)
# arr[,,3] <- t(dissic_mat)
# 
# names(dim(arr)) <- c("depth", "station", "variable")
# 
# dimnames(arr) <- list(
#   depth = NULL,
#   station = NULL,
#   variable = c("nh4", "chl", "dissic")
# )
# 
# # Remove profil with all NA
# na_mask <- apply(arr, 2, function(x) all(is.na(x)))
# 
# # Garder seulement les stations valides
# arr_clean <- arr[, !na_mask, ] %>% as.array()
# 
# # Pi containing the level
# 
# Pi <- st_get_dimension_values(chl_mean, "elevation") %>% as.vector()
# 
# # Apply Bsplin
# library(fda.oce)
# 
# fdobj <- bspl(Pi, arr_clean, nbas = 20, fdn = list("nh4", "chl", "dissic"))
# 
# # Do not work because need complete profil without NA
# 
# fpca_test <- fpca(fdobj)


