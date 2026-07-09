# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# RASTER ANALYSIS AND HARMONIZATION FUNCTIONS
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+

#' Analyze Raster Files
#'
#' Extracts metadata from multiple raster files including dimensions, extent,
#' resolution, and CRS information.
#'
#' @param files Character vector of file paths to raster files
#' @return Data frame with raster metadata for each file
#' @export
analyze_rasters <- function(files) {
  results <- list()
  
  for (i in seq_along(files)) {
    tryCatch({
      r <- rast(files[i])
      
      results[[i]] <- data.frame(
        file = basename(files[i]),
        path = files[i],
        nrow = nrow(r),
        ncol = ncol(r),
        nlyr = nlyr(r),
        xmin = ext(r)[1],
        xmax = ext(r)[2],
        ymin = ext(r)[3],
        ymax = ext(r)[4],
        res_x = res(r)[1],
        res_y = res(r)[2],
        crs = crs(r, proj = TRUE)
      )
    }, error = function(e) {
      results[[i]] <<- data.frame(
        file = basename(files[i]),
        path = files[i],
        error = as.character(e$message)
      )
    })
  }
  
  df_results <- bind_rows(results)
  return(df_results)
}

#' Harmonize Raster Files
#'
#' Harmonizes multiple raster files to have the same resolution and extent.
#' Identifies the most common resolution and extent, then resamples files
#' that don't match.
#'
#' @param files Character vector of file paths to raster files
#' @param output_dir Directory to save harmonized rasters (default: "Datasets/harmonised")
#' @return List containing the harmonized stack, metadata, and common parameters
#' @export
harmonize_rasters <- function(files, output_dir = here::here("Input/Datasets/harmonised")) {
  library(terra)
  library(dplyr)
  
  # 1. Analyze all rasters
  cat("Analyzing rasters...\n")
  raster_info <- analyze_rasters(files)
  
  # 2. Identify most common resolution
  raster_info <- raster_info %>%
    mutate(resolution = paste(round(res_x, 6), round(res_y, 6), sep = "_"))
  
  common_res <- raster_info %>%
    count(resolution, res_x, res_y) %>%
    arrange(desc(n)) %>%
    slice(1)
  
  cat("Most common resolution:", common_res$resolution, 
      "(", common_res$n, "files)\n")
  
  # 3. Identify most common extent
  raster_info <- raster_info %>%
    mutate(extent_id = paste(round(xmin, 4), round(xmax, 4), 
                             round(ymin, 4), round(ymax, 4), sep = "_"))
  
  common_extent <- raster_info %>%
    count(extent_id, xmin, xmax, ymin, ymax) %>%
    arrange(desc(n)) %>%
    slice(1)
  
  cat("Most common extent:", common_extent$extent_id, 
      "(", common_extent$n, "files)\n")
  
  # 4. Create reference raster
  ref_ext <- ext(common_extent$xmin, common_extent$xmax, 
                 common_extent$ymin, common_extent$ymax)
  
  ref_file <- raster_info %>%
    filter(resolution == common_res$resolution,
           extent_id == common_extent$extent_id) %>%
    slice(1) %>%
    pull(path)
  
  ref_raster <- rast(ref_file)
  
  # 5. Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 6. Process each raster
  cat("\nProcessing rasters...\n")
  harmonized_rasters <- list()
  
  for (i in seq_along(files)) {
    file_name <- basename(files[i])
    cat("  -", file_name)
    
    r <- rast(files[i])
    info <- raster_info[i, ]
    
    needs_adjustment <- (info$resolution != common_res$resolution) | 
      (info$extent_id != common_extent$extent_id)
    
    if (needs_adjustment) {
      cat(" [ADJUSTING]")
      r_harmonized <- resample(r, ref_raster, method = "near")
      
      output_path <- file.path(output_dir, file_name)
      writeRaster(r_harmonized, output_path, overwrite = TRUE)
      
      harmonized_rasters[[i]] <- r_harmonized
    } else {
      cat(" [OK]")
      harmonized_rasters[[i]] <- r
    }
    cat("\n")
  }
  
  # 7. Create final stack
  cat("\nCreating final stack...\n")
  final_stack <- rast(harmonized_rasters)
  
  # 8. Summary
  cat("\n=== SUMMARY ===\n")
  cat("Rasters processed:", length(files), "\n")
  cat("Rasters adjusted:", sum(raster_info$resolution != common_res$resolution | 
                                 raster_info$extent_id != common_extent$extent_id), "\n")
  cat("Final resolution:", common_res$res_x, "x", common_res$res_y, "\n")
  cat("Final extent: xmin=", common_extent$xmin, ", xmax=", common_extent$xmax,
      ", ymin=", common_extent$ymin, ", ymax=", common_extent$ymax, "\n")
  
  return(list(
    stack = final_stack,
    info = raster_info,
    common_resolution = common_res,
    common_extent = common_extent
  ))
}

#' Load Layers with Optimization
#'
#' Loads raster layers, automatically using harmonized versions when available
#' and necessary, while using original files when they already match the
#' common resolution and extent.
#'
#' @param main_dir Main directory containing original raster files
#' @param harmonized_dir Directory containing harmonized raster files
#' @param pattern Regex pattern to match raster files (default: NetCDF and TIFF)
#' @return List containing the final stack and file information
#' @export
load_optimized_layers <- function(main_dir = here::here("Datasets"),
                                  harmonized_dir = here::here("Datasets/harmonised"),
                                  pattern = "\\.(nc|tiff?)$") {
  library(terra)
  library(dplyr)
  
  # 1. List all original files
  original_files <- list.files(path = main_dir, 
                               pattern = pattern, 
                               full.names = TRUE,
                               recursive = FALSE)
  
  # 2. Analyze rasters to identify which are OK
  cat("Analyzing original rasters...\n")
  raster_info <- analyze_rasters(original_files)
  
  # 3. Identify most common resolution and extent
  raster_info <- raster_info %>%
    mutate(
      resolution = paste(round(res_x, 6), round(res_y, 6), sep = "_"),
      extent_id = paste(round(xmin, 4), round(xmax, 4), 
                        round(ymin, 4), round(ymax, 4), sep = "_")
    )
  
  common_res <- raster_info %>%
    count(resolution) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    pull(resolution)
  
  common_extent <- raster_info %>%
    count(extent_id) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    pull(extent_id)
  
  # 4. Identify files that are OK (no adjustment needed)
  ok_files <- raster_info %>%
    filter(resolution == common_res, extent_id == common_extent) %>%
    pull(path)
  
  # 5. Identify files to replace with harmonized versions
  files_to_replace <- raster_info %>%
    filter(resolution != common_res | extent_id != common_extent) %>%
    pull(file)
  
  # 6. Build final file list
  final_files <- ok_files
  
  # 7. Add harmonized versions
  if (length(files_to_replace) > 0 && dir.exists(harmonized_dir)) {
    harmonized_files <- file.path(harmonized_dir, files_to_replace)
    
    existing_harmonized <- harmonized_files[file.exists(harmonized_files)]
    final_files <- c(final_files, existing_harmonized)
    
    missing_files <- files_to_replace[!file.exists(harmonized_files)]
    if (length(missing_files) > 0) {
      warning("Missing harmonized files:\n", 
              paste("  -", missing_files, collapse = "\n"))
    }
  }
  
  # 8. Display summary
  cat("\n=== LOADING SUMMARY ===\n")
  cat("Original files OK:", length(ok_files), "\n")
  cat("Harmonized files used:", 
      sum(file.exists(file.path(harmonized_dir, files_to_replace))), "\n")
  cat("Total files to load:", length(final_files), "\n\n")
  
  # 9. Display file provenance
  cat("File provenance:\n")
  for (f in final_files) {
    provenance <- if (dirname(f) == normalizePath(harmonized_dir)) {
      "[HARMONIZED]"
    } else {
      "[ORIGINAL]"
    }
    cat("  ", provenance, basename(f), "\n")
  }
  
  # 10. Load and return stack
  cat("\nLoading rasters...\n")
  final_stack <- rast(final_files)
  
  return(list(
    stack = final_stack,
    files = final_files,
    original_files = ok_files,
    harmonized_files = file.path(harmonized_dir, files_to_replace)
  ))
}

#' Standardize Layer Names
#'
#' Standardizes layer names using a flexible pattern matching dictionary.
#' Handles various naming conventions for common environmental variables.
#'
#' @param layers SpatRaster object with layers to rename
#' @return SpatRaster object with standardized layer names
#' @export
standardize_layer_names <- function(layers) {
  library(stringr)
  
  # Dictionary with multiple patterns for each variable
  dictionary <- list(
    chl_max = c("chl_max", "chl.*max", "chlorophyll.*max"),
    chl_min = c("chl_min", "chl.*min", "chlorophyll.*min"),
    no3_max = c("no3_max", "no3.*max", "nitrate.*max"),
    no3_min = c("no3_min", "no3.*min", "nitrate.*min"),
    par_min = c("par.*min", "radiation.*min"),
    par_max = c("par.*max", "radiation.*max"),
    phyc_max = c("phyc_max", "phyc.*max", "phytoplankton.*max"),
    phyc_min = c("phyc_min", "phyc.*min", "phytoplankton.*min"),
    po4_min = c("po4_min", "po4.*min", "phosphate.*min"),
    po4_max = c("po4_max", "po4.*max", "phosphate.*max"),
    si_min = c("^si_min$", "silicate_min"),
    si_max = c("^si_max$", "silicate_max"),
    ice_conc_min = c("siconc_min", "ice.*conc.*min", "sea.*ice.*min"),
    ice_conc_max = c("siconc_max", "ice.*conc.*max", "sea.*ice.*max"),
    salinity_max = c("so_max", "salinity.*max", "sal.*max"),
    salinity_min = c("so_min", "salinity.*min", "sal.*min"),
    wave_min = c("sws_min", "wave.*min", "significant.*wave.*min"),
    wave_max = c("sws_max", "wave.*max", "significant.*wave.*max"),
    depth_min = c("bathymetry_min", "depth_min", "bathy.*min"),
    depth_max = c("bathymetry_max", "depth_max", "bathy.*max"),
    temp_max = c("thetao_max", "temp.*max", "temperature.*max"),
    temp_min = c("thetao_min", "temp.*min", "temperature.*min"),
    light_max = c(".*[Ll]ight.*[Mm]ax", ".*benthic.*light.*max"),
    light_min = c(".*[Ll]ight.*[Mm]in", ".*benthic.*light.*min")
  )
  
  current_names <- names(layers)
  new_names <- current_names
  
  # Function to find best match
  find_match <- function(name, dictionary) {
    clean_name <- str_replace_all(name, "\\.", "_") %>% str_to_lower()
    
    for (new_name in names(dictionary)) {
      patterns <- dictionary[[new_name]]
      for (pattern in patterns) {
        if (str_detect(clean_name, regex(pattern, ignore_case = TRUE))) {
          return(new_name)
        }
      }
    }
    return(NA)
  }
  
  # Apply matching
  for (i in seq_along(current_names)) {
    match <- find_match(current_names[i], dictionary)
    if (!is.na(match)) {
      new_names[i] <- match
    }
  }
  
  # Check for duplicates
  if (any(duplicated(new_names))) {
    duplicates <- new_names[duplicated(new_names)]
    warning("WARNING: Duplicate names detected:\n",
            paste("  -", unique(duplicates), collapse = "\n"),
            "\nCheck your dictionary!")
  }
  
  # Apply new names
  names(layers) <- new_names
  
  # Display transformation
  cat("=== NAME TRANSFORMATION ===\n")
  mapping <- data.frame(
    Index = 1:length(current_names),
    Original = current_names,
    New = new_names,
    Modified = current_names != new_names
  )
  print(mapping, row.names = FALSE)
  
  cat("\n=== SUMMARY ===\n")
  cat("Layers renamed:", sum(current_names != new_names), "/", length(current_names), "\n")
  cat("Layers unchanged:", sum(current_names == new_names), "\n")
  
  return(layers)
}

#' Extract Study Area
#'
#' Crops raster layers to a specified geographic extent.
#'
#' @param layers SpatRaster object to crop
#' @param lon_min Minimum longitude (default: -110)
#' @param lon_max Maximum longitude (default: 30)
#' @param lat_min Minimum latitude (default: 30)
#' @param lat_max Maximum latitude (default: 79)
#' @return Cropped SpatRaster object
#' @export
extract_study_area <- function(layers, 
                               lon_min = -110, lon_max = 30,
                               lat_min = 30, lat_max = 79) {
  
  cat("\n========================================\n")
  cat("STUDY AREA EXTRACTION\n")
  cat("========================================\n\n")
  
  cat(sprintf("Geographic area:\n"))
  cat(sprintf("  Longitudes: %.1f°W to %.1f°E\n", abs(lon_min), lon_max))
  cat(sprintf("  Latitudes: %.1f°N to %.1f°N\n", lat_min, lat_max))
  
  # Create extent
  study_extent <- ext(lon_min, lon_max, lat_min, lat_max)
  
  # Crop layers
  cat("\nCropping layers...\n")
  layers_cropped <- crop(layers, study_extent)
  
  cat(sprintf("  ✓ Reduced dimensions: %d x %d cells\n", 
              nrow(layers_cropped), ncol(layers_cropped)))
  
  return(layers_cropped)
}

#' Calculate Mean Layer
#'
#' Calculates the mean of minimum and maximum raster layers.
#'
#' @param min_raster SpatRaster with minimum values
#' @param max_raster SpatRaster with maximum values
#' @param param_name Name of the parameter (for display purposes)
#' @return SpatRaster with mean values
#' @export
calculate_mean_layers <- function(min_raster, max_raster, param_name = "parameter") {
  cat(sprintf("  Calculating mean (min + max) / 2 for %s...\n", param_name))
  mean_raster <- (min_raster + max_raster) / 2
  return(mean_raster)
}

#' Calculate Global Thresholds
#'
#' Calculates Q3 (75th percentile) and D9 (90th percentile) thresholds
#' for salinity and wave data within marine areas.
#'
#' @param salinity_mean SpatRaster with mean salinity values
#' @param wave_mean SpatRaster with mean wave height values
#' @param marine_mask SpatRaster with marine area mask (1 = marine, 0 = land)
#' @return List with salinity_q3 and wave_d9 threshold values
#' @export
calculate_thresholds <- function(salinity_mean, wave_mean, marine_mask) {
  cat("\nCalculating global thresholds (Q3 and D9) on all marine data...\n")
  
  # Extract marine values only
  sal_values <- values(salinity_mean)[values(marine_mask) == 1]
  wave_values <- values(wave_mean)[values(marine_mask) == 1]
  
  # Remove NAs
  sal_values <- sal_values[!is.na(sal_values)]
  wave_values <- wave_values[!is.na(wave_values)]
  
  # Calculate Q3 (3rd quartile = 75th percentile)
  salinity_q3 <- quantile(sal_values, probs = 0.75, na.rm = TRUE)
  
  # Calculate D9 (9th decile = 90th percentile)
  wave_d9 <- quantile(wave_values, probs = 0.90, na.rm = TRUE)
  
  cat(sprintf("  Salinity Q3 = %.2f\n", salinity_q3))
  cat(sprintf("  Wave D9 = %.4f m\n", wave_d9))
  
  return(list(
    salinity_q3 = salinity_q3,
    wave_d9 = wave_d9
  ))
}

# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# PELAGIC HABITAT HIERARCHICAL CLASSIFICATION FUNCTIONS
# Based on Beaugrand et al. (2019) - North Atlantic
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+

#' Classify Pelagic Habitats
#'
#' Performs hierarchical classification of pelagic habitats following the
#' Beaugrand et al. (2019) methodology for the North Atlantic.
#' 
#' Classification hierarchy:
#' 1. Bathymetric zones (oceanic, shelf-edge, continental shelf)
#' 2. Sea ice presence
#' 3. Ice-free zones
#' 4. Benthic light availability
#' 5. Wave classification
#' 6. Salinity distinction
#' 7. SST stratification
#'
#' @param layers SpatRaster containing all environmental layers
#' @param extract_area Logical, whether to extract study area (default: TRUE)
#' @param lon_min Minimum longitude for study area (default: -110)
#' @param lon_max Maximum longitude for study area (default: 30)
#' @param lat_min Minimum latitude for study area (default: 30)
#' @param lat_max Maximum latitude for study area (default: 79)
#' @param use_fixed_thresholds Logical, use fixed thresholds instead of calculating from data (default: FALSE)
#' @param salinity_threshold Fixed salinity Q3 threshold (default: 35.23)
#' @param wave_threshold Fixed wave D9 threshold in meters (default: 2.0)
#' @return SpatRaster with habitat classifications (1-15)
#' @export
classify_pelagic_habitats <- function(
    layers,
    extract_area = TRUE,
    lon_min = -110, lon_max = 30,
    lat_min = 30, lat_max = 79,
    use_fixed_thresholds = FALSE,
    salinity_threshold = 35.23,
    wave_threshold = 2.0
) {
  
  cat("\n========================================\n")
  cat("HIERARCHICAL CLASSIFICATION OF PELAGIC HABITATS\n")
  cat("Beaugrand et al. (2019) - North Atlantic\n")
  cat("========================================\n\n")
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # STUDY AREA EXTRACTION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  if (extract_area) {
    layers <- extract_study_area(layers, lon_min, lon_max, lat_min, lat_max)
  }
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LAYER VERIFICATION AND EXTRACTION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nVerifying available layers...\n")
  
  required_vars <- c("depth_min", "depth_max", "ice_conc_min", "ice_conc_max",
                     "light_min", "light_max", "wave_min", "wave_max",
                     "salinity_min", "salinity_max", "temp_min", "temp_max")
  
  available_vars <- names(layers)
  missing_vars <- setdiff(required_vars, available_vars)
  
  if (length(missing_vars) > 0) {
    stop(sprintf("Missing variables: %s", paste(missing_vars, collapse=", ")))
  }
  
  cat("  ✓ All required variables present\n\n")
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # STEP 1: CALCULATE MEANS FOR EACH PARAMETER
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("Step 1: Processing min/max layers for each parameter\n")
  
  bathymetry <- calculate_mean_layers(
    layers[["depth_min"]], 
    layers[["depth_max"]], 
    "bathymetry"
  )
  
  sic <- calculate_mean_layers(
    layers[["ice_conc_min"]], 
    layers[["ice_conc_max"]], 
    "sea ice concentration"
  )
  
  light_seabed <- calculate_mean_layers(
    layers[["light_min"]], 
    layers[["light_max"]], 
    "benthic light"
  )
  
  wave <- calculate_mean_layers(
    layers[["wave_min"]], 
    layers[["wave_max"]], 
    "wave height"
  )
  
  salinity <- calculate_mean_layers(
    layers[["salinity_min"]], 
    layers[["salinity_max"]], 
    "salinity"
  )
  
  sst <- calculate_mean_layers(
    layers[["temp_min"]], 
    layers[["temp_max"]], 
    "sea surface temperature"
  )
  
  cat("  ✓ All means calculated\n")
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # STEP 2: CREATE MARINE MASK
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nStep 2: Creating marine mask...\n")
  
  bathy_range <- global(bathymetry, "range", na.rm=TRUE)
  cat(sprintf("  Bathymetry range: [%.1f, %.1f]\n", bathy_range[1], bathy_range[2]))
  
  # Convert to negative if positive
  if (bathy_range[1] >= 0) {
    cat("  → Converting bathymetry to negative values...\n")
    bathymetry <- -abs(bathymetry)
  }
  
  marine_mask <- bathymetry < 0
  
  n_marine_cells <- global(marine_mask, "sum", na.rm=TRUE)$sum
  cat(sprintf("  ✓ %d marine cells identified\n", n_marine_cells))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # STEP 3: CALCULATE Q3 AND D9 THRESHOLDS
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  if (use_fixed_thresholds) {
    cat("\nStep 3: Using fixed thresholds\n")
    cat(sprintf("  Salinity Q3 = %.2f\n", salinity_threshold))
    cat(sprintf("  Wave D9 = %.2f m\n", wave_threshold))
    thresholds <- list(
      salinity_q3 = salinity_threshold,
      wave_d9 = wave_threshold
    )
  } else {
    cat("\nStep 3: Calculating Q3 and D9 thresholds from data\n")
    thresholds <- calculate_thresholds(salinity, wave, marine_mask)
  }
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # STEP 4: INITIALIZE CLASSIFICATION RASTER
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nStep 4: Initializing classification raster...\n")
  
  habitat_class <- bathymetry
  values(habitat_class) <- NA
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # HIERARCHICAL CLASSIFICATION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nStep 5: Hierarchical classification\n")
  cat("========================================\n\n")
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 1: BATHYMETRIC ZONES
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("LEVEL 1: Bathymetric zone definition\n")
  
  oceanic <- (bathymetry <= -1000) & marine_mask
  shelf_edge <- (bathymetry > -1000) & (bathymetry <= -200) & marine_mask
  continental_shelf <- (bathymetry > -200) & marine_mask
  
  cat(sprintf("  - Oceanic zones (>1000m): %.1f%%\n", 
              100 * global(oceanic, "sum", na.rm=TRUE)$sum / n_marine_cells))
  cat(sprintf("  - Shelf-edge (200-1000m): %.1f%%\n", 
              100 * global(shelf_edge, "sum", na.rm=TRUE)$sum / n_marine_cells))
  cat(sprintf("  - Continental shelf (<200m): %.1f%%\n", 
              100 * global(continental_shelf, "sum", na.rm=TRUE)$sum / n_marine_cells))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 2: SEA ICE PRESENCE
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 2: Sea ice identification\n")
  
  # Convert to percentage if needed
  sic_range <- global(sic, "range", na.rm=TRUE)
  if (sic_range[2] <= 1) {
    cat("  → Converting ice concentration to percentage...\n")
    sic <- sic * 100
  }
  
  ice_present <- (sic > 0) & marine_mask
  
  cat(sprintf("  - Zones with ice: %.1f%%\n", 
              100 * global(ice_present, "sum", na.rm=TRUE)$sum / n_marine_cells))
  
  # HABITAT 1: Oceanic with ice
  habitat_1 <- oceanic & ice_present
  habitat_class[habitat_1] <- 1
  cat(sprintf("    → Habitat 1: %d cells\n", global(habitat_1, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 2: Shelf-edge with ice
  habitat_2 <- shelf_edge & ice_present
  habitat_class[habitat_2] <- 2
  cat(sprintf("    → Habitat 2: %d cells\n", global(habitat_2, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 3: Coastal shelf with ice
  habitat_3 <- (bathymetry > -200) & (bathymetry <= -50) & ice_present & marine_mask
  habitat_class[habitat_3] <- 3
  cat(sprintf("    → Habitat 3: %d cells\n", global(habitat_3, "sum", na.rm=TRUE)$sum))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 3: ICE-FREE ZONES
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 3: Ice-free zone classification\n")
  
  no_ice <- (sic == 0) & marine_mask
  
  cat(sprintf("  - Ice-free zones: %.1f%%\n", 
              100 * global(no_ice, "sum", na.rm=TRUE)$sum / n_marine_cells))
  
  # HABITAT 5: Shelf-edge without ice
  habitat_5 <- shelf_edge & no_ice
  habitat_class[habitat_5] <- 5
  cat(sprintf("    → Habitat 5: %d cells\n", global(habitat_5, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 6: Coastal shelf without ice
  habitat_6 <- (bathymetry > -200) & (bathymetry <= -50) & no_ice & marine_mask
  habitat_class[habitat_6] <- 6
  cat(sprintf("    → Habitat 6: %d cells\n", global(habitat_6, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 7: Shallow coastal zones
  habitat_7 <- (bathymetry > -50) & (bathymetry <= 0) & no_ice & marine_mask
  habitat_class[habitat_7] <- 7
  cat(sprintf("    → Habitat 7: %d cells\n", global(habitat_7, "sum", na.rm=TRUE)$sum))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 4: BENTHIC LIGHT
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 4: Benthic photic zone identification\n")
  
  # HABITAT 8: Benthic photic zones
  habitat_8 <- (bathymetry > -200) & (bathymetry <= 0) & no_ice & 
    (light_seabed > 0) & marine_mask
  habitat_class[habitat_8] <- 8
  cat(sprintf("    → Habitat 8: %d cells\n", global(habitat_8, "sum", na.rm=TRUE)$sum))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 5: WAVE CLASSIFICATION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 5: Wave classification (oceanic zones)\n")
  
  oceanic_no_ice <- oceanic & no_ice
  
  # HABITAT 9: Strong waves
  strong_waves <- (wave >= thresholds$wave_d9)
  habitat_9 <- oceanic_no_ice & strong_waves
  habitat_class[habitat_9] <- 9
  
  cat(sprintf("  - High wave zones: %.1f%%\n", 
              100 * global(habitat_9, "sum", na.rm=TRUE)$sum / max(global(oceanic_no_ice, "sum", na.rm=TRUE)$sum, 1)))
  cat(sprintf("    → Habitat 9: %d cells\n", global(habitat_9, "sum", na.rm=TRUE)$sum))
  
  weak_waves <- (wave < thresholds$wave_d9)
  oceanic_weak <- oceanic_no_ice & weak_waves
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 6: SALINITY DISTINCTION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 6: Salinity distinction\n")
  
  # HABITAT 4: Low salinity
  low_salinity <- (salinity < thresholds$salinity_q3)
  habitat_4 <- oceanic_weak & low_salinity
  habitat_class[habitat_4] <- 4
  
  cat(sprintf("  - Freshened zones: %.1f%%\n", 
              100 * global(habitat_4, "sum", na.rm=TRUE)$sum / max(global(oceanic_weak, "sum", na.rm=TRUE)$sum, 1)))
  cat(sprintf("    → Habitat 4: %d cells\n", global(habitat_4, "sum", na.rm=TRUE)$sum))
  
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  # LEVEL 7: SST STRATIFICATION
  # =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  cat("\nLEVEL 7: SST isotherm stratification\n\n")
  
  high_salinity <- (salinity >= thresholds$salinity_q3)
  oceanic_high_sal <- oceanic_weak & high_salinity
  
  # HABITAT 10: SST 7-10°C
  habitat_10 <- oceanic_high_sal & (sst >= 7) & (sst < 10)
  habitat_class[habitat_10] <- 10
  cat(sprintf("  - SST 7-10°C: %d cells\n", global(habitat_10, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 11: SST 10-13°C
  habitat_11 <- oceanic_high_sal & (sst >= 10) & (sst < 13)
  habitat_class[habitat_11] <- 11
  cat(sprintf("  - SST 10-13°C: %d cells\n", global(habitat_11, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 12: SST 13-16°C
  habitat_12 <- oceanic_high_sal & (sst >= 13) & (sst < 16)
  habitat_class[habitat_12] <- 12
  cat(sprintf("  - SST 13-16°C: %d cells\n", global(habitat_12, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 13: SST 16-19°C
  habitat_13 <- oceanic_high_sal & (sst >= 16) & (sst < 19)
  habitat_class[habitat_13] <- 13
  cat(sprintf("  - SST 16-19°C: %d cells\n", global(habitat_13, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 14: SST 19-22°C
  habitat_14 <- oceanic_high_sal & (sst >= 19) & (sst < 22)
  habitat_class[habitat_14] <- 14
  cat(sprintf("  - SST 19-22°C: %d cells\n", global(habitat_14, "sum", na.rm=TRUE)$sum))
  
  # HABITAT 15: SST 22-25°C
  habitat_15 <- oceanic_high_sal & (sst >= 22) & (sst <= 25)
  habitat_class[habitat_15] <- 15
  cat(sprintf("  - SST 22-25°C: %d cells\n", global(habitat_15, "sum", na.rm=TRUE)$sum))
  
  names(habitat_class) <- "pelagic_habitat"
  
  cat("\n========================================\n")
  cat("Classification complete!\n")
  cat("========================================\n\n")
  
  return(habitat_class)
}

#' Calculate Habitat Statistics
#'
#' Computes detailed statistics for each habitat class including cell counts,
#' percentages, and groupings.
#'
#' @param habitat_raster SpatRaster with habitat classifications
#' @return Data frame with habitat statistics
#' @export
habitat_statistics <- function(habitat_raster) {
  
  cat("\n========================================\n")
  cat("HABITAT STATISTICS\n")
  cat("========================================\n\n")
  
  vals <- values(habitat_raster)
  vals <- vals[!is.na(vals)]
  
  if (length(vals) == 0) {
    cat("WARNING: No habitat values found!\n")
    return(NULL)
  }
  
  habitat_counts <- table(vals)
  
  stats_df <- data.frame(
    habitat = as.numeric(names(habitat_counts)),
    n_cells = as.numeric(habitat_counts),
    percentage = round(as.numeric(habitat_counts) / sum(habitat_counts) * 100, 2)
  )
  
  habitat_descriptions <- c(
    "Oceanic with ice (>1000m, SIC>0)",
    "Shelf-edge with ice (200-1000m, SIC>0)",
    "Coastal shelf with ice (50-200m, SIC>0)",
    "Freshened oceanic (>1000m, Sal<Q3)",
    "Shelf-edge ice-free (200-1000m, SIC=0)",
    "Coastal shelf ice-free (50-200m, SIC=0)",
    "Shallow coastal (0-50m, SIC=0)",
    "Benthic photic zone (0-200m)",
    "Dynamic zone (>1000m, Wave≥D9)",
    "Cold waters (SST 7-10°C)",
    "Cool waters (SST 10-13°C)",
    "Cool temperate waters (SST 13-16°C)",
    "Temperate waters (SST 16-19°C)",
    "Warm waters (SST 19-22°C)",
    "Very warm waters (SST 22-25°C)"
  )
  
  stats_df$description <- habitat_descriptions[stats_df$habitat]
  
  group_map <- c(
    "Ice zones", "Ice zones", "Ice zones",
    "Special oceanic", "Continental shelf", "Continental shelf",
    "Coastal", "Photic zone", "Dynamic zones",
    "SST stratification", "SST stratification", "SST stratification",
    "SST stratification", "SST stratification", "SST stratification"
  )
  
  stats_df$group <- group_map[stats_df$habitat]
  
  print(stats_df[, c("habitat", "group", "description", "n_cells", "percentage")])
  
  cat("\n\nGROUP SUMMARY:\n")
  group_summary <- stats_df %>%
    group_by(group) %>%
    summarise(
      n_habitats = n(),
      total_cells = sum(n_cells),
      percentage = sum(percentage)
    ) %>%
    arrange(desc(percentage))
  
  print(group_summary)
  
  return(stats_df)
}

#' Plot Habitats with Projection
#'
#' Creates a projected map of habitat classifications with land overlay.
#' Supports multiple projection systems optimized for the North Atlantic.
#'
#' @param habitat_raster SpatRaster with habitat classifications
#' @param output_file Optional file path to save the plot (PNG format)
#' @param land_shapefile Optional path to land shapefile (uses Natural Earth if NULL)
#' @param projection Projection type: "LAEA", "LAEA_EU", or "EQDC" (default: "LAEA")
#' @param lon_0 Central longitude for projection (default: -40)
#' @param lat_0 Central latitude for projection (default: 55)
#' @return Projected habitat raster
#' @export
plot_habitats_projected <- function(habitat_raster, 
                                    output_file = NULL,
                                    land_shapefile = NULL,
                                    projection = "LAEA",
                                    lon_0 = -40, lat_0 = 55) {
  
  cat("\n========================================\n")
  cat("PROJECTED MAPPING\n")
  cat("========================================\n\n")
  
  # Define projection adapted for North Atlantic
  if (projection == "LAEA") {
    # Lambert Azimuthal Equal Area - ideal for North Atlantic
    target_crs <- sprintf("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                          lat_0, lon_0)
    proj_name <- "Lambert Azimuthal Equal Area"
  } else if (projection == "LAEA_EU") {
    # LAEA Europe (EPSG:3035)
    target_crs <- "EPSG:3035"
    proj_name <- "LAEA Europe (EPSG:3035)"
  } else {
    # Default: equidistant cylindrical
    target_crs <- "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    proj_name <- "Equidistant Cylindrical"
  }
  
  cat(sprintf("Projection used: %s\n", proj_name))
  
  # Reproject raster
  cat("Reprojecting raster...\n")
  habitat_proj <- project(habitat_raster, target_crs, method="near")
  
  # Load or create land shapefile
  if (is.null(land_shapefile)) {
    cat("Loading continents from Natural Earth...\n")
    land <- ne_countries(scale = 10, returnclass = "sf")
  } else {
    cat(sprintf("Loading shapefile: %s\n", land_shapefile))
    land <- st_read(land_shapefile, quiet = TRUE)
  }
  
  # Reproject continents
  cat("Reprojecting continents...\n")
  land_proj <- st_transform(land, target_crs)
  
  # Color palette
  colors <- c(
    "#B3E5FC", "#81D4FA", "#4FC3F7",  # Blues (ice)
    "#FFF9C4",                         # Yellow (freshened)
    "#C8E6C9", "#A5D6A7", "#81C784",  # Greens (shelf)
    "#FFE082",                         # Light orange (photic)
    "#F44336",                         # Red (dynamic)
    "#E1BEE7", "#CE93D8", "#BA68C8",  # Purples (stratification)
    "#FFB74D", "#FFA726", "#FF9800"   # Oranges (warm)
  )
  
  # Open graphics device
  if (!is.null(output_file)) {
    png(output_file, width = 4000, height = 3000, res = 300)
  }
  
  par(mar = c(2, 2, 3, 6), bg = "white")
  
  # Plot raster
  plot(habitat_proj, 
       main = "Pelagic Habitats - North Atlantic\n(Beaugrand et al., 2019 adapted)",
       col = colors,
       legend = FALSE,
       axes = FALSE,
       box = FALSE,
       mar = c(2, 2, 3, 6))
  
  # Add continents
  plot(st_geometry(land_proj), add = TRUE, col = "grey85", border = "grey50", lwd = 0.5)
  
  # Add legend
  legend("left", 
         legend = c(
           "1: Oceanic ice-influenced pelagic habitat",
           "2: Shelf-edges ice-influenced pelagic habitat",
           "3: Continental shelves ice-influenced pelagic habitat",
           "4: Oceanic subarctic pelagic habitat",
           "5: Shelf-edges pelagic habitat",
           "6: Continental shelves deep (50-200m) pelagic habitat",
           "7: Continental shelves shallow (0-50m) pelagic habitat",  
           "8: Continental shelves (light) pelagic habitat",   
           "9: Gulf Stream pelagic habitat",   
           "10: Oceanic subpolar pelagic habitat",   
           "11: Oceanic cold-temperate pelagic habitat",     
           "12: Oceanic temperate pelagic habitat",    
           "13: Oceanic warm-temperate pelagic habitat",  
           "14: Oceanic subtropical (north) pelagic habitat",   
           "15: Oceanic subtropical (south) pelagic habitat"
         ),
         fill = colors,
         cex = 0.7,
         bg = "white",
         box.col = "grey50",
         title = "Habitats")
  
  if (!is.null(output_file)) {
    dev.off()
    cat(sprintf("\n✓ Map saved: %s\n", output_file))
  }
  
  cat("========================================\n\n")
  
  return(habitat_proj)
}


# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# HABITAT RASTER TO SHAPEFILE CONVERSION FUNCTIONS
# =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+

#' Convert Habitat Raster to Shapefile (Simple)
#'
#' Converts a habitat classification raster to a shapefile with basic attributes.
#' Dissolves adjacent pixels with the same habitat value into single polygons.
#'
#' @param habitat_raster SpatRaster with habitat classifications
#' @param output_file Path for output shapefile (default: "habitats.shp")
#' @return sf object with habitat polygons
#' @export
#' @examples
#' \dontrun{
#' habitat_sf <- convert_habitat_to_shapefile_simple(
#'   habitat_raster = my_habitats,
#'   output_file = "output/habitats.shp"
#' )
#' }
convert_habitat_to_shapefile_simple <- function(habitat_raster, 
                                                output_file = "habitats.shp") {
  
  cat("\n========================================\n")
  cat("SIMPLE HABITAT RASTER TO SHAPEFILE CONVERSION\n")
  cat("========================================\n\n")
  
  # Convert raster to polygons (dissolve by value)
  cat("Converting raster to polygons...\n")
  habitat_vect <- as.polygons(habitat_raster, dissolve = TRUE)
  
  # Convert to sf object for better handling
  cat("Converting to sf object...\n")
  habitat_sf <- st_as_sf(habitat_vect)
  
  # Rename the habitat column
  names(habitat_sf)[1] <- "habitat_id"
  
  # Save as shapefile
  cat("Saving shapefile...\n")
  st_write(habitat_sf, output_file, delete_dsn = TRUE, quiet = TRUE)
  
  cat("\n✓ Shapefile saved to:", output_file, "\n")
  cat("✓ Number of habitat types:", nrow(habitat_sf), "\n")
  cat("========================================\n\n")
  
  return(habitat_sf)
}

#' Convert Habitat Raster to Shapefile (Advanced)
#'
#' Converts a habitat classification raster to a detailed shapefile with
#' habitat names, area calculations, and optional geometry simplification.
#'
#' @param habitat_raster SpatRaster with habitat classifications
#' @param habitat_names Named vector or list of habitat descriptions
#' @param output_file Path for output shapefile (default: "habitats_detailed.shp")
#' @param simplify_tolerance Tolerance for geometry simplification in map units (default: 0 = no simplification)
#' @param calculate_area Logical, whether to calculate area in km² (default: TRUE)
#' @return sf object with detailed habitat polygons
#' @export
#' @examples
#' \dontrun{
#' habitat_names <- c(
#'   "1" = "Oceanic with ice",
#'   "2" = "Shelf-edge with ice",
#'   "3" = "Coastal shelf with ice"
#' )
#' 
#' habitat_sf <- convert_habitat_to_shapefile_advanced(
#'   habitat_raster = my_habitats,
#'   habitat_names = habitat_names,
#'   output_file = "output/habitats_detailed.shp",
#'   simplify_tolerance = 100
#' )
#' }
convert_habitat_to_shapefile_advanced <- function(habitat_raster, 
                                                  habitat_names = NULL,
                                                  output_file = "habitats_detailed.shp",
                                                  simplify_tolerance = 0,
                                                  calculate_area = TRUE) {
  
  cat("\n========================================\n")
  cat("ADVANCED HABITAT RASTER TO SHAPEFILE CONVERSION\n")
  cat("========================================\n\n")
  
  # Convert raster to polygons (dissolve by habitat value)
  cat("Converting raster to polygons...\n")
  habitat_vect <- as.polygons(habitat_raster, dissolve = TRUE)
  
  # Convert to sf object
  cat("Converting to sf object...\n")
  habitat_sf <- st_as_sf(habitat_vect)
  
  # Rename primary column
  names(habitat_sf)[1] <- "habitat_id"
  
  # Add habitat names if provided
  if (!is.null(habitat_names)) {
    cat("Adding habitat names...\n")
    habitat_sf$name <- habitat_names[as.character(habitat_sf$habitat_id)]
  }
  
  # Calculate area (in km²)
  if (calculate_area) {
    cat("Calculating areas...\n")
    habitat_sf$area_km2 <- as.numeric(st_area(habitat_sf)) / 1e6
    habitat_sf$area_pct <- round(habitat_sf$area_km2 / sum(habitat_sf$area_km2, na.rm = TRUE) * 100, 2)
  }
  
  # Simplify geometry if requested (reduces file size)
  if (simplify_tolerance > 0) {
    cat(sprintf("Simplifying geometry (tolerance: %.0f)...\n", simplify_tolerance))
    original_size <- object.size(habitat_sf)
    habitat_sf <- st_simplify(habitat_sf, dTolerance = simplify_tolerance)
    new_size <- object.size(habitat_sf)
    cat(sprintf("  Size reduction: %.1f%%\n", 
                100 * (1 - as.numeric(new_size) / as.numeric(original_size))))
  }
  
  # Save as shapefile
  cat("Saving shapefile...\n")
  st_write(habitat_sf, output_file, delete_dsn = TRUE, quiet = TRUE)
  
  cat("\n✓ Advanced shapefile saved to:", output_file, "\n")
  cat("✓ Total features:", nrow(habitat_sf), "\n")
  if (calculate_area) {
    cat(sprintf("✓ Total area: %.2f km²\n", sum(habitat_sf$area_km2, na.rm = TRUE)))
  }
  cat("========================================\n\n")
  
  return(habitat_sf)
}

#' Export Habitats as Separate Shapefiles
#'
#' Exports each habitat type as an individual shapefile for easier
#' manipulation in GIS software.
#'
#' @param habitat_raster SpatRaster with habitat classifications
#' @param output_dir Directory to save individual habitat shapefiles (default: "habitat_shapefiles")
#' @param habitat_names Named vector or list of habitat descriptions
#' @param calculate_area Logical, whether to calculate area in km² (default: TRUE)
#' @return Named list of sf objects, one per habitat type
#' @export
#' @examples
#' \dontrun{
#' habitat_list <- export_habitats_separately(
#'   habitat_raster = my_habitats,
#'   output_dir = "output/individual_habitats",
#'   habitat_names = habitat_names
#' )
#' }
export_habitats_separately <- function(habitat_raster, 
                                       output_dir = "habitat_shapefiles",
                                       habitat_names = NULL,
                                       calculate_area = TRUE) {
  
  cat("\n========================================\n")
  cat("EXPORTING HABITATS AS SEPARATE SHAPEFILES\n")
  cat("========================================\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n\n")
  }
  
  # Get unique habitat values
  habitat_values <- unique(values(habitat_raster))
  habitat_values <- sort(habitat_values[!is.na(habitat_values)])
  
  cat(sprintf("Found %d unique habitat types\n\n", length(habitat_values)))
  
  results <- list()
  
  for (hab_val in habitat_values) {
    cat(sprintf("Processing habitat %d...\n", hab_val))
    
    # Create mask for this habitat
    hab_mask <- habitat_raster == hab_val
    
    # Convert to polygons
    hab_vect <- as.polygons(hab_mask, dissolve = TRUE)
    hab_vect <- hab_vect[hab_vect[[1]] == 1, ]  # Keep only habitat pixels
    
    if (nrow(hab_vect) > 0) {
      # Convert to sf
      hab_sf <- st_as_sf(hab_vect)
      
      # Add attributes
      hab_sf$habitat_id <- hab_val
      
      # Determine filename
      if (!is.null(habitat_names) && as.character(hab_val) %in% names(habitat_names)) {
        hab_sf$name <- habitat_names[as.character(hab_val)]
        safe_name <- gsub("[^A-Za-z0-9_]", "_", habitat_names[as.character(hab_val)])
        filename <- file.path(output_dir, 
                              sprintf("habitat_%02d_%s.shp", hab_val, safe_name))
      } else {
        filename <- file.path(output_dir, sprintf("habitat_%02d.shp", hab_val))
      }
      
      # Calculate area if requested
      if (calculate_area) {
        hab_sf$area_km2 <- as.numeric(st_area(hab_sf)) / 1e6
      }
      
      # Save
      st_write(hab_sf, filename, delete_dsn = TRUE, quiet = TRUE)
      results[[as.character(hab_val)]] <- hab_sf
      
      if (calculate_area) {
        cat(sprintf("  ✓ Exported: %s (%.2f km²)\n", basename(filename), sum(hab_sf$area_km2)))
      } else {
        cat(sprintf("  ✓ Exported: %s\n", basename(filename)))
      }
    } else {
      cat(sprintf("  ⚠ No polygons found for habitat %d\n", hab_val))
    }
  }
  
  cat(sprintf("\n✓ Successfully exported %d habitat shapefiles\n", length(results)))
  cat(sprintf("✓ Output directory: %s\n", output_dir))
  cat("========================================\n\n")
  
  return(results)
}

#' Create Habitat Summary Table
#'
#' Generates a summary table with statistics for all habitat types.
#'
#' @param habitat_sf sf object with habitat polygons (from conversion functions)
#' @param habitat_names Named vector or list of habitat descriptions
#' @param output_csv Optional path to save summary as CSV
#' @return Data frame with habitat summary statistics
#' @export
create_habitat_summary <- function(habitat_sf, 
                                   habitat_names = NULL,
                                   output_csv = NULL) {
  
  cat("\n========================================\n")
  cat("CREATING HABITAT SUMMARY TABLE\n")
  cat("========================================\n\n")
  
  # Prepare summary data frame
  summary_df <- data.frame(
    habitat_id = habitat_sf$habitat_id,
    n_polygons = 1
  )
  
  # Add names if available
  if (!is.null(habitat_names)) {
    summary_df$name <- habitat_names[as.character(summary_df$habitat_id)]
  } else if ("name" %in% names(habitat_sf)) {
    summary_df$name <- habitat_sf$name
  }
  
  # Add area if available
  if ("area_km2" %in% names(habitat_sf)) {
    summary_df$area_km2 <- habitat_sf$area_km2
    summary_df$area_pct <- round(summary_df$area_km2 / sum(summary_df$area_km2, na.rm = TRUE) * 100, 2)
  }
  
  # Sort by habitat_id
  summary_df <- summary_df[order(summary_df$habitat_id), ]
  
  # Display summary
  print(summary_df, row.names = FALSE)
  
  # Save to CSV if requested
  if (!is.null(output_csv)) {
    write.csv(summary_df, output_csv, row.names = FALSE)
    cat(sprintf("\n✓ Summary saved to: %s\n", output_csv))
  }
  
  cat("========================================\n\n")
  
  return(summary_df)
}

#' Batch Convert Multiple Habitat Rasters
#'
#' Converts multiple habitat rasters to shapefiles in batch mode.
#'
#' @param raster_files Character vector of raster file paths
#' @param output_dir Directory to save shapefiles (default: "batch_output")
#' @param method Conversion method: "simple" or "advanced" (default: "simple")
#' @param habitat_names Named vector or list of habitat descriptions
#' @param simplify_tolerance Tolerance for geometry simplification (only for "advanced" method)
#' @return List of sf objects
#' @export
batch_convert_habitats <- function(raster_files,
                                   output_dir = "batch_output",
                                   method = "simple",
                                   habitat_names = NULL,
                                   simplify_tolerance = 0) {
  
  cat("\n========================================\n")
  cat("BATCH CONVERSION OF HABITAT RASTERS\n")
  cat("========================================\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("Processing %d raster files...\n\n", length(raster_files)))
  
  results <- list()
  
  for (i in seq_along(raster_files)) {
    raster_file <- raster_files[i]
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(raster_files), basename(raster_file)))
    
    tryCatch({
      # Load raster
      hab_raster <- rast(raster_file)
      
      # Generate output filename
      base_name <- tools::file_path_sans_ext(basename(raster_file))
      output_file <- file.path(output_dir, paste0(base_name, ".shp"))
      
      # Convert based on method
      if (method == "advanced") {
        hab_sf <- convert_habitat_to_shapefile_advanced(
          habitat_raster = hab_raster,
          habitat_names = habitat_names,
          output_file = output_file,
          simplify_tolerance = simplify_tolerance
        )
      } else {
        hab_sf <- convert_habitat_to_shapefile_simple(
          habitat_raster = hab_raster,
          output_file = output_file
        )
      }
      
      results[[base_name]] <- hab_sf
      cat(sprintf("  ✓ Success: %d features exported\n\n", nrow(hab_sf)))
      
    }, error = function(e) {
      cat(sprintf("  ✗ Error: %s\n\n", e$message))
    })
  }
  
  cat(sprintf("✓ Batch conversion complete: %d/%d files processed\n", 
              length(results), length(raster_files)))
  cat(sprintf("✓ Output directory: %s\n", output_dir))
  cat("========================================\n\n")
  
  return(results)
}
