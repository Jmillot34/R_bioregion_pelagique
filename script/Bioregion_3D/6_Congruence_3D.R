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
library(usdm)
library(NbClust)
library(clue)
library(maptree)
library(purrr)
library(readr)
library(RColorBrewer)

# ------------------------
# 0. Elements for mapping
# ------------------------

coastline_med_nw <- readRDS(here("input", "srm_elements", "coastline_med_nw.rds"))

ocean_med_nw <- readRDS(here("input", "srm_elements", "ocean_med_nw.rds"))

srm_med_out <- readRDS(here("input", "srm_elements", "srm_med_out.rds"))

# -----------------------
# 0. Choose the scenario
# -----------------------

scen <- "srm_0-200_month"

month_list <- c("jan", "feb", "mar","apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

for (m in 1:12) {

month <- month_list[m]

message("=== MONTH: ", month, " ===")

# Output directory

boundaries_dir <- here("output", "congruence", scen, month, "boundaries")
dir.create(boundaries_dir, recursive = TRUE, showWarnings = FALSE)

map_boundaries_dir <- here("output", "congruence", scen, month, "map_boundaries")
dir.create(map_boundaries_dir, recursive = TRUE, showWarnings = FALSE)

congru_smooth_raster_dir <- here("output", "congruence", scen, month, "congru_smooth_raster")
dir.create(congru_smooth_raster_dir, recursive = TRUE, showWarnings = FALSE)

congru_smooth_map_dir <- here("output", "congruence", scen, month, "congru_smooth_map")
dir.create(congru_smooth_map_dir, recursive = TRUE, showWarnings = FALSE)

# Path to raster directory

rast_dir <- here("output/Congruence", scen, month, "rasters")

# List raster files

rast_files <- list.files(
  rast_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Raster names

rast_names <-  tools::file_path_sans_ext(basename(rast_files))

# -----------------------
# 1. Identify boundaries
# -----------------------

# Compute boundaries layers

for (i in 1:length(rast_files)) {
  
  p <- rast(rast_files[i])
  
  # Neighbors assignment
  
  adj <- adjacent(p, cells = 1:ncell(p), directions = 8, pairs = TRUE)
  
  # Only select the pixel with neighbors associated with several values
  
  boundaries_cells <- unique(adj[ p[adj[,1]] != p[adj[,2]], 1 ])
  
  # Build raster
  
  boundaries <- p 
  
  values(boundaries) <- ifelse(is.na(values(p)), NA, 0)
  
  values(boundaries)[boundaries_cells] <- 1
  
  # Mask on the Med only
  
  boundaries_rast <- boundaries %>% mask(p)
  
  # Save
  
  filename_boundaries <- file.path(boundaries_dir, paste0(rast_names[[i]], "_", month, ".tif"))
  
  # Write raster
  
  writeRaster(boundaries_rast, filename_boundaries, overwrite = TRUE)
  
}

# Function to map boundaries

map_boundaries_function <- function(partition, name) {
  
    #Graticules
    tm_graticules(
      labels.inside.frame = FALSE, 
      lines       = TRUE,
      col         = "grey60",
      lwd         = 0.5,
      labels.size = 2
    ) +
    
    # Main layer: congruence raster
    tm_shape(partition, is.main = TRUE) +
    tm_raster(
      col.scale  = tm_scale_categorical(values = c("black","salmon")),
      col.legend = tm_legend(title = "Boundaries")
    ) +
    
    # Coastline med nw
    tm_shape(coastline_med_nw) + 
    
    tm_polygons(
      fill = "grey94",
      col  = "black"
    ) + 
    
    tm_borders(col = "black", lwd = 1.5) +
    
    # Layout & legend
    tm_layout(
      legend.show             = FALSE,
      legend.outside          = FALSE,
      legend.position = c("right","top"),
      legend.stack = "horizontal",
      legend.title.size       = 2,
      legend.text.size        = 2,
      legend.frame            = FALSE,
      inner.margins           = 0,
      legend.width            = 5,   
      legend.height           = 3.5,  
      
      # Add main title above the map
      main.title              = name,
      main.title.size         = 3,
      main.title.position     = "center"
    )
}

# List boundary files

boundaries_files <- list.files(
  boundaries_dir,
  full.names = TRUE
)

# Read boundary data

boundaries_list <- lapply(boundaries_files, rast)

# Assign names based on file names

names(boundaries_list) <- lapply(boundaries_files, function(x) tools::file_path_sans_ext(basename(x)))

# Vector of the number of clusters k associated to each partition

k_list <- as.numeric(stringr::str_extract(names(boundaries_list), "(?<=_k|km_)\\d+"))

# Boundary maps

map_boundaries_list <- lapply(seq_along(boundaries_list), function(i) {
  map_boundaries_function(
    partition = boundaries_list[[i]],
    name      = names(boundaries_list)[[i]])
})

names(map_boundaries_list) <- names(boundaries_list)

# Saving

for (nm in names(map_boundaries_list)) {
  tmap_save(
    map_boundaries_list[[nm]],
    filename = here("output","congruence", scen, month, "map_boundaries", paste0("map_frontier_", nm, ".png")),
    width = 27,
    height = 15,
    dpi = 600
  )
}

###########################
# 2. Calculate congruence #
###########################

# Weighted each partition by 1/k : boundaries in a partition of 22 should contribute less proportionally 
boundaries_weight_list <- lapply(seq_along(boundaries_list), function(i) {
  boundaries_list[[i]] * (1 / k_list[i])
})

# Max value of congru_wh for normalization
r_max <- global(sum(rast(boundaries_weight_list)), "max", na.rm = TRUE)[1,1] # max value

# Sum with ponderation + normalization between 0 and 1
congru_pond <- sum(rast(boundaries_weight_list))/ r_max

# Sum without ponderation
congru <- sum(rast(boundaries_list))

# Convertir le contour en SpatVector terra (utile pour le masquage raster)
coast_vect <- vect(coastline_med_nw)

####################
# Lissage Gaussien #
####################

# -- Augmentation de résolution pour le lissage --
# Cible : ~0.017° (~1.5 km) — upscale x6 depuis 0.1°

congruence_i_list <- lapply(list(congru, congru_pond), function(r) {
  disagg(r, fact = 6, method = "bilinear")
})

# -- Fonction de lissage gaussien --
# [FIX 4] Normalisation du kernel identique à l'original :
#   soustraction du min puis division par le max (plage [0, 1]).
#   focal() utilise fun = "mean" (et non "sum") comme dans l'original.
gaussian_smooth_terra <- function(x, radius_deg) {
  cell_size <- res(x)[1]
  diameter  <- round(radius_deg * 2 / cell_size)
  if (diameter %% 2 == 0) diameter <- diameter + 1   # doit être impair
  
  # Kernel gaussien normalisé [0, 1]
  dn  <- dnorm(seq(-1, 1, length.out = diameter), mean = 0, sd = 0.8)
  mat <- outer(dn, dn)
  mat <- mat - min(mat)    # [FIX 4] soustraction du min
  mat <- mat / max(mat)    # [FIX 4] division par le max → [0, 1]
  
  focal(x, w = mat, fun = "mean", na.rm = TRUE)   # [FIX 4] fun = "mean"
}

# Lissage au rayon 0.5

r <- 0.5 # c(0.45, 0.50, 0.55)

cat("Smoothing congruence surface at radii:", r, "degrees...\n")

smoothed_congru_list <- lapply(congruence_i_list, function(l) {
  
  sm <- gaussian_smooth_terra(l, r) %>% 
    mask(coast_vect, inverse = TRUE) # Mask land
  
  names(sm) <- "congru"
  
  return(sm)
  
})

names(smoothed_congru_list) <- c("congruence", "congruence_pond")

# Save

writeRaster(smoothed_congru_list[["congruence"]],
            file.path(congru_smooth_raster_dir, "congruence_smooth.tif"),
            overwrite = TRUE)

writeRaster(smoothed_congru_list[["congruence_pond"]],
            file.path(congru_smooth_raster_dir, "congruence_smooth_pond.tif"),
            overwrite = TRUE)

cat("Smoothed rasters saved.\n")

# Read congruence

congru_smooth <- rast(file.path(congru_smooth_raster_dir, "congruence_smooth.tif"))
congru_smooth_pond <- rast(file.path(congru_smooth_raster_dir, "congruence_smooth_pond.tif"))

####################
## Mapping #########
####################

# Elements

grad <- brewer.pal(n = 9, name = "Blues")

pal <- rev(grad[-1])

# -- Map : Surface de congruence lissée --

map_congru_smooth <- 
  
  # Congruence raster
  
  tm_shape(congru_smooth, is.main = TRUE) +
  tm_raster(
    style = "cont",
    palette = pal, 
    col.na = "grey",
    title = "Congruence",
  ) +
  
  # Subregion out 
  tm_shape(srm_med_out) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "grey80",
    alpha = 0.3
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Coastline med nw
  tm_shape(coastline_med_nw) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "black"
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Graticules
  
  tm_graticules(
    labels.inside.frame = FALSE, 
    lines       = FALSE,
    col         = "grey60",
    lwd         = 0.5,
    labels.size = 1.7
  ) +
  
  # Layout & legend
  
  tm_layout(
    legend.show             = TRUE,
    legend.outside          = TRUE,
    legend.outside.position = c("right","top"),
    legend.stack = "horizontal",
    legend.title.size       = 1.7, #2
    legend.text.size        = 1.7, #2
    legend.format = list(fun = function(x) floor(x)),
    legend.frame            = FALSE,
    inner.margins           = 0,
    legend.width            = 6,  #4 
    legend.height           = 20,  #6
    
    # Add main title above the map
    main.title = paste0("Congruence - ", month),
    main.title.size = 1,
    main.title.position     = "center")

map_congru_smooth

tmap_save(map_congru_smooth, file.path(congru_smooth_map_dir, paste0("map_congru_smooth_", month, ".tiff")), width = 14, height = 10, dpi = 600)

# -- Map : Surface de congruence catégorielle --

map_congru_cat <- 
  
  # Congruence raster
  
  tm_shape(congru_smooth, is.main = TRUE) +
  
  tm_raster(
    style = "pretty",
    palette = pal, 
    col.na = "grey",
    title = "Congruence",
  ) +
  
  # Subregion out 
  
  tm_shape(srm_med_out) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "grey80",
    alpha = 0.3
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Coastline med nw
  
  tm_shape(coastline_med_nw) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "black"
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Graticules
  
  tm_graticules(
    labels.inside.frame = FALSE, 
    lines       = FALSE,
    col         = "grey60",
    lwd         = 0.5,
    labels.size = 1.7
  ) +
  
  # Layout & legend
  
  tm_layout(
    legend.show             = TRUE,
    legend.outside          = TRUE,
    legend.outside.position = c("right","top"),
    legend.stack = "horizontal",
    legend.title.size       = 1.7, #2
    legend.text.size        = 1.7, #2
    legend.format = list(fun = function(x) floor(x)),
    legend.frame            = FALSE,
    inner.margins           = 0,
    legend.width            = 6,  #4 
    legend.height           = 20,  #6
    
    # Add main title above the map
    main.title = paste0("Congruence - ", month),
    main.title.size = 1,
    main.title.position     = "center")

map_congru_cat

tmap_save(map_congru_cat, file.path(congru_smooth_map_dir, paste0("map_congru_cat_", month, ".tiff")), width = 14, height = 10, dpi = 600)

# -- Map : Surface de congruence lissée et pondérée --

map_congru_smooth_pond <- 
  
  # Congruence raster
  
  tm_shape(congru_smooth_pond, is.main = TRUE) +
  tm_raster(
    style = "cont",
    palette = pal, 
    col.na = "grey",
    title = "Congruence",
  ) +
  
  # Subregion out 
  tm_shape(srm_med_out) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "grey80",
    alpha = 0.3
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Coastline med nw
  tm_shape(coastline_med_nw) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "black"
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Graticules
  
  tm_graticules(
    labels.inside.frame = FALSE, 
    lines       = FALSE,
    col         = "grey60",
    lwd         = 0.5,
    labels.size = 1.7
  ) +
  
  # Layout & legend
  
  tm_layout(
    legend.show             = TRUE,
    legend.outside          = TRUE,
    legend.outside.position = c("right","top"),
    legend.stack = "horizontal",
    legend.title.size       = 1.7, #2
    legend.text.size        = 1.7, #2
    legend.format = list(fun = function(x) floor(x)),
    legend.frame            = FALSE,
    inner.margins           = 0,
    legend.width            = 6,  #4 
    legend.height           = 20, #6  
    
    # Add main title above the map
    main.title = paste0("Congruence pondérée - ", month),
    main.title.size = 1,
    main.title.position     = "center")

map_congru_smooth_pond

tmap_save(map_congru_smooth_pond, file.path(congru_smooth_map_dir, paste0("map_congru_smooth_pond_", month, ".tiff")), width = 14, height = 10, dpi = 600)

# -- Map : Surface de congruence catégorielle et pondérée --

map_congru_cat_pond <- 

  # Congruence raster
  
  tm_shape(congru_smooth_pond, is.main = TRUE) +
  tm_raster(
    style = "pretty",
    palette = pal, 
    col.na = "grey",
    title = "Congruence",
  ) +
  
  # Subregion out 
  tm_shape(srm_med_out) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "grey80",
    alpha = 0.3
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Coastline med nw
  tm_shape(coastline_med_nw) + 
  
  tm_polygons(
    fill = "grey94",
    col  = "black"
  ) + 
  
  tm_borders(col = "black", lwd = 1) +
  
  # Graticules
  
  tm_graticules(
    labels.inside.frame = FALSE, 
    lines       = FALSE,
    col         = "grey60",
    lwd         = 0.5,
    labels.size = 1.7
  ) +
  
  # Layout & legend
  
  tm_layout(
    main.title = paste0("Congruence pondérée - ", month),
    main.title.size = 1,
    legend.show             = TRUE,
    legend.outside          = TRUE,
    legend.outside.position = c("right","top"),
    legend.stack = "horizontal",
    legend.title.size       = 1.7, #2
    legend.text.size        = 1.7, #2
    legend.format = list(fun = function(x) floor(x)),
    legend.frame            = FALSE,
    inner.margins           = 0,
    legend.width            = 6,  #4 
    legend.height           = 20,  #6
  )

map_congru_cat_pond

tmap_save(map_congru_cat_pond, file.path(congru_smooth_map_dir, paste0("map_congru_cat_pond_", month, ".tiff")), width = 14, height = 10, dpi = 600)

}