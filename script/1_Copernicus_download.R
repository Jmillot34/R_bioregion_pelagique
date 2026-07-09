# Download Copernicus Data

# Website : https://github.com/pepijn-devries/CopernicusMarine

# Libraries

library(terra)
library(dplyr)
library(here)
library(sf)
library(stars)
library(CopernicusMarine)
library(blosc)

source("script/f_copernicus_download.R")

copernicus_dir <- here::here("input/copernicus/2D_med_raw")

# Physical variables
# ------------------

product <- "MEDSEA_MULTIYEAR_PHY_006_004"

details <- cms_product_details(product)

layers <- list(c("cmems_mod_med_phy-cur_my_4.2km_P1M-m"),
               #c("cmems_mod_med_phy-hflux_my_4.2km_P1D-m"), # no heat flux because daily
               c("cmems_mod_med_phy-mld_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-sal_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-ssh_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-temp_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-wflux_my_4.2km_P1M-m"))

variables <- list(c("uo", "vo"), 
                  #c("rsntds", "hfds", "hfls", "hfss", "rlds"), # couche plus lourde que les autres
                  c("mlotst"), 
                  c("so"), 
                  c("zos"), 
                  c("thetao"),
                  c("pr","evs", "friver"))

extent <- c(-6, 30, 36, 46)

elevation <- c(0, -2)

time_range <- c("2000-01-01 UTC", "2024-12-01 UTC")

copernicus_download(product, layers, variables, extent, elevation, time_range, copernicus_dir)

# Biochemical variables : 1999-2023
# ---------------------------------

product <- "MEDSEA_MULTIYEAR_BGC_006_008"

layers <- list(c("med-ogs-car-rean-m"),
               c("med-ogs-co2-rean-m"),
               c("med-ogs-nut-rean-m"),
               c("med-ogs-pft-rean-m"),
               c("med-ogs-bio-rean-m"))

variables <- list(c("dissic", "talk", "ph"),
                  c("fpco2","spco2"),
                  c("nh4", "no3", "po4"),
                  c("chl","phyc"),
                  c("o2", "nppv"))

extent <- c(-6, 30, 36, 46)

elevation <- c(0, -2)

time_range <- c("2000-01-01 UTC", "2024-12-01 UTC")

copernicus_download(product, layers, variables, extent, time_range, copernicus_dir)

# Optic, transparency variables : 2000-2024
# -----------------------------------------

# Download NA for plankton groups + only surface product (prendre le 2D déjà dowloaded)

product <- "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"

layers <- list(c("cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M"), c("cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M"),
               c("cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M"))

variables <- list(c("GREEN", "HAPTO", "MICRO", "PROCHLO", "PROKAR"), c("CDM", "BBP"), c("SPM", "ZSD", "KD490"))

extent <- c(-6, 30, 36, 46)

elevation <- c(0, -200)

years <- 2000:2024

copernicus_download(product, layers, variables, extent, elevation, years, copernicus_dir)

# Wind variables : 1999-2024
# --------------------------

# The function doesn't work = "indice hors limite" + only surface product (prendre le 2D déjà dowloaded)

product <- "WIND_GLO_PHY_CLIMATE_L4_MY_012_003"

layers <- list(c("cmems_obs-wind_glo_phy_my_l4_P1M"))

variables <- list(c("eastward_wind", "northward_wind","wind_stress_magnitude", "wind_speed"))

elevation <- c(0, -200)

years <- 2000:2024

# Download per year : calls the copernicus_download function

copernicus_download(product, layers, variables, extent, elevation, years, copernicus_dir)

# Water turbidity
# ---------------

# Not downloaded : pas dispo sur plus de 5 ans et pas de nom de variables

product <- "OCEANCOLOUR_MED_BGC_HR_L4_NRT_009_211"

layers <- list(c("cmems_obs_oc_med_bgc_geophy_nrt_l4-hr_P1M-m"), c("cmems_obs_oc_med_bgc_optics_nrt_l4-hr_P1M-m"), 
               c("cmems_obs_oc_med_bgc_transp_nrt_l4-hr_P1M-m"))

# Variables which cannot be downloaded automatically
# --------------------------------------------------

## List the files downloaded from the Copernicus platform
wind_files <- list.files(path = here::here("input/copernicus/wind_nc"),
                         pattern = "\\.(nc)$", 
                         full.names = TRUE)

## Read files and set dimensions
wind_list <- lapply(wind_files, read_stars)
wind_str <- do.call(c, wind_list) %>% 
  st_set_dimensions(., which = "x", names = "longitude") %>% 
  st_set_dimensions(., which = "y", names = "latitude") %>% 
  st_set_crs(., "+proj=longlat +datum=WGS84")
names(wind_str) <- gsub("\\.nc", "", names(wind_str))

## Save
wind_filename_rds <- file.path(here::here("input/copernicus/med_wind.rds"))
saveRDS(wind_str, wind_filename_rds)

# Phytoplancton functional group : only NA when apply the function on "OCEANCOLOUR_GLO_BGC_L4_MY_009_104" for "cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M"
# ------------------------------

## List the files downloaded from the Copernicus platform
pft_files <- list.files(path = here::here("Input/copernicus/phyto_ft_nc"),
                        pattern = "\\.(nc)$", 
                        full.names = TRUE)

test <- read_stars(here::here("Input/copernicus/phyto_ft_nc/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_DIATO.nc"))
test[1,,,1]
plot(test[1,,,1])

## Read files and set dimensions
pft_list <- lapply(pft_files, read_stars)
pft_str <- do.call(c, pft_list) %>% 
  st_set_dimensions(., which = "x", names = "longitude") %>% 
  st_set_dimensions(., which = "y", names = "latitude") %>% 
  st_set_crs(., "+proj=longlat +datum=WGS84")

## Save
pft_filename_rds <- file.path(here::here("Input/copernicus/med_wind.rds"))
saveRDS(wind_str, wind_filename_rds)

# Emodnet : bathymetry
# --------------------

# List of files
emodnet_tiles <- list.files(path = here::here("Input/emodnet/emodnet_2024/"), 
                            pattern = "\\.(nc)$", 
                            full.names = TRUE)

# Transformation in raster
emodnet_rast <- lapply(emodnet_tiles, rast)
emodnet_elevation <- lapply(emodnet_rast, function(x) subset(x, "elevation"))

# Merge in one raster
emodnet_merge <- do.call(mosaic, emodnet_elevation)
names(emodnet_merge) <- "bathymetry"

writeRaster(emodnet_merge, here::here("Input/emodnet/emodnet_2024/emodnet_merge.tif"))

# Not complet products
# --------------------

# Biochemical variables: Global Ocean (existe en Med, res 1km mais daily)

product <- "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"
layers <- list(c("cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M"))
variables <- list(c("CDM", "BBP"))
copernicus_download(product, layers, variables)

product <- "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"
layers <- list(c("cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M"))
variables <- list(c("SPM", "ZSD", "KD490"))
copernicus_download(product, layers, variables)

# Biochemical variables : 2022-2025 durée trop courtes

product <- "MEDSEA_ANALYSISFORECAST_BGC_006_014"
layers <- list(c("cmems_mod_med_bgc-car_anfc_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-co2_anfc_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-nut_anfc_4.2km_P1M-m"), # Si new
               c("cmems_mod_med_bgc-optics_anfc_4.2km_P1M-m"), # New
               c("cmems_mod_med_bgc-pft_anfc_4.2km_P1M-m"), # New variables
               c("cmems_mod_med_bgc-bio_anfc_4.2km_P1M-m"))

variables <- list(c("dissic", "talk", "ph"), 
                  c("fpco2","spco2"), 
                  c("nh4", "no3", "po4", "si"),
                  c("kd490"), 
                  c("chl", "diatoChla", "dinoChla", "nanoChla", "picoChla", "diatoC", "dinoC", "nanoC", 
                    "phyc", "picoC", "zooc"), 
                  c("o2", "nppv"))

copernicus_download(product, layers, variables)

# Wave : en hourly seulement
product <- "MEDSEA_MULTIYEAR_WAV_006_012"
layers <- list(c("cmems_mod_med_wav_my_4.2km-climatology_P1M-m"))
variables <- list(c("VTM02", "VHM0"))

copernicus_download(product, layers, variables)


# Plot d'ammonium

# Hiver
test <- readRDS(here("input", "copernicus", "3D_med_raw", "single_year", "nh4_2000_2000.rds"))
r <- test["nh4",,,6,12]

elev <- st_get_dimension_values(test, "elevation")

tm_shape(r) +
  tm_raster(
    style = "quantile",
    n = 5,
    palette = viridis(5),
    title = "NH4 (mmol.m-3)"
  ) +
  tm_layout(
    legend.outside = TRUE
  )

# Eté
plot(test["nh4",,,1,6],
      col = viridis(100),
      main = "NH4 concentration (mmol.m-3)")



