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

source("script/Bioregion_3D/f_copernicus_download.R") 

copernicus_dir <- here::here("input/copernicus/3D_med_raw/single_year")

# Physical variables
# ------------------

# Download the 0-200m for each year : to loud to download the 24 years simultaneously

product <- "MEDSEA_MULTIYEAR_PHY_006_004"

details <- cms_product_details(product)

layers <- list(c("cmems_mod_med_phy-cur_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-cur_my_4.2km_P1M-m"),
               #c("cmems_mod_med_phy-hflux_my_4.2km_P1D-m"), # no heat flux because daily
               c("cmems_mod_med_phy-mld_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-sal_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-ssh_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-temp_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-wflux_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-wflux_my_4.2km_P1M-m"),
               c("cmems_mod_med_phy-wflux_my_4.2km_P1M-m"))

variables <- list(c("uo"), 
                  c("vo"),
                  #c("rsntds", "hfds", "hfls", "hfss", "rlds"), # couche plus lourde que les autres
                  c("mlotst"), 
                  c("so"), 
                  c("zos"), 
                  c("thetao"),
                  c("pr"),
                  c("evs"),
                  c("friver"))

years <- 2000:2024

elevation <- c(0, -200)

extent <- c(-6, 30, 36, 46) # Med 

# Download per year : calls the copernicus_download function

copernicus_year(product, layers, variables, extent, elevation, years, copernicus_dir)

# Biochemical variables : 1999-2023 plot(test[1,,,1,1])
# ---------------------------------

# Download the 0-200m for each year

product <- "MEDSEA_MULTIYEAR_BGC_006_008"

layers <- list(c("cmems_mod_med_bgc-car_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-car_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-car_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-co2_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-co2_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-nut_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-nut_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-nut_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-plankton_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-plankton_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-bio_my_4.2km_P1M-m"),
               c("cmems_mod_med_bgc-bio_my_4.2km_P1M-m"))

variables <- list(c("dissic"),
                  c("talk"),
                  c("ph"),
                  c("fpco2"),
                  c("spco2"),
                  c("nh4"),
                  c("no3"),
                  c("po4"),
                  c("chl"),
                  c("phyc"),
                  c("o2"),
                  c("nppv")
                  )

elevation <- c(0, -200)

years <- 2000:2024

extent <- c(-6, 30, 36, 46)

# Download per year : calls the copernicus_download function

copernicus_year(product, layers, variables, extent, elevation, years, copernicus_dir)

# Optic, transparency variables : 2000-2024
# -----------------------------------------

# Download NA for plankton groups + only surface product (consider the 2D product already downloaded)

product <- "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"

layers <- list(c("cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M"), c("cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M"),
               c("cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M"))

variables <- list(c("GREEN", "HAPTO", "MICRO", "PROCHLO", "PROKAR"), c("CDM", "BBP"), c("SPM", "ZSD", "KD490"))

# Wind variables : 1999-2024
# --------------------------

# Only surface product (consider the 2D product already downloaded)

product <- "WIND_GLO_PHY_CLIMATE_L4_MY_012_003"

layers <- list(c("cmems_obs-wind_glo_phy_my_l4_P1M"))

variables <- list(c("eastward_wind", "northward_wind","wind_stress_magnitude", "wind_speed"))

# Water turbidity
# ---------------

# Not downloaded : pas dispo sur plus de 5 ans et pas de noms de variables

product <- "OCEANCOLOUR_MED_BGC_HR_L4_NRT_009_211"

layers <- list(c("cmems_obs_oc_med_bgc_geophy_nrt_l4-hr_P1M-m"), c("cmems_obs_oc_med_bgc_optics_nrt_l4-hr_P1M-m"), 
               c("cmems_obs_oc_med_bgc_transp_nrt_l4-hr_P1M-m"))

# Phytoplancton functional group : 
# ------------------------------

# only NA when apply the function on "OCEANCOLOUR_GLO_BGC_L4_MY_009_104" for "cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M"

# Emodnet : bathymetry
# --------------------

# Only surface product (consider the 2D product already downloaded)

test <- rast("E:/Post_Doc_MNHN/R/R_bioregion_pelagique/input/copernicus/2D_med_raw/test/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1779802479377.nc")

