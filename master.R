library(reticulate)
library(rgee, lib.loc = "/projappl/project_2003061/Rpackages/")
library(sf)
library(tidyverse)
library(foreach)
library(doParallel)
library(terra)
library(R.utils)
library(Rsagacmd)
library(randomForestSRC, lib.loc = "/projappl/project_2003061/Rpackages/")
library(zoo)
library(qgam)
library(googledrive)

setwd("/projappl/project_2003061/repos/SnowSensing")

# module load r-env-singularity
# start-r
# library(rgee, lib.loc = "/projappl/project_2003061/Rpackages/")
# ee_Initialize("pekkaniittynen", drive = FALSE)
# ee$Authenticate(auth_mode='notebook')
# ee_install_upgrade()


ee_Initialize("pekkaniittynen", drive = T)
# If not working try:
ee$Authenticate(auth_mode='notebook')
ee$Initialize(project='radiant-planet-412111')
ee$String('Hello from the Earth Engine servers!')$getInfo()

invisible(lapply(list.files("R", pattern = ".R$", full.names = T), source))

workers <- min(c(length(future::availableWorkers()),
                 future::availableCores()))
print(workers)

saga <- saga_gis(raster_backend = "terra",
                 cores = workers)

base_landsat_dir <- "/scratch/project_2003061/temp/"
model_path <- "/scratch/project_2003061/landsat_teach/"

site_name <- "halti"
aoi <- list(X = 21.272014, Y = 69.311827)

st <- Sys.time()
system.time({
  image_df <- extract_landsat(aoi = aoi, 
                              site_name = site_name, 
                              aoi_size = 10, 
                              start_date = "2020-01-01", 
                              end_date = Sys.Date(), 
                              months = 1:12,
                              sats = "LC08", # c("LT04","LT05","LE07","LC08","LC09"), 
                              minclouds = 50, 
                              my_timeout = 400, 
                              info_timeout = 60, 
                              base_landsat_dir = base_landsat_dir,
                              workers = workers)
})

if(!"image_df" %in% ls()){
  image_df <- search_image_df(site_name, base_landsat_dir, workers)
}

# Calculate predictors

system.time({
  calc_predictors(image_df, site_name, base_landsat_dir)
})

# Classify imagery
system.time({
  lss <- classify_landsat(image_df, site_name, base_landsat_dir, model_path, workers, saga)
})
gc()

# Calculate snow stuff
system.time({
  snow_vars <- calc_snow_variables(image_df, site_name, base_landsat_dir, workers)
})
# plot(snow_vars$scd)
writeRaster(snow_vars, paste0(base_landsat_dir,"/",site_name,"/","snow_variables.tif"))
Sys.time() - st

layerCor(snow_vars, fun = "pearson", na.rm=TRUE)
names(snow_vars)
