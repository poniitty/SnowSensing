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
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires, lib.loc = "/projappl/project_2003061/Rpackages/")

setwd("/projappl/project_2003061/repos/SnowSensing")

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
results_path <- "/scratch/project_2003061/snow_iceland/"

world <- ne_countries(scale = "large", returnclass = "sf")

ice <- world %>% filter(admin == "Iceland") %>% st_transform(32627)

plot(st_geometry(ice))

icegrid <- ice %>% 
  st_make_grid(cellsize = 20000) %>% 
  st_geometry() %>% 
  st_as_sf

icegrid <- icegrid[ice,]

plot(st_geometry(icegrid))

icegrid <- icegrid %>% 
  mutate(id = paste0("iceland_", 1:nrow(.))) %>% 
  st_buffer(100)

for(site_name in icegrid$id){
  if(!file.exists(paste0(results_path,"/",site_name,".tif")) & 
     !file.exists(paste0(base_landsat_dir, "/", site_name, "/done.csv"))){
    # site_name <- "iceland_135"
    print(site_name)
    aoi <- icegrid %>% filter(id == site_name)
    
    st <- Sys.time()
    image_df <- extract_landsat(aoi = aoi, 
                                site_name = site_name, 
                                aoi_size = 10, 
                                start_date = "2014-01-01", 
                                end_date = "2023-12-31", 
                                months = 1:12,
                                sats = "LC08", # c("LT04","LT05","LE07","LC08","LC09"), 
                                minclouds = 50, 
                                my_timeout = 400, 
                                info_timeout = 60, 
                                base_landsat_dir = base_landsat_dir,
                                workers = workers)
    
    if(!"image_df" %in% ls()){
      image_df <- search_image_df(site_name, base_landsat_dir, workers)
    }
    
    write_csv(tibble(yeah = "done"), 
              paste0(base_landsat_dir, "/", site_name, "/done.csv"))
    
    # Calculate predictors
    calc_predictors(image_df, site_name, base_landsat_dir)
    gc()

    # # Classify imagery
    # lss <- classify_landsat(image_df, site_name, base_landsat_dir, model_path, workers, saga)
    # gc()
    # 
    # # Calculate snow stuff
    # snow_vars <- calc_snow_variables(image_df, site_name, base_landsat_dir, workers)
    # gc()
    # # plot(snow_vars$scd)
    # writeRaster(snow_vars, paste0(results_path,"/",site_name,".tif"))
    Sys.time() - st
    
    # unlink(paste0(base_landsat_dir, "/", site_name, "/"), recursive=TRUE)
  }
}

for(site_name in icegrid$id){
  if(!file.exists(paste0(results_path,"/",site_name,".tif")) & 
     file.exists(paste0(base_landsat_dir, "/", site_name, "/done.csv")) & 
     !file.exists(paste0(base_landsat_dir, "/", site_name, "/started.csv"))){
    # site_name <- "iceland_2"
    print(site_name)
    
    write_csv(tibble(yeah = "started"), 
              paste0(base_landsat_dir, "/", site_name, "/started.csv"))
    
    aoi <- icegrid %>% filter(id == site_name)
    
    st <- Sys.time()
    # image_df <- extract_landsat(aoi = aoi, 
    #                             site_name = site_name, 
    #                             aoi_size = 10, 
    #                             start_date = "2014-01-01", 
    #                             end_date = "2023-12-31", 
    #                             months = 1:12,
    #                             sats = "LC08", # c("LT04","LT05","LE07","LC08","LC09"), 
    #                             minclouds = 50, 
    #                             my_timeout = 400, 
    #                             info_timeout = 60, 
    #                             base_landsat_dir = base_landsat_dir,
    #                             workers = workers)
    
    if(!"image_df" %in% ls()){
      image_df <- search_image_df(site_name, base_landsat_dir, workers)
    }
    
    # # Calculate predictors
    # calc_predictors(image_df, site_name, base_landsat_dir)
    # gc()
    
    # Classify imagery
    lss <- classify_landsat(image_df, site_name, base_landsat_dir, model_path, workers, saga)
    
    # Calculate snow stuff
    snow_vars <- calc_snow_variables(image_df, site_name, base_landsat_dir, workers)
    
    # plot(snow_vars$scd)
    writeRaster(snow_vars, paste0(results_path,"/",site_name,".tif"))
    Sys.time() - st
    
    unlink(paste0(base_landsat_dir, "/", site_name, "/"), recursive=TRUE)
    rm(image_df)
  }
}


