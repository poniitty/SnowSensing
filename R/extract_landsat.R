extract_landsat <- function(aoi, site_name, aoi_size = 5, start_date = "2014-01-01", end_date = Sys.Date(), months = 1:12,
                            sats = "LC08", minclouds = 50, my_timeout = 400, info_timeout = 60, base_landsat_dir, workers = workers){
  
  utmall <- st_read("data/utm_zones.gpkg") %>% 
    filter(ZONE != 0)
  
  if(class(aoi)[1] == "sf"){
    aoi_mid <- aoi %>% 
      st_centroid() %>% 
      st_transform(crs = 4326)
  }
  if(class(aoi)[1] == "list"){
    aoi_mid <- as_tibble(aoi) %>% 
      mutate(name = site_name) %>% 
      st_as_sf(coords = c("X", "Y"), crs = 4326)
    aoi <- as_tibble(aoi) %>% 
      mutate(name = site_name) %>% 
      st_as_sf(coords = c("X", "Y"), crs = 4326)
  }
  
  # WGS84 UTM zones to set the correct projection
  utm <- utmall[aoi_mid,] # Which zone the study points falls in
  
  if(nrow(utm) == 0){
    utm <- utmall[st_nearest_feature(aoi_mid, utmall),]
  }
  
  lat <- st_coordinates(aoi_mid)[,"Y"]
  utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0",utm$ZONE), utm$ZONE)
  epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
  
  # From point to polygon (20km x 20km)
  
  if(st_geometry_type(aoi) == "POLYGON"){
    aoi <- aoi %>% st_transform(crs = epsg) %>% 
      mutate(name = site_name)
  } else {
    aoi <- aoi %>% st_transform(crs = epsg) %>% 
      st_buffer(aoi_size*1000) %>% st_bbox() %>% 
      st_as_sfc() %>% st_as_sf() %>% 
      mutate(name = site_name)
  }
  
  # AOI polygon to GEE format
  aoi_ee <- aoi %>% st_transform(crs = 4326) %>% 
    st_geometry() %>% 
    sf_as_ee()
  
  # Create sub-directory where imagery will be saved if does not exists
  drive_folder <- paste0("TEMP_Landsat_",site_name)
  area_landsat_dir <- paste0(base_landsat_dir,"/",site_name)
  if(!dir.exists(area_landsat_dir)){
    dir.create(area_landsat_dir)
  }
  
  # # List existing files
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  tif_dates <- tibble(satid = unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][2])),
                      date = as.character(ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][5])))))
  if(nrow(tif_dates) == 0){
    tif_dates <- tibble(satid = NA,
                        date = NA)
  }
  
  lss <- extract_landsat_gee(workers = workers, start_date = start_date, end_date = end_date, months = months,
                             aoi_ee = aoi_ee, site_name = site_name, epsg = epsg, excl_dates = tif_dates,
                             sats = sats, my_timeout = my_timeout, info_timeout = info_timeout, 
                             area_landsat_dir = area_landsat_dir, minclouds = minclouds, drive_folder = drive_folder)
  
  if(file.exists(paste0(area_landsat_dir, "/lss.csv"))){
    lss <- bind_rows(read_csv(paste0(area_landsat_dir, "/lss.csv")) %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED)),
                     lss) %>% 
      distinct()
  }
  write_csv(lss, paste0(area_landsat_dir, "/lss.csv"))
  
  # List downloaded images
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  
  lss$file <- lapply(lss$LANDSAT_PRODUCT_ID, function(x){
    xx <- tifs[grepl(x, tifs)]
    xx <- ifelse(length(xx) == 0, NA, xx)
    return(xx)
  }) %>% unlist
  
  lss$area <- site_name
  lss$collection <- gsub("0","C",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
  lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
  lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
  lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),1,3))
  lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),4,6))
  lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
  lss$time <- gsub(".tif","",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
  
  lss <- lss %>% arrange(date)
  
  # Check if all rasters are working. Remove corrupted files.
  img_remove <- unlist(mclapply(lss$file, check_raster, image_dir = area_landsat_dir, mc.cores = workers))
  
  if(length(img_remove) > 0){
    unlink(paste0(area_landsat_dir,"/",img_remove))
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
    tif_dates <- tibble(satid = unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][2])),
                        date = as.character(ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][5])))))
    if(nrow(tif_dates) == 0){
      tif_dates <- tibble(satid = NA,
                          date = NA)
    }
    lss <- extract_landsat_gee(workers = workers, start_date = start_date, end_date = end_date, months = months,
                               aoi_ee = aoi_ee, site_name = site_name, epsg = epsg, excl_dates = tif_dates,
                               sats = sats, my_timeout = my_timeout, info_timeout = info_timeout, 
                               area_landsat_dir = area_landsat_dir, minclouds = minclouds, drive_folder = drive_folder)
    
    if(file.exists(paste0(area_landsat_dir, "/lss.csv"))){
      lss <- bind_rows(read_csv(paste0(area_landsat_dir, "/lss.csv")) %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED)),
                       lss) %>% 
        distinct()
    }
    write_csv(lss, paste0(area_landsat_dir, "/lss.csv"))
  }
  
  # List downloaded images
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  
  lss$file <- lapply(lss$LANDSAT_PRODUCT_ID, function(x){
    xx <- tifs[grepl(x, tifs)]
    xx <- ifelse(length(xx) == 0, NA, xx)
    return(xx)
  }) %>% unlist
  
  lss$area <- site_name
  lss$collection <- gsub("0","C",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
  lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
  lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
  lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),1,3))
  lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),4,6))
  lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
  lss$time <- gsub(".tif","",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
  
  lss <- lss %>% arrange(date) %>% filter(!is.na(file))
  
  # Check if all rasters are working. Remove corrupted files.
  img_remove <- unlist(mclapply(lss$file, check_raster, image_dir = area_landsat_dir, mc.cores = workers))
  
  if(length(img_remove) > 0){
    print(paste0(length(img_remove), " raster(s) not functional. REMOVED!!"))
  }
  
  unlink(paste0(area_landsat_dir,"/",img_remove))
  
  lss <- lss %>% 
    filter(!file %in% img_remove)
  
  lss_d <- lss %>% group_by(date, satid, path) %>% count %>% filter(n > 1) %>% ungroup()
  if(nrow(lss_d) > 0){
    for(ii in seq_len(nrow(lss_d))){
      lss_dd <- lss_d %>% slice(ii)
      
      lss_dd <- right_join(lss, lss_dd, by = join_by(satid, path, date))
      
      rs <- lapply(lss_dd$file, function(x){
        r <- rast(paste0(area_landsat_dir,"/",x))
        r[r == 0] <- NA
        return(r)
      })
      rs <- sprc(rs)
      rs <- mosaic(rs) %>% round
      writeRaster(rs, paste0(area_landsat_dir,"/",lss_dd$file[[1]]),
                  overwrite = TRUE, datatype = "INT2U")
      unlink(paste0(area_landsat_dir,"/",lss_dd$file[[2:nrow(lss_dd)]]))
      lss <- lss %>% filter(!file %in% lss_dd$file[[2:nrow(lss_dd)]])
    }
  }
  
  # First round of selection
  # remove images with very limited clear coverage over the AOI
  lccs <- mclapply(lss$file, calc_coverages, image_dir = area_landsat_dir, mc.cores = workers-1) %>% 
    bind_rows
  
  lss <- full_join(lss, lccs, by = "file")
  
  lss %>% arrange(desc(fill_proportion)) %>% 
    filter(clear_proportion < 0.1) %>% pull(file) -> img_remove
  
  unlink(paste0(area_landsat_dir,"/",img_remove))
  
  lss <- lss %>% 
    filter(!file %in% img_remove)
  
  unlink(list.files(tempdir(), full.names = T))
  write_csv(lss, paste0(area_landsat_dir, "/lss.csv"))
  
  return(lss)
  
}

