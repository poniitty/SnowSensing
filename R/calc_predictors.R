calc_predictors <- function(image_df, site_name, base_landsat_dir){
  
  
  predictor_dir <- paste0(base_landsat_dir,"/",site_name,"/predictors")
  if(!dir.exists(predictor_dir)){
    dir.create(predictor_dir)
  }
  
  r <- rast(paste0(base_landsat_dir,"/",site_name,"/",image_df$file[1]))
  
  epsg <- st_crs(r)$epsg
  
  aoi <- st_bbox(r) %>% 
    st_as_sfc() %>% 
    st_buffer(1000) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    st_as_sf()
  
  # AOI polygon to GEE format
  aoi_ee <- aoi %>% st_transform(crs = 4326) %>% 
    st_geometry() %>% 
    sf_as_ee()
  
  # ALOS DEM
  dataset <- ee$ImageCollection('JAXA/ALOS/AW3D30/V3_2')$filterBounds(geometry = aoi_ee)
  dataset <- dataset$select("DSM")
  
  e <- try({ei <- dataset$getInfo()}, silent = T)
  
  if(class(e) != "try-error"){
    dataset <- dataset$max()
    dataset <- dataset$reproject(crs = paste0("EPSG:",epsg), scale = 30)
    
    task_img <- ee_image_to_drive(
      image = dataset,
      fileFormat = "GEO_TIFF",
      folder = "dems",
      region = aoi_ee,
      scale = 30,
      fileNamePrefix = paste0(site_name,"_ALOSDEM")
    )
    
    task_img$start()
    ee_monitoring(task_img, max_attempts = 150)
    
    flds <- drive_find(q = sprintf("name contains '%s'","dems"), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- lapply(gls$id, function(x){
      # x <- gls$id[1]
      fnm <- paste0("ALOSDEM", ".tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    })
    
    drive_rm(drive_find(q = sprintf("name contains '%s'","dems"), type = "folder"))
    
  }
  
  # ESA LAND COVER
  dataset <- ee$ImageCollection('ESA/WorldCover/v100')$filterBounds(geometry = aoi_ee)
  dataset <- dataset$select("Map")
  
  e <- try({ei <- dataset$getInfo()}, silent = T)
  
  if(class(e) != "try-error"){
    dataset <- dataset$max()
    dataset <- dataset$reproject(crs = paste0("EPSG:",epsg), scale = 10)
    
    task_img <- ee_image_to_drive(
      image = dataset,
      fileFormat = "GEO_TIFF",
      folder = "dems",
      region = aoi_ee,
      scale = 10,
      fileNamePrefix = paste0(site_name,"_ESALC")
    )
    
    task_img$start()
    ee_monitoring(task_img, max_attempts = 150)
    
    flds <- drive_find(q = sprintf("name contains '%s'","dems"), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- lapply(gls$id, function(x){
      # x <- gls$id[1]
      fnm <- paste0("ESALC", ".tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    })
    
    drive_rm(drive_find(q = sprintf("name contains '%s'","dems"), type = "folder"))
  }
  
  # GLIMS
  
  dataset <- ee$FeatureCollection('GLIMS/20230607')$filterBounds(geometry = aoi_ee)
  
  e <- try({ei <- dataset$getInfo()}, silent = T)
  
  if(class(e)[[1]] != "try-error"){
    
    dataset <- dataset$
      filter(ee$Filter$notNull(list("anlys_id")))$
      reduceToImage(
        properties = list("anlys_id"),
        reducer = ee$Reducer$first()
      )
    
    dataset <- dataset$gt(0.2)
    
    dataset <- dataset$reproject(crs = paste0("EPSG:",epsg), scale = 30)
    
    task_img <- ee_image_to_drive(
      image = dataset,
      fileFormat = "GEO_TIFF",
      folder = "dems",
      region = aoi_ee,
      scale = 30,
      fileNamePrefix = paste0(site_name,"_GLIMS")
    )
    
    task_img$start()
    ee_monitoring(task_img, max_attempts = 150)
    
    flds <- drive_find(q = sprintf("name contains '%s'","dems"), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- lapply(gls$id, function(x){
      # x <- gls$id[1]
      fnm <- paste0("GLIMS", ".tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(predictor_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    })
    
    drive_rm(drive_find(q = sprintf("name contains '%s'","dems"), type = "folder"))
  }
  
  
  ############################################################
  # DEM VARIABLES
  
  dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
  
  # SWI
  dem_filled <- saga$ta_preprocessor$fill_sinks_xxl_wang_liu(elev = dem, minslope = 0.01)
  
  swi <- saga$ta_hydrology$saga_wetness_index(dem = dem_filled, area_type = 2, slope_type = 0,)
  crs(swi$twi) <- crs(dem)
  names(swi$twi) <- "swi"
  writeRaster(round(swi$twi*100), paste0(predictor_dir, "/swi.tif"), 
              filetype = "GTiff", overwrite = T, datatype = "INT2S")
  
  # SLOPE
  slp <- saga$ta_morphometry$slope_aspect_curvature(elevation = dem)
  slp$slope <- 180/pi*slp$slope
  crs(slp$slope) <- crs(dem)
  names(slp$slope) <- "slope"
  writeRaster(round(slp$slope*100), paste0(predictor_dir, "/slope.tif"), 
              filetype = "GTiff", overwrite = T, datatype = "INT2U")
  
  # SVF
  lat <- st_coordinates(st_transform(st_centroid(st_as_sf(st_as_sfc(st_bbox(dem)))),4326))[1,"Y"]
  svf <- saga$ta_lighting$sky_view_factor(dem = dem,
                                          ndirs = 8, radius = 5000)
  crs(svf$svf) <- crs(dem)
  names(svf$svf) <- "svf"
  writeRaster(round(svf$svf*1000), 
              paste0(predictor_dir, "/svf.tif"), 
              filetype = "GTiff", overwrite = T, datatype = "INT2U")
  
  saga_remove_tmpfiles()
  
  ############################################################
  # MEDIAN LANDSATS
  
  if(nrow(image_df) > 0){
    
    imagedf2 <- image_df %>% 
      filter(cloud_proportion < 0.8,
             tier == "T1") %>% 
      mutate(month = month(date)) %>% 
      relocate(month, .after = date)
    
    if(nrow(imagedf2) < 50){
      imagedf2 <- image_df %>% 
        filter(cloud_proportion < 0.8) %>% 
        mutate(month = month(date)) %>% 
        relocate(month, .after = date)
    }
    
    for(mo in unique(imagedf2$month) %>% sort()){
      # mo <- 7
      print(mo)
      imagedf3 <- imagedf2 %>% 
        filter(month == mo)
      
      if(nrow(imagedf3) > 50){
        imagedf3 <- imagedf3 %>% 
          arrange(desc(clear_proportion)) %>% 
          slice_head(n = 50)
      }
      
      rr <- lapply(imagedf3$file, function(x){
        # x <- imagedf3$file[1]
        sat_id <- imagedf3 %>% filter(file == x) %>% pull(satid)
        r1 <- rast(paste0(base_landsat_dir,"/",site_name,"/",x))
        
        if(sat_id == "LC08"){
          names(r1) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
        } else {
          rr <- r1[[1]]
          rr[] <- NA
          r1 <- c(rr, r1)
          names(r1) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
        }
        
        cmask <- r1[[1]]
        cmask[] <- unlist(lapply(as.numeric(values(r1[["QA"]])), mask_cloud))
        
        r1 <- r1[[1:7]]
        # plotRGB(r1, r=5, g=4, b=3, stretch = "lin")
        
        r1[cmask == 1] <- NA
        r1[r1 == 0] <- NA
        
        return(r1)
      })
      
      rr <- sprc(rr)
      
      rr <- terra::mosaic(rr, fun = "median")
      # plotRGB(rr, r=5, g=4, b=3, stretch = "lin")
      # plot(rr[["B7"]])
      for(ilayer in names(rr)){
        if(mean(is.na(values(rr[[ilayer]], mat = F))) == 1){
          rr[[ilayer]][] <- 0
        }
      }
      
      if(sum(is.na(values(rr, mat = F))) > 0){
        rr <- focal(rr, 5, "mean", na.rm = T, na.policy = "only")
      }
      
      # plotRGB(rr, r=5, g=4, b=3, stretch = "lin")
      # plotRGB(rr, r=7, g=5, b=4, stretch = "lin")
      
      writeRaster(round(rr), paste0(predictor_dir, "/medianlandsat_", mo, ".tif"),
                  datatype = "INT2U", overwrite = T)
    }
    
    
    r2 <- rast(paste0(predictor_dir, "/medianlandsat_",
                      unique(imagedf2$month) %>% sort(), ".tif"))
    
    all <- rast()
    for(ilayer in unique(names(r2))){
      # print(ilayer)
      r3 <- r2[[which(names(r2) == ilayer)]]
      
      r4 <- approximate(r3, rule = 2, NArule = 2)
      all <- c(all, r4)
    }
    
    nmo <- unique(imagedf2$month) %>% sort() %>% length()
    for(mo in unique(imagedf2$month) %>% sort()){
      # mo <- 1
      wmo <- which(unique(imagedf2$month) %>% sort() == mo)
      rr <- all[[seq(wmo, nmo*7, nmo)]]
      
      for(ilayer in names(rr)){
        if(mean(is.na(values(rr[[ilayer]], mat = F))) == 1){
          rr[[ilayer]][] <- 0
        }
      }
      
      if(sum(is.na(values(rr, mat = F))) > 0){
        repeat{
          rr <- focal(rr, 5, "mean", na.policy="only", na.rm=TRUE)
          if(sum(is.na(values(rr, mat = F))) == 0){ break }
        }
      }
      
      writeRaster(round(rr), paste0(predictor_dir, "/medianlandsat_", mo, ".tif"),
                  datatype = "INT2U", overwrite = T)
    }
  }
  
  # MEDIAN INDICES
  
  imagedf2 <- image_df %>% 
    filter(fill_proportion < 0.8, 
           cloud_proportion < 0.8,
           tier == "T1") %>% 
    mutate(month = month(date)) %>% 
    relocate(month, .after = date) %>% 
    group_by(month) %>% 
    arrange(desc(clear_proportion)) %>% 
    slice_head(n = 3) %>% 
    ungroup()
  
  if(nrow(imagedf2) < 15){
    imagedf2 <- image_df %>% 
      filter(fill_proportion < 0.8, 
             cloud_proportion < 0.8) %>% 
      mutate(month = month(date)) %>% 
      relocate(month, .after = date) %>% 
      group_by(month) %>% 
      arrange(desc(clear_proportion)) %>% 
      slice_head(n = 3) %>% 
      ungroup()
  }
  
  rr <- lapply(imagedf2$file, function(x){
    # x <- "dobbiaco_C2_T1_LT05_192027_19900719_091802GMT.tif"
    # x <- "dobbiaco_C2_T1_LE07_192027_20010215_094823GMT.tif"
    sat_id <- imagedf2 %>% filter(file == x) %>% pull(satid)
    r <- rast(paste0(base_landsat_dir,"/", site_name,"/",x))
    
    if(sat_id == "LC08"){
      names(r) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
    } else {
      rr <- r[[1]]
      rr[] <- NA
      r <- c(rr, r)
      names(r) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
    }
    r[r[["QA"]] == 0] <- NA
    r[[1:7]] <- r[[1:7]]*2.75e-05-0.2
    # plotRGB(r1, r=5, g=4, b=3, stretch = "lin")
    
    cmask <- r[[1]]
    cmask[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cloud))
    smask <- r[[1]]
    smask[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_snow))
    
    ndvi <- (r[["B5"]]-r[["B4"]])/(r[["B5"]]+r[["B4"]])
    ndvi[smask == 1] <- NA
    names(ndvi) <- "ndvi"
    
    ndsi <- (r[["B3"]]-r[["B6"]])/(r[["B3"]]+r[["B6"]])
    names(ndsi) <- "ndsi"
    
    kbri <- (r[["B6"]] - r[["B5"]])/(20*sqrt((r[["B6"]] + r[["B5"]])))
    kbri[smask == 1] <- NA
    names(kbri) <- "kbri"
    
    swm <- (r[["B2"]] + r[["B3"]])/(r[["B5"]] + r[["B6"]])
    swm[smask == 1] <- NA
    names(swm) <- "swm"
    
    bitm <- (((r[["B2"]]**2.0)+(r[["B3"]]**2.0)+(r[["B4"]]**2.0))/3.0)**0.5
    bitm[smask == 1] <- NA
    names(bitm) <- "bitm"
    
    r1 <- c(ndvi, ndsi, kbri, swm, bitm)
    
    return(r1)
  })
  
  ndvi <- lapply(rr, function(x){
    x <- x[["ndvi"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  ndvi <- sprc(ndvi)
  ndvi <- terra::mosaic(ndvi, fun = "mean")
  
  if(sum(is.na(values(ndvi, mat = F))) > 0){
    repeat{
      ndvi <- focal(ndvi, 5, "mean", na.policy="only", na.rm=TRUE)
      if(sum(is.na(values(ndvi, mat = F))) == 0){ break }
    }
  }
  
  kbri <- lapply(rr, function(x){
    x <- x[["kbri"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  kbri <- sprc(kbri)
  kbri <- terra::mosaic(kbri, fun = "mean")
  if(sum(is.na(values(kbri, mat = F))) > 0){
    repeat{
      kbri <- focal(kbri, 5, "mean", na.policy="only", na.rm=TRUE)
      if(sum(is.na(values(kbri, mat = F))) == 0){ break }
    }
  }
  
  ndsi <- lapply(rr, function(x){
    x <- x[["ndsi"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  ndsi <- sprc(ndsi)
  ndsi <- terra::mosaic(ndsi, fun = "mean")
  if(sum(is.na(values(ndsi, mat = F))) > 0){
    repeat{
      ndsi <- focal(ndsi, 5, "mean", na.policy="only", na.rm=TRUE)
      if(sum(is.na(values(ndsi, mat = F))) == 0){ break }
    }
  }
  
  swm <- lapply(rr, function(x){
    x <- x[["swm"]]
    x[x > 10] <- NA
    x[x < -10] <- -NA
    return(x)
  })
  swm <- sprc(swm)
  swm <- terra::mosaic(swm, fun = "mean")
  if(sum(is.na(values(swm, mat = F))) > 0){
    repeat{
      swm <- focal(swm, 5, "mean", na.policy="only", na.rm=TRUE)
      if(sum(is.na(values(swm, mat = F))) == 0){ break }
    }
  }
  
  bitm <- lapply(rr, function(x){
    x <- x[["bitm"]]
    # x[x > 10] <- NA
    # x[x < -10] <- -NA
    return(x)
  })
  bitm <- sprc(bitm)
  bitm <- terra::mosaic(bitm, fun = "median")
  if(sum(is.na(values(bitm, mat = F))) > 0){
    repeat{
      bitm <- focal(bitm, 5, "mean", na.policy="only", na.rm=TRUE)
      if(sum(is.na(values(bitm, mat = F))) == 0){ break }
    }
  }
  
  r1 <- c(ndvi, ndsi, kbri, swm, bitm)
  names(r1) <- c("ndvi", "ndsi", "kbri", "swm", "bitm")
  
  writeRaster(round(r1, 4), paste0(predictor_dir, "/medianindices.tif"),
              overwrite = T)
  
  unlink(list.files(tempdir(), full.names = T, recursive = T))
  
}
