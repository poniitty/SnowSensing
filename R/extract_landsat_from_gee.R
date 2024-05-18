# i <- "2023-05-23"
#
extract_landsat_gee <- function(aoi_ee,
                                epsg,
                                excl_dates,
                                site_name,
                                sats,
                                start_date, 
                                end_date, 
                                months = 1:12,
                                minclouds,
                                my_timeout, 
                                info_timeout, 
                                area_landsat_dir,
                                drive_folder,
                                workers){
  
  idsALL <- tibble()
  
  # LANDSAT TM4
  if("LT04" %in% sats){
    drive_mkdir(drive_folder)
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT04/C02/T1_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids) > 0){
      ids <- ids %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LT04") %>% pull("date")))
    }
    
    ress <- mclapply(ids$id, function(x, drive_folder, my_timeout, epsg, ids){
      # x <- 1
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids$id[(lapply(ids$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    # T2, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT04/C02/T2_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids2 <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids2) > 0){
      ids2 <- ids2 %>% 
        filter(!DATE_ACQUIRED %in% ids$DATE_ACQUIRED) %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LT04") %>% pull("date")))
    }
    
    ress <- mclapply(ids2$id, function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids2$id[(lapply(ids2$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    idsALL <- bind_rows(idsALL, bind_rows(ids, ids2))
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- mclapply(gls$id, function(x){
      # x <- gls$id[1]
      nm <- paste(str_split(gls[gls$id == x,"name"], "_")[[1]][1:7], collapse = "_")
      fnm <- paste0(nm, "_", gsub(":", "", substr(idsALL[grepl(nm, idsALL$LANDSAT_PRODUCT_ID),"SCENE_CENTER_TIME"], 1, 8)), "GMT.tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    }, mc.cores = workers)
    
    drive_rm(drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder"))
    
  }
  
  # LANDSAT TM5
  if("LT05" %in% sats){
    drive_mkdir(drive_folder)
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT05/C02/T1_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids) > 0){
      ids <- ids %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LT05") %>% pull("date")))
    }
    
    ress <- mclapply(ids$id, function(x, drive_folder, my_timeout, epsg, ids){
      # x <- 1
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids$id[(lapply(ids$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    # T2, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT05/C02/T2_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids2 <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids2) > 0){
      ids2 <- ids2 %>% 
        filter(!DATE_ACQUIRED %in% ids$DATE_ACQUIRED) %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LT05") %>% pull("date")))
    }
    
    ress <- mclapply(ids2$id, function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids2$id[(lapply(ids2$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    idsALL <- bind_rows(idsALL, bind_rows(ids, ids2))
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- mclapply(gls$id, function(x){
      # x <- gls$id[1]
      nm <- paste(str_split(gls[gls$id == x,"name"], "_")[[1]][1:7], collapse = "_")
      fnm <- paste0(nm, "_", gsub(":", "", substr(idsALL[grepl(nm, idsALL$LANDSAT_PRODUCT_ID),"SCENE_CENTER_TIME"], 1, 8)), "GMT.tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    }, mc.cores = workers)
    
    drive_rm(drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder"))
    
  }
  
  # LANDSAT ETM7
  if("LE07" %in% sats){
    drive_mkdir(drive_folder)
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LE07/C02/T1_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids) > 0){
      ids <- ids %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LE07") %>% pull("date")))
    }
    
    ress <- mclapply(ids$id, function(x, drive_folder, my_timeout, epsg, ids){
      # x <- 1
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids$id[(lapply(ids$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    # T2, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LE07/C02/T2_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids2 <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY = x$properties$IMAGE_QUALITY,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids2) > 0){
      ids2 <- ids2 %>% 
        filter(!DATE_ACQUIRED %in% ids$DATE_ACQUIRED) %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LE07") %>% pull("date")))
    }
    
    ress <- mclapply(ids2$id, function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids2$id[(lapply(ids2$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    idsALL <- bind_rows(idsALL, bind_rows(ids, ids2))
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- mclapply(gls$id, function(x){
      # x <- gls$id[1]
      nm <- paste(str_split(gls[gls$id == x,"name"], "_")[[1]][1:7], collapse = "_")
      fnm <- paste0(nm, "_", gsub(":", "", substr(idsALL[grepl(nm, idsALL$LANDSAT_PRODUCT_ID),"SCENE_CENTER_TIME"], 1, 8)), "GMT.tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    }, mc.cores = workers)
    
    drive_rm(drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder"))
    
  }
  
  # LANDSAT OLI8
  if("LC08" %in% sats){
    drive_mkdir(drive_folder)
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY_OLI = x$properties$IMAGE_QUALITY_OLI,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids) > 0){
      ids <- ids %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LC08") %>% pull("date")))
    }
    
    ress <- mclapply(ids$id, function(x, drive_folder, my_timeout, epsg, ids){
      # x <- 1
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids$id[(lapply(ids$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    # T2, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC08/C02/T2_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids2 <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY_OLI = x$properties$IMAGE_QUALITY_OLI,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids2) > 0){
      ids2 <- ids2 %>% 
        filter(!DATE_ACQUIRED %in% ids$DATE_ACQUIRED) %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LC08") %>% pull("date")))
    }
    
    ress <- mclapply(ids2$id, function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids2$id[(lapply(ids2$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    idsALL <- bind_rows(idsALL, bind_rows(ids, ids2) %>% rename(IMAGE_QUALITY = IMAGE_QUALITY_OLI))
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- mclapply(gls$id, function(x){
      # x <- gls$id[1]
      nm <- paste(str_split(gls[gls$id == x,"name"], "_")[[1]][1:7], collapse = "_")
      fnm <- paste0(nm, "_", gsub(":", "", substr(idsALL[grepl(nm, idsALL$LANDSAT_PRODUCT_ID),"SCENE_CENTER_TIME"], 1, 8)), "GMT.tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    }, mc.cores = workers)
    
    drive_rm(drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder"))
    
  }
  
  # LANDSAT OLI9
  if("LC09" %in% sats){
    drive_mkdir(drive_folder)
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC09/C02/T1_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY_OLI = x$properties$IMAGE_QUALITY_OLI,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids) > 0){
      ids <- ids %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LC09") %>% pull("date")))
    }
    
    ress <- mclapply(ids$id, function(x, drive_folder, my_timeout, epsg, ids){
      # x <- 1
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids$id[(lapply(ids$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids[ids$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids[ids$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids = ids)
    
    # T2, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC09/C02/T2_L2')$filterDate(as.character(start_date), as.character(end_date))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", minclouds)
    
    e <- try(withTimeout({ei <- dataset$getInfo()}, timeout = info_timeout), silent = T)
    
    ids2 <- lapply(ei$features, function(x){
      xd <- tibble(LANDSAT_PRODUCT_ID = x$properties$LANDSAT_PRODUCT_ID,
                   DATE_ACQUIRED = x$properties$DATE_ACQUIRED,
                   SCENE_CENTER_TIME = x$properties$SCENE_CENTER_TIME,
                   IMAGE_QUALITY_OLI = x$properties$IMAGE_QUALITY_OLI,
                   CLOUD_COVER = x$properties$CLOUD_COVER,
                   SUN_ELEVATION = x$properties$SUN_ELEVATION,
                   GEOMETRIC_RMSE_MODEL = x$properties$GEOMETRIC_RMSE_MODEL)
      return(xd)
    }) %>% bind_rows() %>% 
      rownames_to_column("id")
    if(nrow(ids2) > 0){
      ids2 <- ids2 %>% 
        filter(!DATE_ACQUIRED %in% ids$DATE_ACQUIRED) %>% 
        filter(month(DATE_ACQUIRED) %in% months) %>% 
        filter(!DATE_ACQUIRED %in% (excl_dates %>% filter(satid == "LC09") %>% pull("date")))
    }
    
    ress <- mclapply(ids2$id, function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    ress <- mclapply(ids2$id[(lapply(ids2$LANDSAT_PRODUCT_ID, function(x){ return(any(grepl(x, gls$name)))}) %>% lapply(., isFALSE) %>% unlist)], function(x, drive_folder, my_timeout, epsg, ids2){
      edl <- try({
        dataset2 <- dataset$filterMetadata("LANDSAT_PRODUCT_ID", "equals", as.character(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"]))
        dataset2 <- dataset2$toBands()
        dataset2 <- dataset2$reproject(crs = paste0("EPSG:",epsg), scale = 30)
        dataset2 <- dataset2$toUint16()
        # dataset2$getInfo()
        task_img <- ee_image_to_drive(
          image = dataset2,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(ids2[ids2$id == x,"LANDSAT_PRODUCT_ID"])
        )
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          return(TRUE)
        }
      }, silent = T)
      if(class(edl)[1] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    }, mc.cores = workers-1, drive_folder = drive_folder, my_timeout = my_timeout, epsg = epsg, ids2 = ids2)
    
    idsALL <- bind_rows(idsALL, bind_rows(ids, ids2) %>% rename(IMAGE_QUALITY = IMAGE_QUALITY_OLI))
    
    flds <- drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder")
    gls <- lapply(flds$id, drive_ls) %>% bind_rows()
    
    dl <- mclapply(gls$id, function(x){
      # x <- gls$id[1]
      nm <- paste(str_split(gls[gls$id == x,"name"], "_")[[1]][1:7], collapse = "_")
      fnm <- paste0(nm, "_", gsub(":", "", substr(idsALL[grepl(nm, idsALL$LANDSAT_PRODUCT_ID),"SCENE_CENTER_TIME"], 1, 8)), "GMT.tif")
      
      tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
      if(class(tdl)[1] == "try-error"){
        tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        if(class(tdl)[1] == "try-error"){
          tdl <- try(dl <- drive_download(x, path = paste0(area_landsat_dir,"/",fnm), overwrite = TRUE), silent = T)
        } else {
          return(TRUE)
        }
      } else {
        return(TRUE)
      }
    }, mc.cores = workers)
    
    drive_rm(drive_find(q = sprintf("name contains '%s'",drive_folder), type = "folder"))
    
  }
  
  return(idsALL)
}
