mask_cirrus <- function(x) {as.numeric(intToBits(x)[3])}
mask_cloud <- function(x) {as.numeric(intToBits(x)[4])}
mask_cshadow <- function(x) {as.numeric(intToBits(x)[5])}
mask_snow <- function(x) {as.numeric(intToBits(x)[6])}
mask_clear <- function(x) {as.numeric(intToBits(x)[7])}
mask_water <- function(x) {as.numeric(intToBits(x)[8])}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

classify_landsat <- function(image_df, site_name, base_landsat_dir, model_path, workers, saga){
  
  predictor_dir <- paste0(base_landsat_dir,"/",site_name,"/predictors")
  class_landsat_dir <- paste0(base_landsat_dir,"/",site_name,"/classifications")
  if(!dir.exists(class_landsat_dir)){
    dir.create(class_landsat_dir)
  }
  
  mod5 <- readRDS(paste0(model_path, "RF_model_TM05_3.rds"))
  mod7 <- readRDS(paste0(model_path, "RF_model_ETM07_3.rds"))
  mod8 <- readRDS(paste0(model_path, "RF_model_OLI08_3.rds"))
  
  r <- rast(paste0(base_landsat_dir, "/", site_name, "/", image_df$file[1]), 1)
  
  slope <- rast(paste0(predictor_dir, "/slope.tif"))/100
  slope <- project(slope, r)
  
  swi <- rast(paste0(predictor_dir, "/swi.tif"))/100
  swi <- project(swi, r)
  
  esalc <- rast(paste0(predictor_dir, "/ESALC.tif"))/10
  esalc <- aggregate(esalc, 3, getmode)
  esalc <- project(esalc, r, method="near")
  
  esaw <- rast(paste0(predictor_dir, "/ESALC.tif"))
  esaw[esaw == 80] <- 1
  esaw[esaw != 1] <- 0
  esaw <- aggregate(esaw, 3, mean)
  esaw <- project(esaw, r)
  if(file.exists(paste0(predictor_dir,"/OSM_WATERS.tif"))){
    wt <- rast(paste0(predictor_dir,"/OSM_WATERS.tif"))
    wt <- project(wt, esaw)
    esaw <- sum(wt, esaw)
    esaw[esaw > 1] <- 1
  }
  slope[esaw > 0.2] <- 0
  
  if(file.exists(paste0(predictor_dir, "/GLIMS.tif"))){
    glaciers <- rast(paste0(predictor_dir, "/GLIMS.tif"))
    glaciers <- project(glaciers, r, method="near")
  } else {
    glaciers <- esaw
    glaciers[] <- 0
  }
  
  dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
  dem <- project(dem, r)
  
  svf <- rast(paste0(predictor_dir, "/svf.tif"))/1000
  svf <- project(svf, r)
  
  suppressWarnings(suppressMessages(lss <- lapply(image_df$file, classifying_function,
                                                  image_df = image_df, site_name = site_name, saga = saga,
                                                  mod7 = mod7, mod5 = mod5, mod8 = mod8,
                                                  predictor_dir = predictor_dir, class_landsat_dir = class_landsat_dir,
                                                  swi = swi, slope = slope, esalc = esalc, esaw = esaw, 
                                                  glaciers = glaciers, dem = dem, svf = svf)))
  
  return(unlist(lss))
  
}

classifying_function <- function(imageid, image_df, predictor_dir, class_landsat_dir, site_name, saga, 
                                 mod8, mod7, mod5, swi, slope, esalc, esaw, glaciers, dem, svf){
  # imageid <- "LC08_L2SP_196013_20220504_20220518_02_T1_101656GMT.tif"
  if(!file.exists(paste0(class_landsat_dir, "/", imageid))){
    
    # print(imageid)
    satid <- image_df %>% filter(file == imageid) %>% pull(satid)
    
    if(satid %in% c("LT05","LT04")){
      mod_vars <- mod5$xvar.names
    }
    if(satid == "LE07"){
      mod_vars <- mod7$xvar.names
    }
    if(satid == "LC08"){
      mod_vars <- mod8$xvar.names
    }
    
    r <- rast(paste0(base_landsat_dir, "/", site_name, "/", imageid))
    
    if(satid == "LC08"){
      names(r) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
    } else {
      rr <- r[[1]]
      rr[] <- 0
      r <- c(rr, r)
      names(r) <- c("B1","B2","B3","B4","B5","B6","B7","B10","QA","RADSAT")
    }
    r[r[["QA"]] == 0] <- NA
    
    if("swi" %in% mod_vars){
      r[["swi"]] <- swi
    }
    if("esalc" %in% mod_vars){
      r[["esalc"]] <- esalc
    }
    if("esaw" %in% mod_vars){
      r[["esaw"]] <- esaw
    }
    if("slope" %in% mod_vars){
      r[["slope"]] <- slope
    }
    if("GLIMS" %in% mod_vars){
      r[["GLIMS"]] <- glaciers
    }
    
    if(sum(!is.na(values(r[[1]]))) > 100){
      
      # Pisr
      
      idate <- as.character(ymd(str_split(imageid, "_")[[1]][4]))
      if("pisr" %in% mod_vars){
        
        lat <- st_coordinates(st_transform(st_centroid(st_as_sf(st_as_sfc(st_bbox(raster::raster(dem))))),4326))[1,"Y"]
        lon <- st_coordinates(st_transform(st_centroid(st_as_sf(st_as_sfc(st_bbox(raster::raster(dem))))),4326))[1,"X"]
        
        time_offset <- lon * 24 / 360
        gmt_time_dec <- as.numeric(paste0(substr(tail(str_split(imageid, "_")[[1]], 1), 1, 2),
                                          ".",
                                          round(as.numeric(substr(tail(str_split(imageid, "_")[[1]], 1), 3, 4))/60*100)))
        
        moment <- gmt_time_dec + time_offset
        
        if(moment > 24){
          moment <- gmt_time_dec + time_offset - 24
        }
        if(moment < 0){
          moment <- gmt_time_dec + time_offset + 24
        }
        pisr <- saga$ta_lighting$potential_incoming_solar_radiation(grd_dem = dem,
                                                                    grd_svf = svf,
                                                                    latitude = lat,
                                                                    period = 0,
                                                                    day = idate,
                                                                    moment = moment, .verbose = F)$grd_total
        
        if(is.null(pisr)){
          pisr <- saga$ta_lighting$potential_incoming_solar_radiation(grd_dem = dem,
                                                                      grd_svf = svf,
                                                                      latitude = lat,
                                                                      period = 0,
                                                                      day = idate,
                                                                      moment = moment+0.10, .verbose = F)$grd_total
        }
        if(is.null(pisr)){
          pisr <- saga$ta_lighting$potential_incoming_solar_radiation(grd_dem = dem,
                                                                      grd_svf = svf,
                                                                      latitude = lat,
                                                                      period = 0,
                                                                      day = idate,
                                                                      moment = moment-0.10, .verbose = F)$grd_total
        }
        crs(pisr) <- crs(r)
        pisr[esaw > 0.2] <- median(values(pisr, mat = F), na.rm = T)
        r[["pisr"]] <- pisr
      }
      
      # Median landsat
      
      mo <- month(idate)
      e <- try({
        mr <- rast(paste0(predictor_dir, "/medianlandsat_",mo,".tif"))
      })
      if(class(e) == "try-error"){
        e <- try({
          mr <- rast(paste0(predictor_dir, "/medianlandsat_",mo+1,".tif"))
        })
      }
      if(class(e) == "try-error"){
        e <- try({
          mr <- rast(paste0(predictor_dir, "/medianlandsat_",mo-1,".tif"))
        })
      }
      
      mr_diff <- mr - r[[1:7]]
      names(mr_diff) <- paste0("mr_diff_",names(mr_diff))
      if(any(grepl("mr_diff_B5_std5", mod_vars))){mr_diff[["mr_diff_B5_std5"]] <- focal(r[["B5"]], w = 5, fun = "sd")-focal(mr[["B5"]], w = 5, fun = "sd")}
      if(any(grepl("mr_diff_B7_std5", mod_vars))){mr_diff[["mr_diff_B7_std5"]] <- focal(r[["B7"]], w = 5, fun = "sd")-focal(mr[["B7"]], w = 5, fun = "sd")}
      
      # median indices
      mi <- rast(paste0(predictor_dir, "/medianindices.tif"))
      
      # Focal stuff with the focal imagery
      
      filter <- matrix(1, nrow=5, ncol=5) 
      filter[ceiling(length(filter)/2)] <- 0
      
      if(any(grepl("B2_tpi", mod_vars))){r[["B2_tpi"]] <- r[["B2"]] - focal(r[["B2"]], w = filter, fun = mean,
                                                                            na.policy = "omit", na.rm = T)}
      if(any(grepl("B5_tpi", mod_vars))){r[["B5_tpi"]] <- r[["B5"]] - focal(r[["B5"]], w = filter, fun = mean,
                                                                            na.policy = "omit", na.rm = T)}
      if(any(grepl("B7_tpi", mod_vars))){r[["B7_tpi"]] <- r[["B7"]] - focal(r[["B7"]], w = filter, fun = mean,
                                                                            na.policy = "omit", na.rm = T)}
      
      if(any(grepl("B2_min5", mod_vars))){r[["B2_min5"]] <- focal(r[["B2"]], w = 5, fun = "min", na.policy = "omit", na.rm = T)}
      if(any(grepl("B5_min5", mod_vars))){r[["B5_min5"]] <- focal(r[["B5"]], w = 5, fun = "min", na.policy = "omit", na.rm = T)}
      if(any(grepl("B7_min5", mod_vars))){r[["B7_min5"]] <- focal(r[["B7"]], w = 5, fun = "min", na.policy = "omit", na.rm = T)}
      
      if(any(grepl("B2_max5", mod_vars))){r[["B2_max5"]] <- focal(r[["B2"]], w = 5, fun = "max", na.policy = "omit", na.rm = T)}
      if(any(grepl("B5_max5", mod_vars))){r[["B5_max5"]] <- focal(r[["B5"]], w = 5, fun = "max", na.policy = "omit", na.rm = T)}
      if(any(grepl("B7_max5", mod_vars))){r[["B7_max5"]] <- focal(r[["B7"]], w = 5, fun = "max", na.policy = "omit", na.rm = T)}
      
      if(any(grepl("B2_std5", mod_vars))){r[["B2_std5"]] <- focal(r[["B2"]], w = 5, fun = "sd", na.policy = "omit", na.rm = T)}
      if(any(grepl("B5_std5", mod_vars))){r[["B5_std5"]] <- focal(r[["B5"]], w = 5, fun = "sd", na.policy = "omit", na.rm = T)}
      if(any(grepl("B7_std5", mod_vars))){r[["B7_std5"]] <- focal(r[["B7"]], w = 5, fun = "sd", na.policy = "omit", na.rm = T)}
      
      # Previous classification
      rrrr <- rast(ncols=180, nrows=180, xmin=0)
      fm <- focalMat(rrrr, 5, "circle")
      fm[fm > 0] <- 1
      
      if("cirrus_probs5" %in% mod_vars){
        cirrus_probs5 <- r[[1]]
        cirrus_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cirrus))
        r[["cirrus_probs5"]] <- focal(cirrus_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("cloud_probs5" %in% mod_vars){
        cloud_probs5 <- r[[1]]
        cloud_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cloud))
        r[["cloud_probs5"]] <- focal(cloud_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("cshadow_probs5" %in% mod_vars){
        cshadow_probs5 <- r[[1]]
        cshadow_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cshadow))
        r[["cshadow_probs5"]] <- focal(cshadow_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("snow_probs5" %in% mod_vars){
        snow_probs5 <- r[[1]]
        snow_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_snow))
        r[["snow_probs5"]] <- focal(snow_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("water_probs5" %in% mod_vars){
        water_probs5 <- r[[1]]
        water_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_water))
        r[["water_probs5"]] <- focal(water_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("clear_probs5" %in% mod_vars){
        clear_probs5 <- r[[1]]
        clear_probs5[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_clear))
        r[["clear_probs5"]] <- focal(clear_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      if("cloud_probs25" %in% mod_vars){
        fm <- focalMat(rrrr, 25, "circle")
        fm[fm > 0] <- 1
        
        cloud_probs25 <- r[[1]]
        cloud_probs25[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cloud))
        r[["cloud_probs25"]] <- focal(cloud_probs25, w = fm, fun = "mean", na.policy = "omit", na.rm = T)
      }
      
      r <- c(r, mr_diff, mi)
      r[["QA"]][is.na(r[["B1"]])] <- NA
      r[["QA"]][is.na(r[["B2"]])] <- NA
      r[["QA"]][is.na(r[["B3"]])] <- NA
      r[["QA"]][is.na(r[["B4"]])] <- NA
      r[["QA"]][is.na(r[["B5"]])] <- NA
      r[["QA"]][is.na(r[["B6"]])] <- NA
      r[["QA"]][is.na(r[["B7"]])] <- NA
      
      d <- as.data.frame(r, na.rm=FALSE) %>% 
        filter(!is.na(QA)) %>% 
        mutate(across(ends_with("tpi"), ~ifelse(is.na(.x), 0, .x)))
      
      if(nrow(d) == 0){
        
        rr <- r[["QA"]]
        rr[["class"]] <- NA
        rr[["land"]] <- NA
        rr[["water"]] <- NA
        rr[["snow"]] <- NA
        rr[["cloud"]] <- NA
        rr[["artif"]] <- NA
        
        rr <- rr[[-1]]
        # plot(rr)
        writeRaster(rr, paste0(class_landsat_dir, "/", imageid),
                    datatype = "INT1U", overwrite = T)
        
        saga_remove_tmpfiles()
        
        return(imageid)
        
      } else {
        layers <- c("B1","B2","B3","B4","B5","B6","B7",
                    names(d %>% select(ends_with("_min5"), ends_with("_max5"), ends_with("_std5")) %>% 
                            select(-starts_with("B2"))))
        for(i in layers){
          # print(i)
          layers <- layers[-1]
          for(ii in layers){
            d[,paste(i, ii, sep = "_")] <- (d[,i]-d[,ii])/(d[,i]+d[,ii])
          }
        }
        
        d <- d %>% 
          mutate(satid = ifelse(satid == "LT04", "LT05", satid))
        
        d <- d %>% 
          mutate(across(everything(), ~ifelse(is.nan(.x), 0, .x))) #%>% 
        # mutate(across(all_of(c("esalc","glcfcs")), factor))
        
        std_names <- names(d)[grepl("_std", names(d))]
        d <- d %>% 
          mutate(across(all_of(std_names), ~ifelse(is.na(.x), 0, .x)))
        
        # PREDICTIONS
        if(satid %in% c("LT05","LT04")){
          preds <- predict(mod5, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`5`, ~round(.x*100)))
        }
        if(satid == "LE07"){
          preds <- predict(mod7, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`4`, ~round(.x*100)))
        }
        if(satid == "LC08"){
          preds <- predict(mod8, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`4`, ~round(.x*100)))
        }
        # d %>% select(-RADSAT) %>% filter(!complete.cases(.))
        rm(d)
        
        if(nrow(preds) == length(na.omit(values(r[["QA"]], mat = F)))){
          rr <- r[["QA"]]
          rr[["class"]] <- 0
          rr[["class"]][!is.na(r[["QA"]])] <- as.numeric(preds$pred)
          rr[["land"]] <- 0
          rr[["land"]][!is.na(r[["QA"]])] <- as.numeric(preds$`1`)
          rr[["water"]] <- 0
          rr[["water"]][!is.na(r[["QA"]])] <- as.numeric(preds$`2`)
          rr[["snow"]] <- 0
          rr[["snow"]][!is.na(r[["QA"]])] <- as.numeric(preds$`3`)
          rr[["cloud"]] <- 0
          rr[["cloud"]][!is.na(r[["QA"]])] <- as.numeric(preds$`4`)
          rr[["artif"]] <- 0
          if(satid %in% c("LT05","LT04")){
            rr[["artif"]][!is.na(r[["QA"]])] <- as.numeric(preds$`5`)
          }
          
          rr[is.na(r[["QA"]])] <- NA
          rr <- rr[[-1]]
          # plot(rr)
          # plotRGB(r, r=7, g=5, b=2, stretch = "lin")
          writeRaster(rr, paste0(class_landsat_dir, "/", imageid),
                      datatype = "INT1U", overwrite = T)
        } else {
          stop("FUUUUUUUUUUUUUUUUU")
        }
        
        saga_remove_tmpfiles()
        
        return(imageid)
      }
    }
  }
}
