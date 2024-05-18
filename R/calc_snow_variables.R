calc_snow_variables <- function(image_df, site_name, base_landsat_dir, workers){
  
  class_landsat_dir <- paste0(base_landsat_dir,"/",site_name,"/classifications")
  predictor_dir <- paste0(base_landsat_dir,"/",site_name,"/predictors")
  
  class_tifs <- list.files(class_landsat_dir, pattern = "GMT.tif$")
  
  image_df <- image_df %>% 
    filter(file %in% class_tifs)
  
  #########################################################################
  # Calc class stats
  
  if(nrow(image_df) > 0){
    
    image_df <- image_df %>%
      mutate(cloud_proportion_own = as.numeric(NA),
             clear_proportion_own = as.numeric(NA),
             land_proportion_own = as.numeric(NA),
             snow_proportion_own = as.numeric(NA))
    
    for(imageid in (image_df %>% pull(file))){
      # imageid <- "LC08_L2SP_217013_20230611_20230614_02_T1_122620GMT.tif"
      
      r <- rast(paste0(class_landsat_dir, "/", imageid))
      vr2 <- values(r[[1]], mat = F)
      
      image_df[image_df$file == imageid, "cloud_proportion_own"] <- round(mean(vr2 == 4, na.rm =  T), 3)
      image_df[image_df$file == imageid, "clear_proportion_own"] <- round(mean(vr2 < 4, na.rm = T), 3)
      image_df[image_df$file == imageid, "land_proportion_own"] <- round(mean(vr2 == 1, na.rm = T), 3)
      image_df[image_df$file == imageid, "snow_proportion_own"] <- round(mean(vr2 == 3, na.rm = T), 3)
    }
  }
  
  # Filter fully cloudy
  image_df <- image_df %>% 
    filter(cloud_proportion_own < 0.80)
  
  
  if(nrow(image_df) >= 50){
    
    image_df %>% 
      filter(clear_proportion_own > 0.8 & fill_proportion < 0.2) %>% 
      mutate(doy = yday(date)) %>% 
      gam(snow_proportion_own ~ s(doy, bs = "cc"), data = .,
          knots=list(doy=c(1,365)), family = "binomial") -> gmod
    # summary(gmod)
    preds <- tibble(doy = 1:365)
    preds$preds <- predict(gmod, preds, type = "response")
    
    maxsnowdoy <- preds$doy[which.max(preds$preds)]
    maxsnowdoy <- ifelse(maxsnowdoy > 250, 30, maxsnowdoy)
    maxsnowdoy <- ifelse(maxsnowdoy < (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min)+14,
                         (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min)+14, maxsnowdoy)
    minsnowdoy <- preds$doy[which.min(preds$preds)]
    gc()
    
    # image_df %>%
    #   filter(clear_proportion_own > 0.8 & fill_proportion < 0.2) %>%
    #   mutate(doy = yday(date)) %>%
    #   ggplot(aes(y = snow_proportion_own, x = doy)) +
    #   geom_point() +
    #   geom_line(mapping = aes(y = preds, x = doy), data = preds, color = "gray", size = 1) +
    #   geom_point(mapping = aes(y = preds, x = doy), data = preds %>% filter(doy == maxsnowdoy), color = "red", size = 3) +
    #   geom_point(mapping = aes(y = preds, x = doy), data = preds %>% filter(doy == minsnowdoy), color = "blue", size = 3) +
    #   ylab("Image snow proportion")
    
    rs <- read_classifications(image_df, class_landsat_dir)
    
    # Raster list to data frame
    mm <- rs %>%
      map(nosnow_values) %>% 
      reduce(bind_rows)
    
    rtemp <- rs[[1]][[1]]
    
    # Doys and years of the scenes to the data frame
    mm$doy <- rep(yday(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                  each = ncell(rs[[1]]))
    mm$year <- rep(year(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                   each = ncell(rs[[1]]))
    mm$week <- rep(week(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                   each = ncell(rs[[1]]))
    mm$month <- rep(month(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                    each = ncell(rs[[1]]))
    
    mm <- mm %>% 
      mutate(week = ifelse(week > 52, 52, week))
   
    middoys <- tibble(date = seq.Date(as_date("2022-01-01"),as_date("2022-12-31"),1)) %>% 
      mutate(month = month(date),
             doy = yday(date)) %>% 
      group_by(month) %>% 
      summarise(ndays = n(),
                middoy = mean(doy))
    
    mcs <- image_df %>% 
      mutate(month = month(date)) %>% 
      group_by(month) %>% 
      count
    
    if(nrow(mcs) != 12){
      extradf <- middoys %>% 
        filter(!month %in% mcs$month) %>% 
        mutate(doy = round(middoy),
               snow = 1,
               prop = 80) %>% 
        select(snow, prop, doy, month) %>% 
        sample_n(size = mean(mcs$n)*nrow(.), replace = T)
    } else {
      extradf <- NULL
    }
    
    rm(rs)
    gc()
    
    # Function to calculate snow melting doy and trend
    cal_scd <- function(d_all, extradf){
      # cellid <- 7802
      # cellid <- 3961
      # d_all <- mm %>% mutate(cell2 = cell) %>% nest(data = -cell2) %>% slice(66) %>% pull(data)
      # d_all <- d_all[[1]]
      # d_all <- mm %>% filter(cell == 72)
      cellid <- d_all$cell[1]
      
      d_all <- bind_rows(d_all %>% select(nosnow, nosnow_prop, doy, year, month, week) %>% drop_na() %>% 
                           rename(snow = nosnow, prop = nosnow_prop),
                         d_all %>% select(snow, snow_prop, doy, year, month, week) %>% drop_na() %>% 
                           rename(prop = snow_prop))
      
      
      results <- tibble(cell = cellid,
                        nobs = d_all %>% select(doy, year) %>% distinct() %>% nrow,
                        nyears = length(unique(d_all$year)),
                        gamr2 = as.numeric(NA), 
                        max_snow_doy = as.numeric(NA), min_snow_doy = as.numeric(NA), 
                        max_snow_prop = as.numeric(NA), min_snow_prop = as.numeric(NA),
                        scd_raw = as.numeric(NA), snow_days = as.numeric(NA),
                        scd = as.numeric(NA), melt = as.numeric(NA), ns = as.numeric(NA))
      
      if(!is.null(extradf)){
        d_all <- bind_rows(extradf, d_all)
      }
      
      if(nrow(d_all) > 20){
        
        results$scd_raw <- round(weighted.mean(d_all$snow, d_all$prop)*365,1)
        
        if((results$scd_raw/365) > 0.01 | (results$scd_raw/365) < 0.99 | !is.nan(results$scd_raw)){
          
          # pred_grid <- expand.grid(doy = 1:365,
          #                          year = unique(d_all$year))
          pred_grid <- expand.grid(doy = 1:365)
          
          # d_all %>%
          #   ggplot(aes(y = snow, x = doy)) +
          #   geom_point()
          
          e <- try({
            # d_all %>% 
            #   mutate(year = as_factor(year)) %>% 
            #   gam(snow ~  s(doy, bs = "cc") + s(doy, by = year), 
            #       data = ., family = "binomial", weights = prop, method = "REML") -> gammod
            d_all %>% 
              mutate(year = as_factor(year)) %>% 
              gam(snow ~  s(doy, bs = "cc", k = 5), 
                  knots=list(doy=c(1,365)),
                  data = ., family = "binomial", weights = prop, method = "REML") -> gammod
            
            # d_all %>% 
            #   qgam(snow ~  s(doy, bs = "cc"), qu = 0.2, data = .,
            #        argGam = list(knots=list(doy=c(1,365)), 
            #                      weights = d_all$prop, method = "REML")) -> gammod_l
            # 
            # d_all %>% 
            #   qgam(snow ~  s(doy, bs = "cc"), qu = 0.8, data = .,
            #        argGam = list(knots=list(doy=c(1,365)), 
            #                      weights = d_all$prop, method = "REML")) -> gammod_h
            
          }, silent = T)
          
          if(!class(e)[1] == "try-error"){
            gamsum <- summary(gammod)
            preds <- predict(gammod, pred_grid, type = "response", se.fit = T)
            pred_grid$pred <- preds$fit
            pred_grid$ci <- preds$se.fit*1.96
            # pred_grid$pred_l <- (preds$fit - preds$se.fit*1.96)
            # pred_grid$pred_h <- (preds$fit + preds$se.fit*1.96)
            
            # pred_grid$pred <- predict(gammod, pred_grid, type = "response", exclude = paste0("s(doy):year",2014:2023))
            results$gamr2 <- gamsum$r.sq
            
            results$scd <- round(sum(pred_grid$pred),1)
            # results$scd_ci <- round(sum(pred_grid$ci),3)
            
            results$snow_days <- round(sum(pred_grid$pred >= 0.5),1)
            
            pred_grid <- bind_rows(pred_grid, pred_grid) %>% 
              mutate(doy = row_number())
            
            pred_grid <- pred_grid %>% 
              mutate(mv_mean = rollmean(pred, k = 30, align = "center", na.pad = T))
            
            max_doy_orig <- which.max(pred_grid$mv_mean)
            max_doy_orig <- ifelse(max_doy_orig > 365, max_doy_orig - 365, max_doy_orig)
            
            min_doy_orig <- which.min(pred_grid$mv_mean)
            min_doy_orig <- ifelse(min_doy_orig > 365, min_doy_orig - 365, min_doy_orig)
            
            results$max_snow_doy <- max_doy_orig
            results$min_snow_doy <- min_doy_orig
            
            results$max_snow_prop <- round(max(pred_grid$mv_mean, na.rm = T),3)
            results$min_snow_prop <- round(min(pred_grid$mv_mean, na.rm = T),3)
            
            # pred_grid %>%
            #   ggplot(aes(x = doy)) +
            #   geom_line(aes(y = pred), color= "black") +
            #   geom_line(aes(y = pred + ci), color= "red") +
            #   geom_line(aes(y = pred - ci), color= "blue")
            
            if(min(pred_grid$pred) < 0.5 & max(pred_grid$pred) > 0.5){
              
              pred_grid <- pred_grid %>% 
                slice(max_doy_orig:nrow(.)) %>% 
                mutate(doy = row_number())
              
              max_doy <- which.max(pred_grid$mv_mean)
              min_doy <- which.min(pred_grid$mv_mean)
              
              results$melt <- which.max(pred_grid$pred[max_doy:730] < 0.5) + max_doy-1 + max_doy_orig-1
              # results$melt_ci <- pred_grid$ci[which.max(pred_grid$pred[max_doy:730] < 0.5) + max_doy-1]
              
              results$ns <- which.max(pred_grid$pred[min_doy:730] > 0.5) + min_doy-1 + max_doy_orig-1
              # results$ns_ci <- pred_grid$ci[which.max(pred_grid$pred[min_doy:730] > 0.5) + min_doy-1]
              
            }
          }
        }
      }
      return(results)
    }
    
    # results <- parallel::mclapply(mm %>% mutate(cell2 = cell) %>% nest(data = -cell2) %>% pull(data), cal_scd, mc.cores = future::availableCores())
    system.time({
      options(show.error.messages = FALSE)
      results <- mclapply(mm %>% 
                            mutate(cell2 = cell) %>% 
                            nest(data = -cell2) %>% 
                            pull(data), 
                          cal_scd, mc.cores = workers, extradf = extradf
      )
      options(show.error.messages = TRUE)
      results <- bind_rows(results) %>% 
        mutate(melt = ifelse(melt > 365, melt - 365, melt),
               ns = ifelse(ns > 365, ns - 365, ns))
    })
    
    gc()
    
    rs <- lapply(names(results)[-1], function(x){
      rtemp[] <- results %>% pull(x)
      names(rtemp) <- x
      # plot(rtemp, main = x)
      return(rtemp)
    })
    
    rs <- rast(rs)
    # plot(rs[["scd"]])
    # 
    esalc <- rast(paste0(predictor_dir, "/ESALC.tif"))/10
    esalc <- aggregate(esalc, 3, getmode)
    rs[["esalc"]] <- project(esalc, rs, method="near")
    
    if(file.exists(paste0(predictor_dir, "/GLIMS.tif"))){
      glaciers <- rast(paste0(predictor_dir, "/GLIMS.tif"))
      glaciers <- project(glaciers, rs, method="near")
    } else {
      glaciers <- rs[["esalc"]]
      glaciers[] <- 0
    }
    rs[["GLIMS"]] <- glaciers
    
    dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
    rs[["elevation"]] <- project(dem, rs)
    
    mi <- rast(paste0(predictor_dir, "/medianindices.tif"))
    rs[["ndvi"]] <- project(mi$ndvi, rs)
    
  } else {
    print("Less than 50 individual images! Snow variables not calculated.")
  }
  
  return(rs)
}

read_classifications <- function(imagedf, basedir){
  
  rs <- lapply(imagedf$file, function(x){
    
    r <- rast(paste0(basedir,"/",x))
    
    cloudm <- r[["class"]]
    cloudm[r[["class"]] > 3] <- NA
    cloudm <- focal(cloudm, 3, min, na.rm = F, fillvalue=1, expand=F)
    snow <- r[["snow"]]
    nosnow <- r[["land"]] + r[["water"]]
    snow <- mask(snow, cloudm)
    nosnow <- mask(nosnow, cloudm)
    snow[snow < 20] <- NA
    nosnow[nosnow < 20] <- NA
    
    rs <- c(nosnow, snow)
    names(rs) <- c("nosnow","snow")
    
    return(rs)
  })
  
  names(rs) <- imagedf$file
  return(rs)
  
}

nosnow_values <- function(x){
  r <- x
  x[["nosnow"]][!is.na(x[["nosnow"]])] <- 0
  x[["snow"]][!is.na(x[["snow"]])] <- 1
  m <- values(c(x,r), dataframe = T) %>% 
    set_names(c("nosnow","snow","nosnow.1","snow.1")) %>% 
    rownames_to_column("cell") %>% 
    mutate(cell = as.integer(cell)) %>% 
    rename(nosnow_prop = nosnow.1,
           snow_prop = snow.1)
  return(m)
}

# function(d_all){
#   # cellid <- 7802
#   # cellid <- 3961
#   # d_all <- mm %>% mutate(cell2 = cell) %>% nest(data = -cell2) %>% slice(102) %>% pull(data)
#   # d_all <- d_all[[1]]
#   # d_all <- mm %>% filter(cell == cellid)
#   cellid <- d_all$cell[1]
#   
#   d_all <- bind_rows(d_all %>% select(nosnow, nosnow_prop, doy, year, month, week) %>% drop_na() %>% 
#                        rename(snow = nosnow, prop = nosnow_prop),
#                      d_all %>% select(snow, snow_prop, doy, year, month, week) %>% drop_na() %>% 
#                        rename(prop = snow_prop))
#   
#   
#   results <- tibble(cell = cellid,
#                     nobs = d_all %>% select(doy, year) %>% distinct() %>% nrow,
#                     nyears = length(unique(d_all$year)),
#                     gamr2 = as.numeric(NA), scd = as.numeric(NA),
#                     mean_melt = as.numeric(NA), early_melt = as.numeric(NA), late_melt = as.numeric(NA),
#                     mean_ns = as.numeric(NA), late_ns = as.numeric(NA), early_ns = as.numeric(NA))
#   
#   if(nrow(d_all) > 20){
#     
#     wm <- weighted.mean(d_all$snow, d_all$prop)
#     
#     if(wm < 0.01 | wm > 0.99 | is.nan(wm)){
#       
#       results$scd <- wm
#       
#     } else {
#       # pred_grid <- expand.grid(doy = 1:365,
#       #                          year = unique(d_all$year))
#       pred_grid <- expand.grid(doy = 1:365)
#       
#       e <- try({
#         # d_all %>% 
#         #   mutate(year = as_factor(year)) %>% 
#         #   gam(snow ~  s(doy, bs = "cc") + s(doy, by = year), 
#         #       data = ., family = "binomial", weights = prop, method = "REML") -> gammod
#         d_all %>% 
#           mutate(year = as_factor(year)) %>% 
#           gam(snow ~  s(doy, bs = "cc"), 
#               knots=list(doy=c(1,365)),
#               data = ., family = "binomial", weights = prop, method = "REML") -> gammod
#         
#         d_all %>% 
#           qgam(snow ~  s(doy, bs = "cc"), qu = 0.2, data = .,
#                argGam = list(knots=list(doy=c(1,365)), 
#                              weights = d_all$prop, method = "REML")) -> gammod_l
#         
#         d_all %>% 
#           qgam(snow ~  s(doy, bs = "cc"), qu = 0.8, data = .,
#                argGam = list(knots=list(doy=c(1,365)), 
#                              weights = d_all$prop, method = "REML")) -> gammod_h
#         
#       }, silent = T)
#       
#       if(!class(e)[1] == "try-error"){
#         gamsum <- summary(gammod)
#         pred_grid$pred <- predict(gammod, pred_grid, type = "response")
#         pred_grid$pred_l <- predict(gammod_l, pred_grid, type = "response")
#         pred_grid$pred_h <- predict(gammod_h, pred_grid, type = "response")
#         # pred_grid$pred <- predict(gammod, pred_grid, type = "response", exclude = paste0("s(doy):year",2014:2023))
#         
#         scd <- mean(pred_grid$pred)
#         
#         # pred_grid %>% 
#         #   mutate(pred_l = ifelse(pred_l < 0, 0, pred_l),
#         #          pred_h = ifelse(pred_h > 1, 1, pred_h)) %>% 
#         #   ggplot(aes(x = doy)) +
#         #   geom_line(aes(y = pred), color= "black") +
#         #   geom_line(aes(y = pred_l), color= "red") +
#         #   geom_line(aes(y = pred_h), color= "blue")
#         
#         pred_grid <- pred_grid %>% 
#           mutate(mv_mean = rollmean(pred, k = 20, align = "right", na.pad = T)) %>% 
#           mutate(mv_mean_l = rollmean(pred_l, k = 20, align = "right", na.pad = T)) %>% 
#           mutate(mv_mean_h = rollmean(pred_h, k = 20, align = "right", na.pad = T))
#         
#         pred_grid <- bind_rows(pred_grid, pred_grid) %>% 
#           mutate(doy = row_number())
#         
#         max_doy_orig <- which.max(pred_grid$mv_mean)
#         pred_grid <- pred_grid %>% 
#           slice(max_doy_orig:nrow(.)) %>% 
#           mutate(doy = row_number())
#         
#         max_doy <- which.max(pred_grid$mv_mean)
#         min_doy <- which.min(pred_grid$mv_mean)
#         max_doy_l <- which.max(pred_grid$mv_mean_l)
#         min_doy_l <- which.min(pred_grid$mv_mean_l)
#         max_doy_h <- which.max(pred_grid$mv_mean_h)
#         min_doy_h <- which.min(pred_grid$mv_mean_h)
#         
#         mean_melt <- which.max(pred_grid$pred[max_doy:730] < 0.5) + max_doy-1 + max_doy_orig-1
#         early_melt <- which.max(pred_grid$pred_l[max_doy_l:730] < 0.5) + max_doy_l-1 + max_doy_orig-1
#         late_melt <- which.max(pred_grid$pred_h[max_doy_h:730] < 0.5) + max_doy_h-1 + max_doy_orig-1
#         
#         mean_ns <- which.max(pred_grid$pred[min_doy:730] > 0.5) + min_doy-1 + max_doy_orig-1
#         late_ns <- which.max(pred_grid$pred_l[min_doy_l:730] > 0.5) + min_doy_l-1 + max_doy_orig-1
#         early_ns <- which.max(pred_grid$pred_h[min_doy_h:730] > 0.5) + min_doy_h-1 + max_doy_orig-1
#         
#         results <- tibble(cell = cellid,
#                           nobs = d_all %>% select(doy, year) %>% distinct() %>% nrow,
#                           nyears = length(unique(d_all$year)),
#                           gamr2 = gamsum$r.sq,
#                           scd,
#                           mean_melt, early_melt, late_melt,
#                           mean_ns, late_ns, early_ns)
#         
#       }
#     }
#   }
#   return(results)
# }
