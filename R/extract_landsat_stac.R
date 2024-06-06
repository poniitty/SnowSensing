
make_vsicurl_url <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=landsat-c2-l2",
    "&url=",
    base_url
  )
}

extract_landsat_stac <- function(aoi,
                                epsg,
                                excl_dates,
                                site_name,
                                sats = c("LT04","LT05","LE07","LC08","LC09"),
                                start_date, 
                                end_date, 
                                months = 1:12,
                                minclouds,
                                area_landsat_dir,
                                workers){
  
  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  
  sats2 <- sats %>% 
    ifelse(. == "LT04", "landsat-4", .) %>% 
    ifelse(. == "LT05", "landsat-5", .) %>% 
    ifelse(. == "LE07", "landsat-7", .) %>% 
    ifelse(. == "LC08", "landsat-8", .) %>% 
    ifelse(. == "LC09", "landsat-9", .)
  
  it_obj <- s_obj %>% 
    stac_search(collections = "landsat-c2-l2",
                bbox = st_bbox(aoi %>% st_transform(4326)),
                datetime = paste0(start_date,"/",end_date),
                limit = 1000) %>%
    get_request()
  
  if(length(it_obj$features) == 1000){
    dr <- tibble(date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day")) %>% 
      mutate(gr = cut_number(row_number(.), 10)) %>% 
      group_by(gr) %>% 
      group_split()
    
    it_obj <- s_obj %>% 
      stac_search(collections = "landsat-c2-l2",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  datetime = paste0(min(dr[[1]]$date),"/",max(dr[[1]]$date)),
                  limit = 1000) %>%
      get_request()
    
    for(i in 2:length(dr)){
      # i <- 2
      it_obj2 <- s_obj %>% 
        stac_search(collections = "landsat-c2-l2",
                    bbox = st_bbox(aoi %>% st_transform(4326)),
                    datetime = paste0(min(dr[[i]]$date),"/",max(dr[[i]]$date)),
                    limit = 1000) %>%
        get_request()
      
      it_obj$features <- c(it_obj$features, it_obj2$features)
      
    }
    
  }
  
  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) {x$properties$`eo:cloud_cover` < minclouds}) %>%
    items_filter(filter_fn = function(x) {x$properties$platform %in% sats2}) %>% 
    assets_select(asset_names = c("coastal","blue","green","red","nir08","swir16","swir22","lwir","lwir11","qa_pixel","qa_radsat")) %>%
    items_filter(filter_fn = function(x) {month(ymd_hms(x$properties$datetime)) %in% months}) %>%
    items_filter(filter_fn = function(x) {!as_date(ymd_hms(x$properties$datetime)) %in% excl_dates$date})
  
  if(length(it_obj$features) > 0){
    
    itst <- items_as_sf(it_obj)
    itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
    
    itst$area <- as.numeric(st_area(itst))/(1000*1000)
    
    tt <- itst %>% 
      mutate(date = as_date(ymd_hms(datetime))) %>% 
      group_by(date, `landsat:wrs_path`, .add = TRUE) %>% 
      mutate(n = n()) %>% 
      arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()
    
    tt <- lapply(tt, function(x){
      # x <- tt[[1]]
      if(nrow(x) > 1){
        if(diff(x %>% pull(area)) == 0){
          x <- x %>% slice(1)
        } else {
          stint <- st_intersection(x) %>% slice(2) %>% st_area %>% as.numeric()
          x <- x %>% filter(area > stint/(1000*1000)+1)
        }
      }
      return(x)
    }) %>% bind_rows()
    
    itst <- itst %>% 
      filter(`landsat:scene_id` %in% (tt %>% pull(`landsat:scene_id`)))
    
    it_obj <- it_obj %>%
      items_filter(filter_fn = function(x) {x$properties$`landsat:scene_id` %in% (itst %>% pull(`landsat:scene_id`))})
    
    juuh <- mclapply(it_obj$features, function(ft){
      # ft <- it_obj$features[[1]]
      
      nm <- paste0(paste(str_split(ft$id, "_")[[1]][1:4], collapse = "_"), "_",
                   gsub("-","",substr(ft$properties$created, 1, 10)), "_",
                   paste(str_split(ft$id, "_")[[1]][5:6], collapse = "_"), "_",
                   gsub(":","",substr(ft$properties$datetime, 12, 19)), "GMT.tif")
      
      if(!file.exists(paste0(area_landsat_dir,"/", nm))){
        # print(ft$id)
        full_url <- make_vsicurl_url(assets_url(ft) %>% sort)
        file_names <- gsub("TIF$","tif",basename(full_url))
        
        juuh <- lapply(seq_len(length(full_url)), function(nr){
          e <- try({
            gdal_utils(
              "warp",
              source = full_url[[nr]],
              destination = paste0(tempdir(),"/",file_names[[nr]]),
              options = c(
                "-t_srs", st_crs(aoi)$wkt,
                "-te", st_bbox(aoi),
                "-tr", c(30, 30)
              )
            )
          }, silent = TRUE)
          if(class(e)[[1]] == "try-error"){
            return(FALSE)
          } else {
            return(TRUE)
          }
        })
        
        err <- file_names[!unlist(juuh)]
        if(length(err) > 0){
          ll <- lapply(err, function(xx){
            r <- rast(paste0(tempdir(),"/",file_names[unlist(juuh)][1]))
            r[] <- NA
            writeRaster(r, paste0(tempdir(),"/",xx), datatype = "INT2U", overwrite = TRUE)
          })
        }
        
        r <- rast(c(paste0(tempdir(),"/",file_names[grepl("_SR_", file_names)]),
                    paste0(tempdir(),"/",file_names[grepl("_ST_", file_names)]),
                    paste0(tempdir(),"/",file_names[grepl("_QA_", file_names)])))
        names(r) <- lapply(names(r), function(x) paste(str_split(x, "_")[[1]][8:9], collapse = "_")) %>% unlist
        
        writeRaster(r, paste0(area_landsat_dir,"/", nm), datatype = "INT2U", overwrite = TRUE)
        
        unlink(c(paste0(tempdir(),"/",file_names)))
      }
    }, mc.cores = workers)
    
    itst$platform <- itst$platform %>% 
      ifelse(. == "landsat-4", "LT04", .) %>% 
      ifelse(. == "landsat-5", "LT05", .) %>% 
      ifelse(. == "landsat-7", "LE07", .) %>% 
      ifelse(. == "landsat-8", "LC08", .) %>% 
      ifelse(. == "landsat-9", "LC09", .)
    
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
    
    idsALL <- itst %>% 
      mutate(date = gsub("-","", as_date(ymd_hms(datetime)))) %>% 
      mutate(date2 = gsub("-","", as_date(ymd_hms(created)))) %>% 
      mutate(time = paste0(gsub(":","", substr(datetime, 12, 19)), "GMT.tif")) %>% 
      mutate(file = paste0(platform, "_", `landsat:correction`, "_",
                           `landsat:wrs_path`, `landsat:wrs_row`, "_",
                           date, "_", date2, "_",
                           `landsat:collection_number`, "_", `landsat:collection_category`, "_",
                           time)) %>% 
      select(file, date, `eo:cloud_cover`, `view:sun_elevation`) %>% 
      rename(DATE_ACQUIRED = date,
             CLOUD_COVER = `eo:cloud_cover`, 
             SUN_ELEVATION = `view:sun_elevation`) %>% 
      st_drop_geometry() %>% 
      filter(file %in% tifs) %>% 
      mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>% 
      rownames_to_column("id")
  } else {
    print("No new Landsat scenes downloaded!")
    idsALL <- NULL
  }
  
  return(idsALL)
}
