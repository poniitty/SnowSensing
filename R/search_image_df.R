search_image_df <- function(site_name, base_landsat_dir, workers){
  
  area_landsat_dir <- paste0(base_landsat_dir,"/",site_name)
  
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  
  lss <- tibble(area = site_name,
                file = tifs)
  lss$collection <- gsub("0","C",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
  lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
  lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
  lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),1,3))
  lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),4,6))
  lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
  lss$time <- gsub(".tif","",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
  
  lccs <- mclapply(lss$file, calc_coverages, image_dir = area_landsat_dir, mc.cores = workers) %>% 
    bind_rows
  
  lss <- full_join(lss,
                        lccs, by = "file")
  
  lss <- lss %>% arrange(date)
  
  return(lss)
}
