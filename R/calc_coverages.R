mask_fill <- function(x) {as.numeric(intToBits(x)[1])}
mask_dilcloud <- function(x) {as.numeric(intToBits(x)[2])}
mask_cloud <- function(x) {as.numeric(intToBits(x)[4])}
mask_cshadow <- function(x) {as.numeric(intToBits(x)[5])}
mask_snow <- function(x) {as.numeric(intToBits(x)[6])}
mask_clear <- function(x) {as.numeric(intToBits(x)[7])}
mask_water <- function(x) {as.numeric(intToBits(x)[8])}
conf_cloud <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[9:10]), collapse = ""))}
conf_cshadow <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[11:12]), collapse = ""))}
conf_snow <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[13:14]), collapse = ""))}
conf_cirrus <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[15:16]), collapse = ""))}

# image <- lss$file[1]
calc_coverages <- function(image, image_dir){
  # image <- lss$file[[5]]
  require(terra)
  
  rs <- rast(paste0(image_dir,"/",image))
  rsn <- names(rs)
  # plot(rs)
  rs[[4]][is.na(rs[[4]])] <- 0
  
  cmask <- rs[[1]]
  cmask[] <- unlist(lapply(as.numeric(values(rs[[rsn[grepl("_pixel",rsn, ignore.case = T)]]])), mask_cloud))
  
  fill_prop <- mean(values(rs[[4]]) == 0)
  rs[[4]][cmask == 1] <- NA
  cloud_prop <- mean(is.na(values(rs[[4]])))
  
  df <- tibble(file = image,
               fill_proportion = fill_prop,
               cloud_proportion = cloud_prop,
               clear_proportion = 1-(fill_prop+cloud_prop))
  
  return(df)
  
}
