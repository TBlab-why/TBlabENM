
#' Environment variable processing
#' @description
#' 将不同来源的环境变量进行格式转换、重新投影到指定的研究区.
#'
#' @param radir 要处理的环境变量所在的路径
#' @param ref SpatRaster或SpatVector,作为参考研究区，如果是SpatRaster裁剪后的变量范围和分辨率与
#'     SpatRaster一致,如果是SpatVector,裁剪后的变量范围与SpatVector一致，分辨率由res参数指定。
#' @param proname 投影时期的名字，必须是radir下的子文件夹的名字
#' @param factors 分类变量名，要与radir下的变量名称一致
#' @param method 连续变量的重采样方法，对于分类变量内部设置为near，可选bilinear,cubic,near等
#' @param format 输出变量格式
#' @param outdir 输出文件路径
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param parallel 是否并行
#' @param ncpu 并行的cpu数
#' @param res numeric. Can be used to set the resolution of the output raster if ref is a SpatVector
#'
#' @importFrom purrr map_chr map_int map
#' @importFrom dplyr mutate
#' @return 栅格文件
#' @export
#'
#' @examples
#' ENMunityenv(radir = "F:/example/env",
#'             ref = rast("F:/var2/eblf_proj/present/tif/bio1.tif"),
#'             proname = c("present", "acc2030ssp126","acc2030ssp245"),
#'             outdir = "F:/example/env",
#'             overwrite = T,
#'             parallel = T,
#'             ncpu = 2)
#'
ENMunityenv <- function(radir, ref, proname = NULL, factors = NULL, method = "bilinear", res,
                 format = "tif", outdir = NULL, overwrite = F, parallel = F, ncpu = 2){
  dir.create(paste0(outdir, "/env/"), recursive = TRUE, showWarnings = FALSE)
  if(is.null(outdir)){outdir <- "."}

  #创建栅格列表
  if(is.null(proname)){
    ralist <- list.files(radir, full.names = TRUE, pattern = ".asc$|.tif$|.TIF$|.ASC$")} else {
      ralist <- c()
      for (i in seq_along(proname)) {
          ralist1 <- list.files(paste0(radir, "/", proname[i]),
                               full.names = T,  pattern = ".asc$|.tif$|.TIF$|.ASC$")
          ralist <- c(ralist, ralist1)
    }}
   #构建栅格列表数据框
    radf <- as.data.frame(ralist) %>%
      mutate(name = map_chr(.x = ralist, .f = function(x){
        a <- stringr::str_split_1(x, "/")[length(stringr::str_split_1(x, "/"))]
        a <- stringr::str_split_1(a, ".asc$|.tif$|.TIF$|.ASC$")[1]
      })) %>%
      mutate(method = map_chr(.x = name, .f = function(x){
        if (x %in% factors) {method <- "near"} else {method <- method}
      })) %>%
      mutate(factor = map_int(.x = name, .f = function(x){
        if (x %in% factors) {factor <- 1} else {method <- 0}
      })) %>%
      mutate(proname = map_chr(.x = ralist, .f = function(x){
        if (is.null(proname)) {
          a <- NA
        }else{a <- stringr::str_split_1(x, "/")[length(stringr::str_split_1(x, "/"))-1]}

      }))

    #检查要处理的变量是否含有factors
    if (is.null(factors) == FALSE){
      if (sum(radf$factor) == 0){
        stop("The specified categorical variable was not found. Please check parameter factors!/n")}}


    fun1 <- function(x){
      ra <- terra::rast(x)
      if (terra::crs(ref, proj = TRUE) == terra::crs(ra, proj = TRUE)) { #投影相同的情况
        if (class(ref) == "SpatRaster") {
          ra_r <- terra::crop(ra, ref, mask = T) %>% terra::mask(ref)} else {
            ra_r <- terra::crop(ra, ref)
          }
        if (terra::ext(ra_r) != terra::ext(ref)) {
          if (class(ref) == "SpatRaster") {
            ra_r <- terra::resample(ra_r, ref, "near") }
                                                    }
         } else{ #投影不同的情况

          if (class(ref) == "SpatRaster") {
            ra <- terra::project(ra, ref, method = radf$method[which(x == radf[1])]) } else {
            ra <- terra::project(ra, terra::crs(ref), method = radf$method[[which(x == radf[1])]], res = res) }

          ra_r <- terra::crop(ra, ref, mask = T) %>% terra::mask(ref)
          if (terra::ext(ra_r) != terra::ext(ref)) {
            ra_r <- terra::resample(ra_r, ref, "near")
          }

          }


      if(is.null(proname)){
        terra::writeRaster(
          ra_r, paste0(outdir, "/env/",radf$name[which(x==radf[1])], ".", format),
                        NAflag = -9999, overwrite = overwrite) } else {
          dir.create(paste0(outdir, "/env/", radf$proname[which(x==radf[1])]),
                     recursive = TRUE, showWarnings = FALSE)
        terra::writeRaster(ra_r, paste0(outdir, "/env/",radf$proname[which(x==radf[1])], "/",
                                        radf$name[which(x==radf[1])], ".", format),
                                         NAflag = -9999, overwrite = overwrite) }
    }
  if (parallel == TRUE) {
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
   # snowfall::sfExport("star_time")
    snowfall::sfLibrary(purrr)
    k <- snowfall::sfLapply(radf$ralist, fun1)
    snowfall::sfStop()  # 关闭集群

  } else {

    for (i in 1:nrow(radf)) {
      fun1(radf[i,1])

    }
  }

  }


