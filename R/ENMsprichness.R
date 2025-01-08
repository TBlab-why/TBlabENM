
#' @title Species Richness Analysis
#' @description 叠加不同物种的存在/不存在栅格，计算物种丰富度.
#' @param d 数据框或字符向量, 提供要计算的物种名称. 为数据框时, 第一列为物种名. (物种名称要与booleandir下的表示物种存在/不存在栅格文件名保持一致)
#' @param x 数值向量, 表示要纳入分析的物种序号, 与d中的物种序号一致. 默认使用全部物种.
#' @param booleandir 物种的存在/不存在栅格文件所在路径. 文件名要包含d中指定的物种.
#' @param key 列表, 用于提取booleandir中对应的栅格进行分析. 例如分别计算不同时期的物种丰富度.
#' @param bioindex 丰富度指数, 可选sr, we, cwe中的任意一个. 如果为cwe, 同时计算sr和we.
#' @param overwrite logical. If TRUE, filename is overwritten.
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行时cpu个数
#'
#' @importFrom terra writeRaster
#' @importFrom dplyr mutate
#' @importFrom purrr map map2 %>%
#' @importFrom stringr str_detect
#' @export
#'
#' @return 物种丰富度图，.tif格式。
#' @examples
#' ENMsprichness(d = read.csv("F:/eblf/TBlabENM/spdata.csv"),
#'               x = c(3,5,7,9),
#'               booleandir = "F:/eblf/TBlabENM/BooleanmapT_MTSStr",
#'               key = list("2050RCP8.5"= "2050RCP8.5",
#'                          "2070RCP8.5"= "2070RCP8.5",
#'                          "Present"= "Present"),
#'               bioindex = "cwe",
#'               overwrite = T,
#'               outdir = "F:/4",
#'               parallel = T,
#'               ncpu = 2)
ENMsprichness <- function(d, x = NULL, booleandir, key = NULL,
                          bioindex = "sr", overwrite = FALSE, outdir = NULL,
                          parallel = F, ncpu = 2) {

  start_time <- Sys.time()
  #参数检验
  if (is.character(d)) {d <- as.data.frame(d)}
  names(d)[1] <- "species"
  if (is.null(x)) {x <- 1:nrow(d)}
  x <- unique(x)
  if (max(x) > nrow(d)) {
    stop("x provides an invalid species number.")
  }
  if (length(x) == 1) {stop("It has to be more than 2 species, otherwise it doesn't make sense.")}
  if (file.exists(booleandir) == FALSE){
    stop(paste0("Dir '", booleandir, "' not exists."))}
  if (length(bioindex) != 1) {stop("length bioindex must be 1.")}
  if (bioindex %in% c("sr", "we", "cwe") == FALSE) {
      stop("bioindex must be 'sr', 'we' or 'cwe'.")
  }
  if (is.null(key)) {
    warning("All files in the booleandir will be treated as the same period, please confirm that you want to do this.")
    key <- list("all" = ".tif")
    }
  if (is.null(names(key))) {
    names(key) <- key
  }

  #创建保存路径
  if (is.null(outdir)) {outdir = "."}

  ##读取所有的二值图路径
  rasterlist <- list.files(booleandir, pattern = ".tif$|.TIF$", all.files = F, full.names = T)
  if (length(rasterlist) == 0) {stop("No .tif files were found under booleandir.")}
  #如果指定key，则提取想要的时期
  if (is.null(key) == FALSE) {
    rasterlist <- unlist(key) %>% as.data.frame() %>%
      #ra为不同时期的所有路径
      mutate(ra  = map(.x = ., .f = function(x){
        as.data.frame(rasterlist[str_detect(rasterlist, x)])
      }))
    bb <- row.names(rasterlist)
    row.names(rasterlist) <- rasterlist[[1]]
    b <- rasterlist$.
  }
  #计算SR y=b[1]
  sr <- function(y) {

    if (is.null(key) == FALSE) {
      rasterlist = rasterlist[[y,2]]
      if (nrow(rasterlist) == 0) {
        stop(paste0("No file name contains the string '", y, "', please check the key parameter."))}
    }

    names(rasterlist) <- "path"
    rasterlist <- rasterlist$path
    raster <- c()
    for (i in x) {
      raster1 <- rasterlist[str_detect(rasterlist, d[i,1])]
      raster <- c(raster, raster1)
    }

    ra <- sum(terra::rast(raster))
    terra::writeRaster(ra, paste0(outdir, "/biodiversity/sr/", bb[b %in% y], ".tif"), overwrite = overwrite)
  }
#计算we

  we <- function(z){
    #临时文件夹
    random_num <- sample(1:100000, 1)
    dir.create(paste0(outdir, "/TBlabENMtemp", random_num, "/", bb[b %in% z], "/"),
               showWarnings = FALSE, recursive = TRUE)
    if (is.null(key) == FALSE) {
    rasterlist = rasterlist[[z, 2]]
    names(rasterlist) <- "path"
    rasterlist <- rasterlist$path }
    raster <- c()
    for (i in x) {
      raster1 <- rasterlist[str_detect(rasterlist, d[i,1])]
      raster <- c(raster, raster1)
    }
    df <- raster %>%
      as.data.frame() %>%
      mutate(path = 1:length(raster)) %>%
      mutate(ra = map(.x = ., .f = function(x){
        terra::rast(x)
      })) %>%
     mutate(we = map(.x = ra, .f = function(x){
       x*(1/(terra::freq(x)[2,3]))
     })) %>%
      mutate(map2(.x = we, .y = path, .f = function(x,y){
        terra::writeRaster(x, paste0(outdir, "/TBlabENMtemp",random_num, "/",  bb[b %in% z], "/", y, ".tif"), overwrite = overwrite)
      }))

    raster <- list.files(paste0(outdir, "/TBlabENMtemp",random_num, "/",  bb[b %in% z], "/"), full.names = TRUE)
    ra <- sum(terra::rast(raster))
    terra::writeRaster(ra, paste0(outdir, "/biodiversity/we/",  bb[b %in% z],".tif"), overwrite = overwrite)
    unlink(paste0(outdir, "/TBlabENMtemp",random_num), recursive = TRUE)
  }


cwe <- function(y){

  cwe <- terra::rast(paste0(outdir, "/biodiversity/we/", bb[b %in% y],".tif"))/terra::rast(paste0(outdir, "/biodiversity/sr/",bb[b %in% y], ".tif"))
  terra::writeRaster(cwe, paste0(outdir, "/biodiversity/cwe/", bb[b %in% y],".tif"), overwrite = overwrite)
}


  if (parallel == TRUE){
    # 开启集成
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
    snowfall::sfExport("b")
    snowfall::sfExport("bb")
    snowfall::sfExport("d")
    snowfall::sfExport("x")
    snowfall::sfExport("overwrite")
    snowfall::sfExport("outdir")
    snowfall::sfExport("booleandir")
    snowfall::sfExport("key")
    snowfall::sfExport("bioindex")
    snowfall::sfExport("sr")
    snowfall::sfExport("we")
    snowfall::sfExport("cwe")
    snowfall::sfLibrary(tidyverse)

    if (bioindex == "sr") {
      dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
    snowfall::sfLapply(b, sr) }
    if (bioindex == "we") {
      dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
      snowfall::sfLapply(b, we) }
    if (bioindex == "cwe") {
  dir.create(paste0(outdir, "/biodiversity/cwe"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
  snowfall::sfLapply(b, sr)
  snowfall::sfLapply(b, we)
  snowfall::sfLapply(b, cwe) }
    snowfall::sfStop()  # 关闭集群
} else{
  if (bioindex == "sr") {
    dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
    for (i in b) {sr(i)} }
  if (bioindex == "we") {
    dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
    for (i in b) {we(i)} }
  if (bioindex == "cwe") {
    for (i in b) {
      dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
      dir.create(paste0(outdir, "/biodiversity/cwe"), showWarnings = FALSE, recursive = TRUE)
      sr(i)
      we(i)
      cwe(i)} }
}
end_time <- Sys.time()
print(end_time - start_time)

}
