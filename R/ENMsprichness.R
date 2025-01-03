
#' @title Species Richness Analysis
#' @description 根据布尔图进行叠加，计算物种丰富度。
#'
#' @param x 数值向量。表示要进行转换的物种的序号。这个序号来自ENMspdata()生成的标准文件。
#' @param spdata spdata 由ENMspdata()生成的标准文件。也可手动创建或修改，但是不建议。
#' 当maxent结果文件出现变动后，强烈建议重新运行ENMspdata()。
#' @param proname 想要计算的物种丰富度的时期名称.
#' @param booleandir 物种单个时期的布尔图所在路径。
#' @param path 输出结果路径。
#' @param bioindex 物种丰富度指标，有sr（物种丰富度），we（加权物种丰富度），cwe（特有加权物种丰富度）
#'
#' @return 物种丰富度图，.tif格式。
#' @importFrom terra writeRaster
#' @importFrom dplyr mutate
#' @importFrom purrr map map2
#' @importFrom stringr str_detect

#'
#' @param parameters 数据框，第一列为物种名(要与resultdir下的名称保持一致),其他列位置不做要求，但是必须有一列或多列
#'     为划分适生区和非适生区的阈值。
#' @param x 数值向量,表示物种序号,与parameters中的物种序号一致。
#' @param booleandir 适生区和非适生区结果文件路径。
#' @param key 要提取的关键词，一般对应投影时期
#' @param bioindex 丰富度指数，可选sr,we,cwe
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行时cpu个数
#'
#' @export
#' @examples
#' ENMsprichness(parameters = read.csv("F:/eblf/TBlabENM/spdata.csv"),
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
ENMsprichness <- function(parameters, x = NULL, booleandir, key = NULL,
                          bioindex = "sr", overwrite = FALSE, outdir = NULL,
                          parallel = F, ncpu = 2) {
  start_time <- Sys.time()
  #unlink(paste0(outdir, "/TBlabENM/TBlabENMtemp3"), recursive = TRUE)
  if (file.exists(booleandir)==FALSE){
    stop("booleandir not find.")}
  #创建保存路径
  if(is.null(outdir)){outdir = "."}
  if(is.null(x)){x <- 1:nrow(parameters)}
  ##读取所有的二值图路径
  rasterlist <- list.files(booleandir, pattern = ".tif$|.TIF$", all.files = F, full.names = T)
  if(length(rasterlist) == 0) {rasterlist <- booleandir}
  #提取想要的时期
  if(is.null(key)==FALSE){
    rasterlist <- unlist(key) |> as.data.frame() |>
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

    if(is.null(key)==FALSE){ rasterlist = rasterlist[[y,2]]}

    names(rasterlist) <- "path"
    rasterlist <- rasterlist$path
    raster <- c()
    for (i in x) {
      raster1 <- rasterlist[str_detect(rasterlist, parameters[i,1])]
      raster <- c(raster, raster1)
    }

    ra <- sum(terra::rast(raster))
    terra::writeRaster(ra, paste0(outdir, "/biodiversity/sr/",bb[b%in%y], ".tif"), overwrite = overwrite)
  }
#计算we

  we <- function(z){

    dir.create(paste0(outdir, "/TBlabENMtemp3/", bb[b%in%z], "/"), showWarnings = FALSE, recursive = TRUE)
    if(is.null(key)==FALSE){
    rasterlist = rasterlist[[z,2]]
    names(rasterlist) <- "path"
    rasterlist <- rasterlist$path }
    raster <- c()
    for (i in x) {
      raster1 <- rasterlist[str_detect(rasterlist, parameters[i,1])]
      raster <- c(raster, raster1)
    }
    df <- raster |>
      as.data.frame() |>
      mutate(path = 1:length(raster)) |>
      mutate(ra = map(.x = ., .f = function(x){
        terra::rast(x)
      })) |>
     mutate(we = map(.x = ra, .f = function(x){
       x*(1/(terra::freq(x)[2,3]))
     })) |>
      mutate(map2(.x = we, .y = path, .f = function(x,y){
        terra::writeRaster(x, paste0(outdir, "/TBlabENMtemp3/",  bb[b%in%z], "/", y, ".tif"), overwrite = overwrite)
      }))

    raster <- list.files(paste0(outdir, "/TBlabENMtemp3/",  bb[b%in%z], "/"), full.names = TRUE)
    ra <- sum(terra::rast(raster))
    terra::writeRaster(ra, paste0(outdir, "/biodiversity/we/",  bb[b%in%z],".tif"), overwrite = overwrite)
    unlink(paste0(outdir, "/TBlabENMtemp3"), recursive = TRUE)
  }


cwe <- function(y){

  cwe <- terra::rast(paste0(outdir, "/biodiversity/we/", bb[b%in%y],".tif"))/terra::rast(paste0(outdir, "/biodiversity/sr/",bb[b%in%y], ".tif"))
  terra::writeRaster(cwe, paste0(outdir, "/biodiversity/cwe/", bb[b%in%y],".tif"), overwrite = overwrite)
}


  if(parallel == TRUE){
    # 开启集成
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
    snowfall::sfExport("b")
    snowfall::sfExport("bb")
    snowfall::sfExport("parameters")
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

    if(bioindex == "sr"){
      dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
    snowfall::sfLapply(b, sr) }
    if(bioindex == "we"){
      dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
      snowfall::sfLapply(b, we) }
    if(bioindex == "cwe"){
  dir.create(paste0(outdir, "/biodiversity/cwe"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
  snowfall::sfLapply(b, sr)
  snowfall::sfLapply(b, we)
  snowfall::sfLapply(b, cwe) }
    snowfall::sfStop()  # 关闭集群
} else{
  if(bioindex == "sr"){
    dir.create(paste0(outdir, "/biodiversity/sr"), showWarnings = FALSE, recursive = TRUE)
    for(i in b) {sr(i)} }
  if(bioindex == "we"){
    dir.create(paste0(outdir, "/biodiversity/we"), showWarnings = FALSE, recursive = TRUE)
    for(i in b) {we(i)} }
  if(bioindex == "cwe"){
    for(i in b) {
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
