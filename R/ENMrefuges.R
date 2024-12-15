
#' @title Individual Species Refuge
#' @description
#'    根据适生区/非适生区栅格计算物种的避难所。
#'
#' @param x 数值向量,表示物种序号,与parameters中的物种序号一致。
#' @param booleandir 适生区和非适生区结果文件路径。
#' @param parameters 数据框，第一列为物种名(要与resultdir下的名称保持一致)
#' @param key 要提取的关键词，一般对应投影时期
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行时cpu个数
#' 本文件夹，请自己指定。
#' @return 返回表示避难所区域的Boolean图

#' @export
#' @importFrom stringr str_detect
#' @examples
#' ENMrefuges(parameters = read.csv("F:/eblf/TBlabENM/spdata.csv"),
#'            x = c(3,5,7,9),
#'            booleandir = "F:/eblf/TBlabENM/BooleanmapT_MTSStr",
#'            key = list("future"= c("Present", "2050RCP8.5", "2070RCP8.5"),
#'                       "past"= c("Present", "Mid-Holocene", "Last-Glacial-Maximum")),
#'            overwrite = T,
#'            outdir = "F:/4",
#'            parallel = T,
#'            ncpu = 2)

ENMrefuges <- function(parameters, x = NULL, booleandir, key = NULL,
                      overwrite = FALSE, outdir = NULL,
                      parallel = F, ncpu = 2) {
  #unlink(paste0(outdir, "/TBlabENM/TBlabENMtemp1"), recursive = TRUE)
 #计算单个物种的避难所
  fun1 <- function(x){
    ra <- b[stringr::str_detect(b, parameters[x,1])]
    raster <- sum(terra::rast(ra))
    raster[raster<length(ra)] <- 0
    raster[raster==length(ra)] <- 1
    terra::writeRaster(raster, paste0(outdir, "/TBlabENMtemp1/", parameters[x,1], ".tif")) }

  star_time <- Sys.time() ## 记录程序开始时间
  #读取单个物种单个时期的二值图结果文件夹,该文件夹下不要有其他的无关tif
  if (file.exists(booleandir)==FALSE){
    stop("booleandir not find.")}

  booleanlist <- list.files(booleandir,  full.names = T, pattern = "tif$|TIF$")
  if(length(booleanlist) == 0) {booleanlist <- booleandir}
  #创建保存路径
  if(is.null(outdir)){outdir = "."}
  if(is.null(x)){x <- 1:nrow(parameters)}
  dir.create(paste0(outdir, "/refuges"), showWarnings = FALSE, recursive = TRUE)

  ##选择指定的时期的所有二值图
  if(is.null(key)){b <- rasterlist} else{
    for (i in 1:length(key)) {
      dir.create(paste0(outdir, "/TBlabENMtemp1"), showWarnings = FALSE, recursive = TRUE)
      a <- key[[i]]
     #提取所有a中的栅格
      b <- c()
      for (j in a) {
        b1 <- booleanlist[stringr::str_detect(booleanlist, j)]
        b <- c(b,b1)}
     #按物种名提取
      if(parallel == TRUE){
        # 开启集成
        snowfall::sfInit(parallel = TRUE, cpus = ncpu)
        # 注册每个环境变量
        #snowfall::sfExportAll()
        #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
        snowfall::sfExport("x")
        snowfall::sfLapply(x, fun1)
        snowfall::sfStop()  # 关闭集群
      } else{
        for (i in x) {
          fun1(i)
        }
      }

      ra <- list.files(paste0(outdir, "/TBlabENMtemp1"), full.names = TRUE)
      raster <- sum(terra::rast(ra))
      terra::writeRaster(raster, paste0(outdir, "/refuges/", names(key)[i], ".tif"), overwrite = overwrite)
      unlink(paste0(outdir, "/TBlabENMtemp1"), recursive = TRUE)
    }

    }
    end_time <- Sys.time()
    print(end_time-star_time)
}

