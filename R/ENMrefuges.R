
#' @title Species Refuge Richness
#' @description
#'    根据适生区/非适生区栅格计算一组物种的避难所丰富度.
#'
#' @param x 数值向量, 表示要纳入分析的物种序号, 与d中的物种序号一致. 默认使用全部物种.
#' @param booleandir 物种的存在/不存在栅格文件所在路径. 文件名要包含d中指定的物种.
#' @param d 数据框或字符向量, 提供要计算的物种名称. 为数据框时, 第一列为物种名. (物种名称要与booleandir下的表示物种存在/不存在栅格文件名保持一致)
#' @param key 列表, 用于提取booleandir中对应的栅格进行分析. 例如分别计算不同时期的物种丰富度.
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行时cpu个数

#' @return 返回表示避难所区域的Boolean图

#' @export
#' @importFrom stringr str_detect
#' @examples
#' ENMrefuges(d = read.csv("F:/eblf/TBlabENM/spdata.csv"),
#'            x = c(3,5,7,9),
#'            booleandir = "F:/eblf/TBlabENM/BooleanmapT_MTSStr",
#'            key = list("future"= c("Present", "2050RCP8.5", "2070RCP8.5"),
#'                       "past"= c("Present", "Mid-Holocene", "Last-Glacial-Maximum")),
#'            overwrite = T,
#'            outdir = "F:/4",
#'            parallel = T,
#'            ncpu = 2)

ENMrefuges <- function(d, x = NULL, booleandir, key = NULL,
                       overwrite = FALSE, outdir = NULL,
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
  if (file.exists(booleandir) == FALSE) {
    stop(paste0("Dir '", booleandir, "' not exists."))}
  #创建保存路径
  if (is.null(outdir)) {outdir = "."}
  dir.create(paste0(outdir, "/refuges"), showWarnings = FALSE, recursive = TRUE)


  #计算单个物种的避难所
  refuges <- function(x){ #x为物种序号
    ra <- b[stringr::str_detect(b, d[x,1])]
    raster <- sum(terra::rast(ra))
    raster[raster < length(ra)] <- 0
    raster[raster == length(ra)] <- 1
    terra::writeRaster(raster, paste0(outdir, "/TBlabENMtemp", random_num, "/", d[x,1], ".tif") ) }

  #读取单个物种单个时期的二值图结果文件夹,该文件夹下不要有其他的无关tif
  ##读取所有文件路径
  booleanlist <- list.files(booleandir, full.names = TRUE, pattern = "tif$|TIF$")
  if (length(booleanlist) == 0) {stop("No tif files in the booleandir directory.")}

  ##选择指定的时期的所有二值图

  if (is.null(key)) {b <- booleanlist} else {
    for (m in 1:length(key)) {#i
      #创建临时文件夹
      random_num <- sample(1:100000, 1)
      dir.create(paste0(outdir, "/TBlabENMtemp", random_num), showWarnings = FALSE, recursive = TRUE)

      a <- key[[m]] #a为不同时期的文件，比如present-2030-2050
      #提取所有a中的栅格
      b <- c()
      for (j in a) {
        b1 <- booleanlist[stringr::str_detect(booleanlist, j)]
        b <- c(b,b1)} #b为包含key[[i]]的所有时期的栅格
      #按物种名提取
      if (parallel == TRUE) {
        # 开启集成
        snowfall::sfInit(parallel = TRUE, cpus = ncpu)
        # 注册每个环境变量
        #snowfall::sfExportAll()
        #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
        snowfall::sfExport("x")
        snowfall::sfLapply(x, refuges)
        snowfall::sfStop()  # 关闭集群
      } else {
        for (i in x) {
          refuges(i)
        }
      }

      #每个物种的避难所相加得到避难所丰富度
      ra <- list.files(paste0(outdir, "/TBlabENMtemp", random_num), full.names = TRUE)
      raster <- sum(terra::rast(ra))
      print(paste0(outdir, "/refuges/", names(key)[m], ".tif"))
      terra::writeRaster(raster, paste0(outdir, "/refuges/", names(key)[m], ".tif"), overwrite = overwrite)
      unlink(paste0(outdir, "/TBlabENMtemp", random_num), recursive = TRUE)
    }

  }
  end_time <- Sys.time()
  print(end_time - start_time)
}
