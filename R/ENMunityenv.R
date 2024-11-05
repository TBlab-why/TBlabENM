
#' Environment variable processing
#' @description
#' 将不同来源的环境变量转换为格式、分辨率、范围统一的变量，用于生态位建模。目前仅
#' 支持tif和asc格式变量的处理。
#' @param envdir 要处理的原始环境变量集所在路径，该路径下包括代表不同时期或地理空间的子文件夹，
#' 每个子文件夹内存放对应的环境变量。
#' @param range SpatRaster, SpatVector
#' @param proname 不同时期或地理空间的名称，与envdir下的子文件夹名称对应
#' @param factors 分类变量名称
#' @param method 连续变量重采样的方法
#' @param crop 是否将原始变量裁剪至研究区范围
#' @param res 数值型，当range为SpatVector时输出栅格的分辨率
#' @param outformat 输出文件格式,可选asc,tif或二者组合
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行使用的cpu数
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map map2 map_chr
#' @importFrom terra crs project rast crop mask writeRaster
#' @importFrom stringr str_split_1
#' @importFrom tibble as_tibble
#'
#' @return 返回tif或asc文件
#' @export
#'
#' @examples
#' ENMunityenv(envdir = "F:/example/TBlabENM/env1",
#'             proname <- c("present", "lgm"),
#'             factors = NULL,
#'           # range <- rast("C:/Users/why/TBlabENM/inst/extdata/envar/tif/bio1.tif"),
#'             range <- vect("D:/Desktop/ENMdata/Chinamap/2022年省界/sheng2022.shp") ,
#'             crop = T,
#'             method  = "bilinear",
#'             res = 10000,
#'             outformat = "tif",
#'             outdir = "F:/example",
#'             parallel = F,
#'             ncpu = 2)
ENMunityenv <- function(envdir,
                        proname,
                        factors = NULL,
                        range,
                        crop = T,
                        method = "bilinear",
                        res = NULL,
                        outformat = "tif",
                        overwrite = F,
                        outdir = NULL,
                        parallel = F,
                        ncpu = 2){

  if(is.null(outdir)){outdir <- "."}
#当输入为栅格
  fun1 <- function(j){ #j为单个变量,j=2
    #对单行进行分析
       df <- df1[j,]
      # 读取栅格数据
        ra <- rast(df$value)
        names(ra) <- df$envname
      #将ra的crs转为range
        pro <- terra::project(range, ra, method = "near", align = T, threads = T)
        pro1 <- terra::project(ra, crs(range), method = df$method,res =1000, align = T, threads = T)
      #如果crop == T，则裁剪，否则不裁剪，保存至raster
        if(crop == T){
          if((ext(pro)[1] < ext(pro)[1] & ext(pro)[2] > ext(range)[2]&
              ext(pro)[3] < ext(pro)[3]& ext(pro)[4] > ext(range)[4])==FALSE){
            stop("The range extends beyond the grid to be processed")}
          ra <- crop(ra, pro, mask = T) %>% mask(., pro)
          }
        #将投影转为range投影
        pro <- project(ra, range, method = 'near', align = T, threads = T)
      # 保存栅格数据到同一文件夹
        for (i in 1:length(outformat)) {
          writeRaster(pro, paste0(df$path, ".",outformat[i]), overwrite=overwrite)}
}

  for (i in 1:length(proname)) {#i=1
    #创建对应时期的文件夹
    dir.create(paste0(outdir, "/TBlabENM/env/",proname[i]), recursive = TRUE, showWarnings = FALSE)

    #提取每个时期的原始变量
    df1 <- list.files(paste0(envdir, "/", proname[i]), recursive = T, full.names = T, pattern = "asc$|tif$") %>%
      as_tibble() %>%
      #变量名
      mutate(envname = map_chr(.x = value, .f = function(x){
        name <- str_split_1(x, "/")[length(str_split_1(x, "/"))]
        name1 <-  str_split_1(name, "\\.")[length(str_split_1(name, "\\."))-1]
      })) %>%
      #保存路径
      mutate(path = map_chr(.x = envname, .f = function(x){
        paste0(outdir, "/TBlabENM/env/",proname[i],"/", x)
        })) %>%
      #method
      mutate(method = map_chr(.x = envname, .f = function(x){
           if(x %in% factors){"near"} else {method}
        }))

    if(parallel == TRUE){
      ncpu = ncpu
      # 开启集成
      snowfall::sfInit(parallel = TRUE, cpus = ncpu)
      # 注册每个环境变量
      #snowfall::sfExportAll()
      #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
      snowfall::sfExport("radf")
      snowfall::sfLibrary(tidyverse)
      snowfall::sfLapply(1:length(df), fun1)
      snowfall::sfStop()  # 关闭集群
    } else { for(j in 1:nrow(df1)) { fun1(j) } }


}
  }

