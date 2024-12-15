
#' @title The suitable area was reclassified and the area was calculated
#' @description
#'     依据选定的阈值从maxent的结果文件中将连续型的分布图转为Boolean图
#'    （用1和0表示适生区和非适生区）。自动生成paste0("./TBlabENM/Booleanmap", threshold)文件保存结果。
#' @details
#' rcl:rcl参数提供了栅格数据重分类的参数，如果提供的是列表，重分类后各类的值分别为0,1,2...,
#' 该顺序代表的适生区类型与列表中的顺序一致，例如当rcl = list("非适生区" = c(0,0.2), "低适生区" = c(0.2,0.4),
#' "中适生区" = c(0.4,0.6), "高适生区" = c(0.6,0.1))时，重分类后值为0代表非适生区，值为1代表低适生区...,值为
#' 3代表高适生区.当列表中的顺序不是从小到大的顺序时，对应重分类后的值所对应的适生区类型也随之改变。
#'     当为数值型或字符型时，则进行二元划分，重分类后的值为0和1，0代表非适生区，1代表适生区.当为数值型
#' 时，对所有栅格都使用该值,常用于经验阈值，例如0.5.当为字符型时，使用parameters参数提供的数据框中与之
#' 对应的列的值作为划分依据,常用于固定阈值的情况，如MTSS阈值.
#'
#'
#' @param parameters 数据框，第一列为物种名(要与resultdir下的名称保持一致),其他列位置不做要求，但是必须有一列或多列
#'     为划分适生区和非适生区的阈值。
#' @param x 数值向量,表示物种序号,与parameters中的物种序号一致。
#' @param resultdir maxent结果文件路径。
#' @param prefix 要转换的栅格文件的前缀
#' @param suffix 要转换的栅格文件的后缀
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param outdir 输出文件夹
#' @param parallel 是否并行
#' @param ncpu 并行时cpu个数
#' @param rcl 列表,数值型或字符型,为重分类设置参数.当为列表时,分别按顺序指定非适生区
#'     和适生区(例如低适生区、中适生区、高适生区等)和对应的值,数值型和字符型
#'     是列表的特例，即二元划分：只分为适生区和非适生区。当为数值型，所有小于该值的
#'     都视为非适生区，重分类后值为0，所有大于该值的都视为适生区，重分类后值为1.当为
#'     字符型时，则在parameters参数提供的数据框中要有相同名称的列，指定以哪个值进行
#'     二元重分类。详见detail.
#'
#' @return 生成每个栅格重分类后的tif文件，并统计面积
#' @export
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom purrr map_dbl
#' @importFrom purrr map map_chr pmap
#' @importFrom purrr map2
#' @importFrom purrr map2_chr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate relocate
#'
#'
#' @examples
#' ENMclassify(parameters = read.csv("F:/eblf/TBlabENM/spdata.csv"),
#'            x = c(3,5,7,9),
#'            resultdir="F:/eblf/TBlabENM/result",
#'            prefix = NULL,
#'            suffix = "tif",
#'            rcl = list("非适生区" = c(0,0.2),
#'                       "低适生区" = c(0.2,0.4),
#'                       "中适生区" = c(0.4,0.6),
#'                       "高适生区" = c(0.6,1)),
#'            overwrite = T,
#'            outdir = "F:/4",
#'            parallel = T,
#'            ncpu = 2)

ENMclassify <- function(parameters, x = NULL, resultdir,
                        prefix = NULL, suffix, rcl,
                        overwrite = FALSE, outdir = NULL, parallel = F, ncpu = 2) {
  star_time <- Sys.time() ## 记录程序开始时间

  #根据threshold的类型分析
  threshold <- rcl
  if(is.list(threshold)){
    a <- c()
    for (i in 1:length(rcl)) {
      a1 <- c(rcl[[i]], i-1)
      a <- c(a,a1)}
    parameters$reclass =  list(matrix(a, ncol=3, byrow=TRUE))
  }
  # threshold = 0.2
  if(is.numeric(threshold)){
    a <- c(0, threshold, 0, threshold, 1000, 1)
    parameters$reclass =  list(matrix(a, ncol=3, byrow=TRUE))
  }
  # threshold = "T_MTSSte"
  if(is.character(threshold)){
    parameters$reclass <- NA
    for (i in 1:nrow(parameters)) {
      a <- list(
        matrix(c(0, parameters[i,threshold],0, parameters[i,threshold], 1000, 1),
               ncol=3, byrow=TRUE))
      parameters$reclass[i] <- a
    }
  }


  #读取模拟结果列表(只读文件夹)（以物种为单位）
  if (file.exists(resultdir)==FALSE){
    stop("Resultdir not find.")}
  if(is.null(outdir)){outdir = "."}
  if(is.null(x)){x <- 1:nrow(parameters)}
  spdata <- parameters
  #临时文件夹
  star_time <- sample(1:100000,1)
  dir.create(paste0(outdir, "/TBlabENM", star_time), recursive = T, showWarnings = FALSE)
  #创建文件夹用于保存二值图
  dir.create(paste0(outdir, "/reclass"), recursive = T, showWarnings = FALSE)

  #################################################函数
  fun1 <- function(x){
    spname <- spdata[x,1]

    #读取栅格数据
    ra_df <- list.files(paste0(resultdir, "/", spdata[x,1]), full.names = TRUE) %>%
      list.files(., full.names = TRUE, pattern = paste0("^", prefix, ".*", suffix, "$")) %>%
      as.data.frame() %>%
      mutate(reclass = spdata[x, "reclass"]) %>%
      #生成保存文件名
      mutate(path = map_chr(.x = ., .f = function(x){
        a <- stringr::str_split_1(x, "/")
        paste0(a[length(a)-1], ".tif")
      })) %>%
      #投影名
      mutate(pro = map(.x = path, .f = function(x){
        str_split_1(x, ".tif")[1]
      }))

    ###############################

    radf <- ra_df %>%
      mutate(reclass = map2(.x = ., .y = reclass, .f = function(x,y){
        terra::classify(terra::rast(x), y)
      })) %>%
      #单个像元面积
      mutate(unit = map(.x = reclass, .f = function(x){
        prod(terra::res(x))
      })) %>%
      #统计各个分类值的像元个数*单个像元面积
      mutate(area = map2(.x = reclass, .y = unit, .f = function(x,y){
        count <- terra::freq(x)
        count$count <- count$count*y
        tidyr::spread(count, key = "value",
                      value = "count")
      })) %>%
      mutate(area1 = map2(.x = area, .y = pro, .f = function(x,y){
        x$pro = y
        x$species = spname
        x
      })) %>%
      mutate(ss = map2(.x = reclass, .y = path, .f = function(x,y){
        terra::writeRaster(x,
                           paste0(outdir, "/reclass/", spname, "_", y ),overwrite = overwrite)
      }))

    #新建数据框保存单个物种的结果
    data <- data.frame(matrix(NA, nrow = 0, ncol = ncol(radf[,7][[1]])))
    names(data) <- names(radf[,7][[1]])

    for (i in 1:length(radf[,7])) {
      b1 <- radf[,7][[i]]
      data <- rbind(data,b1)
    }



    data <- data[-1]%>%
      relocate(pro)%>%
      relocate(species)
    write.csv(data, paste0(outdir, "/TBlabENM", star_time, "/", spname,".csv"), row.names = FALSE)

  }

  #并行


  if(parallel == TRUE){
    ncpu = ncpu
    # 开启集成
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    # 注册每个环境变量
    #snowfall::sfExportAll()
    #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
    snowfall::sfExport("resultdir")
    snowfall::sfExport("parameters")
    snowfall::sfExport("x")
    snowfall::sfExport("prefix")
    snowfall::sfExport("suffix")
    snowfall::sfExport("rcl")
    snowfall::sfExport("overwrite")
    snowfall::sfExport("outdir")
    snowfall::sfExport("fun1")

    snowfall::sfLibrary(tidyverse)

    snowfall::sfLapply(x, fun1)
    snowfall::sfStop()  # 关闭集群
  } else {
    ## 第一个位置：新建一个起始进度条
    pb <- utils::txtProgressBar(style=3)
    ##
    for (i in x) {
      fun1(i)
      ##第二个位置：实时显示进度
      utils::setTxtProgressBar(pb, which(i == x)/length(x))
      print(paste0(outdir, "/reclass/", spdata[i,1]))
    }
  }

  #组合每个物种的面积
  splist <- list.files(paste0(outdir, "/TBlabENM" , star_time), full.names = TRUE)
  v <- read.csv(splist[1])
  v <- v[0,]
  for (i in 1:length(splist)) {
    a <- read.csv(splist[i])
    v <- rbind(v, a)
  }
  if(is.list(threshold)){names(v) <- c("species", "pro", names(rcl))} else{
    names(v) <- c("species", "pro", "USA", "SA")
  }
  write.csv(v, paste0(outdir, "/suitable_area.csv"), row.names = FALSE, fileEncoding = "GB18030")
  unlink(paste0(outdir, "/TBlabENM", star_time), recursive = TRUE)
  end_time <- Sys.time()  ## 记录程序结束时间
  ## 第三个位置关闭进度条
  if(exists("pb")){close(pb)}
  print(end_time - star_time)

}
