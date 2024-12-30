
#' @title Perform the MaxEnt model for a single species
#' @description 本函数为单个物种构建MaxEnt模型.
#' @details
#' 如果已经对MaxEet模型参数进行了优化(例如使用\code{\link[ENMeval]{ENMevaluate}})，
#' 用户可以使用参数args在一定程度上修改模型参数，来获得质量更高的模型。可以使用evlist
#' 参数来选择需要的环境变量集。当提供参数prodir时，则将模型投影至指定的地理空间，否则仅拟合模型。
#'
#' @param x 要模拟的物种文件的路径，包括文件名（.csv）。
#' @param evdir 环境变量路径,格式为.asc
#' @param evlist 用于建模的环境变量，是环境变量路径下所有变量的子集。用下标表示，
#' 数字代表要选用的变量，默认使用所有变量。
#' @param factors 字符型向量。指定哪些变量是分类变量。未指定时，则默认全为连续变量。
#' @param nbg 随机背景点数量，当指定参数mybgfile时忽略。
#' @param bgwidth 数值型(单位m), 以发生点为中心, 以bgwidth为半径创建一个缓冲区, 用来选择背景点.如果为NULL则在整个环境层范围内选择背景点.
#' @param mybgfile 自定义的背景点。包含两列（经度、纬度）
#' 未指定时，随机从环境层生成。
#' @param args 自定义的MaxEnt模型参数，详见\code{\link[TBlabENM]{maxent_args}}。
#' @param prodir 用列表储存的投影文件路径，包含了投影时期的名称和与之对应的环境变量路径
#' @param outdir 输出文件路径，若未指定则在当前目录下生成./TBlabENM/maxent文件夹保存输出结果
#' @param parallel 是否并行计算
#' @param ncpu 如果并行，使用的cpu数

#' @return 包含maxent模型结果的一系列文件
#' @export
#' @importFrom utils read.csv
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows
#' @importFrom predicts MaxEnt
#' @importFrom magrittr %>%
#' @examples
#'
#' maxent_single(
#'   x = system.file("extdata", "species/Phoebe sheareri.csv", package = "TBlabENM"),
#'   evdir = system.file("extdata", "envar/asc", package = "TBlabENM"),
#'   evlist = 1:7,
#'   factors = c("dr", "fao90"),
#'   nbg = 1000,
#'   bgwidth = 500000,
#'   args = maxent_args(),
#'   prodir = list(present = system.file("extdata", "envar/asc", package = "TBlabENM"),
#'               present1 = system.file("extdata", "envar/asc", package = "TBlabENM")
#'   ),
#'   parallel = T,
#'   ncpu =2)

maxent_single <- function(x,
                          evdir,
                          evlist = NULL,
                          factors = NULL,
                          mybgfile = NULL,
                          nbg = 10000,
                          bgwidth = NULL,
                          args = maxent_args(),
                          prodir = NULL,
                          outdir = NULL,
                          parallel = F,
                          ncpu = 2) {
  star_time <- Sys.time() ## 记录程序开始时间
  #检查prodir参数是否包含名字，否则报错
  if (is.null(prodir) == FALSE) {
    if (is.null(names(prodir))) {
      stop("The projection path does not have a valid name, please check parameter prodir")
    }
  }
  #获取物种名 对路径拆分并取倒数第一个字符串
  spname1 <- stringr::str_split_1(x, "/")[length(stringr::str_split_1(x, "/"))]
  sp_name <- stringr::str_split_1(spname1, ".csv$")[1]
  cat(paste("Buildding MaxEnt models for", sp_name, "\n"))

  #读取环境变量
  biolist <- list.files(evdir, pattern = ".asc$|.tif$", full.names = TRUE)
  if (is.null(evlist) == FALSE) {
    biolist <- biolist[evlist]
  }
  biostack <- terra::rast(biolist)

  #读取、提取存在点的环境值并转化为数据框，生成存在环境数据
  occ <- utils::read.csv(x) #读取物种坐标数据
  occ <- occ[c(2, 3)]
  names(occ) <- c("x", "y")

  if (is.null(bgwidth) == FALSE & is.null(mybgfile)) {
    #将发生点转为sf对象
    occs.sf <- sf::st_as_sf(occ, coords = c("x","y"), crs = terra::crs(biostack))
    #转为等面积投影
    eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)
    #创建缓冲区
    occs.buf <- sf::st_buffer(occs.sf, dist = bgwidth) |>
      sf::st_union() |>
      sf::st_sf() |>
      sf::st_transform(crs = terra::crs(biostack)) #再转为经纬度投影
    #将环境变量裁剪至缓冲区
    biostack <- terra::crop(biostack, occs.buf)
    biostack <- terra::mask(biostack, occs.buf)}

  occdata <- terra::extract(biostack, occ, ID = FALSE)

  #读取、提取背景点的环境值并转化为数据框，生成环境背景数据
  if (is.null(mybgfile)) {
    mybg0 <- terra::spatSample(biostack, nbg, na.rm = T, xy = T)
    mybg <- mybg0[1:2]
    mybgdata <- mybg0[-(1:2)]
  } else {
    mybg <- mybgfile
    names(mybg) <- c("x", "y")
    mybgdata <- terra::extract(biostack, mybg, ID = FALSE)
  }

  #将背景点赋值为0
  mybg$'p/b' <- rep(0, times = nrow(mybg))



  #将存在点赋值为1
  occ$'p/b' <- rep(1, times = nrow(occ))

  ####全变量模拟####
  #组合存在环境数据与环境背景数据构成环境变量数据框
  evdata <- dplyr::bind_rows(occdata, mybgdata)
  #指定分类变量
  #获取变量名
  bio_name <- c()
  for (i in seq_along(biolist)) {
    bioname1 <- stringr::str_split_1(biolist[i], "/")[length(stringr::str_split_1(biolist[i], "/"))]
    bio_name0 <- stringr::str_split_1(bioname1, ".asc|.tif")[1]
    bio_name <- c(bio_name, bio_name0)
  }
  #
  if (is.null(factors) == FALSE) {
    site <- bio_name[bio_name %in% factors]
    evdata <- evdata %>% dplyr::mutate_at(site, as.factor)
  }

  #组合存在点和背景点构成坐标点数据框（包含0或1的）（要与上面的顺序一致）
  xydata <- dplyr::bind_rows(occ, mybg)
  args2 <- args
  #当坐标点少于25个时，将replicates=设置为坐标点数。
  if (nrow(occ) < 25) {
    args2[1] <- paste0("replicates=", nrow(occdata))
    n_na <-  nrow(occ) - nrow(occdata)}

  # if (nrow(occ) >= 25) {
  #   args2[1] <- "replicates=10"
  #   n_na <-  nrow(occ) - nrow(occdata)}

  if (n_na > 0) {
    warning(paste0(n_na, "occurs points has a missing value on the environment grid, removed by default."))
    }

  #模型模拟, xydata$`p/b`
  if (is.null(outdir)) {
    outdir <- "./maxent/"
  } else{
    outdir <- paste0(outdir, "/maxent/")
  }

  if (is.null(prodir) == TRUE) {
    me1 <- predicts::MaxEnt(evdata, xydata$`p/b` , #环境变量、坐标数据和背景点
                         args2, paste0(outdir, sp_name))
  }#输出路径

  ##当提供prodir时，对不同时期循环
  if (is.null(prodir) == FALSE) {
    cat(paste("Perform the projection for", sp_name, "\n"))

    if (parallel == T) {
      # library(snowfall)

      # 开启集成
      snowfall::sfInit(parallel = TRUE, cpus = ncpu)

      #加载需要用到的变量或函数 因为下面函数fff中要用到prodir参数
      snowfall::sfExport("prodir")
      #构建函数fff
      fff <- function(y) {
        predicts::MaxEnt(
          evdata,
          xydata$`p/b` ,
          #环境变量、坐标数据和背景点
          append(args2, paste0("projectionlayers=", prodir[[y]])),
          #新建文件夹保存模拟结果
          path = paste0(outdir, sp_name, "/", names(prodir)[y])
        )

        #将结果文件的asc格式转为tif格式以节约内存
        df <- list.files(
          paste0(outdir, sp_name, "/", names(prodir)[y]),
          pattern = "asc$",full.names = TRUE) %>%
          as.data.frame()
        names(df) <- "file"
        df1 <- df %>%
          mutate(path = map_chr(
            .x = file,
            .f = function(x) {
              stringr::str_replace(x, pattern = ".asc", ".tif")
            }
          )) %>%
          dplyr::mutate(ra = purrr::map(
            .x = file,
            .f = function(x) {
              terra::rast(x)
            }
          )) %>%
          mutate(purrr::map2(
            .x = ra,
            .y = path,
            .f = function(x, y) {
              terra::writeRaster(x, y, overwrite = TRUE)
            }
          ))

      }
      snowfall::sfLapply(seq_along(prodir), fff)
      snowfall::sfStop()  # 关闭集群

      #删除asc格式
      gg <- file.remove(
        list.files(paste0(outdir, sp_name), full.names = TRUE) %>%
          list.files(., pattern = "asc$", full.names = TRUE)
      )

      end_time <- Sys.time()
      print(end_time - star_time)

    } else {
      for (b in seq_along(prodir)) {
        args3 <- append(args2, paste0("projectionlayers=", prodir[[b]]))
        cat(paste0("Project into "), prodir[[b]])
        #模型模拟, xydata$`p/b`
        me1 <- predicts::MaxEnt(evdata,
                             xydata$`p/b` ,
                             #环境变量、坐标数据和背景点
                             args3,
                             #新建文件夹保存模拟结果
                             path = paste0(outdir, sp_name, "/", names(prodir)[b]))

      }
      #将结果文件的asc格式转为tif格式以节约内存
      df <- list.files(
        paste0(outdir, sp_name, "/", names(prodir)[b]),
        pattern = "asc$", full.names = TRUE) %>%
        as.data.frame()
      names(df) <- "file"
      df1 <- df %>%
        mutate(path = map_chr(
          .x = file,
          .f = function(x) {
            stringr::str_replace(x, pattern = ".asc", ".tif")
          }
        )) %>%
        mutate(purrr::map2(
          .x = file,
          .y = path,
          .f = function(x, y) {
            terra::writeRaster(terra::rast(x), y, overwrite = TRUE)
          }
        ))

      #删除asc格式
      # gg <- file.remove(list.files(paste0(outdir, sp_name), full.names = TRUE) %>%
      #                     list.files(., pattern = "asc$", full.names = TRUE))

      end_time <- Sys.time()
      print(end_time - star_time)
    }
  #
  }
  #提取参数
  df <- data.frame(matrix(NA, 1, 10))
  names(df) <- c(
    "species",
    "number",
    "env",
    "fc",
    "rm",
    "replicates",
    "AUCtrain",
    "AUCtest",
    "MTSStrain",
    "MTSStest"
  )
  df$species <- sp_name

  df$env <- paste0(sort(names(occdata)), collapse = ",")
  fc <- args2[3:7][stringr::str_detect(args2[3:7], "TRUE")]
  fc1 <- c()
  for (i in fc) {
    fc2 <- stringr::str_split_1(i, "")[1]
    fc1 <- c(fc1, fc2)
  }
  df$fc <- toupper(paste0(fc1, collapse = ""))
  df$rm <- stringr::str_split_1(args2[2], "=")[2]
  df$replicates <- stringr::str_split_1(args2[1], "=")[2]
  if (is.null(prodir)) {
    data <- read.csv(paste0(outdir, sp_name, "/maxentResults.csv"))
  } else {
    data <- read.csv(paste0(outdir, sp_name, "/", names(prodir)[1], "/maxentResults.csv"))
  }

  df$AUCtest <- data$Test.AUC[length(data$Test.AUC)]
  df$AUCtrain <- data$Training.AUC[length(data$Training.AUC)]
  df$MTSStrain <- data$Maximum.training.sensitivity.plus.specificity.Logistic.threshold[length(data$Maximum.training.sensitivity.plus.specificity.Logistic.threshold)]
  df$MTSStest <- data$Maximum.test.sensitivity.plus.specificity.Logistic.threshold[length(data$Maximum.test.sensitivity.plus.specificity.Logistic.threshold)]
  df$number <- data$X.Training.samples[1] + data$X.Test.samples[1]
  cat(paste0("The results are saved ", outdir, sp_name, "\n"))
  write.csv(df, paste0(outdir, sp_name, "/parameters.csv"), row.names = FALSE)


  return(df)
}
