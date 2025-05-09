#' @title Calculate the area of the reclassified raster
#' @description 计算重分类栅格面积.
#' @param reclassdir 重分类栅格所在文件夹路径, 不包括物种名.
#' @param crs 给缺乏投影的栅格定义投影. 例如"epsg:4326". 栅格的投影应与环境变量保持一致, 错误的投影会得到错误的结果.
#' @param prefix 字符型. 要重分类的栅格文件名的前缀, 用于筛选栅格.
#' @param suffix 字符型. 要重分类的栅格文件名的后缀, 用于筛选栅格.
#' @param overwrite logical. If TRUE, filename is overwritten.
#' @param outdir 输出文件夹路径.
#'
#' @return
#'    列表. reclass_area_path为计算的面积文件存储路径; area为重分类后的面积(单位KM²).
#' @export
#'
#' @examples
#' ENMreclass_area(
#'   reclassdir = "/ifs1/User/wuhaiyang/sdm/reclass25部分",
#'   crs = "epsg:4326",
#'   prefix = "wz0001",
#'   suffix = "tif",
#'   outdir = "/ifs1/User/wuhaiyang/sdm"
#' )
#'
ENMreclass_area <- function(reclassdir, crs,
                            prefix = NULL, suffix,
                            overwrite = FALSE, outdir = NULL) {
  star_time <- Sys.time() ## 记录程序开始时间

  # 读取模拟结果列表(只读文件夹)（以物种为单位）
  if (file.exists(reclassdir) == FALSE) {
    stop("reclassdir not find.")
  }
  if (is.null(outdir)) {
    outdir <- "."
  }
  ################################################# 函数
  # 读取栅格数据
  ra_df <- list.files(reclassdir,
    full.names = TRUE, recursive = TRUE,
    pattern = paste0("^", prefix, ".*", suffix, "$")
  ) %>%
    as.data.frame() %>%
    # 获取文件名
    dplyr::mutate(name = purrr::map(.x = ., .f = function(x) {
      str_split_1(basename(x), ".tif")[1]
    }))

  ###############################

  radf <- ra_df %>%
    dplyr::mutate(ra = map(.x = ., .f = function(x) {
      terra::rast(x)
    })) %>%
    dplyr::mutate(area = purrr::map(.x = ra, .f = function(x) {
      if (terra::crs(x) == "") {
        terra::crs(x) <- crs
      }
      terra::expanse(x, unit = "km", byValue = TRUE, wide = TRUE) # 计算唯一值对应的面积
    })) %>% # 计算面积
    dplyr::mutate(area1 = purrr::map2(.x = area, .y = name, .f = function(x, y) {
      x$name <- y
      x
    }))

  # 新建数据框保存单个物种的结果
  data <- data.frame(matrix(NA, nrow = 0, ncol = ncol(radf[, 5][[1]])))
  names(data) <- names(radf[, 5][[1]])

  for (i in 1:length(radf[, 5])) {
    b1 <- radf[, 5][[i]]
    data <- merge(data, b1, all = TRUE)
  }

  data <- data[-1] %>%
    dplyr::relocate(name)
  utils::write.csv(data, paste0(outdir, "/", "reclass_area.csv"),
    row.names = FALSE
  )


  end_time <- Sys.time() ## 记录程序结束时间
  print(end_time - star_time)
  # 原始数据, 重分类数据, 面积表, 存储路径
  s3 <- list(reclass_area_path = paste0(outdir, "/", "reclass_area.csv"), area = data)
  return(s3)
}
