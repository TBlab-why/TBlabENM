
#' @title Occurrence Data Check
#' @description 检查物种发生数据是否存在缺失值以及是否在指定的海拔范围内。
#'
#' @param evdir 环境变量路径。格式为 .asc 或者 .tif。
#' @param elevation 海拔范围，单位m。字符型向量。
#' @param elrange 数值型向量。允许的海拔误差，单位m。默认500m。
#' @param spdir 要模拟的物种文件的路径，包括文件名（.csv）。
#' @param deal_NA 字符型,当发生点所在的栅格像元为缺失值时的处理方法.可选择"mv","rm","none".
#' 当为"mv"时,对物种坐标进行微调,移动至最近的非NA值像元上，当为"rm"时,直接删除该点,当为"none"时,
#' 不做任何处理.
#' @param width 数值型,指定缓冲区的宽度,当deal_NA选择"mv"时,在给定宽度的缓冲区内如果栅格都为缺失值则
#' 直接删除该点.单位通常为m.该值应大于像元的宽度，否则没有意义.
#' @param outdir 当deal_NA选择"mv"或"rm"时,处理后的数据的保存路径.
#' @param evselectev 字符型；用于检查数据的变量，第一个要为海拔变量。一种类型的变量选择一个即可。
#'
#' @return 返回可能存在问题的数据，建议对这些数据进行检查。
#' @export
#'
#' @importFrom stringr str_which
#' @importFrom dplyr mutate filter arrange
#' @importFrom purrr map map2 map2_dbl
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @author 吴海洋 和 田斌
#'
#' @examples
#' ENMspcheck(spdir = "./xx.csv",
#'            evdir = paste0(system.file(package="TBlabENM"), "/extdata/present"),
#'            #用于检查的变量，一种类型的变量选择一个即可
#'            evselectev = c("elev30", "bio01", "SU_CODE", "uvb1"),
#'            #海拔范围从植物志查询
#'            elevation = c(1000,1200),
#'            elrange = 500,
#'            deal_NA = "none")
ENMspcheck <- function(spdir, evdir, evselectev, elevation,
                       elrange = 500, deal_NA = "none", width,
                       outdir) {

  #获取物种名 对路径拆分并取倒数第一个字符串
    spname1 <- stringr::str_split_1(spdir, "/")[length(stringr::str_split_1(spdir, "/"))]
    sp_name <- stringr::str_split_1(spname1, ".csv$")[1]

  #读取栅格数据列表

    biolist <- list.files(evdir, pattern = "(tif|asc)$", full.names = TRUE)
    #选择要检查的变量，并stack
    evselectevdir <- c()
    for (i in seq_along(evselectev)) {
      evselectevdir1 <- biolist[stringr::str_which(biolist, paste0(evselectev[i],"\\."))]
      evselectevdir <- c(evselectevdir, evselectevdir1)
    }
    biostack_t <- terra::rast(evselectevdir)
  #读取物种坐标数据并提取出经纬度列
    occ1 <- utils::read.csv(spdir)
    occ <- occ1[c(2,3)]
    names(occ) <- c("x", "y")

  #提取每个坐标点对应的环境值并转为数据框
  occdata_t <- as.data.frame(terra::extract(biostack_t, occ, ID = FALSE))
  #检查缺失值所在行
  row_na <- which(rowSums(is.na(occdata_t)) != 0)
  #合并occ和occdata_t
  occdata <- cbind(occ1[1:3], occdata_t)
  #判断坐标点是否在海拔范围之内,这里的F为符合要求的数据

  occ_e <- ifelse(occdata_t[[1]] <= (max(elevation) + elrange)
                  & occdata_t[[1]] >= (min(elevation) - elrange), F, T)
  occ_e[which(is.na(occ_e))] <- F

  #sum(occ_e)==0和length(row_na)==0说明数据都符合要求，否则发出警告
  if (sum(occ_e) == 0 & length(row_na) == 0) {
    cat("The data set passes the altitude range and ENM envionment variables check.")
  } else {
    warning("The above coordinate points have some problem in elevation or ENM envionment variables, please check now.")
  }


#相应处理
  if (deal_NA == "mv" & length(row_na) > 0) {
    de_n <- c()
  #提取含有缺失值的点所在栅格的中心坐标
  point_err <- occ[row_na, ] #有缺失值的点
  cell_extract <- terra::extract(biostack_t, point_err[1:2], xy = T, ID = FALSE)
  for (i in seq_along(row_na)) {
  point <- occ[row_na[i], ]
  n <- as.numeric(rownames(point)) #有缺失值的点的行号
  cellsite_ext <- cell_extract[i,]
  cellsite_xy <- cellsite_ext[(length(cellsite_ext) - 1):length(cellsite_ext)]
  if (is.nan(cellsite_xy[1,1])) {de_n <- c(de_n, n)
  next}
  #将栅格裁剪至有问题的点周围，这里设置为1°应该可以满足大部分需求
  ex <- terra::ext(cellsite_xy$x-1, cellsite_xy$x+1, cellsite_xy$y-1, cellsite_xy$y+1)
  biostack_t1 <- terra::crop(biostack_t, ex)
  #以其中一个栅格为模版创建空栅格并对site赋值
  new_ra <- terra::rast(ncols = ncol(biostack_t1), nrows  = nrow(biostack_t1),
                     crs = terra::crs(biostack_t1), extent = terra::ext(biostack_t1))
  new_ra[terra::cellFromXY(new_ra, cellsite_xy)] <- 1
  #构建缓冲区
  new_ra_buffer <- terra::buffer(new_ra, width = width)

  #计算缓冲区的每个栅格与point之间的距离
  df <- which(new_ra_buffer[] == TRUE) %>%  #提取缓冲区内每个栅格的索引
    as.data.frame() %>%
    mutate(value = map(.x = ., .f = function(x){
      biostack_t1[x]})) %>%  #提取栅格值
    mutate(sum = map(.x = value, .f = function(x){
      sum(x)})) %>%     #求和，当值为NA时说明至少在一个图层中有缺失，下面移除这些缺失的栅格
    filter(., !grepl('NA', sum))

    if (nrow(df) > 0) {
      df <- df %>%
      mutate(xy = map(.x = ., .f = function(x){
        as.data.frame(terra::xyFromCell(new_ra_buffer, x))
        })) %>%  #获取栅格中心坐标
      mutate(point = map(.x = ., .f = function(x){
        point
      })) %>%
      mutate(distant = map2_dbl(.x = xy, .y = point, .f = function(x, y){
         stats::dist(rbind(x,y))
      })) %>%
      dplyr::arrange(., distant)}
  #将原数据替换为距离最近的栅格的中心坐标
    if (nrow(df) > 0) {
    occ1[n, 2:3] <- df$xy[[1]]
    } else {
      de_n <- c(de_n, n)
    }

  }
  if (length(de_n)>0) { occ1 <- occ1[-de_n, ]}
  if (exists("outdir")) {
    write.csv(occ1, paste0(outdir, "/", sp_name, "_check.csv"), row.names = FALSE)}
  }
  if (deal_NA == "rm" & length(row_na) > 0) {
    occ1 <- occ1[-n, ]
    if (exists("outdir")) {
      write.csv(occ1, paste0(outdir, "/", sp_name, "_check.csv"), row.names = FALSE)}
    }
#对有问题的3种情况进行相应的输出
if (sum(occ_e) != 0 & length(row_na) == 0) {return(occdata[occ_e,])}
if (sum(occ_e) == 0 & length(row_na) != 0) {return(occdata[row_na,])}
if (sum(occ_e) != 0 & length(row_na) != 0) {return(rbind(occdata[occ_e,], occdata[row_na,]))

}

}


