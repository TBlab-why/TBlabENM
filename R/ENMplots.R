
#' Draw pictures for ENM results
#' @description
#' 主要用于绘制适生区的图, 例如maxent生成的连续型适宜性栅格, 重分类后的分类适宜性栅格等. 预设一套背景模版, 可以在该模版上替换不同的栅格和矢量绘制常用的地图. 也可以结合ggplot2及其扩展绘制更个性化的地图.
#'
#' @param spatraster spatraster对象, 通常为物种分布的适宜性栅格.
#' @param spatvector spatvector对象, 通常为研究区的矢量地图.
#' @param rastercolors 字符型, 创建渐变色(continuous scale)或分类色(discrete scale)用于设置适宜性栅格的颜色.
#' @param crs The coordinate reference system (CRS) into which all data should be projected before plotting. If not specified, will use the CRS defined in the first sf layer of the plot.
#' @param zones NULL or SpatRaster with the same geometry identifying zones in spatraster.
#' @param landcolor 除研究区外的其他陆地颜色. 用来做背景色.
#' @param watercolor 缺失值的颜色, 通常代表水体颜色.
#' @param arrange 包含3个元素的列表, 分别为nrow, ncol, dir. 当spatraster为多层时, 将使用分面绘制每层栅格, nrow和ncol 设置一张图片上分图的行列数, dir设置分图的排列方向, dir = "v"则按水平方向排列, dir = "h"则按垂直方向排列. 如果层数大于nrow* ncol, 则会生成多张图片, 直到绘制完所有层.
#' @param categories 逻辑值. 将spatraster视为类别型还是连续型栅格.
#' @param maxcell positive integer. Maximum number of cells to use for the plot.
#' @param expand 数值向量, 大于等于1, 长度为1或4. 按比例扩展绘图范围. 当长度为4时分别对应左右下上.
#' @param outdir 图片保存路径.
#' @param filename 图片名称.
#' @param width 图片宽度(cm).
#' @param height 图片高度(cm).
#'
#' @return
#' ggplot2对象
#' @export
#'
#' @examples
#'
#' # 单层栅格绘图
#' ##下载澳大利亚行政边界数据
#' aus_map <- geodata::gadm(country = "AUS", level = 1, path = tempdir())
#' ##下载全球12个月的平均温度栅格数据
#' tavg <- geodata::worldclim_global("tavg", res = 10, path = tempdir(), version = "2.1")
#' ##裁剪至澳大利亚范围
#' tavg_aus <- terra::crop(tavg, aus_map, mask = T)
#' ##提取1月的数据
#' tavg_aus_Jan <- tavg_aus[[1]]
#' ##绘图
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map)
#' ##更改陆地的颜色和水体的颜色
#' ###如果不需要陆地和水体颜色可以将颜色设置为"transparent", 透明色.
#' ENMplots(spatraster = tavg_aus_Jan_class,
#'          spatvector = aus_map,
#'          landcolor = "grey60",
#'          watercolor = "#BEE8FF")
#' ##自定义温度栅格的颜色
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map,
#'          rastercolors = c("yellow", "green"))
#' ##移除指北针和比例尺
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map,
#'          annotation_north_arrow = F,
#'          annotation_scale = F)
#' ##更改投影坐标系
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map,
#'          crs = "epsg:4544")
#' ##仅绘制区域地图
#' ###区域地图需要提前准备一个与spatraster具有相同地理参考的栅格数据. 如果只想展示研究区的一部分结果可能很有用, 例如某个省.
#' zone <- terra::crop(tavg_aus_Jan, terra::ext(tavg_aus_Jan) + c(-27,-5,0,-30))
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map,
#'          zones = zone)
#' ##保存图片
#' ###可以在函数中使用参数保存图片, 请在filename中指定图片格式. 支持的格式同ggplot2::ggsave().
#' ENMplots(spatraster = tavg_aus_Jan,
#'          spatvector = aus_map,
#'          outdir = "your path",
#'          filename = "plot.jpg",
#'          width = 10,
#'          height = 10)
#'
#' ##添加其他ggplot对象
#' ### ENMplots()返回的是ggplot对象，可以与其他ggplot扩展结合，拥有更大的灵活性
#' ###使用tidyverse包的调色盘并更换主题
#' ENMplots(spatraster = tavg_aus_Jan, spatvector = aus_map) +
#'   scale_fill_whitebox_c(
#'     palette = "high_relief",
#'     labels = scales::label_number(suffix = "º"),
#'     n.breaks = 12,
#'     guide = guide_legend(reverse = TRUE)
#'   )  +
#'   theme_bw()
#' ###如果添加了额外的ggplot对象，则需要使用ggplot2::ggsave()在外部保存图片.
#'
#' ##此外还可以绘制简单的样本分布点
#' ###随机选择100个点当做样本分布点
#' point <- terra::spatSample(tavg_aus_Jan, 100, xy = T, na.rm = T)
#' ENMplots(spatraster = NULL,  #如果不需要温度栅格就设置为NULL
#'          spatvector = aus_map) +
#'   ggplot2::geom_point(data = point, aes(x = x, y = y))

#' ##单层分类栅格绘图
#' ##将澳大利亚1月的平均温重分类为高中低(3, 2, 1)三个类别
#' rcl <- matrix(c(0,20,1,20,30,2,30,35,3), ncol = 3, byrow = TRUE)
#' tavg_aus_Jan_class <- terra::classify(tavg_aus_Jan, rcl = rcl)
#' ##绘图
#' ###如果想要对分类栅格绘图，需要设置categories=T. 同样, 如果想要对连续栅格绘图需要设置categories=F.
#' ENMplots(spatraster = tavg_aus_Jan_class,
#'          spatvector = aus_map,
#'          categories = T)
#' # 多层栅格绘图
#' ##多层栅格绘图主要用于批量绘图并对图片排版
#' ###首先需要创建含有多层的spatraster, tavg_aus包含了12个层(12个月的平均温度).
#' ENMplots(spatraster = tavg_aus,
#'          spatvector = aus_map,
#'          categories = F,
#'          arrange = list(2, 3, "v"), #使用arrange参数设置每张图片上各个分图的数量及排列方式.
#'          outdir = "your path",
#'          device = "pdf",
#'          width = 20,
#'          height = 20)
#' ####注意: 当多层绘图时不用对生成的图片命名, 内部会自动对每张图进行编号. 需要通过device参数设置图片格式.

ENMplots <- function(spatraster = NULL,
                     spatvector = NULL,
                     zones = NULL,
                     expand = 1.05,
                     rastercolors = NULL,
                     landcolor = "white",
                     watercolor = "#BEE8FF",
                     maxcell = 5e+05,
                     categories = F,
                     crs = NULL,
                     annotation_north_arrow = T,
                     annotation_scale = T,
                     arrange = list(nrow = 2, ncol = 2, dir = "v"),
                     width,
                     height,
                     device = "jpg",
                     filename,
                     outdir
                     ) {
  #参数检查
  if (sum(length(expand) == c(1, 4)) == 0) {stop("The length of expand must be 1 or 4.")}
  for (i in length(expand)) {
    if (expand[i] < 1) {stop("Expand must be greater than or equal to 1.")}
  }

  if (is.null(names(arrange))) {
    names(arrange) <- c("nrow", "ncol", "dir")
  }
  if ("" %in% names(arrange)) {
    names(arrange)[which(names(arrange) %in% "")] <- c("nrow", "ncol", "dir")[which(names(arrange) %in% "")]
  }
  if (missing(width)) {width = NA}
  if (missing(height)) {height = NA}
  #分类/数值型栅格
  if (categories == TRUE) { #分类栅格
    if (terra::is.factor(spatraster) == FALSE) {
      spatraster <- terra::as.factor(spatraster)
    }
    ##设置颜色
    if (is.null(rastercolors)) {
      p_color <- tidyterra::scale_fill_whitebox_d(
        palette = "muted",
        guide = ggplot2::guide_legend(reverse = TRUE)
      )
    } else {
      p_color <- #自定义离散颜色
        ggplot2::scale_fill_manual(values = rastercolors, na.value = "transparent")
    }

  } else { #数值栅格
    if (terra::is.factor(spatraster)) {
      spatraster <- terra::as.numeric(spatraster)
    }
    #设置颜色
    if (is.null(rastercolors)) {
      p_color <- tidyterra::scale_fill_whitebox_c(
        palette = "muted",
        n.breaks = 10,
        guide = ggplot2::guide_legend(reverse = TRUE)
      )
    } else {
      p_color <- #自定义离散颜色
        ggplot2::scale_fill_gradientn(colors = rastercolors, na.value = "transparent")
    }
  }

  #添加指北针和比例尺
  if (annotation_north_arrow == FALSE) {zbz <- NULL} else {
    zbz <- ggspatial::annotation_north_arrow(
      location = "tl",  # 指北针位置
      style = ggspatial::north_arrow_fancy_orienteering,  # 更改为 minimal 样式
      which_north = "true",  # 使用真实北方向
      pad_x = ggplot2::unit(0.5, "cm"),  # 调整指北针的水平边距
      pad_y = ggplot2::unit(0.5, "cm")   # 调整指北针的垂直边距
    )
  }
  if (annotation_scale == FALSE) {blc <- NULL} else {
  blc <- ggspatial::annotation_scale(
    location = "bl",  # 比例尺位置
    style = "ticks"
    )
  }

  #转换投影到指定的crs
  # if (is.null(crs) == FALSE) {
  #   if (is.null(spatraster) == FALSE) {
  #     spatraster <- terra::project(spatraster, crs)
  #   }
  # }
  #统一spatvector和spatraster投影
  if (is.null(spatvector) == FALSE) {
    spatvector <- terra::project(spatvector, terra::crs(spatraster))}
  #读取内置世界地图矢量数据
  if (is.null(landcolor) == FALSE) {
    #word_vect <- terra::vect("C:/Users/why/Documents/ArcGIS/中国地图/世界地图/世界国家.shp")
    word_vect <- terra::vect(paste0(system.file(package = "TBlabENM"), "/extdata/land/continent.shp")) %>%
      terra::project(terra::crs(spatraster))
  }

  ##裁剪世界地图至研究区域
  if (is.null(landcolor)) {
    word_vect_crop <- NULL
  } else {
    word_vect_crop <- terra::crop(word_vect, terra::ext(spatvector)*expand) |>
      terra::aggregate()
  }
  #裁剪地图至指定区域
  if (is.null(zones) == FALSE) {
    spatraster <- terra::crop(spatraster, zones, mask = TRUE)
  }


  #没有栅格的情况
  if (is.null(spatraster)) {
    p1 <- ggplot2::ggplot() +
          tidyterra::geom_spatvector(data = word_vect_crop, fill = landcolor) +
          tidyterra::geom_spatvector(data = spatvector, fill = NA) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::scale_x_continuous(expand = c(0, 0))

    p2 <- p1 + p_color +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(color = 'black', fill = watercolor),
        panel.border = ggplot2::element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.25,
          linetype = 'dashed',
          colour = "gray40"
        )
      ) +
      ggplot2::coord_sf(crs = crs) + zbz + blc
    if (missing(outdir) == FALSE) {
      ggplot2::ggsave(
        filename = filename,
        plot = p2,
        device = NULL,
        path = outdir,
        scale = 1,
        width = width,
        height = height,
        units = c("cm"),
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
        create.dir = TRUE
      )}
    return(p2)
  }  else {
    #仅有一个栅格的情况
    if (terra::nlyr(spatraster) == 1) {
        p1 <- ggplot2::ggplot() +
          tidyterra::geom_spatvector(data = word_vect_crop, fill = landcolor) +
          tidyterra::geom_spatraster(data = spatraster, maxcell = maxcell) +
          tidyterra::geom_spatvector(data = spatvector, fill = NA) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::scale_x_continuous(expand = c(0, 0))

      p2 <- p1 + p_color +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(color = 'black', fill = watercolor),
          panel.border = ggplot2::element_rect(color = "black", size = 1, fill = NA),
          panel.grid.major = ggplot2::element_line(
            linewidth = 0.25,
            linetype = 'dashed',
            colour = "gray40"
          )
        ) +
        ggplot2::coord_sf(crs = crs) + zbz + blc
      if (missing(outdir) == FALSE) {
        ggplot2::ggsave(
          filename = filename,
          plot = p2,
          device = NULL,
          path = outdir,
          scale = 1,
          width = width,
          height = height,
          units = c("cm"),
          dpi = 300,
          limitsize = TRUE,
          bg = NULL,
          create.dir = TRUE
        )}
      return(p2)
    }
  }
  #有多个栅格的情况
  page <- arrange$ncol*arrange$nrow
  if (terra::nlyr(spatraster) > 1) {
    #对每页的图进行循环 i = 2
    for (i in 1:ceiling(terra::nlyr(spatraster)/page)) {
      spatraster_i <- spatraster[[(i * page - page + 1):(i*page)]]
      p1 <- ggplot2::ggplot() +
        tidyterra::geom_spatvector(data = word_vect_crop, fill = landcolor) +
        tidyterra::geom_spatraster(data = spatraster_i, maxcell = maxcell) +
        tidyterra::geom_spatvector(data = spatvector, fill = NA) +
        ggplot2::facet_wrap(~ lyr, ncol = arrange$ncol, nrow = arrange$nrow,
                            dir = arrange$dir) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_x_continuous(expand = c(0, 0))

    p2 <- p1 + p_color +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(color = 'black', fill = watercolor),
        panel.border = ggplot2::element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.25,
          linetype = 'dashed',
          colour = "gray40"
        )
      ) + zbz + blc

      ggplot2::coord_sf(crs = crs)
    if (missing(outdir) == FALSE) {
      ggplot2::ggsave(
        filename = paste0("plots_", i, ".", device),
        plot = p2,
        path = outdir,
        scale = 1,
        width = width,
        height = height,
        units = c("cm"),
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
        create.dir = TRUE
      )}
    }
    return(p2)
  }

}
