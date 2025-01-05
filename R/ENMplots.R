# 单层连续栅格绘图
##下载美国行政边界数据
aus_map <- geodata::gadm(country = "AUS", level = 1, path = tempdir())
aus_map1 <- aus_map
aus_map <- aus_map1
aus_map <- geodata::gadm(country = "PRK", level = 1, path = tempdir())
plot(aus_map)
##栅格数据
tavg <- geodata::worldclim_global("tavg", res = 10, path = tempdir(), version = "2.1")
tavg_aus <- terra::crop(tavg, aus_map, mask = T)
plot(tavg_aus)
##1月
tavg_aus_Jan <- tavg_aus[[1]]
terra::plot(tavg_aus_Jan)
##绘图
ENMplots(spatraster = tavg_aus_Jan,
         spatvector = aus_map)
###更改陆地的颜色和水体的颜色
ENMplots(spatraster = tavg_aus_Jan_class,
         spatvector = aus_map,
         landcolor = "grey60",
         watercolor = "#BEE8FF")
###自定义颜色
ENMplots(spatraster = tavg_aus_Jan,
         spatvector = aus_map,
         rastercolors = c("yellow", "green"))
###更改投影坐标系
ENMplots(spatraster = tavg_aus_Jan,
         spatvector = aus_map,
         rastercolors = NULL,
         crs = "epsg:4544")
###仅绘制一个省份
zone <- crop(tavg_aus_Jan, ext(tavg_aus_Jan)+ c(-27,-5,0,-30))
ENMplots(spatraster = tavg_aus_Jan,
         spatvector = aus_map,
         zones = zone)
###保存图片

ENMplots(spatraster = tavg_aus_Jan,
         spatvector = aus_map,
         outdir = "your path",
         filename = "plot.jpg",
         width = 10,
         height = 10)

### 添加其他ggplot对象
### ENMplots()返回的是ggplot对象，可以与其他ggplot扩展结合，拥有更大的灵活性
p <- ENMplots(spatraster = tavg_aus_Jan,
              spatvector = aus_map)
####使用调色盘颜色
p + scale_fill_viridis_c()


#单层分类栅格绘图

##将澳大利亚1月的平均温重分类为高中低三个类别
rcl <- matrix(c(0,20,1,20,30,2,30,35,3), ncol = 3, byrow = TRUE)
tavg_aus_Jan_class <- terra::classify(tavg_aus_Jan, rcl = rcl)
plot(tavg_aus_Jan_class)
##绘图
p <- ENMplots(spatraster = tavg_aus_Jan_class,
              spatvector = aus_map,
              rastercolors = NULL,
              crs = "epsg:4324",
              landcolor = "white",
              #  zones = ra_crop,
              categories = T,
              maxcell = 5e+05,
              watercolor = "#BEE8FF",
              arrange = list(2, 3, "v"),
              outdir = "D:/Desktop",
              filename = "2.jpg",
              width = 10,
              height = 10)
p
##指定坐标系
p <- ENMplots(spatraster = tavg_aus_Jan_class,
              spatvector = aus_map,
              rastercolors = NULL,
              crs = "epsg:4324",
              landcolor = "white",
              #  zones = ra_crop,
              categories = T,
              maxcell = 5e+05,
              watercolor = "#BEE8FF",
              arrange = list(2, 3, "v"),
              outdir = "D:/Desktop",
              filename = "2.jpg",
              width = 10,
              height = 10)
# 多层连续栅格绘图



xy <- read.csv("F:/eblf/species/thin1km/allspecies_tidy/EAcheckoccthin/occthin1km/wz0001_thin1.csv")
ra_crop <- terra::crop(ra1, china_vect, mask = T)

ra1 <- terra::rast("F:/example/ra.tif")


ralist <- list.files("F:/example", pattern = "tif$", full.names = T)
ra <- terra::rast(ralist)
names(ra) <- c("sp_2030", "sp_2050", "sp_2070", "sp_2090")
 +
  ggplot2::geom_point(data = xy, aes(x = longitude, y = latitude))


#' Draw pictures for ENM
#' @description
#' 绘制适生区的图
#'
#' @param spatraster spatraster, 适宜性栅格
#' @param spatvector spatvector, 矢量地图
#' @param rastercolors 适宜性栅格的颜色
#' @param crs 投影
#' @param zones 区域
#' @param landcolor 大陆颜色
#' @param watercolor 水体颜色
#' @param arrange 图形排列方式
#' @param categories 是否为类别型栅格
#' @param maxcell 最大栅格数
#' @param outdir 图片保存路径
#' @param filename 图片名称
#' @param width 图片宽度
#' @param height 图片高度
#'
#' @return
#' ggplot对象和图
#' @export
#'
#' @examples

ENMplots <- function(spatraster,
                     spatvector,
                     rastercolors = c("green2", "#FF4500"),
                     crs = NULL,
                     zones = NULL,
                     landcolor = "white",
                     watercolor = "#BEE8FF",
                     arrange = list(nrow = 2, ncol = 2, dir = "v"),
                     maxcell = 5e+03,
                     categories = F,
                     outdir,
                     filename,
                     width,
                     height
                     ) {
  #参数检查
  if (is.null(names(arrange))) {
    names(arrange) <- c("nrow", "ncol", "dir")
  }
  if ("" %in% names(arrange)) {
    names(arrange)[which(names(arrange) %in% "")] <- c("ncol", "nrow", "dir")[which(names(arrange) %in% "")]
  }

  #分类/数值型栅格
  if (categories == TRUE) { #分类栅格
    if (terra::is.factor(ra2) == FALSE) {
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
        ggplot2::scale_fill_manual(colors = rastercolors, na.value = "transparent")
    }

  } else { #数值栅格
    if (terra::is.factor(ra2)) {
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

  #读取内置世界地图矢量数据
  word_vect <- terra::vect("C:/Users/why/Documents/ArcGIS/中国地图/世界地图/世界国家.shp")

  ##裁剪世界地图至研究区域
  if (is.null(landcolor)) {
    word_vect_crop <- NULL
  } else {
    word_vect_crop <- terra::crop(word_vect, terra::ext(spatvector)) |>
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
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.25,
          linetype = 'dashed',
          colour = "gray40"
        )
      ) +
      ggplot2::coord_sf(crs = crs)
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
          panel.grid.major = ggplot2::element_line(
            linewidth = 0.25,
            linetype = 'dashed',
            colour = "gray40"
          )
        ) +
        ggplot2::coord_sf(crs = crs)
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
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.25,
          linetype = 'dashed',
          colour = "gray40"
        )
      ) +
      ggplot2::coord_sf(crs = crs)
    if (missing(outdir) == FALSE) {
      ggplot2::ggsave(
        filename = paste0("plots_", i, ".jpg"),
        plot = p2,
        #device = "jpg",
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
