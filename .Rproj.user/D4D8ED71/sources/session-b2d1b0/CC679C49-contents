
#封装函数
#' @title Map Suitable Areas
#' @description 根据阈值绘制适生区图，阈值以下代表非适生区。可以采用固定阈值，也可以自由指定。
#' @param x 数值向量。表示要进行转换的物种的序号。这个序号来自ENMspdata()生成的标准文件。
#' @param spdata 由ENMspdata()生成的标准文件。也可手动创建或修改，但是不建议。
#' 当maxent结果文件出现变动后，强烈建议重新运行ENMspdata()。
#' @param proname 想要绘图的时期名称。
#' @param threshold 划分适生区的阈值。可以为"T_MTSSte", "T_MTSStr"（依据spdata）
#' 或者为一个0到1的数值型，比如0.2。
#' @param wordmap 世界地图
#' @param regionmap 区域地图，比如研究区域的地图
#' @param color 指定不同适生区等级的颜色。5个等级，需要指定5种颜色，分别代表适生等级从低到高。
#' @param legend 逻辑值，是否需要图例，默认添加图例
#' @param resultdir maxent结果文件路径。例如paste0(getwd(), "/TBlabENM/result")
#' @param path 结果保存路径
#' @param format 结果保存格式,(e.g. png), or one of "eps", "ps", "tex" (pictex),
#'  "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#'
#' @return 按照阈值进行划分的适生区图
#' @export
#' @importFrom sf st_read
#' @importFrom utils txtProgressBar
#' @importFrom methods as
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom utils setTxtProgressBar
#' @importFrom grid rectGrob
#' @importFrom grid gpar
#' @importFrom grid unit
#' @importFrom ggplot2 annotate
#' @importFrom ggspatial annotation_scale
#' @import tidyverse
#'
#' @examples
#' #读取.spdata.csv文件
#' spdata <- read.csv(paste0(system.file(package="TBlabENM"),"/extdata/.spdata.csv"))
#' #绘图
#' ENMPic(x = 1:10,
#'       spdata = spdata,
#'       proname = "present",
#'       threshold = 0.2,
#'       resultdir = "./TBlabENM/result",
#'       wordmap = "D:/Desktop/ENMdata/worldmap/带投影坐标系/world_pro.shp",
#'       regionmap = "D:/Desktop/ENMdata/Chinamap/2022年省界/sheng2022.shp",
#'       #color = c("", "", "", "", ""),
#'       legend = TRUE,
#'       format = "jpeg",
#'       path= "./TBlabENM/sppicture")
ENMPic <- function(parameters, x, key, threshold = "T_MTSSte",
                   resultdir = "./TBlabENM/result",
                   wordmap = NULL,
                   regionmap = NULL,
                   color = NULL,
                   legend = TRUE,
                   format = "jpeg",
                   path) {
  ####环境准备
 # sysfonts::font_add('TNR','C:/Windows/Fonts/times.ttf')
 # showtext::showtext_auto()

  ####加载地图
  if (is.null(wordmap)==FALSE) {
    wordmap1 <- sf::st_read(wordmap)
  }
  if (is.null(regionmap)==FALSE) {
    regionmap1 <- sf::st_read(regionmap)
  }


  ###读取模拟结果文件夹（以物种为单位）
  if (file.exists(resultdir)==FALSE){
    stop("Picdir not find.")}

  ## 第一个位置：新建一个起始进度条
  pb <- utils::txtProgressBar(style=3)
  star_time <- Sys.time() ## 记录程序开始时间

  #设置默认颜色
  if (is.null(color)) {
    color <- c("1" = rgb(242,242,242,maxColorValue = 255),
               "2" = rgb(61,230,0,maxColorValue = 255),
               "3" = rgb(16,191,0,maxColorValue = 255),
               "4" = rgb(0,152,5,maxColorValue = 255),
               "5" = rgb(0,128,0,maxColorValue = 255))}
  for (i in x) {
    for (j in seq_along(proname)) {
      #获取物种名
      spname <- spdata[i,1]
      #读取栅格数据

      name <- list.files(paste0(resultdir, "/", spname, "/", proname[j]),
                         pattern = "_avg.asc$", full.names = TRUE)
      if (length(name)==0) {
        name <- list.files(paste0(resultdir, "/", spname, "/", proname[j]),
                           pattern = "_avg.tif$", full.names = TRUE)
      }
      if (length(name)==0) {
        stop("Raster does not exist")
      }
      pro <- raster::raster(name)

      #定义投影
      #raster::crs(pro)="+proj=aea +lat_0=0 +lon_0=105 +lat_1=25 +lat_2=47 +x_0=0 +y_0=0 +ellps=krass +units=m +no_defs"
      #很关键的一步：将栅格数据转成可供ggplot2绘制的数据格式
      pro_spdf <- methods::as(pro, "SpatialPixelsDataFrame")
      pro_df <- as.data.frame(pro_spdf)
      #赋值列名
      colnames(pro_df) <- c("value", "x", "y")


      if(is.character(threshold)) {
        col <- which(names(spdata) == threshold)

        pro_df$class <- as.factor(
          ifelse(pro_df$value <= spdata[i,col], 1,
                 ifelse(pro_df$value <= 0.4, 2,
                        ifelse(pro_df$value <= 0.6, 3,
                               ifelse(pro_df$value <= 0.8, 4, 5)))))}

      if(is.numeric(threshold)) {

        pro_df$class <- as.factor(
          ifelse(pro_df$value <= threshold, 1,
                 ifelse(pro_df$value <= 0.4, 2,
                        ifelse(pro_df$value <= 0.6, 3,
                               ifelse(pro_df$value <= 0.8, 4, 5)))))}
      #计算xy取值范围和xy值域比例
      xrange <- (max(pro_df$x)-min(pro_df$x))
      yrange <- (max(pro_df$y)-min(pro_df$y))
      ratio <- (max(pro_df$x)-min(pro_df$x))/(max(pro_df$y)-min(pro_df$y))


      ###修改ggplot2函数，使图例键之间的距离调大
      # function to increase vertical spacing between legend keys
      # @clauswilke
      draw_key_polygon3 <- function (data, params, size) {
        grid::rectGrob(width = grid::unit(1.4, "npc"),
                       height = grid::unit(0.8, "npc"),
                       gp = grid::gpar(col = NA, fill = alpha(data$fill %||%
                                                                data$colour %||% "grey20", data$alpha), lty = data$linetype %||%
                                         1))}

      # register new key drawing function,
      # the effect is global & persistent throughout the R session
      GeomRaster$draw_key = draw_key_polygon3


      ####绘制分类适生区图####
      #设置图形基本参数
      #设置主题
      bgtm <- theme(#清空默认的背景网格设置
        # panel.grid.major = element_blank(),
        #设置背景色
        panel.background = element_rect(fill = "#d6ecf0"),
        #设置绘图区域页边距
        plot.margin=unit(c(0.14,0.10,0.45,0.42),'cm'),
        #设置边框
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        #设置坐标轴标签（删除标签及其所留下的空白）
        axis.title = element_blank(),
        #坐标刻度线位置及长度（设置内部显示，负值代表内部）
        axis.ticks.length.x = unit(-0.20,"cm"),
        axis.ticks.length.y = unit(-0.20,"cm"),
        #设置刻度文本距离绘图边界的位置（设置为内部显示）
        axis.text.x = element_text(margin = unit(c(t=-0.75, r=0, b=0, l=0), 'cm')),
        axis.text.y = element_text(margin = unit(c(0, -1.4, 0, 0), 'cm')),
        #设置坐标轴刻度标签文本字体格式
        axis.text = element_text(family = "TNR", size = 45,
                                 face = 'bold', colour = "black") )
      # p_d1为没有图例的参数
      if (is.null(regionmap)&is.null(wordmap)&legend==FALSE) {
        p_d1 <- ggplot2::ggplot() +

          #自定义离散颜色
          scale_fill_manual(
            # 这里使用colorRampPalette()创建5等分颜色，使颜色接近渐变
            # values = colorRampPalette(c("Snow1",	"ForestGreen"))(4),
            #这里自己指定颜色，"1"设置无颜色
            values = color,
            guide = "none" ) +

          #设置主题
          theme(#清空默认的背景网格设置
            panel.grid.major = element_blank(),
            #设置背景色
            panel.background = element_rect(fill = "#d6ecf0"),
            #设置绘图区域页边距
            plot.margin=unit(c(0.15,0.15,0.15,0.15),'cm'),
            #设置边框
            panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
            #设置坐标轴标签（删除标签及其所留下的空白）
            axis.title = element_blank(),
            #坐标刻度线位置及长度（设置内部显示，负值代表内部）
            axis.ticks.length.x = unit(-0.10,"cm"),
            axis.ticks.length.y = unit(-0.10,"cm"),
            #不显示坐标轴刻度文本
            axis.text = element_blank()
          )
      } else {

        if (legend==FALSE) {


          p_d1 <- ggplot2::ggplot() +

            #自定义离散颜色
            scale_fill_manual(
              # 这里使用colorRampPalette()创建5等分颜色，使颜色接近渐变
              # values = colorRampPalette(c("Snow1",	"ForestGreen"))(4),
              #这里自己指定颜色，"1"设置无颜色
              values = color,
              guide = "none" ) +

            #设置主题
            bgtm
        } }

      ##添加图例
      lg <- theme( #图例位置
        legend.position = c(0.9,0.21),
        # legend.direction = "vertical",
        #  legend.margin = margin(t = 0, r = 2, b = 5, l = 0, unit = "cm"),
        #图例文本
        legend.text = element_text(family = "TNR", size = 40, margin = margin(l= -15)),
        #图例标题
        legend.title = element_text(family = "TNR", size = 50, margin = margin(b=-14, l=-8)),

        #图例背景
        legend.background = element_blank(),
        #图例键
        legend.key = element_rect(color = NA, fill = NA),
        #图例键高度和宽度
        #legend.key.height = unit(0.8, "cm"),
        # legend.key.width = unit(10.5, "cm"),
        legend.key.size = unit(0.8, "cm"))
      # p_d1为有图例的参数
      if (is.null(regionmap)&is.null(wordmap)&legend==TRUE) {
        p_d1 <- ggplot2::ggplot() +

          #自定义离散颜色
          scale_fill_manual(
            labels = c("0-MTSS","MTSS-0.4","0.4-0.6","0.6-0.8","0.8-1"),
            #反转图例顺序
            guide = guide_legend(reverse = TRUE),
            # 这里使用colorRampPalette()创建5等分颜色，使颜色接近渐变
            # values = colorRampPalette(c("Snow1",	"ForestGreen"))(4),
            #这里自己指定颜色，"1"设置无颜色
            values = color) +

          #设置主题
          theme(#清空默认的背景网格设置
            panel.grid.major = element_blank(),
            #设置背景色
            panel.background = element_rect(fill = "#d6ecf0"),
            #设置绘图区域页边距
            plot.margin=unit(c(0.15,0.15,0.15,0.15),'cm'),
            #设置边框
            panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
            #设置坐标轴标签（删除标签及其所留下的空白）
            axis.title = element_blank(),
            #坐标刻度线位置及长度（设置内部显示，负值代表内部）
            axis.ticks.length.x = unit(-0.10,"cm"),
            axis.ticks.length.y = unit(-0.10,"cm"),
            #不显示坐标轴刻度文本
            axis.text = element_blank() ) + lg


      } else {

        if (legend==TRUE) {

          p_d1 <- ggplot2::ggplot() +

            #自定义离散颜色
            scale_fill_manual(name="Value",
                              labels = c("0-MTSS","MTSS-0.4","0.4-0.6","0.6-0.8","0.8-1"),
                              #反转图例顺序
                              guide = guide_legend(reverse = TRUE),
                              # 这里使用colorRampPalette()创建5等分颜色，使颜色接近渐变
                              # values = colorRampPalette(c("Snow1",	"ForestGreen"))(4),
                              #这里自己指定颜色，"1"设置无颜色
                              values = color) +

            #设置主题
            bgtm + lg

        } }

      #绘制栅格（分布图）
      graster <- geom_raster(data = pro_df, aes(x = x, y = y, fill = class), alpha=0.8)
      #绘制区域地图
      if (is.null(regionmap)==FALSE) {
        gregionmap1 <- geom_sf(data = regionmap1, aes( geometry = `geometry`), size = 0.2,
                               color = rgb(110,110,110,maxColorValue = 255), fill = NA)
        gregionmap2 <- geom_sf(data = regionmap1, aes( geometry = `geometry`), size = 0.2,
                               color = rgb(110,110,110,maxColorValue = 255), fill = NA,
                               #中国地图线形设为虚线
                               linetype="dotted")
      }

      #设置为虚线的区域地图
      #绘制世界地图
      if (is.null(wordmap)==FALSE) {
        gwordmap1 <- geom_sf(data = wordmap1, aes(geometry = `geometry`), size = 0.2,
                             color = rgb(79, 79, 79, maxColorValue = 255), fill = NA)}
      #放大地图，只绘制想要的区域
      cf <- coord_sf(xlim = c(min(pro_df$x)+xrange*0.042, max(pro_df$x)-xrange*0.042),
                     ylim = c(min(pro_df$y)+yrange*0.042, max(pro_df$y)-yrange*0.042))
      #添加文字


      at <- ggplot2::annotate("text", x = (min(pro_df$x)+xrange*0.27), y = (max(pro_df$y)-yrange*0.12), label = proname[j],
                              size = 30, family = "TNR", fontface = "bold", color="black")

      #比例尺
      as <- ggspatial::annotation_scale(width_hint = 0.4,
                                        text_cex = 3,
                                        text_family = "TNR",
                                        pad_x = unit(1.5, "cm"),
                                        pad_y = unit(1.1, "cm"))
      if (is.null(wordmap)){
        if (is.null(regionmap)) {
          #当wordmap和regionmap都没有时只绘制栅格图
          if (legend == TRUE) {
            p_d <-
              p_d1 +  #底图设置
              graster + #绘制栅格（分布图）
              cf + #放大地图，只绘制想要的区域
              at + #添加文本
              as #比例尺
          } else {p_d <- p_d1 +  graster + cf}

        } else {
          ###绘制栅格图和区域地图
          if (legend == TRUE) {
            p_d <- p_d1 + graster + gregionmap1 + cf + at +as
          } else {p_d <- p_d1 + graster + gregionmap1 + cf}


        } } else {
          if (is.null(regionmap)) {
            ###绘制栅格图和世界地图
            if (legend == TRUE) {
              p_d <- p_d1 + graster + gwordmap1 + cf + at +as
            } else {p_d <- p_d1 + graster + gwordmap1 + cf}


          } else {
            if (legend == TRUE) {
              ###绘制栅格图和区域地图和世界地图
              p_d <- p_d1 + graster + gregionmap2 + gwordmap1 +cf + at + as
            } else {
              p_d <- p_d1 + graster + gregionmap2 + gwordmap1 +cf}
          }}


      ##保存
      if (format=="ps") {ggplot2::ggsave(plot = p_d,
                                         filename = paste0(paste(spdata[i,1], proname[j], sep = "-"), ".psd" ),
                                         path = path,
                                         device = format, scale=1, dpi = "print" ,units = "cm",
                                         width = 18,
                                         height = 18/ratio )} else {
                                           ggplot2::ggsave(plot = p_d, filename = paste0(
                                             paste(spdata[i,1], proname[j], sep = "-"),
                                             paste0(".", format) ),
                                             path = path,
                                             device = format, scale=1, dpi = "print" ,units = "cm",
                                             width = 18,
                                             height = 18/ratio ) }
      ##第二个位置：实时显示进度
      utils::setTxtProgressBar(pb,
                               ((which(x==i)-1)*length(proname)+j)/(length(proname)*length(x)))
      print(name)
    }

  }

  end_time <- Sys.time()  ## 记录程序结束时间
  ## 第三个位置关闭进度条
  close(pb)
  print(end_time - star_time)
}
