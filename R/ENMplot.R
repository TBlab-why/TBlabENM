
#' @title Draw pictures of ENM related results
#' @description 对相关结果进行分类的单张图片绘制，比如适生区图、避难所图、丰富度图等。
#'
#' @param wordmap 世界地图
#' @param regionmap 区域地图，比如研究区域的地图
#' @param color 指定不同适生区等级的颜色。默认5个等级。
#' @param legend 逻辑值，是否需要图例，默认添加图例
#' @param path 结果保存路径
#' @param format 结果保存格式,(e.g. png), or one of "eps", "ps", "tex" (pictex),
#'  "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param raster 栅格文件完整路径
#' @param rcl 向量，用于指定重分类的区间，
#' @param Legend_title 图例标题
#' @param filename 保存的文件名称
#' @param label 标注文字，例如时期的标注
#'
#' @return picture
#' @export
#' @importFrom sf st_read
#' @importFrom terra classify
#' @importFrom terra as.data.frame
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom grid rectGrob
#' @importFrom grid gpar
#' @importFrom grid unit
#' @importFrom ggplot2 annotate
#' @importFrom ggspatial annotation_scale
#' @import tidyverse
#'
#' @examples
#' #绘图
#' ENMplot(raster = "G:/Pittosporum/TBlabENM/sr_refuges_future-past/spp-refuges.tif",
#'         rcl = c(0, 5, 10, 15, 20, 24),
#'         color = NULL,
#'         wordmap = "D:/Desktop/ENMdata/worldmap/带投影坐标系/world_pro.shp",
#'         regionmap = "D:/Desktop/ENMdata/Chinamap/2022年省界/sheng2022.shp",
#'         label = NULL,
#'         legend = TRUE,
#'         Legend_title = "Value",
#'         format = "jpeg",
#'         filename = "2.jpeg",
#'         path = "G:/Pittosporum")



ENMplot <- function(raster,
                    rcl = NULL,
                    color = NULL,
                    wordmap = NULL,
                    regionmap = NULL,
                    label = NULL,
                    legend = TRUE,
                    Legend_title = "Value",
                    format = "jpeg",
                    filename = "myplot.jpeg",
                    path) {
  ####环境准备
  ##windows环境R如果无法识别中文路径时运行,最好不要有中文，也不要运行这行
  #Sys.setlocale(category = "LC_ALL", locale="German")
  #加载times new Roman 字体
  #sysfonts::font_add('TNR','C:/Windows/Fonts/times.ttf')
  #showtext::showtext_auto()

  ####加载地图
  if (is.null(wordmap)==FALSE) {
    wordmap1 <- sf::st_read(wordmap)
  }
  if (is.null(regionmap)==FALSE) {
    regionmap1 <- sf::st_read(regionmap)
  }


  ###读取模拟结果文件夹（以物种为单位）
  #设置默认颜色
  if (is.null(color)) {
    color <- c("1" = rgb(242,242,242,maxColorValue = 255),
               "2" = rgb(61,230,0,maxColorValue = 255),
               "3" = rgb(16,191,0,maxColorValue = 255),
               "4" = rgb(0,152,5,maxColorValue = 255),
               "5" = rgb(0,128,0,maxColorValue = 255))}

  #读取栅格数
  raster1 <- terra::rast(raster)
  #重分类
  raster1_classify <- terra::classify(raster1, rcl, include.lowest=TRUE, brackets=TRUE)
  #很关键的一步：将栅格数据转成可供ggplot2绘制的数据格式
  pro_df <- terra::as.data.frame(raster1_classify, xy = TRUE)
  #赋值列名
  colnames(pro_df) <- c("x", "y", "level")
  #重命名因子水平
  pro_df$class <- as.factor(as.numeric(pro_df$level))

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
    legend.position = c(0.88,0.21),
    # legend.direction = "vertical",
    #  legend.margin = margin(t = 0, r = 2, b = 5, l = 0, unit = "cm"),
    #图例文本
    legend.text = element_text(family = "TNR", size = 40, margin = margin(l= 10)),
    #图例标题
    legend.title = element_text(family = "TNR", size = 50, margin = margin(b=11, l=-8)),

    #图例背景
    legend.background = element_blank(),
    #图例键
    legend.key = element_rect(color = NA, fill = NA),
    #图例键高度和宽度
    #legend.key.height = unit(0.8, "cm"),
    # legend.key.width = unit(10.5, "cm"),
    legend.key.size = unit(0.8, "cm"))
  # p_d1为有图例的参数
  #生成图例标签
  lab <- c()
  for (i in 1:length(rcl)) {
    lab1 <- paste0(rcl[i]+1, "-", rcl[i+1])
    lab <- c(lab, lab1)
  }
  lab <- lab[1:length(lab) - 1]
  if(lab[1] == "1-1") {lab[1] = "0"
                       lab[2] = paste0(rcl[2], "-", rcl[3])}
  if (is.null(regionmap)&is.null(wordmap)&legend==TRUE) {
    p_d1 <- ggplot2::ggplot() +

      #自定义离散颜色
      scale_fill_manual(
        name = Legend_title,
        labels = lab,
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
        scale_fill_manual(name = Legend_title,
                          labels = lab,
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

  at <- ggplot2::annotate("text", x = (min(pro_df$x)+xrange*0.27), y = (max(pro_df$y)-yrange*0.12), label = label,
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
                                     filename = filename,
                                     path = path,
                                     device = format, scale=1, dpi = "print" ,units = "cm",
                                     width = 18,
                                     height = 18/ratio )} else {
                                       ggplot2::ggsave(plot = p_d, filename = filename,
                                                       path = path,
                                                       device = format, scale=1, dpi = "print" ,units = "cm",
                                                       width = 18,
                                                       height = 18/ratio ) }

}
