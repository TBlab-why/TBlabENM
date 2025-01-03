
#' @title Statistical Variables According to Contribution Size and Plot
#' @description 本函数根据环境变量对模型的贡献大小进行排序，并统计对模型贡献最大的前三个变量的个数
#' 绘制堆叠条形统计图。
#' @param parameters 数据框，第一列为物种名(要与resultdir下的名称保持一致),第二列为坐标点数,第三列为使用的环境变量
#' @param x 数值向量,表示物种序号,与parameters中的物种序号一致。
#' @param wordsize 用于调整图的标签字体大小，数值型。
#' @param resultdir maxent结果文件路径。
#' @param device 输出图片格式，or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param outdir 输出文件路径
#'
#' @return 生成"./Variable_contribution"文件夹，该文件夹包括一个按贡献大小排列变量的.csv文件
#' 和.jpeg的最重要的三个变量贡献堆叠条形统计图。
#' @export
#' @importFrom utils read.csv write.csv
#' @importFrom dplyr arrange mutate
#' @importFrom stringr str_split str_count
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot ggsave geom_bar
#' @importFrom purrr map_dbl map map2 map2_chr
#' @importFrom ggsci scale_fill_npg
#' @importFrom stats aggregate
#'
#' @examples
#' maxent_envcb(parameters = read.csv("F:/4/TBlabENM/maxent/parameters.csv"),
#'              x = 1:10,
#'              resultdir = "F:/4/TBlabENM/maxent",
#'              wordsize = 1,
#'              device = "jpeg",
#'              outdir = "F:/4")
maxent_envcb <- function(parameters,
                         x = NULL,
                         resultdir,
                         wordsize = 1,
                         device = "jpeg",
                         outdir = NULL) {
  #检查输入变量
  if (is.null(x)) {
    x <- 1:nrow(parameters)
  }
  if (file.exists(resultdir) == FALSE) {
    stop("Resultdir not find.")
  }
  #检查x值是否正确
  if (max(x) > nrow(parameters)) {
    stop("x exceeds the total number of species")
  }
  if (is.null(outdir)) {
    outdir = "."
  }
  #选择要统计作图的物种
  spdata <- parameters[x, 1:3]
  names(spdata) <- c("species", "number", "env")
  #新建空数据框保存结果
  maxselectev <- data.frame(matrix(ncol = 5, nrow = length(x)))
  names(maxselectev) <- c(
    "species",
    "occurrence",
    "Variable contribution size",
    "Variable contribution size1",
    "Variable contribution size2"
  )

  ss <- spdata |>
    mutate(nenv = map_dbl(
      .x = env,
      .f = function(x) {
        str_count(x, ",") + 1
      }
    )) |>
    mutate(path = map_chr(
      .x = species,
      .f = function(x) {
        list.files(
          paste0(resultdir, "/", x),
          full.names = TRUE,
          pattern = "maxentResults.csv$" ,
          recursive = TRUE
        )[1]
      }
    )) |>
    mutate(par = map(
      .x = path,
      .f = function(x) {
        read.csv(x)
      }
    )) |>
    #提取变量重要性
    mutate(ev_cb = map2(
      .x = par,
      .y = nenv,
      .f = function(x, y) {
        dt <- as.data.frame(t(x[nrow(x), 12:(11 + y), drop = FALSE]))
        names(dt) <- "value"
        dplyr::arrange(dt, dplyr::desc(value))
      }
    )) |>
    mutate(var = map(
      .x = ev_cb,
      .f = function(x) {
        site <- unlist(stringr::str_split(rownames(x), "\\."))
        site <- site[seq(1, length(site), 2)]
      }
    )) |>
    mutate(v1 = map2_chr(
      .x = var,
      .y = ev_cb,
      .f = function(x, y) {
        paste(paste0(x, "(", y$value, "%)", collapse = ","))
      }
    )) |>
    mutate(v2 = map_chr(
      .x = var,
      .f = function(x) {
        paste(paste0(x, collapse = ","))
      }
    )) |>
    mutate(v3 = map2_chr(
      .x = var,
      .y = ev_cb,
      .f = function(x, y) {
        paste(paste0(x, " ", y$value, collapse = ","))
      }
    ))

  maxselectev$species <- ss$species
  maxselectev$occurrence <- ss$number
  maxselectev$`Variable contribution size` <- ss$v1
  maxselectev$`Variable contribution size1` <- ss$v2
  maxselectev$`Variable contribution size2` <- ss$v3

  #保存本地
  dir.create(paste0(outdir, "/Variable_contribution"),
             showWarnings = FALSE)
  utils::write.csv(
    maxselectev[1:3],
    row.names = FALSE,
    paste0(outdir, "/Variable_contribution/Variable.csv")
  )

  ####变量贡献条形统计图之按贡献最大的前三个变量排序
  #取出前三重要变量
  maxselectev1 <- data.frame(matrix(ncol = 3, nrow = nrow(maxselectev)))
  names(maxselectev1) <- c("First", "Second", "Thirdly")
  for (i in 1:nrow(maxselectev)) {
    maxselectev1[i, ] <- stringr::str_split(maxselectev[, 4], ",")[[i]][1:3]
  }

  ##
  data1 <- tidyr::pivot_longer(
    maxselectev1,
    #需要转置的数据
    c(First, Second, Thirdly),
    #需要转置的列，这里是除了X列其他都转置
    names_to = "Contribution",
    #转置后原本的列名会变成观测生成新列，指定这列的列名(列名包含了分类信息)
    values_to = "Variate"
  )  #转置后原本的各列的值会合并成一个新列，指定这列的列名

  #统计每个变量出现次数，按照次数指定因子水平
  times <- table(data1$Variate)

  data1$Variate <- factor(data1$Variate, levels =  rev(names(times[order(times)])))
  p <- ggplot2::ggplot(data1) +
    geom_bar(aes(Variate, fill = Contribution), #反转条形的堆积顺序
             position = position_stack(reverse = TRUE)) +
    #反转图例顺序
    guides(fill = guide_legend(
      reverse = TRUE,
      title.theme = element_text(family = "serif", size = wordsize * 15)
    )) +
    theme_bw() + ylab("Frequency") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        family = "serif"
      ),
      axis.text = element_text(family = "serif", size = wordsize * 15),
      axis.title = element_text(family = "serif", size = wordsize *
                                  20),
      legend.text = element_text(family = "serif", size = wordsize *
                                   12)
    ) +
    ggsci::scale_fill_npg()

  ggplot2::ggsave(
    plot = p,
    filename = paste0("Variable contribution statistics.", device),
    path = paste0(outdir, "/Variable_contribution"),
    scale = 1,
    dpi = "print",
    width = 10,
    height = 7
  )



  ####变量贡献条形统计图之按变量总贡献排序
  #
  maxselectev2 <- sort(unlist(stringr::str_split(maxselectev[, 5], ",")))
  maxselectev3 <- unlist(stringr::str_split(maxselectev2, " "))
  ee <- data.frame(matrix(
    data = maxselectev3,
    ncol = 2,
    nrow = length(maxselectev3) / 2,
    byrow = TRUE
  ))
  names(ee) <- c("Variate", "Contribution value")
  ee[2] <- as.numeric(ee[, 2])
  ee1 <- stats::aggregate(ee[, 2], list(Variate = ee$Variate), sum)
  ee1[, 3] <- ee1$x / (100 * length(x))
  names(ee1) <- c("Variate", "Contribution value", "Percent contribution")
  write.csv(
    ee1,
    paste0(
      outdir,
      "/Variable_contribution/Variable Percent contribution.csv"
    ),
    row.names = FALSE
  )
  p <- ggplot2::ggplot(ee1) +
    geom_col(aes(
      reorder(Variate, `Percent contribution`, decreasing = TRUE),
      `Percent contribution`,
      fill = Variate
    )) +
    scale_y_continuous(labels = scales::percent) + xlab("Variate") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        family = "serif"
      ),
      axis.text = element_text(family = "serif", size = wordsize * 15),
      axis.title = element_text(family = "serif", size = wordsize *
                                  20)
    )

  ggplot2::ggsave(
    plot = p,
    filename = paste0("Variable Percent contribution.", device),
    path = paste0(outdir, "/Variable_contribution"),
    scale = 1,
    dpi = "print",
    width = 10,
    height = 7
  )
}
