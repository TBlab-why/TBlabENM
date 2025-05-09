#' @title MaxEnt Models Parameters
#' @description 本函数生成MaxEnt模型的默认参数(TBlab)，您也可以自己指定。常与其他
#'    建模函数连用.
#' @param replicates 模型重复次数。默认10次。
#' @param betamultiplier 正则化乘数。默认为1。
#' @param pictures 是否生成模拟图。默认为true。
#' @param outputgrids 是否为每次模拟生成栅格结果。默认为false。
#' @param l 是否在模型中使用特征函数"linear"。默认为true。
#' @param q 是否在模型中使用特征函数"quadratic"。默认为true。
#' @param p 是否在模型中使用特征函数"product"。默认为false。
#' @param t 是否在模型中使用特征函数"threshold"。默认为false。
#' @param h 是否在模型中使用特征函数"hinge"。默认为false。
#' @param responsecurves 是否生成响应曲线图。默认为true。
#' @param replicatetype 指定重复类型。默认为"crossvalidate"。
#' @param jackknife 是否进行折刀分析。默认为true。
#' @param outputformat 选择输出格式。默认为"logistic"。

#' @return 包含MaxEnt建模参数的字符型向量。
#' @export
#' @author 吴海洋
#' @examples
#' # 使用默认参数(TBlab)
#' args <- maxent_args()
#' args
#' # 自定义参数。
#' args <- maxent_args(
#'   replicates = 5,
#'   betamultiplier = 2
#' )
#' args
maxent_args <- function(replicates = 10,
                        betamultiplier = 1,
                        l = TRUE,
                        q = TRUE,
                        p = FALSE,
                        t = FALSE,
                        h = FALSE,
                        replicatetype = "crossvalidate",
                        responsecurves = TRUE,
                        jackknife = TRUE,
                        pictures = TRUE,
                        outputgrids = FALSE,
                        outputformat = "logistic") {
  c(
    paste0("replicates=", replicates),
    paste0("betamultiplier=", betamultiplier), # 重复次数和正则化乘数
    paste0("linear=", l), # 5种特征函数
    paste0("quadratic=", q),
    paste0("product=", p),
    paste0("threshold=", t),
    paste0("hinge=", h),
    paste0("replicatetype=", replicatetype), # 重复类型
    paste0("responsecurves=", responsecurves), # 响应曲线
    paste0("jackknife=", jackknife), # 折刀分析
    paste0("pictures=", pictures),
    paste0("outputgrids=", outputgrids),
    paste0("outputformat=", outputformat) # 输出文件格式
  )
}
